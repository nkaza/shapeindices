## -- 1. constrained Delaunay triangulation --------------------------------

#' Triangulate one simple POLYGON (may contain holes) via sfdct.
#' @param part a single POLYGON sfg (one part of a possibly-multi geometry)
#' @param crs the CRS to attach to the output (sfdct doesn't preserve it)
#' @return sfc of triangle POLYGONs, or NULL if triangulation produced none
#' @noRd
.cdt_part <- function(part, crs) {
    part_sf <- st_sf(geometry = st_sfc(part, crs = crs))
    tri <- tryCatch(sfdct::ct_triangulate(part_sf), error = function(e) NULL)
    if (is.null(tri) || nrow(tri) == 0) return(NULL)
    tri <- st_collection_extract(st_geometry(tri), "POLYGON")

    # keep only triangles genuinely inside the source ring (guards against
    # leakage across holes, numerical edge cases)
    cen  <- suppressWarnings(st_centroid(tri))
    keep <- lengths(st_intersects(cen, part_sf)) > 0
    tri[keep]
}

# Splits a (MULTI)POLYGON sfc of length 1 into its POLYGON parts without
# st_cast(), which is unreliable here (post-st_make_valid() geometry can be
# tagged "GEOMETRY"/"GEOMETRYCOLLECTION" even when the content is a plain
# POLYGON, and st_cast() dispatches on that tag). POLYGON/MULTIPOLYGON sfg
# objects are plain nested lists, so parts are pulled out directly instead.
#' @param poly a single (length-1) sfc (MULTI)POLYGON, normally already
#'   st_make_valid()'d
#' @return an sfc of POLYGON parts, or NULL if nothing polygonal is left
#'   after cleanup
#' @noRd
.split_polygon_parts <- function(poly) {
    if (length(poly) == 0 || any(st_is_empty(poly))) return(NULL)

    crs  <- st_crs(poly)
    g    <- poly[[1]]
    type <- class(g)[2]

    if (type == "GEOMETRYCOLLECTION") {
        poly <- st_collection_extract(poly, "POLYGON", warn = FALSE)
        if (length(poly) == 0 || any(st_is_empty(poly))) return(NULL)
        g    <- poly[[1]]
        type <- class(g)[2]
    }

    if (type == "POLYGON") {
        return(st_sfc(g, crs = crs))
    }
    if (type == "MULTIPOLYGON") {
        return(st_sfc(lapply(unclass(g), sf::st_polygon), crs = crs))
    }
    NULL  # POINT/LINESTRING/etc - no polygonal area left
}

#' Constrained Delaunay triangulation of a (multi)polygon
#'
#' Triangulates each part separately so disjoint parts never bridge, and
#' hole boundaries are always respected as constraints.
#' @param poly a single sfg/sfc (MULTI)POLYGON (length 1); does not run
#'   st_make_valid() itself - callers are expected to have cleaned it up
#' @return an sf data frame with one row per triangle (tri_id, area,
#'   geometry), or NULL if the polygon has no triangulatable area
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' tri <- cdt_triangles(nc[nc$NAME == "Wake", ])
#' nrow(tri)
#' plot(sf::st_geometry(tri), border = "grey40")
#' @export
cdt_triangles <- function(poly) {
    poly <- st_geometry(poly)
    stopifnot(length(poly) == 1)

    parts <- .split_polygon_parts(poly)
    if (is.null(parts)) return(NULL)

    crs   <- st_crs(parts)
    tri_list <- lapply(parts, .cdt_part, crs = crs)
    tri_list <- Filter(Negate(is.null), tri_list)
    if (length(tri_list) == 0) return(NULL)
    tri <- do.call(c, tri_list)
    st_sf(tri_id = seq_along(tri), area = as.numeric(st_area(tri)),
          geometry = tri, crs = crs)
}

#' Validate and triangulate a polygon once
#'
#' Builds the `prep` object every mesh-based index accepts, so one
#' triangulation can be shared across several index calls on the same
#' polygon (see shape_indices(), which does this automatically).
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param simplify_tolerance optional st_simplify() tolerance. NULL
#'   (default) skips simplification entirely. Applied after cleaning,
#'   before triangulating - trades a small, bounded amount of boundary
#'   detail for a smaller mesh (fewer CDT triangles and fewer total
#'   boundary edges, both of which reduce every downstream index's cost
#'   and memory use). Always in the geometry's OWN CRS unit at the point
#'   simplification runs, NOT necessarily metres: geographic (lon/lat)
#'   input gets auto-projected to a metric CRS first (always metres), but
#'   input already in a projected CRS keeps whatever linear unit that CRS
#'   uses - many US State Plane CRSs are US survey feet, for example.
#'   Warns if that unit isn't metres, since a silently-wrong unit makes
#'   the same tolerance value mean a very different amount of
#'   simplification than intended.
#' @return list(poly, tri) - cleaned, planar geometry and its triangle mesh
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' prep <- prepare_polygon(nc[nc$NAME == "Wake", ])
#' nrow(prep$tri)
#' @export
prepare_polygon <- function(poly, simplify_tolerance = NULL) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    if (!is.null(simplify_tolerance)) {
        .warn_if_nonmetric_tolerance(poly, simplify_tolerance)
        poly <- st_simplify(poly, dTolerance = simplify_tolerance, preserveTopology = TRUE)
        poly <- .make_valid_warn(poly)
    }
    list(poly = poly, tri = cdt_triangles(poly))
}
