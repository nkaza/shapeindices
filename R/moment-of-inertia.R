## =========================================================================
## 4. Polar moment-of-inertia index (uses the same CDT triangulation)
## =========================================================================
##
## J = Ixx + Iyy, the polar second moment of area about the shape's own
## mass centroid (rho-weighted average of triangle centroids - see
## .moment_of_inertia_core() for why that's the correct reference point
## for any density, not just the area centroid), summed triangle-by-
## triangle via the standard shoelace formula. index = A^2 / (2*pi*J).
##
## Rests on a real theorem, not a heuristic: for fixed area, a disk
## centred at its own centroid uniquely minimises J among all measurable
## planar sets (bathtub principle - J weights area by squared distance
## from the centroid). index in (0, 1], = 1 iff (almost everywhere) a
## disk. Penalises dispersal more sharply than convexity_index(), but
## conflates "not round" with "not convex" (a thin convex rectangle
## scores low too).
##
## WEIGHTED VERSION: a per-triangle weight makes density rho_i = W_i/A_i
## piecewise-constant instead of uniform, so the disk-minimises-J theorem
## no longer applies - the minimiser for a given density profile is
## concentric rings, densest at the centre (same rearrangement-inequality
## family). For parcels sorted by descending density, with S_k =
## cumulative area of the k densest:
##
##   J_ref = (1 / 2*pi) * sum_k rho_(k) * (S_k^2 - S_(k-1)^2)
##
## index = J_ref / J_actual, again in (0, 1]. weight = NULL (uniform
## density) collapses this exactly to the original A^2/(2*pi). `weight`
## is normalised to sum to 1 before use; `total_weight` in the return
## value reports the caller's original, un-normalised sum instead.

#' Polar moment-of-inertia compactness/dispersal index of a (multi)polygon.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, used as that triangle's mass instead of its own area.
#'   NULL (default) reproduces the unweighted index exactly. See the
#'   concentric-rings note above for why the reference changes once
#'   density is non-uniform.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return list(index, J, Ixx, Iyy, Ixy, J_ref, area, total_weight,
#'   centroid, triangles). `centroid` is the plain geometric centroid when
#'   weight = NULL, or the true mass centroid when weighted - J/Ixx/Iyy/Ixy
#'   are always computed about this point (`Ixy`, the product of inertia,
#'   isn't used by this index at all - it's here for
#'   `moment_isotropy_index()`, which shares this same triangulation and mass
#'   centroid). `total_weight` is sum(weight) as supplied (before
#'   normalisation), or area when weight = NULL.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' moment_of_inertia_index(wake)$index
#' # Camden: a long, narrow county - convex, but far from a disk
#' moment_of_inertia_index(nc[nc$NAME == "Camden", ])$index
#'
#' # weight: substitutes for each triangle's own mass - weighting by the
#' # triangle's own area exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' moment_of_inertia_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
moment_of_inertia_index <- function(poly, prep = NULL, weight = NULL, simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)
    poly <- prep$poly
    tri  <- prep$tri
    n <- if (is.null(tri)) 0 else nrow(tri)
    if (n == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, J = NA_real_, Ixx = NA_real_, Iyy = NA_real_, Ixy = NA_real_,
                    J_ref = NA_real_, area = NA_real_, total_weight = NA_real_,
                    centroid = NULL, triangles = tri))
    }
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }

    tri_area  <- tri$area
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri_area else .normalize_weight(weight)
    rho       <- w / tri_area   # density: 1 everywhere iff weight = NULL

    res <- .moment_of_inertia_core(st_geometry(tri), tri_area, rho, poly)
    if (!is.null(raw_total)) res$total_weight <- raw_total
    res$triangles <- tri
    res
}

#' Shared core, also used by shape_indices_sf(byrow = FALSE) on a combined
#' multi-row triangle mesh. The reference point `G` is the mass centroid
#' (mass-weighted average of triangle centroids), not the plain geometric
#' one - computing Ixx/Iyy about the wrong point inflates J by
#' mass * distance^2 (parallel axis theorem) whenever density is
#' non-uniform. For uniform density this is identical to
#' st_centroid(st_union(poly)) (verified to floating-point precision), so
#' unweighted behaviour is unchanged.
#' @param tri_geom sfc of triangle POLYGON geometries
#' @param tri_area numeric vector, each triangle's own physical area
#' @param rho numeric vector, each triangle's density (weight / tri_area);
#'   uniformly 1 for the unweighted case
#' @param poly the (multi)polygon `tri_geom` triangulates, used only to
#'   attach a CRS to the returned centroid point
#' @return list(index, J, Ixx, Iyy, Ixy, J_ref, area, total_weight,
#'   centroid, triangles) - triangles is always NULL here, the caller
#'   attaches its own. `Ixy` (product of inertia, about the same mass
#'   centroid as Ixx/Iyy) is computed for moment_isotropy_index()'s benefit -
#'   moment_of_inertia_index() itself only uses Ixx/Iyy.
#' @noRd
.moment_of_inertia_core <- function(tri_geom, tri_area, rho, poly) {
    mass <- rho * tri_area
    G <- .mass_centroid(tri_geom, mass)

    moms <- vapply(tri_geom, function(t) {
        v <- st_coordinates(t)[1:3, 1:2, drop = FALSE]
        v[, 1] <- v[, 1] - G[1]
        v[, 2] <- v[, 2] - G[2]

        # ensure CCW winding for a consistent shoelace sign
        signed2A <- sum(v[, 1] * v[c(2, 3, 1), 2] - v[c(2, 3, 1), 1] * v[, 2])
        if (signed2A < 0) v <- v[c(1, 3, 2), ]

        x <- v[, 1]; y <- v[, 2]
        xn <- x[c(2, 3, 1)]; yn <- y[c(2, 3, 1)]
        cross <- x * yn - xn * y
        c(Ixx = sum((y^2 + y * yn + yn^2) * cross) / 12,
          Iyy = sum((x^2 + x * xn + xn^2) * cross) / 12,
          # standard polygon product-of-inertia formula (e.g. Eberly,
          # "Polygon Mass Properties") - same cross/winding setup as
          # Ixx/Iyy above, just the xy cross term instead of x^2/y^2
          Ixy = sum((x * yn + 2 * x * y + 2 * xn * yn + xn * y) * cross) / 24)
    }, numeric(3))

    Ixx <- sum(rho * moms["Ixx", ])
    Iyy <- sum(rho * moms["Iyy", ])
    Ixy <- sum(rho * moms["Ixy", ])
    J   <- Ixx + Iyy
    A   <- sum(tri_area)
    W   <- sum(rho * tri_area)

    # concentric-rings reference: sort by density descending, cumulative
    # area gives each triangle's annulus, closed-form J of that arrangement
    ord    <- order(rho, decreasing = TRUE)
    S      <- cumsum(tri_area[ord])
    S_prev <- c(0, S[-length(S)])
    J_ref  <- sum(rho[ord] * (S^2 - S_prev^2)) / (2 * pi)

    index <- if (J > 0) J_ref / J else NA_real_

    list(index = index, J = J, Ixx = Ixx, Iyy = Iyy, Ixy = Ixy, J_ref = J_ref,
         area = A, total_weight = W,
         centroid = st_sfc(st_point(G), crs = st_crs(poly)), triangles = NULL)
}

#' Vectorised wrapper for an sf data frame - each row indexed independently.
#' @param x an sf data frame
#' @param ... passed to moment_of_inertia_index() for every row
#' @return `x` with one new column, moi_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- moment_of_inertia_index_sf(nc[1:5, ])
#' res$moi_index
#' @export
moment_of_inertia_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    # geoms[i], not lapply(geoms, .) directly - the latter strips to a bare
    # sfg with no CRS of its own
    res <- lapply(seq_along(geoms), function(i) moment_of_inertia_index(geoms[i], ...))
    x$moi_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
