## -- 1c. adaptive mesh refinement (opposite direction from convex_decompose) -
##
## Medial subdivision of every CDT triangle, area-adaptively: the mesh's
## own largest triangle gets `max_depth` levels (4^max_depth children);
## every smaller triangle gets only as many levels as it needs to reach
## that same final sub-triangle size (see .adaptive_tri_depth() in
## utils.R, shared with radial_concentration_index()'s own internal
## quadrature refinement). More, smaller triangles, not fewer, larger
## convex pieces - the opposite of convex_decompose().

#' Adaptive subdivision of a CDT mesh into a finer triangle mesh.
#'
#' Standalone utility, not currently called by any of this package's own
#' index functions - `radial_concentration_index()` does the same
#' area-adaptive subdivision internally, but collapses each sub-triangle
#' straight to a weighted centroid rather than keeping triangle geometries
#' around, so it doesn't go through this function. See `convex_decompose()`
#' for the opposite direction (fewer, larger convex pieces rather than
#' more, smaller triangles).
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating
#' @param max_depth subdivision depth ceiling, applied to the mesh's own
#'   largest triangle only (4^max_depth children for that triangle; fewer
#'   for smaller ones, adaptively - see `.adaptive_tri_depth()`)
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return an sf data frame shaped like `cdt_triangles()`'s output (tri_id,
#'   area, geometry), or NULL if the polygon triangulates to no pieces
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' tri  <- cdt_triangles(nc[nc$NAME == "Dare", ])
#' fine <- subdivide_mesh(nc[nc$NAME == "Dare", ])
#' nrow(tri)    # coarse CDT triangles
#' nrow(fine)   # more, smaller triangles
#' @export
subdivide_mesh <- function(poly, prep = NULL, max_depth = 4, simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)
    tri <- prep$tri
    n_tri <- if (is.null(tri)) 0 else nrow(tri)
    if (n_tri == 0) return(NULL)
    crs <- st_crs(tri)

    if (max_depth == 0) {
        return(st_sf(tri_id = seq_len(n_tri), area = tri$area,
                      geometry = st_geometry(tri), crs = crs))
    }

    tri_area <- tri$area
    corner_mat <- vapply(st_geometry(tri), function(g) {
        v <- st_coordinates(g)[1:3, 1:2, drop = FALSE]
        c(v[1, ], v[2, ], v[3, ])
    }, numeric(6))
    A_all <- t(corner_mat[1:2, , drop = FALSE])
    B_all <- t(corner_mat[3:4, , drop = FALSE])
    C_all <- t(corner_mat[5:6, , drop = FALSE])

    depth_i <- .adaptive_tri_depth(tri_area, max_depth)

    A_list <- vector("list", max_depth + 1)
    B_list <- vector("list", max_depth + 1)
    C_list <- vector("list", max_depth + 1)
    for (d in 0:max_depth) {
        grp <- which(depth_i == d)
        if (length(grp) == 0) next
        A <- A_all[grp, , drop = FALSE]
        B <- B_all[grp, , drop = FALSE]
        C <- C_all[grp, , drop = FALSE]
        for (k in seq_len(d)) {
            sub <- .subdivide_tri_batch(A, B, C)
            A <- sub$A; B <- sub$B; C <- sub$C
        }
        A_list[[d + 1]] <- A; B_list[[d + 1]] <- B; C_list[[d + 1]] <- C
    }
    A <- do.call(rbind, A_list); B <- do.call(rbind, B_list); C <- do.call(rbind, C_list)
    n_sub <- nrow(A)

    if (requireNamespace("sfheaders", quietly = TRUE)) {
        poly_df <- data.frame(
            polygon_id = rep(seq_len(n_sub), each = 4),
            x = as.vector(t(cbind(A[, 1], B[, 1], C[, 1], A[, 1]))),
            y = as.vector(t(cbind(A[, 2], B[, 2], C[, 2], A[, 2])))
        )
        geoms <- sfheaders::sfc_polygon(poly_df, x = "x", y = "y", polygon_id = "polygon_id")
        st_crs(geoms) <- crs
    } else {
        geoms <- st_sfc(lapply(seq_len(n_sub), function(i) {
            st_polygon(list(rbind(A[i, ], B[i, ], C[i, ], A[i, ])))
        }), crs = crs)
    }
    st_sf(tri_id = seq_len(n_sub), area = as.numeric(st_area(geoms)), geometry = geoms)
}
