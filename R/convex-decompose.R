## -- 1b. convex decomposition (approximate alternative to CDT) ------------
##
## Hertel-Mehlhorn: merges adjacent CDT triangles across a shared diagonal
## whenever both endpoints stay convex. Convexity is tracked as a running
## sum of interior angles per (piece, vertex) via union-find - no GEOS
## calls until pieces are dissolved at the end - so this is O(n) plus
## O(n log n) union-find, not the O(k^2)+ cost of literally re-testing
## st_union/st_convex_hull per merge.
##
## Not guaranteed-minimal (Hertel-Mehlhorn bounds piece count at 4x
## optimal, not optimal itself), and can hide reflex vertices inside one
## piece's own boundary the same way CDT can (see convexity_index()'s
## known limitation) - fewer pairs to test, not a more correct index.

#' Convex decomposition of a (multi)polygon (Hertel-Mehlhorn).
#'
#' Not currently called by any of this package's own index functions -
#' every index works directly off `cdt_triangles()`'s triangles (optionally
#' refined finer by `subdivide_mesh()`), never off merged convex pieces.
#' This is a standalone utility for a caller who wants an approximate
#' convex decomposition for their own purposes; see `subdivide_mesh()` for
#' the opposite direction (more, smaller triangles rather than fewer,
#' larger convex pieces).
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return an sf data frame shaped like cdt_triangles()'s output (piece_id,
#'   area, geometry), or NULL if the polygon triangulates to no pieces
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' tri    <- cdt_triangles(nc[nc$NAME == "Dare", ])
#' pieces <- convex_decompose(nc[nc$NAME == "Dare", ])
#' nrow(tri)     # many small triangles
#' nrow(pieces)  # fewer, larger convex pieces
#' @export
convex_decompose <- function(poly, prep = NULL, simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)
    poly <- prep$poly
    tri  <- prep$tri
    n_tri <- if (is.null(tri)) 0 else nrow(tri)
    if (n_tri == 0) return(NULL)
    crs <- st_crs(poly)

    if (n_tri == 1) {
        return(st_sf(piece_id = 1L, area = tri$area, geometry = st_geometry(tri), crs = crs))
    }

    # corner coordinates -> deduplicated integer vertex ids (CDT corners
    # come from the polygon's own vertex set, so exact-match works)
    coords_list <- lapply(st_geometry(tri), function(g) st_coordinates(g)[1:3, 1:2, drop = FALSE])
    all_coords  <- do.call(rbind, coords_list)
    key      <- sprintf("%.9f_%.9f", all_coords[, 1], all_coords[, 2])
    uniq_key <- unique(key)
    vid      <- match(key, uniq_key)
    vid_mat  <- matrix(vid, ncol = 3, byrow = TRUE)          # n_tri x 3 corner ids
    vcoord   <- all_coords[match(uniq_key, key), , drop = FALSE]

    # interior angle of every triangle at each of its 3 corners
    ang <- matrix(0, n_tri, 3)
    for (t in seq_len(n_tri)) {
        p <- vcoord[vid_mat[t, ], , drop = FALSE]
        for (k in 1:3) {
            a  <- p[k, ]; b <- p[(k %% 3) + 1, ]; cc <- p[((k + 1) %% 3) + 1, ]
            u  <- b - a;  v  <- cc - a
            lu <- sqrt(sum(u^2)); lv <- sqrt(sum(v^2))
            ang[t, k] <- if (lu < 1e-12 || lv < 1e-12) 0 else
                acos(max(-1, min(1, sum(u * v) / (lu * lv))))
        }
    }

    # edge -> triangle-pair map (internal diagonals are shared by exactly 2
    # triangles; boundary edges by exactly 1)
    edge_key <- function(a, b) paste(pmin(a, b), pmax(a, b))
    all_edges   <- c(edge_key(vid_mat[, 1], vid_mat[, 2]),
                      edge_key(vid_mat[, 2], vid_mat[, 3]),
                      edge_key(vid_mat[, 3], vid_mat[, 1]))
    all_tri_idx <- rep(seq_len(n_tri), 3)
    all_vpair   <- rbind(vid_mat[, 1:2], vid_mat[, 2:3], vid_mat[, c(3, 1)])
    internal    <- Filter(function(idx) length(idx) == 2, split(seq_along(all_edges), all_edges))

    # union-find over triangles, tracking a running interior-angle total per
    # (current-root, vertex) - merge test for diagonal (u,v): does the
    # angle sum at u, and at v, stay <= pi? O(1) per diagonal.
    parent <- seq_len(n_tri)
    find <- function(x) {
        while (parent[x] != x) { parent[x] <<- parent[parent[x]]; x <- parent[x] }
        x
    }
    angle_store <- new.env(parent = emptyenv())
    for (t in seq_len(n_tri)) for (k in 1:3) {
        assign(paste0(t, "_", vid_mat[t, k]), ang[t, k], envir = angle_store)
    }
    get_angle <- function(root, v) {
        ky <- paste0(root, "_", v)
        if (exists(ky, envir = angle_store, inherits = FALSE)) get(ky, envir = angle_store) else NA_real_
    }

    for (e in internal) {
        t1 <- all_tri_idx[e[1]]; t2 <- all_tri_idx[e[2]]
        u  <- all_vpair[e[1], 1]; v <- all_vpair[e[1], 2]
        r1 <- find(t1); r2 <- find(t2)
        if (r1 == r2) next   # already joined via another path (e.g. around a hole)
        a1u <- get_angle(r1, u); a2u <- get_angle(r2, u)
        a1v <- get_angle(r1, v); a2v <- get_angle(r2, v)
        if (anyNA(c(a1u, a2u, a1v, a2v))) next
        new_u <- a1u + a2u
        new_v <- a1v + a2v
        if (new_u <= pi + 1e-9 && new_v <= pi + 1e-9) {
            parent[r2] <- r1
            assign(paste0(r1, "_", u), new_u, envir = angle_store)
            assign(paste0(r1, "_", v), new_v, envir = angle_store)
        }
    }

    # dissolve accepted groups into real geometries, once
    group      <- vapply(seq_len(n_tri), find, integer(1))
    grp_ids    <- sort(unique(group))
    piece_geom <- do.call(c, lapply(grp_ids, function(g) st_union(st_geometry(tri)[group == g])))
    areas      <- as.numeric(st_area(piece_geom))
    st_sf(piece_id = seq_along(grp_ids), area = areas, geometry = piece_geom, crs = crs)
}
