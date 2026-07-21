## =========================================================================
## Exact weighted mesh for shape_indices_sf(byrow = FALSE), via constrained
## Delaunay triangulation with row boundaries as constraint segments
## =========================================================================
##
## .weighted_mesh() (shape-indices.R) triangulates poly_u alone via sfdct,
## then area-weight-averages each row's density onto whichever triangles
## it overlaps via st_intersection() - a triangle straddling several rows
## gets their AVERAGE density baked in uniformly, smearing concentration
## and making the weighted index a function of the arbitrary triangulation,
## not just of (shape, weights). Verified on real Census block data: this
## overestimates compactness (moment_of_inertia_index()) by ~0.05-0.10
## index points on realistic urban-area-scale inputs.
##
## FIX: supply every kept row's own boundary as a CONSTRAINT SEGMENT to the
## triangulation itself (RTriangle::pslg()/triangulate() - the constrained-
## triangulation engine sfdct already wraps via ct_triangulate(), but
## doesn't expose segments through its own API; sfdct's own path-extraction
## internals are unexported, so this calls RTriangle directly rather than
## relying on them). A constrained Delaunay triangulation never produces a
## triangle crossing a constraint edge, so every triangle sits entirely
## inside exactly one row by construction - weight allocation is exact,
## not approximated, with no depth/resolution parameter to tune. Verified:
## weight recovery is exactly 100.00000% (to displayed precision) on three
## real Urban Area files, including a 38,175-row/1,979,609-triangle case,
## cross-checked against an independent implementation to the last digit.
##
## SCOPE: only fixes shape_indices_sf(byrow = FALSE)'s WEIGHTED mesh -
## never touches byrow = TRUE or unweighted byrow = FALSE, neither of
## which builds a per-row density field at all. .weighted_mesh() itself is
## unchanged and stays in the codebase as the automatic fallback whenever
## this path can't be built safely (overlapping input rows, or any
## RTriangle failure) - see .constrained_weighted_mesh()'s own comments
## for exactly when that triggers.
##
## RAW ARRAYS, NOT SF PIECES: the returned mesh is RTriangle's own output
## format directly (P: point matrix, T: 3-column triangle-vertex-index
## matrix) - never converted to sf POLYGON objects. Building ~2 million sf
## polygons (one per triangle, at real Urban-Area scale) measured at ~500
## seconds by itself, dominated by R's own per-object construction
## overhead, not triangulation cost (the triangulation itself takes ~2s at
## that scale). Every index that actually runs on this mesh - moment of
## inertia, moment isotropy, radial concentration, directional balance,
## and convexity/span's Monte Carlo sampling - only ever needs per-triangle
## centroids, areas, and corner coordinates, all read directly off P/T
## with no sf involved at all (see the array-native helpers below).
## Convexity/span's DETERMINISTIC mode is the one consumer that genuinely
## needs real polygon pieces (for the same-piece diagonal-term shortcut) -
## and that mode is exactly the one forced off above the memory-derived
## safety ceiling (see shape-indices.R's dispatch), so it never actually
## runs on this mesh in practice.

#' Vectorised sf-geometry -> Planar Straight Line Graph (PSLG) extraction:
#' every ring (exterior AND holes, across every part of a MULTIPOLYGON) of
#' every feature in `geoms` becomes a closed loop of constraint segments,
#' sharing vertices between touching features via coordinate-rounded
#' deduplication. One st_coordinates() call for the whole input, one
#' vectorised match() for vertex dedup, and a run-length-encoding-based
#' "close each ring" shift - no per-feature R loop. A naive per-point
#' environment-lookup version of this same dedup showed clear O(n^2)-ish
#' blowup at scale (successive 5k-feature chunks took 0.1, 0.2, 0.4, 0.9
#' minutes); this vectorised version processes a real 38,175-feature
#' input's worth of coordinates in ~5 seconds.
#'
#' st_coordinates() labels every row with ring/part/feature indices in
#' columns named "L1", "L2", ... - for a plain POLYGON input the last such
#' column is the feature index and any earlier ones are the ring within
#' that feature (1 = exterior, >1 = a hole); for MULTIPOLYGON there's one
#' more column, ring-within-part then part-within-feature then feature.
#' Grouping by every "L" column except the last (the feature index) gives
#' exactly one group per ring instance, whatever the geometry type - this
#' is what makes holes (an easy thing to silently drop by filtering only
#' `L1 == 1`, a real mistake caught mid-development via a real block with
#' a hole in real test data - area conservation immediately exposed the
#' dropped ring) and multi-part features both work correctly with the same
#' code, rather than needing separate handling.
#' @param geoms an sfc of POLYGON/MULTIPOLYGON geometries
#' @param snap_digits decimal places to round coordinates to before
#'   deduplicating vertices - touching features' shared boundary points
#'   must round to the same key to be treated as one shared vertex
#' @return list(P, S, n_rings) - P an Nx2 point matrix, S an Mx2 matrix of
#'   1-based indices into P (each row a constraint segment), n_rings the
#'   number of closed loops found (features + holes, informational only)
#' @noRd
.extract_pslg <- function(geoms, snap_digits = 6L) {
    co <- st_coordinates(geoms)
    lcols <- grep("^L", colnames(co), value = TRUE)
    feat_col  <- lcols[length(lcols)]
    ring_cols <- setdiff(lcols, feat_col)
    ring_key <- if (length(ring_cols) == 0) co[, feat_col] else
        interaction(co[, feat_col], do.call(paste, as.data.frame(co[, ring_cols, drop = FALSE])), drop = TRUE)
    grp_id <- as.integer(ring_key)

    key <- sprintf("%.*f_%.*f", snap_digits, co[, "X"], snap_digits, co[, "Y"])
    uk  <- unique(key)
    vid <- match(key, uk)
    ux  <- co[match(uk, key), "X"]
    uy  <- co[match(uk, key), "Y"]

    r <- rle(grp_id)
    end_idx   <- cumsum(r$lengths)
    start_idx <- end_idx - r$lengths + 1L
    is_close_dup <- seq_along(grp_id) %in% end_idx   # sf always repeats a ring's first point as its last

    from_rows <- which(!is_close_dup)
    next_row  <- from_rows + 1L
    grp_of_from  <- grp_id[from_rows]
    is_last_kept <- is_close_dup[next_row]
    start_by_grp <- integer(max(grp_id))
    start_by_grp[r$values] <- start_idx
    next_row[is_last_kept] <- start_by_grp[grp_of_from[is_last_kept]]

    list(P = cbind(ux, uy), S = cbind(vid[from_rows], vid[next_row]), n_rings = length(r$lengths))
}

#' Builds the exact constrained-triangulation weighted mesh for
#' shape_indices_sf(byrow = FALSE) - see this file's own header for why,
#' and .weighted_mesh() (shape-indices.R) for the coarse approach this
#' replaces (and falls back to - see `ok` below).
#'
#' Two conditions make constrained triangulation unsafe to attempt, and
#' both are checked before calling into RTriangle rather than left to
#' surface as an opaque failure: kept rows that genuinely overlap (not
#' just touch) violate the assumption that row boundaries form a valid
#' planar partition, since a constraint segment can't sensibly represent
#' two rows claiming the same territory; and any other RTriangle failure
#' (e.g. a genuinely degenerate/self-intersecting row boundary slipping
#' through st_make_valid()). Either signals `ok = FALSE`, and the caller
#' (.shape_indices_sf_grouped()) falls back to .weighted_mesh() with a
#' warning - never a hard error, since the coarse approach has always
#' tolerated exactly this kind of messy input.
#' @param x an sf data frame, already projected
#' @param weights NULL, a column-name string, or a numeric vector - see
#'   .resolve_row_weights()
#' @param simplify_tolerance see .resolve_union() - applied to poly_u's
#'   outer boundary only. Unlike the coarse .weighted_mesh() path, this no
#'   longer controls overall mesh size (that's now driven by the rows'
#'   own boundary detail, which is never simplified, since simplifying it
#'   would reintroduce exactly the smearing this mesh exists to avoid) -
#'   only the union's own outer perimeter detail.
#' @return list(ok, P, T, tri_area, tri_weight, poly_u, raw_total, crs).
#'   When `ok = FALSE`, only `ok` and `reason` are meaningful - the caller
#'   ignores every other field and calls .weighted_mesh() instead. P/T
#'   are RTriangle::triangulate()'s own output format directly (a point
#'   matrix and a 3-column triangle-vertex-index matrix) - see this
#'   file's own header for why these are never converted to sf polygons.
#' @noRd
.constrained_weighted_mesh <- function(x, weights, simplify_tolerance = NULL) {
    ru    <- .resolve_union(x, weights, simplify_tolerance)
    geoms <- ru$geoms
    w_row <- ru$w_row
    crs   <- st_crs(x)

    if (length(geoms) > 1 && any(lengths(suppressWarnings(st_overlaps(geoms))) > 0)) {
        return(list(ok = FALSE, reason = paste(
            "kept rows overlap (not just touch), which violates the",
            "planar-partition assumption constrained triangulation needs")))
    }

    pslg_res <- tryCatch(.extract_pslg(geoms), error = function(e) e)
    if (inherits(pslg_res, "error")) {
        return(list(ok = FALSE, reason = paste("boundary extraction failed:", conditionMessage(pslg_res))))
    }

    tr <- tryCatch(
        RTriangle::triangulate(RTriangle::pslg(P = pslg_res$P, S = pslg_res$S), Y = TRUE),
        error = function(e) e
    )
    if (inherits(tr, "error")) {
        return(list(ok = FALSE, reason = paste("constrained triangulation failed:", conditionMessage(tr))))
    }
    if (nrow(tr$T) == 0) {
        return(list(ok = FALSE, reason = "constrained triangulation produced no triangles"))
    }

    # assign each triangle to its row via centroid-in-polygon - exact,
    # since a triangle can never straddle a row boundary by construction;
    # triangles matching no row (holes' interiors, correctly triangulated
    # then discarded here - the same pattern sfdct's own ct_triangulate()
    # uses for holes) get dropped from the returned mesh entirely
    cx <- rowMeans(matrix(tr$P[t(tr$T), 1], ncol = 3, byrow = TRUE))
    cy <- rowMeans(matrix(tr$P[t(tr$T), 2], ncol = 3, byrow = TRUE))
    cent_sfc <- st_sfc(lapply(seq_along(cx), function(i) st_point(c(cx[i], cy[i]))), crs = crs)
    hits <- st_intersects(cent_sfc, geoms)
    row_idx <- vapply(hits, function(h) if (length(h)) h[1L] else NA_integer_, integer(1))

    keep_tri <- !is.na(row_idx)
    T_keep   <- tr$T[keep_tri, , drop = FALSE]
    row_idx  <- row_idx[keep_tri]

    tri_area <- vapply(seq_len(nrow(T_keep)), function(i) {
        v <- tr$P[T_keep[i, ], ]
        abs((v[2, 1] - v[1, 1]) * (v[3, 2] - v[1, 2]) - (v[3, 1] - v[1, 1]) * (v[2, 2] - v[1, 2])) / 2
    }, numeric(1))

    # density from RAW (pre-normalisation) row weight/area - the same
    # invariant .resolve_row_weights()/.resolve_union() already keep for
    # `raw_total`, kept explicit here since it's what makes each
    # triangle's weight exactly reproduce its own row's true density,
    # regardless of how many other rows/triangles exist
    row_area    <- as.numeric(st_area(geoms))
    row_density <- w_row / row_area
    tri_weight  <- tri_area * row_density[row_idx]

    list(ok = TRUE, P = tr$P, T = T_keep, tri_area = tri_area, tri_weight = tri_weight,
         poly_u = ru$poly_u, raw_total = ru$raw_total, crs = crs)
}

## -- array-native computation, used only by the mesh above --------------
##
## Every index that actually runs on this mesh only needs per-triangle
## centroids, areas, and corner coordinates - never a real sf POLYGON
## object (see this file's own header for the measured cost of building
## one per triangle at real scale). These three helpers mirror the exact
## same formulas the sf-based functions they stand in for already use
## (moment-of-inertia.R's .moment_of_inertia_core(), convexity-index.R's
## .sample_weighted_points()), just reading coordinates directly out of a
## triangulation's own P (point matrix) / T (triangle-vertex-index
## matrix) instead of round-tripping through sf geometries.

#' Array-native mirror of .moment_of_inertia_core() (moment-of-inertia.R) -
#' same mass centroid / Ixx/Iyy/Ixy / concentric-rings-reference formulas,
#' vectorised over P/T directly. Skips the per-triangle CCW-winding
#' correction .moment_of_inertia_core() carries (needed there because CDT
#' triangulation doesn't guarantee consistent winding) - verified
#' empirically across several real and synthetic inputs that
#' RTriangle::triangulate()'s own output is always consistently
#' CCW-wound (every triangle's signed area positive, zero exceptions), so
#' the shoelace formulas below can rely on that directly.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix - a
#'   .constrained_weighted_mesh() result's own P/T
#' @param tri_area numeric vector, length nrow(T)
#' @param rho numeric vector, length nrow(T) - each triangle's density
#'   (weight / area)
#' @param crs the mesh's CRS, attached to the returned centroid point
#' @return same shape as .moment_of_inertia_core(): list(index, J, Ixx,
#'   Iyy, Ixy, J_ref, area, total_weight, centroid, triangles = NULL)
#' @noRd
.moment_of_inertia_core_array <- function(P, T, tri_area, rho, crs) {
    mass <- rho * tri_area
    tc <- (P[T[, 1], , drop = FALSE] + P[T[, 2], , drop = FALSE] + P[T[, 3], , drop = FALSE]) / 3
    G <- c(sum(mass * tc[, 1]), sum(mass * tc[, 2])) / sum(mass)

    vx1 <- P[T[, 1], 1] - G[1]; vy1 <- P[T[, 1], 2] - G[2]
    vx2 <- P[T[, 2], 1] - G[1]; vy2 <- P[T[, 2], 2] - G[2]
    vx3 <- P[T[, 3], 1] - G[1]; vy3 <- P[T[, 3], 2] - G[2]

    cross1 <- vx1 * vy2 - vx2 * vy1
    cross2 <- vx2 * vy3 - vx3 * vy2
    cross3 <- vx3 * vy1 - vx1 * vy3

    Ixx_tri <- ((vy1^2 + vy1 * vy2 + vy2^2) * cross1 +
                (vy2^2 + vy2 * vy3 + vy3^2) * cross2 +
                (vy3^2 + vy3 * vy1 + vy1^2) * cross3) / 12
    Iyy_tri <- ((vx1^2 + vx1 * vx2 + vx2^2) * cross1 +
                (vx2^2 + vx2 * vx3 + vx3^2) * cross2 +
                (vx3^2 + vx3 * vx1 + vx1^2) * cross3) / 12
    # standard polygon product-of-inertia formula, same cross/winding
    # setup as Ixx_tri/Iyy_tri above, just the xy cross term
    Ixy_tri <- ((vx1 * vy2 + 2 * vx1 * vy1 + 2 * vx2 * vy2 + vx2 * vy1) * cross1 +
                (vx2 * vy3 + 2 * vx2 * vy2 + 2 * vx3 * vy3 + vx3 * vy2) * cross2 +
                (vx3 * vy1 + 2 * vx3 * vy3 + 2 * vx1 * vy1 + vx1 * vy3) * cross3) / 24

    Ixx <- sum(rho * Ixx_tri)
    Iyy <- sum(rho * Iyy_tri)
    Ixy <- sum(rho * Ixy_tri)
    J   <- Ixx + Iyy
    A   <- sum(tri_area)
    W   <- sum(rho * tri_area)

    ord    <- order(rho, decreasing = TRUE)
    S      <- cumsum(tri_area[ord])
    S_prev <- c(0, S[-length(S)])
    J_ref  <- sum(rho[ord] * (S^2 - S_prev^2)) / (2 * pi)

    index <- if (J > 0) J_ref / J else NA_real_

    list(index = index, J = J, Ixx = Ixx, Iyy = Iyy, Ixy = Ixy, J_ref = J_ref,
         area = A, total_weight = W,
         centroid = st_sfc(st_point(G), crs = crs), triangles = NULL)
}

#' Every triangle's own centroid, vectorised over P/T at once - used
#' directly for row-assignment in .constrained_weighted_mesh(), and as
#' the mass-centroid building block for the directional-balance/radial-
#' concentration array mirrors below. NOT used as their actual point
#' cloud on its own - a constrained triangulation only adds vertices
#' where ROW BOUNDARIES need them, so a large, simple, uniform-density
#' row (a big rural block, say) can still triangulate into just a
#' handful of large interior triangles, and collapsing one of those to
#' its bare centroid would reintroduce exactly the centroid-collapse
#' bias .adaptive_tri_depth()'s own docs warn against - see
#' .radial_point_cloud_array() below for the point cloud those two
#' indices actually use.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @return Mx2 matrix of triangle centroids
#' @noRd
.tri_centroids_array <- function(P, T) {
    (P[T[, 1], , drop = FALSE] + P[T[, 2], , drop = FALSE] + P[T[, 3], , drop = FALSE]) / 3
}

#' Array-native mirror of .radial_point_cloud() (radial-concentration-
#' index.R) - identical area-adaptive medial-subdivision math (reusing
#' .subdivide_tri_batch()/.adaptive_tri_depth() unchanged, since both
#' already operate on plain coordinate matrices, not sf), just reading
#' each triangle's own corners directly out of P/T instead of an sf
#' st_coordinates() loop. See .tri_centroids_array()'s own comment for
#' why subdivision still matters on this mesh, not just the coarse one.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param w numeric vector, length nrow(T) - each triangle's own mass
#' @param max_depth subdivision depth ceiling
#' @return list(p, w) - p an Nx2 coordinate matrix, w length N
#' @noRd
.radial_point_cloud_array <- function(P, T, tri_area, w, max_depth = 4) {
    A_all <- P[T[, 1], , drop = FALSE]
    B_all <- P[T[, 2], , drop = FALSE]
    C_all <- P[T[, 3], , drop = FALSE]

    depth_i <- .adaptive_tri_depth(tri_area, max_depth)

    p_list <- vector("list", max_depth + 1)
    w_list <- vector("list", max_depth + 1)
    for (d in 0:max_depth) {
        grp <- which(depth_i == d)
        if (length(grp) == 0) next
        A <- A_all[grp, , drop = FALSE]
        B <- B_all[grp, , drop = FALSE]
        C <- C_all[grp, , drop = FALSE]
        orig_idx <- seq_along(grp)
        for (k in seq_len(d)) {
            sub <- .subdivide_tri_batch(A, B, C)
            A <- sub$A; B <- sub$B; C <- sub$C
            orig_idx <- rep(orig_idx, times = 4)
        }
        n_sub <- 4^d
        p_list[[d + 1]] <- (A + B + C) / 3
        w_list[[d + 1]] <- (w[grp] / n_sub)[orig_idx]
    }
    list(p = do.call(rbind, p_list), w = unlist(w_list))
}

#' Array-native mirror of .mesh_directional_balance_index() (directional-
#' balance-index.R) - same mass-centroid/resultant-length math, point
#' cloud from .radial_point_cloud_array() instead of the sf-based one.
#' Only ever called with a real weight vector - this mesh only exists
#' for shape_indices_sf(byrow = FALSE, weights = ...).
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param weight numeric vector, length nrow(T)
#' @param crs the mesh's CRS, attached to the returned centroid point
#' @return same shape as .mesh_directional_balance_index(): list(index, R,
#'   mean_angle, area, total_weight, centroid, triangles = NULL)
#' @noRd
.mesh_directional_balance_index_array <- function(P, T, tri_area, weight, crs) {
    area <- sum(tri_area)
    w    <- .normalize_weight(weight)
    tc   <- .tri_centroids_array(P, T)
    G    <- c(sum(w * tc[, 1]), sum(w * tc[, 2])) / sum(w)

    cloud <- .radial_point_cloud_array(P, T, tri_area, w)
    res   <- .resultant_from_cloud(cloud$p, cloud$w, G)

    list(index = 1 - res$R, R = res$R, mean_angle = res$mean_angle, area = area,
         total_weight = sum(weight), centroid = st_sfc(st_point(G), crs = crs), triangles = NULL)
}

#' Array-native mirror of .mesh_radial_concentration_index() (radial-
#' concentration-index.R) - same geometric-median/annulus-reference math,
#' point cloud from .radial_point_cloud_array() instead of the sf-based
#' one. Only ever called with a real weight vector, same reason as
#' .mesh_directional_balance_index_array() above.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param weight numeric vector, length nrow(T)
#' @param crs the mesh's CRS, attached to the returned center point
#' @return same shape as .mesh_radial_concentration_index(): list(index,
#'   D1, D1_ref, area, total_weight, center, triangles = NULL)
#' @noRd
.mesh_radial_concentration_index_array <- function(P, T, tri_area, weight, crs) {
    area <- sum(tri_area)
    w    <- .normalize_weight(weight)

    cloud <- .radial_point_cloud_array(P, T, tri_area, w)
    gm    <- .geometric_median(cloud$p, cloud$w)

    D1_ref <- .annulus_reference_D1(tri_area, weight)
    index  <- D1_ref / gm$D1

    list(index = index, D1 = gm$D1, D1_ref = D1_ref, area = area,
         total_weight = sum(weight), center = st_sfc(st_point(gm$center), crs = crs), triangles = NULL)
}

#' Array-native mirror of .random_point_directional_balance_index()'s
#' points-supplied branch (directional-balance-index.R) - that function
#' needs the mass centroid (via st_geometry(tri)) unconditionally, even
#' once points are pre-drawn, so a lightweight tri stand-in isn't enough;
#' this recomputes G directly off P/T instead.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param weight numeric vector, length nrow(T)
#' @param points pre-drawn n_lines x 2 coordinate matrix
#' @param crs the mesh's CRS, attached to the returned centroid point
#' @return same shape as .random_point_directional_balance_index()
#' @noRd
.random_point_directional_balance_index_array <- function(P, T, tri_area, weight, points, crs) {
    area <- sum(tri_area)
    w    <- .normalize_weight(weight)
    tc   <- .tri_centroids_array(P, T)
    G    <- c(sum(w * tc[, 1]), sum(w * tc[, 2])) / sum(w)

    res <- .resultant_from_cloud(points, rep(1, nrow(points)), G)

    list(index = 1 - res$R, R = res$R, mean_angle = res$mean_angle, area = area,
         total_weight = sum(weight), centroid = st_sfc(st_point(G), crs = crs), triangles = NULL)
}

#' Array-native mirror of .random_point_radial_index()'s points-supplied
#' branch (radial-concentration-index.R) - same geometric-median/annulus-
#' reference math over the pre-drawn Monte Carlo sample.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param weight numeric vector, length nrow(T)
#' @param points pre-drawn n_lines x 2 coordinate matrix
#' @param crs the mesh's CRS, attached to the returned center point
#' @return same shape as .random_point_radial_index()
#' @noRd
.random_point_radial_index_array <- function(P, T, tri_area, weight, points, crs) {
    area <- sum(tri_area)
    gm   <- .geometric_median(points, rep(1, nrow(points)))

    D1_ref <- .annulus_reference_D1(tri_area, weight)
    index  <- D1_ref / gm$D1

    list(index = index, D1 = gm$D1, D1_ref = D1_ref, area = area,
         total_weight = sum(weight), center = st_sfc(st_point(gm$center), crs = crs), triangles = NULL)
}

#' Array-native mirror of .mesh_depth_index() (depth-index.R) - same
#' point-cloud/annulus-reference-identity math, point cloud from
#' .radial_point_cloud_array() instead of the sf-based one. Unlike the
#' other array mirrors above, this needs `poly_u` itself (not just P/T) -
#' depth's own field is distance to the ACTUAL polygon boundary, not
#' anything derivable from triangle centroids/areas alone (see
#' depth-index.R's own file header). Only ever called with a real weight
#' vector, same reason as the other array mirrors above.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param weight numeric vector, length nrow(T)
#' @param poly_u the (multi)polygon this mesh triangulates
#' @return same shape as .mesh_depth_index(): list(index, mean_depth,
#'   ref_depth, area, total_weight, triangles = NULL)
#' @noRd
.mesh_depth_index_array <- function(P, T, tri_area, weight, poly_u) {
    area <- sum(tri_area)
    w    <- .normalize_weight(weight)

    cloud <- .radial_point_cloud_array(P, T, tri_area, w)
    bnd   <- st_boundary(poly_u)
    d     <- as.numeric(st_distance(.coords_to_points(cloud$p, st_crs(poly_u)), bnd))
    mean_depth <- sum(cloud$w * d) / sum(cloud$w)

    ref_depth <- .annulus_reference_depth(tri_area, weight)
    index <- mean_depth / ref_depth

    list(index = index, mean_depth = mean_depth, ref_depth = ref_depth, area = area,
         total_weight = sum(weight), triangles = NULL)
}

#' Array-native mirror of .random_point_depth_index()'s points-supplied
#' branch (depth-index.R) - same boundary-distance math over the
#' pre-drawn Monte Carlo sample.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param weight numeric vector, length nrow(T)
#' @param points pre-drawn n_lines x 2 coordinate matrix
#' @param poly_u the (multi)polygon this mesh triangulates
#' @return same shape as .random_point_depth_index()
#' @noRd
.random_point_depth_index_array <- function(P, T, tri_area, weight, points, poly_u) {
    area <- sum(tri_area)
    bnd  <- st_boundary(poly_u)
    d    <- as.numeric(st_distance(.coords_to_points(points, st_crs(poly_u)), bnd))
    mean_depth <- mean(d)

    ref_depth <- .annulus_reference_depth(tri_area, weight)
    index <- mean_depth / ref_depth

    list(index = index, mean_depth = mean_depth, ref_depth = ref_depth, area = area,
         total_weight = sum(weight), triangles = NULL)
}

#' Lazily converts a slice of this mesh's own P/T arrays into real sf
#' TRIANGLE pieces - the one thing convexity_index()'s/span_index()'s
#' deterministic mode genuinely needs (boundary-clipped candidate lines
#' and the same-piece diagonal term for convexity; real per-triangle
#' geometry for span's quadrature/self-distance terms), and exactly the
#' cost this file's own header describes as prohibitive at full mesh
#' scale. Deterministic mode is only ever reached here below
#' .safe_deterministic_tri_ceiling() (see shape-indices.R's dispatch), so
#' n_tri is always small when this actually runs - the ~500s measured
#' cost was at ~2 million triangles, not the few hundred/thousand this
#' is bounded to.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param tri_area numeric vector, length nrow(T)
#' @param crs the mesh's CRS
#' @return sf data frame with an `area` column and POLYGON geometry, same
#'   shape .weighted_mesh()'s own `pieces` has
#' @noRd
.array_mesh_to_sf_pieces <- function(P, T, tri_area, crs) {
    n <- nrow(T)
    geom <- st_sfc(lapply(seq_len(n), function(i) {
        v <- P[T[i, ], , drop = FALSE]
        st_polygon(list(rbind(v, v[1, , drop = FALSE])))
    }), crs = crs)
    st_sf(area = tri_area, geometry = geom)
}

#' Array-native mirror of .sample_weighted_points() (convexity-index.R) -
#' same triangle-then-barycentric-point sampling, reading corner
#' coordinates directly out of P/T instead of a per-triangle
#' st_coordinates() loop, and only for the SAMPLED triangles rather than
#' building a corner-coordinate table for every triangle in the mesh
#' regardless of how many points are drawn.
#' @param P Nx2 point matrix, T Mx3 triangle-vertex-index matrix
#' @param weight numeric vector, length nrow(T), need not sum to 1
#' @param n number of points to draw
#' @return n x 2 coordinate matrix
#' @noRd
.sample_weighted_points_array <- function(P, T, weight, n) {
    tri_idx <- sample.int(nrow(T), size = n, replace = TRUE, prob = weight)
    A <- P[T[tri_idx, 1], , drop = FALSE]
    B <- P[T[tri_idx, 2], , drop = FALSE]
    C <- P[T[tri_idx, 3], , drop = FALSE]

    r1 <- runif(n); r2 <- runif(n)
    flip     <- (r1 + r2) > 1
    r1[flip] <- 1 - r1[flip]
    r2[flip] <- 1 - r2[flip]

    cbind(x = A[, 1] + r1 * (B[, 1] - A[, 1]) + r2 * (C[, 1] - A[, 1]),
          y = A[, 2] + r1 * (B[, 2] - A[, 2]) + r2 * (C[, 2] - A[, 2]))
}

## -- memory-derived safety ceiling for O(n^2) deterministic modes -------
##
## A .constrained_weighted_mesh() can be orders of magnitude larger than
## any CDT mesh this package has ever built (real Urban Area data: 1.98
## million triangles, vs. a few tens of thousands from the coarse
## .weighted_mesh() path) - large enough that convexity_index()'s/
## span_index()'s DETERMINISTIC mode, O(n^2) in triangle-pair count, goes
## from "slow" to genuinely unsafe: at that scale, n_candidate_lines
## exceeds .mesh_convexity_index()'s own existing "large mesh" WARNING
## threshold (20,000) by roughly 98 million times. A fixed hard-coded
## ceiling would be the wrong shape for this - it can't know how much
## memory is actually available on whatever machine is running the
## computation. Instead, reuses this package's own existing
## .available_memory_mb()/.choose_line_chunk_size() pattern (utils.R,
## built for exactly this kind of decision already, for convexity's line-
## clipping chunking) to derive a triangle-count ceiling from real
## available memory, inverting each index's own already-existing cost
## formula rather than inventing a new one.

#' The largest triangle count safe to run convexity_index()'s or
#' span_index()'s deterministic (O(n^2)-in-triangle-pairs) mode with,
#' given a memory budget - inverts whichever index's own existing "large
#' mesh" warning formula (.mesh_convexity_index()'s n_candidate_lines =
#' choose(n,2)*n_quad^2; .mesh_span_index()'s n_points^2 = (n*n_quad)^2)
#' against .available_memory_mb(). This is a genuinely different
#' threshold from those functions' own existing (fixed, much smaller)
#' warning thresholds - those are calibrated for "this might take a
#' noticeable amount of time", a low bar; this one is calibrated for "this
#' would exhaust available memory", a much higher one - so it's expected,
#' not a bug, that this ceiling allows far more triangles than the
#' existing warning threshold would pass without complaint.
#' @param n_quad quadrature points per triangle (1 or 3)
#' @param formula which index's own cost formula to invert
#' @param bytes_per_pair conservative bytes needed per candidate-line/
#'   quadrature-point-pair - same estimate .choose_line_chunk_size() uses
#' @param mem_fraction fraction of estimated available memory to budget
#' @return integer >= 1, the largest triangle count considered safe
#' @noRd
.safe_deterministic_tri_ceiling <- function(n_quad, formula = c("convexity", "span"),
                                             bytes_per_pair = 200, mem_fraction = 0.2) {
    formula <- match.arg(formula)
    budget_bytes <- .available_memory_mb() * 1024^2 * mem_fraction
    max_units <- budget_bytes / bytes_per_pair
    n_max <- if (formula == "convexity") {
        # choose(n,2)*n_quad^2 <= max_units, i.e. n(n-1) <= 2*max_units/n_quad^2
        0.5 + sqrt(0.25 + 2 * max_units / n_quad^2)
    } else {
        # (n*n_quad)^2 <= max_units
        sqrt(max_units) / n_quad
    }
    max(1L, as.integer(floor(n_max)))
}
