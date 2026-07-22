## =========================================================================
## 7. Radial concentration index: mean distance to the geometric median vs
## an equal-area circle
## =========================================================================
##
## D1(rho, c) = (1/W) integral rho(x) |x-c| dx, the mean distance from a
## rho/W-weighted random point to a candidate centre c. D1(rho) =
## min_c D1(rho,c), achieved at the geometric median (Fermat-Weber point) -
## NOT the centroid, which only minimises the SQUARED distance (that's
## moment_of_inertia_index()'s job: variance is minimised at the mean, but
## mean absolute deviation is minimised at the median). index =
## D1_ref/D1(rho), in (0, 1].
##
## REFERENCE SHAPES: a bathtub-principle argument, simpler than
## span_index()'s Riesz argument since this is a single integral, not a
## double one over pairs. For any FIXED centre c, the disk of area A
## centred at c minimises integral_Omega |x-c| dx among all measurable
## sets of that area (|x-c| is radially increasing about c, so the
## bathtub principle applies directly). Evaluating this at c = Omega's
## own geometric median (translated to the origin) shows the disk is the
## exact global minimiser of D1 over every measurable shape of a given
## area - holes, multiple parts, non-convexity: none of that affects the
## proof. Weighted: the same per-level-set argument gives concentric
## rings sorted densest-to-centre, the same construction
## moment_of_inertia_index()'s and span_index()'s references use.
##
## CLOSED FORMS: distance-to-a-fixed-point depends only on radius, not
## angle, so (unlike span_index()) neither reference needs elliptic
## integrals - plain 1D radial integrals suffice:
##   D1_circle(A) = (2/3) * sqrt(A/pi)
##   ring mean radius, inner/outer a, b: (2/3) * (a^2+a*b+b^2)/(a+b)
##
## ACTUAL VALUE (deterministic = TRUE): the geometric median has no closed
## form - found via Weiszfeld's algorithm on a point cloud, then D1 is the
## weighted mean distance from that cloud to the found centre. There's no
## self-term issue the way span_index() had one (point-to-fixed-point
## distance never degenerates the way same-triangle point-pairs did), but
## a DIFFERENT refinement issue takes its place: |x-c| is convex in x, so
## by Jensen's inequality collapsing a triangle to its own centroid
## systematically UNDERESTIMATES its contribution, while tri_quad_points()'s
## fixed 3-point rule (built for accurately estimating AREA integrals, not
## distance-to-an-external-point) overshoots it instead - verified against
## a direct Monte Carlo estimate to be a genuinely non-negligible bias
## (single digits of percent) at the mesh sizes this package actually
## sees, not a rounding-level nicety. The fix is the same one span_index()
## uses for its self-term: recursively split each triangle by medial
## subdivision and use the sub-triangle centroids as the point cloud,
## which converges geometrically as subdivision deepens
## (.subdivide_tri_batch(), a vectorised version of span-index.R's
## .subdivide_tri() that subdivides every triangle in one pass instead of
## looping triangle-by-triangle). Every triangle needs this here (not
## just a same-piece special case), so the point cloud is bigger than
## span's quadrature cloud: n_triangles * 4^4 = 256x points, still only
## O(n) per Weiszfeld iteration rather than O(n^2) like
## convexity_index()/span_index()'s deterministic mode - but that 256x
## constant factor is large enough on real-world meshes (tens of
## thousands of CDT triangles from a complex polygon) that it can
## dominate the whole computation's memory and time - "still only O(n)"
## isn't cheap enough on its own once that 256x constant factor is large,
## which is exactly why deterministic = FALSE exists below.
##
## ACTUAL VALUE (deterministic = FALSE): draws `n_lines` points directly
## from the weighted density (.random_point_radial_index(), reusing
## convexity_index()'s .sample_weighted_points()) and runs Weiszfeld on
## that sample instead - decouples cost from mesh size entirely, the same
## fix convexity_index()/span_index() already have. When called from
## shape_indices()/shape_indices_sf() alongside convexity and/or span in
## Monte Carlo mode, all three reuse ONE shared draw instead of sampling
## independently three times over (see shape-indices.R's own comments on
## its point-sharing mechanism, the Monte Carlo analogue of sharing one
## `prep` triangulation across indices).
##
## WEISZFELD CONVERGENCE: stops on the OBJECTIVE (mean distance) changing
## by less than a tolerance, not on the centre's position. For a
## symmetric multi-part shape (e.g. two identical separated blobs) the
## minimiser is a whole segment, not a single point (the classical
## two-equal-point-masses degeneracy), so the position can keep moving
## along that segment even after the value has converged. A small floor
## on distances guards against division by zero if an iterate lands
## exactly on a point.

#' Weiszfeld's algorithm for the weighted geometric median of a point
#' cloud - the point minimising the weighted mean distance to it.
#' @param p an Nx2 coordinate matrix
#' @param w numeric vector, length N, weights (need not sum to 1)
#' @param max_iter maximum iterations
#' @param tol relative convergence tolerance on the objective (mean
#'   weighted distance) - see file header for why not on position
#' @return list(center, D1) - D1 is the converged mean weighted distance
#' @noRd
.geometric_median <- function(p, w, max_iter = 200, tol = 1e-10) {
    mean_dist <- function(c) sum(w * sqrt((p[, 1] - c[1])^2 + (p[, 2] - c[2])^2)) / sum(w)

    c_t <- colSums(p * w) / sum(w)   # start from the centroid
    obj <- mean_dist(c_t)
    for (i in seq_len(max_iter)) {
        d  <- pmax(sqrt((p[, 1] - c_t[1])^2 + (p[, 2] - c_t[2])^2), 1e-12)
        wd <- w / d
        c_t <- colSums(p * wd) / sum(wd)
        obj_new <- mean_dist(c_t)
        if (abs(obj - obj_new) < tol * max(1, obj_new)) { obj <- obj_new; break }
        obj <- obj_new
    }
    list(center = c_t, D1 = obj)
}

#' One medial-subdivision step, applied to ALL current triangles at once -
#' the vectorised counterpart of span-index.R's .subdivide_tri(), which
#' does the same thing for a single triangle via nested lists. Matches
#' .subdivide_tri()'s exact child ordering (v1,m_ab,m_ca) /
#' (v2,m_bc,m_ab) / (v3,m_ca,m_bc) / (m_ab,m_bc,m_ca), just stacked as
#' four Nx2 blocks instead of returned as a length-4 list, so N triangles
#' in -> 4N triangles out in one vectorised pass instead of N separate
#' recursive calls.
#' @param A,B,C Nx2 matrices, the three vertices of N triangles
#' @return list(A, B, C), each 4Nx2 - the vertices of the 4N children
#' @noRd
.subdivide_tri_batch <- function(A, B, C) {
    AB <- (A + B) / 2
    BC <- (B + C) / 2
    CA <- (C + A) / 2
    list(A = rbind(A, B, C, AB), B = rbind(AB, BC, CA, BC), C = rbind(CA, AB, BC, CA))
}

#' Point cloud used for both the geometric-median search and the D1
#' estimate: every CDT triangle recursively split by medial subdivision
#' into equal-area sub-triangles, one point (its centroid) per
#' sub-triangle - see file header for why tri_quad_points()'s coarser
#' 1/3-point rule isn't accurate enough here.
#'
#' DEPTH IS AREA-ADAPTIVE, NOT FIXED: a fixed depth for every triangle
#' means every triangle pays the same 4^depth cost regardless of its own
#' size - wasteful on a real CDT mesh, where a complex boundary produces
#' many small triangles alongside a handful of large ones, and it's the
#' large triangles specifically whose centroid-collapse bias is worst
#' (Jensen's inequality - see file header). Subdivision halves a
#' triangle's own edge length (quarters its area) per level, so matching
#' a fixed FINAL sub-triangle size across the whole mesh - rather than a
#' fixed depth - means small triangles need less of it: depth_i =
#' ceil(log4(area_i / target_area)), clamped to `[0, max_depth]`. Anchoring
#' target_area to the mesh's OWN largest triangle at max_depth keeps the
#' worst-case triangle's accuracy identical to the old fixed-depth
#' behaviour (no regression there); every smaller triangle needs fewer
#' points to reach that same final resolution, cutting the total point
#' count for any mesh with real size variation - which real-world CDT
#' meshes always have.
#'
#' Triangles are grouped by their own depth_i (at most max_depth + 1
#' distinct groups) and each group is subdivided uniformly via
#' .subdivide_tri_batch(), reusing the exact same vectorised, per-level
#' batch logic - just scoped to one depth-homogeneous group at a time
#' instead of the whole mesh, and concatenated at the end. The batched
#' rbind() stacking interleaves children in triangle-index blocks rather
#' than grouping each triangle's own descendants contiguously, so
#' `orig_idx` (rebuilt fresh per group) tracks, for every row of that
#' group's point cloud, which of the group's OWN triangles it descended
#' from - needed to assign the right weight fraction to each point,
#' since row order alone no longer says so.
#' @param tri sf data frame of CDT triangles, with an `area` column
#' @param w numeric vector, each triangle's own mass (already resolved -
#'   physical area if unweighted, normalised weight if weighted)
#' @param max_depth subdivision depth ceiling, applied to the mesh's own
#'   largest triangle only (4^max_depth points for that triangle; fewer
#'   for smaller ones, adaptively)
#' @return list(p, w) - p an Nx2 coordinate matrix, w length N
#' @noRd
.radial_point_cloud <- function(tri, w, max_depth = 4) {
    n <- nrow(tri)
    tri_area <- tri$area
    corner_mat <- vapply(st_geometry(tri), function(g) {
        v <- st_coordinates(g)[1:3, 1:2, drop = FALSE]
        c(v[1, ], v[2, ], v[3, ])
    }, numeric(6))
    A_all <- t(corner_mat[1:2, , drop = FALSE])
    B_all <- t(corner_mat[3:4, , drop = FALSE])
    C_all <- t(corner_mat[5:6, , drop = FALSE])

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

#' @noRd
.disk_reference_D1 <- function(area) (2 / 3) * sqrt(area / pi)

#' Weighted reference: sort triangles by density descending (same
#' construction as moment_of_inertia_index()'s/span_index()'s rings),
#' then average each ring's own mean radius, weighted by ring mass -
#' closed form throughout, no quadrature needed (see file header).
#' @param tri_area numeric vector, each triangle's own area
#' @param weight numeric vector, one entry per triangle, need not sum to 1
#' @return D1_ref
#' @noRd
.annulus_reference_D1 <- function(tri_area, weight) {
    W    <- .normalize_weight(weight)
    rho  <- W / tri_area
    ord  <- order(rho, decreasing = TRUE)
    r_hi <- sqrt(cumsum(tri_area[ord]) / pi)
    r_lo <- c(0, r_hi[-length(r_hi)])
    W    <- W[ord]

    mean_r <- (2 / 3) * (r_lo^2 + r_lo * r_hi + r_hi^2) / (r_lo + r_hi)
    sum(W * mean_r)
}

#' Deterministic mean-distance-to-geometric-median over a fixed mesh -
#' shared engine behind radial_concentration_index().
#' @param tri sf data frame of CDT triangles, with an `area` column
#' @param weight optional numeric vector, length nrow(tri)
#' @return list(index, D1, D1_ref, area, total_weight, center, triangles)
#' @noRd
.mesh_radial_concentration_index <- function(tri, weight = NULL) {
    n <- if (is.null(tri)) 0 else nrow(tri)
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }
    if (n == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, D1 = NA_real_, D1_ref = NA_real_, area = NA_real_,
                    total_weight = NA_real_, center = NULL, triangles = tri))
    }

    area      <- sum(tri$area)
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri$area else .normalize_weight(weight)

    cloud <- .radial_point_cloud(tri, w)
    gm    <- .geometric_median(cloud$p, cloud$w)

    D1_ref <- if (is.null(weight)) .disk_reference_D1(area) else .annulus_reference_D1(tri$area, weight)
    index  <- D1_ref / gm$D1

    list(index = index, D1 = gm$D1, D1_ref = D1_ref, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total,
         center = st_sfc(st_point(gm$center), crs = st_crs(tri)), triangles = tri)
}

## Stochastic (deterministic = FALSE): draws `n_lines` points directly from
## the weighted density (reusing convexity_index()'s .sample_weighted_points()
## - proportional-by-area sampling when unweighted, which is an exact,
## GEOS-free way to sample uniformly over the polygon, not an approximation)
## and runs Weiszfeld on that sample instead of the depth^4-per-triangle
## subdivision cloud .mesh_radial_concentration_index() always builds. This
## is what decouples cost from mesh size the way convexity_index()/
## span_index() already do - see file header for why this index never had
## that escape valve until now. Named `n_lines`, not `n_points`, for the
## same reason span_index() shares the name with convexity_index(): one
## argument sets every mesh index's Monte Carlo sample count via
## shape_indices()'s `...`, even though no line or point-pair is built here.

#' Monte Carlo (deterministic = FALSE) mean-distance-to-geometric-median -
#' internal engine for radial_concentration_index()'s deterministic =
#' FALSE. See that function's own roxygen for the user-facing parameter
#' reference; parameters here are the same, without defaults.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param n_lines number of random points to sample, or a function(n_tri)
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon()
#' @param seed optional RNG seed
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`; NULL samples proportional to triangle area
#' @param points optional pre-drawn Nx2 coordinate matrix, from
#'   shape_indices()'s point-sharing mechanism (see its own comments) -
#'   when supplied, sampling is skipped and every row is used as one point
#'   in the cloud; N need not equal n_lines (more points only improves
#'   precision)
#' @return list(index, D1, D1_ref, area, total_weight, center, triangles)
#' @noRd
.random_point_radial_index <- function(poly, n_lines, prep, seed, weight = NULL, points = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly)
    tri   <- prep$tri
    n_tri <- if (is.null(tri)) 0 else nrow(tri)

    if (!is.null(weight)) {
        if (n_tri == 0) {
            stop("`weight` needs a triangle mesh to sample from, but this polygon ",
                 "triangulated to no triangles.")
        }
        if (length(weight) != n_tri) {
            stop("`weight` must have one entry per triangle (", n_tri, "), got ", length(weight), ".")
        }
    }
    if (n_tri == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, D1 = NA_real_, D1_ref = NA_real_, area = NA_real_,
                    total_weight = NA_real_, center = NULL, triangles = tri))
    }

    if (is.function(n_lines)) n_lines <- max(1L, round(n_lines(n_tri)))

    area      <- sum(tri$area)
    raw_total <- if (is.null(weight)) NULL else sum(weight)

    if (!is.null(points)) {
        coords <- points
    } else {
        if (!is.null(seed)) set.seed(seed)
        coords <- .sample_weighted_points(tri, weight %||% tri$area, n_lines)
    }

    gm <- .geometric_median(coords, rep(1, nrow(coords)))

    D1_ref <- if (is.null(weight)) .disk_reference_D1(area) else .annulus_reference_D1(tri$area, weight)
    index  <- D1_ref / gm$D1

    list(index = index, D1 = gm$D1, D1_ref = D1_ref, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total,
         center = st_sfc(st_point(gm$center), crs = st_crs(tri)), triangles = tri)
}

#' Radial concentration index: mean distance to the geometric median vs
#' an equal-area circle
#'
#' D1_ref/D1, in (0, 1], where D1 is the mean distance from random
#' interior points to the shape's own geometric median (the point
#' minimising that mean distance - not the centroid) and D1_ref is the
#' same quantity for the reference shape (a circle, unweighted; a
#' concentric annulus, weighted) - both provable minimisers of D1. See
#' `vignette("f-understanding-radial-concentration-index")` for the
#' derivation and the proofs.
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param deterministic if TRUE (default), compute over a fixed
#'   depth-4-subdivision point cloud (4^4 = 256 points per CDT triangle) -
#'   O(n) in triangle count, but with a large constant factor from that
#'   256x expansion. If FALSE, Monte Carlo over `n_lines` points sampled
#'   directly from the weighted density instead - decouples cost from
#'   mesh size, useful when the mesh is large enough that the 256x
#'   subdivision cloud itself becomes the bottleneck.
#' @param n_lines number of random points to sample when
#'   deterministic = FALSE (default 3000), or a function(n_tri) - same
#'   argument name as convexity_index()'s/span_index()'s Monte Carlo mode,
#'   so one value passed through shape_indices() sets every mesh index's
#'   sample count at once, even though this draws single points, not
#'   lines or point-pairs.
#' @param seed optional RNG seed, only used when deterministic = FALSE.
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating when also calling other indices
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, substituting for that triangle's own area/mass throughout.
#'   NULL (default) compares against a circle; otherwise the reference
#'   becomes a concentric annulus matching this weight's own histogram.
#' @param points optional pre-drawn Nx2 coordinate matrix, only meaningful
#'   when deterministic = FALSE - primarily for shape_indices()'s internal
#'   point-sharing mechanism (draw one sample and reuse it across
#'   convexity/span/radial_concentration instead of each drawing
#'   independently), not typically supplied directly.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return list(index, D1, D1_ref, area, total_weight, center, triangles).
#'   `center` is the geometric median found by Weiszfeld's algorithm. For
#'   symmetric shapes the geometric median can be non-unique (the
#'   minimising "point" is a whole segment - e.g. between two identical
#'   separated blobs), so `center` may land anywhere along it, including
#'   inside a hole or the gap between multi-part pieces; the index value
#'   itself is unaffected.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' radial_concentration_index(wake)$index
#' # Camden: a long, narrow county - far from a circle
#' radial_concentration_index(nc[nc$NAME == "Camden", ])$index
#'
#' # deterministic = FALSE: Monte Carlo estimate instead of the full
#' # subdivision cloud - a seed makes it reproducible
#' radial_concentration_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#'
#' # weight: substitutes for each triangle's own area - weighting by the
#' # triangle's own area exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' radial_concentration_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
radial_concentration_index <- function(poly, deterministic = TRUE, n_lines = 3000, seed = NULL,
                                        prep = NULL, weight = NULL, points = NULL,
                                        simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)

    if (!deterministic) {
        return(.random_point_radial_index(poly, n_lines = n_lines, prep = prep, seed = seed,
                                           weight = weight, points = points))
    }
    if (!is.null(points)) {
        stop("`points` (a pre-drawn sample) has no meaning for deterministic = TRUE, which ",
             "computes over the full subdivision point cloud, not a random sample. Drop ",
             "`points`, or set deterministic = FALSE to use it.")
    }

    .mesh_radial_concentration_index(prep$tri, weight = weight)
}

#' radial_concentration_index() for every row of an sf data frame
#'
#' Each row indexed independently and unweighted.
#' @param x an sf data frame
#' @param ... passed to radial_concentration_index() for every row
#' @return `x` with one new column, radial_concentration_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- radial_concentration_index_sf(nc[1:5, ])
#' res$radial_concentration_index
#' @export
radial_concentration_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) radial_concentration_index(geoms[i], ...))
    x$radial_concentration_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
