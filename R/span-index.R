## =========================================================================
## 5. Span index: mean pairwise interior distance relative to a reference shape
## =========================================================================
##
## D(rho) = (1/W^2) integral-integral rho(x) rho(y) |x-y| dx dy - the mean
## distance between two independent points drawn from density rho/W.
## index = D_ref/D(rho), in (0, 1], = 1 iff rho already matches its
## reference shape.
##
## Distinct from moment_of_inertia_index(): E|X-Y|^2 = 2J/W exactly (a
## Var(X-Y) identity), so a squared-distance version of this index would
## just restate J_ref/J_actual, not add information. Plain (non-squared)
## distance has no such shortcut - E|X-Y| < sqrt(E|X-Y|^2) strictly
## (Jensen), by a shape-dependent gap, which is what this index measures.
##
## REFERENCE SHAPES: by Riesz's rearrangement inequality (applied to the
## kernel C - |x-y|, C large enough to keep it nonnegative - equivalent to
## minimising |x-y| itself since rearrangement preserves total mass), D is
## minimised when rho is arranged symmetric-decreasing about a common
## centre: a disk for uniform density (weight = NULL), concentric rings
## sorted densest-to-centre for a piecewise-constant one (weight supplied).
## These are true global minimisers over every rearrangement of the given
## weight/area histogram, not just plausible candidates.
##
## UNWEIGHTED REFERENCE: mean distance between 2 uniform points in a disk
## of radius R is the classical 128R/(45*pi); with R = sqrt(A/pi),
##   D_circle(A) = 128*sqrt(A) / (45 * pi^1.5)
##
## WEIGHTED REFERENCE: no closed form exists ring-to-ring (the |x-y| kernel
## has no polynomial antiderivative in r the way moment_of_inertia_index()'s
## r^2 does), but averaging over the angle between two points at radii
## r1, r2 does:
##   g(r1, r2) = (2/pi)*(r1+r2)*E(k),  k = 2*sqrt(r1*r2)/(r1+r2)
## where E is the complete elliptic integral of the second kind (computed
## via the AGM algorithm below - no new package dependency). This reduces
## the reference to a 2D numerical integral per pair of rings (same sorted
## construction as moment_of_inertia_index()'s), via Gauss-Legendre
## quadrature on each ring's own radial marginal.
##
## ACTUAL VALUE: the same triangulated-mesh quadrature-cloud idea as
## convexity_index(), but simpler - no line-clipping needed, just Euclidean
## distances between quadrature points (tri_quad_points()), summed over all
## pairs including same-triangle ones (which approximate D's self-term; the
## true zero-measure diagonal contributes exactly 0, not an approximation
## error).

## -- special functions and quadrature, local to this file -----------------

#' Complete elliptic integral of the second kind via the AGM algorithm.
#' @param k numeric vector in `[0, 1]`
#' @return E(k), same length as k
#' @noRd
.ellip_E <- function(k) {
    k <- pmin(pmax(k, 0), 1)
    out <- numeric(length(k))
    near1 <- k > 1 - 1e-9
    out[near1] <- 1   # E(1) = 1 exactly; AGM below is 0/Inf indeterminate there
    idx <- which(!near1)
    if (length(idx) == 0) return(out)

    a <- rep(1, length(idx)); b <- sqrt(1 - k[idx]^2); c <- k[idx]
    sum_term <- 0.5 * c^2
    pow2 <- 0.5
    for (i in seq_len(30)) {
        a_new <- (a + b) / 2
        b_new <- sqrt(a * b)
        c <- (a - b) / 2
        a <- a_new; b <- b_new
        pow2 <- pow2 * 2
        sum_term <- sum_term + pow2 * c^2
        if (max(c) < 1e-16) break
    }
    out[idx] <- (pi / (2 * a)) * (1 - sum_term)
    out
}

#' Mean distance between 2 points at radii r1, r2 with uniform angle
#' between them - the angular part of the annulus reference integral.
#' @param r1,r2 numeric vectors (recycled against each other)
#' @return same length as the longer of r1, r2
#' @noRd
.mean_chord <- function(r1, r2) {
    s <- r1 + r2
    k <- ifelse(s > 0, pmin(2 * sqrt(r1 * r2) / s, 1), 0)
    (2 / pi) * s * .ellip_E(k)
}

#' m-point Gauss-Legendre nodes/weights on `[-1, 1]` via the Golub-Welsch
#' eigendecomposition of the Jacobi matrix.
#' @param m number of nodes
#' @return list(x, w)
#' @noRd
.gauss_legendre <- function(m) {
    if (m == 1) return(list(x = 0, w = 2))
    k <- seq_len(m - 1)
    beta <- k / sqrt(4 * k^2 - 1)
    J <- matrix(0, m, m)
    J[cbind(seq_len(m - 1), 2:m)] <- beta
    J[cbind(2:m, seq_len(m - 1))] <- beta
    eig <- eigen(J, symmetric = TRUE)
    ord <- order(eig$values)
    list(x = eig$values[ord], w = 2 * eig$vectors[1, ord]^2)
}

#' Gauss-Legendre nodes/weights transformed to radii in `[r_lo, r_hi]`,
#' weighted by that ring's own radial density f(r) = 2r/(r_hi^2 - r_lo^2)
#' (the radius marginal of a uniform annulus).
#' @param r_lo,r_hi ring bounds
#' @param gl a .gauss_legendre() result, reused across rings
#' @return list(r, p) - p sums to 1
#' @noRd
.radial_nodes <- function(r_lo, r_hi, gl) {
    half <- (r_hi - r_lo) / 2
    mid  <- (r_hi + r_lo) / 2
    r    <- half * gl$x + mid
    dens <- 2 * r / (r_hi^2 - r_lo^2)
    p    <- gl$w * half * dens
    list(r = r, p = p / sum(p))
}

## -- same-triangle ("self") term ------------------------------------------
##
## The cross term between two DIFFERENT triangles is well approximated by
## a handful of quadrature points each (tri_quad_points()), since it's
## averaging over two whole separate regions. The SAME-triangle term
## (mean distance between 2 points both inside one triangle) is not: with
## only 1 or 3 fixed points, most or all same-triangle "pairs" coincide,
## understating it by 20-100% for typically-sized CDT triangles - not a
## negligible correction at the mesh sizes this package actually uses.
## There's no closed form for a general (scalene) triangle, so instead:
## recursively split the triangle via medial subdivision (connect edge
## midpoints -> 4 similar sub-triangles, each 1/4 the area) and take the
## mean pairwise distance among the sub-triangle centroids. Verified by
## Monte Carlo to converge geometrically (~4x smaller error per level);
## depth 4 (256 centroids) gets well under 1% relative error.
##
## DEPTH IS AREA-ADAPTIVE, NOT FIXED, same idea and same
## .adaptive_tri_depth() helper (utils.R) as
## radial_concentration_index()'s own point-cloud subdivision. Unlike
## that O(4^depth) point cloud, this self-term is O((4^depth)^2) PER
## TRIANGLE - .tri_self_mean_distance() computes the full pairwise
## distance matrix among a triangle's own centroids - so giving every
## triangle the mesh's worst-case depth wastes even more here: quartering
## a triangle's depth quarters its point count but cuts its own
## self-distance cost sixteenfold. .mesh_span_index() computes each
## triangle's own depth once (anchored so the mesh's largest triangle
## matches the old fixed-depth-4 behaviour exactly, no regression there)
## and passes it down instead of hardcoding depth = 4 for every triangle.

#' @param v 3x2 matrix of a triangle's vertex coordinates
#' @return list of 4 3x2 vertex matrices (medial subdivision)
#' @noRd
.subdivide_tri <- function(v) {
    m_ab <- (v[1, ] + v[2, ]) / 2
    m_bc <- (v[2, ] + v[3, ]) / 2
    m_ca <- (v[3, ] + v[1, ]) / 2
    list(rbind(v[1, ], m_ab, m_ca), rbind(v[2, ], m_bc, m_ab),
         rbind(v[3, ], m_ca, m_bc), rbind(m_ab, m_bc, m_ca))
}

#' Mean distance between 2 random points within one triangle, approximated
#' via the centroids of its depth-k medial subdivision (see file header).
#' @param v 3x2 matrix of a triangle's vertex coordinates
#' @param depth subdivision depth (4^depth centroids) - .mesh_span_index()
#'   passes this triangle's own area-adaptive depth, not a fixed value
#' @return the mean pairwise distance, in the same length units as v
#' @noRd
.tri_self_mean_distance <- function(v, depth = 4) {
    tris <- list(v)
    for (i in seq_len(depth)) tris <- unlist(lapply(tris, .subdivide_tri), recursive = FALSE)
    cen <- t(vapply(tris, colMeans, numeric(2)))
    mean(as.matrix(dist(cen)))
}

## -- reference D -----------------------------------------------------------

#' @noRd
.disk_reference_D <- function(area) {
    (128 / (45 * pi)) * sqrt(area / pi)
}

#' Weighted reference: sort triangles by density descending (same
#' construction as moment_of_inertia_index()'s rings), group them into
#' coarse rings bounded by an internal density-ratio cap rather than one
#' ring per triangle (see below for why), then average .mean_chord() over
#' every pair of rings' own radial marginals - flattened into one radial
#' "point cloud" so it's a single vectorised sum, the same trick
#' .mesh_span_index() uses for the actual value.
#'
#' RING COUNT: the reference is a property of the weight-vs-radius
#' profile, not of how many CDT triangles the mesh happened to produce -
#' unlike the polygon's own actual D, which genuinely needs one term per
#' triangle-pair (that's exactly the cost deterministic = FALSE exists to
#' avoid on a big mesh). One ring PER TRIANGLE instead would give
#' gl_order * n_triangles quadrature nodes and an O((gl_order *
#' n_triangles)^2) pairwise sum for the reference alone - fine for a
#' small mesh, catastrophic on a real one: a couple thousand triangles
#' already means tens of millions of matrix cells; tens of thousands of
#' triangles is billions, an out-of-memory crash regardless of
#' `deterministic`, since this reference is computed unconditionally
#' whenever `weight` is supplied (verified: one ring per triangle crashes
#' past 18GB on a 2169-triangle real-world mesh, computing nothing else).
#'
#' WHY RATIO-BOUNDED, NOT RANK/AREA/WEIGHT BINNING: a binning rule that
#' groups a FIXED number of triangles per ring (by rank, by area, or by
#' weight share, all tried and rejected here) can be defeated by an
#' adversarial input - a triangle whose density is wildly higher than its
#' immediate neighbours' gets merged in with them anyway once its
#' area/weight/rank share is exhausted, diluting its concentrated mass
#' across a ring far wider than it should occupy and inflating D_ref
#' (verified: a single 1e6x-weighted triangle inflated D_ref enough to
#' push index above 1, under both area-based and pure rank-based
#' binning). Since input is already density-sorted, the fix bounds what's
#' ALLOWED into one ring directly: walk the sorted triangles and close
#' the current ring - starting a new one - whenever EITHER its own top
#' density is more than `max_ratio` times the next triangle's, OR it's
#' used up its fair share of however many triangles remain over however
#' many rings remain. An outlier, however extreme, closes its own ring
#' immediately and cheaply via the ratio test; the ring budget that would
#' otherwise have gone to it is then simply available for the rest,
#' since the "fair share" is recomputed from what's actually left after
#' every close, not fixed once up front from the raw triangle count -
#' this is what keeps the well-behaved bulk of realistic data from being
#' degraded by however many outliers happen to exist. The
#' remaining_rings > 1L guard caps the total ring count at max_rings
#' exactly: the final ring always absorbs whatever's left, however it
#' compares in density, since there's no budget remaining to split it
#' further. Below max_rings triangles, every triangle keeps its own ring
#' exactly as before - no coarsening, no behaviour change.
#' @param tri_area numeric vector, each triangle's own area
#' @param weight numeric vector, one entry per triangle, need not sum to 1
#' @param gl_order Gauss-Legendre nodes per ring (fixed, not user-facing)
#' @param max_rings maximum number of coarse rings, regardless of
#'   triangle count (fixed, not user-facing)
#' @param max_ratio a ring's own top triangle's density may be at most
#'   this many times its bottom triangle's before a new ring starts,
#'   budget permitting (fixed, not user-facing)
#' @return D_ref
#' @noRd
.annulus_reference_D <- function(tri_area, weight, gl_order = 8, max_rings = 100, max_ratio = 4) {
    W    <- .normalize_weight(weight)
    rho  <- W / tri_area
    ord  <- order(rho, decreasing = TRUE)
    tri_area <- tri_area[ord]
    W        <- W[ord]
    rho      <- rho[ord]
    n <- length(tri_area)

    if (n > max_rings) {
        # greedy ring-building with a dynamically re-budgeted target size:
        # a ring closes when EITHER its internal density ratio would
        # exceed max_ratio OR it's used up its fair share of the
        # remaining triangles over the remaining ring budget, whichever
        # comes first. Re-deriving the fair share after every close (from
        # however many triangles and rings are actually left) is what
        # keeps this from degrading the well-behaved bulk of the data:
        # outliers close their own ring early and cheaply, and the
        # budget that would otherwise have been "wasted" on them is
        # simply redistributed across the rest, rather than warping a
        # FIXED equal-count split computed once up front. remaining_rings
        # > 1L guards the total ring count at max_rings exactly - the
        # last ring always absorbs whatever's left, however it compares
        # in density, since there's no budget left to split it further.
        bin_id    <- integer(n)
        cur_bin   <- 1L
        bin_start <- 1L
        top_rho   <- rho[1]
        remaining_rings <- max_rings
        for (i in 2:n) {
            remaining_tri <- n - bin_start + 1L
            target_size   <- ceiling(remaining_tri / remaining_rings)
            if (remaining_rings > 1L &&
                (top_rho / rho[i] > max_ratio || (i - bin_start) >= target_size)) {
                cur_bin <- cur_bin + 1L
                remaining_rings <- remaining_rings - 1L
                bin_start <- i
                top_rho <- rho[i]
            }
            bin_id[i] <- cur_bin
        }
        tri_area <- as.numeric(tapply(tri_area, bin_id, sum))
        W        <- as.numeric(tapply(W, bin_id, sum))
    }

    r_hi <- sqrt(cumsum(tri_area) / pi)
    r_lo <- c(0, r_hi[-length(r_hi)])

    gl    <- .gauss_legendre(gl_order)
    nodes <- Map(.radial_nodes, r_lo, r_hi, MoreArgs = list(gl = gl))
    r_all <- unlist(lapply(nodes, `[[`, "r"))
    p_all <- unlist(Map(function(nd, w_i) nd$p * w_i, nodes, W))

    sum(outer(p_all, p_all) * outer(r_all, r_all, .mean_chord))
}

## -- the span index itself ---------------------------------------------

#' Deterministic (deterministic = TRUE) mean pairwise distance over a fixed
#' mesh - shared engine behind span_index().
#' @param tri sf data frame of CDT triangles, with an `area` column
#' @param n_quad quadrature points per triangle (1 or 3)
#' @param weight optional numeric vector, length nrow(tri)
#' @return list(index, D, D_ref, area, total_weight, triangles)
#' @noRd
.mesh_span_index <- function(tri, n_quad, weight = NULL) {
    n <- if (is.null(tri)) 0 else nrow(tri)
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }
    if (n == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_, area = NA_real_,
                    total_weight = NA_real_, triangles = tri))
    }

    q <- n_quad
    n_points <- n * q
    if (n_points^2 > 4e6) {
        warning(
            "Large mesh: ", n, " triangles, ", n_points, " quadrature points at n_quad = ",
            q, " (deterministic = TRUE is O((n*n_quad)^2), and weight (if supplied) adds an ",
            "O(n^2) reference computation) - this can be slow. Try n_quad = 1, ",
            "deterministic = FALSE, or prep = prepare_polygon() with a coarser mesh."
        )
    }

    area      <- sum(tri$area)
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri$area else .normalize_weight(weight)

    qmat <- if (q == 1) {
        st_coordinates(st_centroid(st_geometry(tri)))[, 1:2, drop = FALSE]
    } else {
        do.call(rbind, lapply(st_geometry(tri), function(g) {
            tri_quad_points(st_coordinates(g)[1:3, 1:2, drop = FALSE])
        }))
    }
    pt_w  <- rep(w / q, each = q)
    same  <- outer(rep(seq_len(n), each = q), rep(seq_len(n), each = q), "==")
    cross_sum <- sum((outer(pt_w, pt_w) * as.matrix(dist(qmat)))[!same])

    # depth adapts to each triangle's own CONTRIBUTION to self_sum
    # (~ w[i]^2 * sqrt(area_i), since self-distance scales with a
    # triangle's own linear size), not to its raw area - weighting can
    # make a physically tiny triangle carry nearly all the mass, and area
    # alone would then wrongly collapse its self-distance to 0 exactly
    # (mean pairwise distance among a single depth-0 centroid), the same
    # adversarial failure mode .annulus_reference_D() above had to guard
    # against for the same reason (see its own comments)
    depth_i <- .adaptive_tri_depth(w^2 * sqrt(tri$area), max_depth = 4L)
    self_sum <- sum(vapply(seq_len(n), function(i) {
        v <- st_coordinates(st_geometry(tri)[i])[1:3, 1:2, drop = FALSE]
        w[i]^2 * .tri_self_mean_distance(v, depth = depth_i[i])
    }, numeric(1)))

    D <- (cross_sum + self_sum) / sum(w)^2

    D_ref <- if (is.null(weight)) .disk_reference_D(area) else .annulus_reference_D(tri$area, weight)
    index <- D_ref / D

    list(index = index, D = D, D_ref = D_ref, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total, triangles = tri)
}

## Stochastic (deterministic = FALSE): draws `n_lines` independent point pairs and
## averages their distance directly - no line construction/clipping needed
## (unlike convexity_index()'s Monte Carlo mode), since span_index() only
## ever needs the distance itself. Named `n_lines`, not `n_pairs`, to match
## convexity_index()'s argument: conceptually a line's length is exactly
## the distance between its two endpoints, even though no LINESTRING
## geometry is ever actually built here - and sharing the name means
## shape_indices()'s `...` sets both functions' Monte Carlo sample count
## with one argument instead of two.

#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param n_lines number of random point pairs, or a function(n_tri)
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon()
#' @param seed optional RNG seed
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`; NULL samples uniformly
#' @param points optional pre-drawn (2*n_lines)x2 coordinate matrix, from
#'   shape_indices()'s point-sharing mechanism (see its own comments) -
#'   when supplied, sampling is skipped entirely
#' @return list(index, D, D_ref, area, total_weight, triangles)
#' @noRd
.random_pair_span_index <- function(poly, n_lines, prep, seed, weight = NULL, points = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly)
    poly_geom <- prep$poly
    tri       <- prep$tri
    n_tri     <- if (is.null(tri)) 0 else nrow(tri)

    if (!is.null(weight)) {
        if (n_tri == 0) {
            stop("`weight` needs a triangle mesh to sample from, but this polygon ",
                 "triangulated to no triangles.")
        }
        if (length(weight) != n_tri) {
            stop("`weight` must have one entry per triangle (", n_tri, "), got ", length(weight), ".")
        }
    }

    if (is.function(n_lines)) n_lines <- max(1L, round(n_lines(n_tri)))

    if (n_tri >= 2) {
        n_pairs_deterministic <- choose(n_tri, 2)
        if (n_lines >= n_pairs_deterministic / 2) {
            warning(sprintf(
                paste("n_lines (%d) is not substantially lower than the %d triangle-pairs",
                      "that deterministic = TRUE (%d triangles) would evaluate for this same",
                      "polygon; deterministic = FALSE is meant as a cheaper approximation for",
                      "meshes too large to enumerate exhaustively - consider",
                      "deterministic = TRUE instead, or a smaller n_lines."),
                n_lines, n_pairs_deterministic, n_tri))
        }
    }

    poly_u <- st_union(poly_geom)
    area   <- sum(tri$area)
    needed <- 2 * n_lines

    if (!is.null(points)) {
        if (nrow(points) < needed) {
            stop("`points` has ", nrow(points), " rows but ", needed, " are needed for n_lines = ", n_lines, ".")
        }
        coords <- points[seq_len(needed), , drop = FALSE]
    } else if (!is.null(weight)) {
        if (!is.null(seed)) set.seed(seed)
        coords <- .sample_weighted_points(tri, weight, needed)
    } else {
        if (!is.null(seed)) set.seed(seed)
        pts <- st_sample(poly_u, size = needed, type = "random", exact = TRUE)
        tries <- 0
        while (length(pts) < needed && tries < 10) {
            pts <- c(pts, st_sample(poly_u, size = needed, type = "random", exact = TRUE))
            tries <- tries + 1
        }
        if (length(pts) < needed) {
            warning("Could not sample enough interior points; index is not defined.")
            return(list(index = NA_real_, D = NA_real_, D_ref = NA_real_, area = area,
                        total_weight = area, triangles = tri))
        }
        pts    <- pts[seq_len(needed)]
        coords <- st_coordinates(pts)[, 1:2, drop = FALSE]
    }

    x1 <- coords[seq(1, needed, 2), , drop = FALSE]
    x2 <- coords[seq(2, needed, 2), , drop = FALSE]
    D  <- mean(sqrt((x1[, 1] - x2[, 1])^2 + (x1[, 2] - x2[, 2])^2))

    raw_total <- if (is.null(weight)) NULL else sum(weight)
    D_ref <- if (is.null(weight)) .disk_reference_D(area) else .annulus_reference_D(tri$area, weight)
    index <- D_ref / D

    list(index = index, D = D, D_ref = D_ref, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total, triangles = tri)
}

#' Span index: mean pairwise interior distance vs. an equal-area circle
#'
#' D_ref/D, in (0, 1], where D is the mean distance between two random
#' interior points and D_ref is the same quantity for the reference shape
#' (a circle, unweighted; a concentric annulus, weighted) - both are
#' provable minimisers of D. Distinct from moment_of_inertia_index(),
#' which a squared-distance version of this index would just collapse to.
#' See `vignette("d-understanding-span-index")` for the derivation and
#' both proofs.
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param deterministic if TRUE (default), compute over the polygon's own
#'   fixed CDT quadrature mesh - O(n^2). If FALSE, Monte Carlo over
#'   `n_lines` random interior point pairs instead - O(n_lines), useful
#'   when the mesh is too large. `n_quad` is unused (and an error if
#'   passed explicitly) when deterministic = FALSE.
#' @param n_quad quadrature points per triangle when deterministic = TRUE
#'   (1 or 3), used for the cross-triangle term only (the same-triangle
#'   term is always computed separately, by subdividing each triangle).
#'   3 (default) is more accurate; drop to 1 on meshes of more than a few
#'   hundred triangles if that matters.
#' @param n_lines number of random point pairs to sample when
#'   deterministic = FALSE (default 3000), or a function(n_tri) - same
#'   argument name as convexity_index()'s Monte Carlo mode, so one value
#'   passed through shape_indices() sets every mesh index's sample count
#'   at once, though no line geometry is actually built here. Warns if
#'   not substantially lower than what deterministic = TRUE would need
#'   for the same polygon.
#' @param seed optional RNG seed, only used when deterministic = FALSE.
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating when also calling other indices
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, substituting for that triangle's own area/mass throughout.
#'   NULL (default) compares against a circle; otherwise the reference
#'   becomes a concentric annulus matching this weight's own histogram.
#' @param points optional pre-drawn (2*n_lines)x2 coordinate matrix, only
#'   meaningful when deterministic = FALSE - primarily for
#'   shape_indices()'s internal point-sharing mechanism (draw one sample
#'   and reuse it across convexity/span/radial_concentration instead of
#'   each drawing independently), not typically supplied directly.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return list(index, D, D_ref, area, total_weight, triangles).
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' span_index(wake)$index
#' # Camden: a long, narrow county - far from a circle
#' span_index(nc[nc$NAME == "Camden", ])$index
#'
#' # deterministic = FALSE: Monte Carlo estimate instead of the exhaustive mesh
#' span_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#'
#' # weight: substitutes for each triangle's own area - weighting by the
#' # triangle's own area exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' span_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
span_index <- function(poly, deterministic = TRUE, n_quad = 3, n_lines = 3000, seed = NULL,
                        prep = NULL, weight = NULL, points = NULL, simplify_tolerance = NULL) {
    n_quad_given <- !missing(n_quad)
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)

    if (!deterministic) {
        if (n_quad_given) {
            stop("`n_quad` selects the quadrature refinement used by deterministic = TRUE; ",
                 "it has no meaning for deterministic = FALSE, which samples fresh random ",
                 "points rather than quadrature points on a fixed mesh. Drop the ",
                 "`n_quad` argument, or set deterministic = TRUE to use it.")
        }
        return(.random_pair_span_index(poly, n_lines = n_lines, prep = prep, seed = seed,
                                        weight = weight, points = points))
    }
    if (!is.null(points)) {
        stop("`points` (a pre-drawn sample) has no meaning for deterministic = TRUE, which ",
             "computes over the full quadrature mesh, not a random sample. Drop `points`, ",
             "or set deterministic = FALSE to use it.")
    }

    if (!n_quad %in% c(1, 3)) {
        stop("n_quad must be 1 (centroid only) or 3 (Hammer-Stroud rule); other ",
             "quadrature orders aren't implemented.")
    }

    .mesh_span_index(prep$tri, n_quad = n_quad, weight = weight)
}

#' span_index() for every row of an sf data frame
#'
#' Each row indexed independently and unweighted.
#' @param x an sf data frame
#' @param ... passed to span_index() for every row
#' @return `x` with one new column, span_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- span_index_sf(nc[1:5, ])
#' res$span_index
#' @export
span_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) span_index(geoms[i], ...))
    x$span_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
