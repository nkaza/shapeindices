## =========================================================================
## 10. Directional balance index: mean resultant length of the mass
## distribution's own angular spread about its centroid
## =========================================================================
##
## R = |(1/W) integral rho(x) * exp(i*theta(x)) dx|, where theta(x) =
## atan2(y-Gy, x-Gx) is the BEARING of point x as seen from the shape's own
## mass centroid G - not its distance, which is discarded entirely. G is
## already the mass centroid by definition, so integral rho(x)*(x-G) dx = 0
## trivially always - this index is NOT that. e^{i*theta} keeps only the
## bearing, discarding how far out the mass extends, so it measures whether
## the shape's own AREA is evenly spread across bearings or concentrated
## toward some of them. index = 1 - R, in [0, 1] (see below for why this
## bound needs no rearrangement inequality or PSD argument - simpler than
## either of this package's other two "why is it bounded" techniques).
##
## NAME: deliberately not "moment" anything, unlike moment_of_inertia_index()/
## moment_isotropy_index() - those integrate a genuine SECOND moment
## (position squared) and earn that word honestly. This index only uses
## bearing, not position magnitude at all, so calling it a moment would be
## the same category of naming mistake this package already made once with
## "eccentricity_index" (see moment-isotropy-index.R's own file header) -
## overclaiming a technical meaning the formula doesn't actually have.
## "Balance", not "uniformity", for the same reason: R = 0 does NOT mean
## the mass is spread uniformly across every bearing - it only means the
## directional pulls CANCEL. A symmetric dumbbell (equal mass due north and
## due south of G, joined by a thin neck) scores R = 0 exactly, index = 1,
## identical to a disk - see BLIND SPOT below.
##
## BOUND: |R| <= 1 by the triangle inequality on a complex/vector-valued
## expectation - |E[e^{i*theta}]| <= E|e^{i*theta}| = E[1] = 1 - holds for
## ANY probability measure over bearings, continuous or a finite weighted
## point cloud alike, with no further argument needed (unlike
## convexity_index()/span_index()/radial_concentration_index()'s
## rearrangement-inequality reference-shape proofs, or
## moment_isotropy_index()'s PSD-matrix argument). Because this holds for
## ANY finite weighted sample, not just in a continuum limit, the DISCRETE
## point-cloud/Monte-Carlo estimate below is exactly bounded in [0, 1] too -
## no clamping is ever needed to enforce the bound, only ordinary
## quadrature/sampling error affects how close the estimate lands to the
## true continuous value.
##
## BLIND SPOT (a feature, not purely a limitation): any shape with 2-fold+
## rotational symmetry about its own centroid scores R = 0 exactly (index =
## 1) - broader than "= 1 iff a disk". A symmetric dumbbell/two-lobed shape
## scores index = 1 here while scoring low on moment_of_inertia_index()/
## span_index()/radial_concentration_index()/moment_isotropy_index() (all
## elongated/dispersed along the shared axis). Paired together, the two
## families distinguish "symmetric two-sided reach" (fine here, bad on the
## others) from "one-sided tentacle/appendage" (bad on BOTH) - plausibly
## relevant to redistricting analysis, where a district reaching out in one
## direction only reads differently from one reaching out symmetrically in
## two.
##
## GENUINELY NEW INFORMATION vs the package's other nine indices: a
## first-harmonic/directional-PREFERENCE measure - does the mass lean
## toward one bearing - structurally different from moment_isotropy_index()'s
## second-moment/axis-preference measure - does the mass prefer an axis,
## agnostic to which end of it. A shape "pointing north" and its mirror
## image "pointing south" score identically on moment_isotropy_index() but
## oppositely (as measured by `mean_angle`, not `index` itself, which only
## reports magnitude) on this one.
##
## NO CLOSED FORM: unlike Ixx/Iyy/Ixy, integral e^{i*theta(x)} dA over a
## triangle isn't polynomial in the vertex coordinates - dividing by radius
## breaks the polynomial structure that makes the moment-of-inertia family
## exact. So both deterministic = TRUE and deterministic = FALSE below are
## genuinely approximate, the same situation convexity_index()/span_index()/
## radial_concentration_index() are in (not moment_of_inertia_index()/
## moment_isotropy_index(), which are exact regardless of mode).
##
## DETERMINISTIC MODE reuses radial_concentration_index()'s own
## area-adaptive subdivision point cloud (.radial_point_cloud(), see its
## own header in radial-concentration-index.R) unchanged - a single point
## per triangle (its centroid) under-resolves e^{i*theta(x)} across a large
## triangle spanning a wide range of bearings as seen from G, the same
## generic single-point-quadrature-is-too-coarse issue
## radial_concentration_index() already has a fix for, even though the
## MECHANISM there (Jensen's inequality on a convex |x-c|) doesn't directly
## apply to e^{i*theta}, which isn't simply convex or concave in x. The fix
## is the same regardless: finer subdivision converges to the true
## integral, verified empirically in this package's tests rather than
## assumed to transfer.
##
## MONTE CARLO MODE (deterministic = FALSE) draws points directly from the
## weighted density (reusing convexity_index()'s .sample_weighted_points(),
## the same mechanism radial_concentration_index()'s own Monte Carlo mode
## uses) and computes theta relative to the EXACT centroid G (found in
## closed form regardless of mode - never randomised itself), then
## R_hat = |mean(e^{i*theta_i})|. KNOWN BIAS: R_hat is a magnitude of a
## sample-mean vector, a nonlinear transform, so by Jensen's inequality
## E[R_hat] >= |E[sample mean]| = R - for a truly balanced shape (R = 0
## exactly), R_hat is essentially never exactly 0 in a finite sample, so
## index_MC is systematically slightly BELOW the true index, shrinking as
## n_lines grows (~1/sqrt(n_lines) at worst; measured up to -0.008 at
## n_lines = 3000 on real near-balanced Urban Area data).
##
## WHY NO CORRECTION IS APPLIED (a decision, made on evidence, 2026-07):
## no unbiased estimator of 1 - R exists at all - provable via the
## two-point bearing family (mass p at angle 0, 1-p at pi gives
## R = |2p-1|; any estimator's expectation is a polynomial in p, |2p-1|
## has a kink at p = 1/2), the same obstruction that makes sample
## standard deviation biased while sample variance isn't. The kink sits
## at R = 0, exactly where the bias is worst, and asymptotic corrections
## carry a 1/R term that diverges there too. What IS exactly estimable is
## R^2: since E[e^{i(theta_j - theta_k)}] = R^2 for independent draws,
## (n*R_hat^2 - 1)/(n - 1) is exactly unbiased for R^2, any distribution,
## any n >= 2 - verified on 4 real UA files (bias <= 4e-4 across 8
## weighting cases, 500 reps at n = 3000, vs -0.001..-0.008 for the naive
## index). The index nevertheless stays defined as 1 - R, not 1 - R^2,
## because on real urban areas (R in 0.01..0.21) the squared scale crams
## everything into [0.95, 1] - measured spans: 1-R over 0.79..0.99 vs
## 1-R^2 over 0.956..0.9999 on the same shapes - and an unbiased-at-the-
## boundary estimator MUST sometimes exceed 1 (measured ~50% of runs at
## exact balance), which would clash with the hard [0,1] bound this
## package promises. So: index = 1 - R (standard circular variance,
## readable scale, exact bound), bias documented in the roxygen, and the
## man page teaches users the one-line unbiased 1 - R^2 formula to apply
## to the returned R themselves when batch comparability matters.

#' Deterministic (deterministic = TRUE) directional balance over a fixed
#' mesh - shared engine behind directional_balance_index() and
#' shape_indices_sf(byrow = FALSE)'s combined mesh.
#' @param tri sf data frame of CDT triangles, with an `area` column
#' @param weight optional numeric vector, length nrow(tri)
#' @return list(index, R, mean_angle, area, total_weight, centroid, triangles)
#' @noRd
.mesh_directional_balance_index <- function(tri, weight = NULL) {
    n <- if (is.null(tri)) 0 else nrow(tri)
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }
    if (n == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, R = NA_real_, mean_angle = NA_real_, area = NA_real_,
                    total_weight = NA_real_, centroid = NULL, triangles = tri))
    }

    area      <- sum(tri$area)
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri$area else .normalize_weight(weight)

    G     <- .mass_centroid(st_geometry(tri), w)
    cloud <- .radial_point_cloud(tri, w)
    res   <- .resultant_from_cloud(cloud$p, cloud$w, G)

    list(index = 1 - res$R, R = res$R, mean_angle = res$mean_angle, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total,
         centroid = st_sfc(st_point(G), crs = st_crs(tri)), triangles = tri)
}

#' Weighted mean resultant length of a point cloud's bearing from a fixed
#' centre - the discrete estimator both deterministic = TRUE (subdivision
#' cloud) and deterministic = FALSE (Monte Carlo sample) reduce to; see file
#' header for the bound that makes this exact for ANY finite weighted
#' cloud, not just in the continuum limit.
#' @param p Nx2 coordinate matrix
#' @param w numeric vector, length N, weights (need not sum to 1)
#' @param G numeric(2), the centre bearings are measured from
#' @return list(R, mean_angle)
#' @noRd
.resultant_from_cloud <- function(p, w, G) {
    theta <- atan2(p[, 2] - G[2], p[, 1] - G[1])
    Wtot  <- sum(w)
    Rx <- sum(w * cos(theta)) / Wtot
    Ry <- sum(w * sin(theta)) / Wtot
    list(R = sqrt(Rx^2 + Ry^2), mean_angle = atan2(Ry, Rx))
}

#' Monte Carlo (deterministic = FALSE) directional balance - internal
#' engine for directional_balance_index()'s deterministic = FALSE. See that
#' function's own roxygen for the user-facing parameter reference;
#' parameters here are the same, without defaults.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param n_lines number of random points to sample, or a function(n_tri)
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon()
#' @param seed optional RNG seed
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`; NULL samples proportional to triangle area
#' @param points optional pre-drawn Nx2 coordinate matrix, from
#'   shape_indices()'s point-sharing mechanism - when supplied, sampling is
#'   skipped and every row is used as one point in the cloud
#' @return list(index, R, mean_angle, area, total_weight, centroid, triangles)
#' @noRd
.random_point_directional_balance_index <- function(poly, n_lines, prep, seed, weight = NULL, points = NULL) {
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
        return(list(index = NA_real_, R = NA_real_, mean_angle = NA_real_, area = NA_real_,
                    total_weight = NA_real_, centroid = NULL, triangles = tri))
    }

    if (is.function(n_lines)) n_lines <- max(1L, round(n_lines(n_tri)))

    area      <- sum(tri$area)
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri$area else .normalize_weight(weight)

    G <- .mass_centroid(st_geometry(tri), w)

    if (!is.null(points)) {
        coords <- points
    } else {
        if (!is.null(seed)) set.seed(seed)
        coords <- .sample_weighted_points(tri, weight %||% tri$area, n_lines)
    }

    res <- .resultant_from_cloud(coords, rep(1, nrow(coords)), G)

    list(index = 1 - res$R, R = res$R, mean_angle = res$mean_angle, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total,
         centroid = st_sfc(st_point(G), crs = st_crs(tri)), triangles = tri)
}

#' Directional balance index of a (multi)polygon's mass distribution
#'
#' `1 - R`, where `R` is the mean resultant length of the (weighted)
#' BEARING of the shape's own interior mass, viewed from its own mass
#' centroid, in `[0, 1]`. `1` means the mass is directionally balanced -
#' no bearing pulls harder than the rest; lower values mean the mass
#' leans toward one direction (a one-sided appendage, an off-centre
#' lobe). One caveat to know before using it: balanced does not mean
#' round - a symmetric dumbbell (equal lobes on opposite sides) scores 1
#' here, same as a disk, because opposite pulls cancel exactly. See
#' `vignette("g-understanding-directional-balance-index")` for the
#' derivation, the bound proof, and that blind spot worked through.
#' @section Monte Carlo bias, and computing an unbiased estimate yourself:
#'   With `deterministic = FALSE`, the returned `index` is a biased
#'   estimate of `1 - R`: a finite sample's resultant length is never
#'   exactly zero, so it systematically overestimates a small true `R`,
#'   and near-balanced shapes score slightly LOW - by up to roughly
#'   `1/sqrt(n_lines)` for an exactly balanced shape (about 0.018 at the
#'   default `n_lines = 3000`), shrinking as the true `R` grows or as
#'   `n_lines` increases. This cannot be corrected away on the `1 - R`
#'   scale: no unbiased estimator of `1 - R` exists for any finite
#'   sample, for the same structural reason a sample standard deviation
#'   cannot be unbiased even though a sample variance can.
#'
#'   What does admit an exactly unbiased estimate is the closely related
#'   `1 - R^2` (the same quantity on a squared scale: same `[0, 1]`
#'   range, same direction, same shapes scoring exactly 1). Compute it
#'   yourself from the returned `R`:
#'   `1 - (n_lines * R^2 - 1) / (n_lines - 1)`.
#'   This is exactly unbiased for `1 - R^2` - any weight distribution,
#'   any `n_lines >= 2` - and is the right quantity to average or compare
#'   across many near-balanced shapes in a batch. It can land slightly
#'   above 1 for near-balanced shapes; that is a necessary feature, not
#'   an error (no estimator confined to `[0, 1]` can be unbiased at the
#'   boundary), so do not truncate it back to 1 if unbiasedness is the
#'   point.
#'
#'   `deterministic = TRUE` carries no statistical bias at all - only
#'   deterministic quadrature error, which shrinks with mesh resolution.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param deterministic if TRUE (default), compute over the same
#'   area-adaptive subdivision point cloud `radial_concentration_index()`
#'   uses. If FALSE, Monte Carlo over `n_lines` points sampled directly
#'   from the weighted density instead - see the section below on the
#'   small downward bias this carries for near-balanced shapes, and how
#'   to compute an unbiased estimate of the squared-scale variant from
#'   the returned `R`.
#' @param n_lines number of random points to sample when
#'   deterministic = FALSE (default 3000), or a function(n_tri) - same
#'   argument name as the package's other Monte Carlo mesh indices, so
#'   one value passed through shape_indices() sets every sample count at
#'   once.
#' @param seed optional RNG seed, only used when deterministic = FALSE.
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating when also calling other indices
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, substituting for that triangle's own area/mass throughout.
#'   NULL (default) reproduces the unweighted (area-based) index exactly.
#' @param points optional pre-drawn Nx2 coordinate matrix, only meaningful
#'   when deterministic = FALSE - primarily for shape_indices()'s internal
#'   point-sharing mechanism, not typically supplied directly.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return list(index, R, mean_angle, area, total_weight, centroid,
#'   triangles). `R` is the mean resultant length the index is built
#'   from (`index = 1 - R`) - kept in the output so the unbiased
#'   squared-scale estimate described above can be computed from it.
#'   `mean_angle` (radians, `atan2()` convention) is the bearing of the
#'   net directional pull - which way the mass leans, if it leans at
#'   all; not meaningful when `R` is near 0 (any direction is roughly as
#'   good as any other). `centroid` is the exact mass centroid, always
#'   found in closed form regardless of `deterministic`.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' directional_balance_index(wake)$index
#'
#' # deterministic = FALSE: Monte Carlo estimate - a seed makes it reproducible
#' directional_balance_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#'
#' # weight: substitutes for each triangle's own area - weighting by the
#' # triangle's own area exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' directional_balance_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
directional_balance_index <- function(poly, deterministic = TRUE, n_lines = 3000, seed = NULL,
                                       prep = NULL, weight = NULL, points = NULL,
                                       simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)

    if (!deterministic) {
        return(.random_point_directional_balance_index(poly, n_lines = n_lines, prep = prep,
                                                         seed = seed, weight = weight, points = points))
    }
    if (!is.null(points)) {
        stop("`points` (a pre-drawn sample) has no meaning for deterministic = TRUE, which ",
             "computes over the full subdivision point cloud, not a random sample. Drop ",
             "`points`, or set deterministic = FALSE to use it.")
    }

    .mesh_directional_balance_index(prep$tri, weight = weight)
}

#' directional_balance_index() for every row of an sf data frame
#'
#' Each row indexed independently and unweighted.
#' @param x an sf data frame
#' @param ... passed to directional_balance_index() for every row
#' @return `x` with one new column, directional_balance_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- directional_balance_index_sf(nc[1:5, ])
#' res$directional_balance_index
#' @export
directional_balance_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) directional_balance_index(geoms[i], ...))
    x$directional_balance_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
