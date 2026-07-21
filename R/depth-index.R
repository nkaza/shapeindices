## =========================================================================
## 8. Depth index: mean distance to the boundary vs an equal-area circle
## =========================================================================
##
## d(s) = min_{b in boundary(P)} |s-b|, the distance transform - how far an
## interior point sits from the nearest point on the shape's OWN boundary
## (unlike radial_concentration_index()'s D1, which is distance to a
## SINGLE interior reference point, the geometric median). depth_bar(rho) =
## (1/W) integral rho(s)*d(s) ds, the rho-weighted mean of that field.
## index = depth_bar(rho) / depth_bar_ref, in (0, 1] - note the ratio is
## the OPPOSITE way up from every other index in this package (ref/actual
## there; actual/ref here), because the reference shape here MAXIMISES the
## quantity rather than minimising it (see below), so actual <= ref always.
##
## THE DEPTH PROPOSITION (Angel, Parent & Civco 2010, stated without proof;
## proved here via Brunn-Minkowski - see
## vignette("h-understanding-depth-index") for the full argument):
## among all measurable shapes of a given area, the disk maximises
## unweighted mean depth. Sketch: integral_P d(s) ds = integral_0^inf
## |erosion_t(P)| dt (layer-cake identity on the distance-transform field),
## and Brunn-Minkowski bounds each erosion's area by the disk's own erosion
## area at every t simultaneously - so bounding holds by construction for
## holes and multiple parts too, with no separate argument needed the way
## some of this package's other indices require one.
##
## CLOSED-FORM UNWEIGHTED REFERENCE: for a disk of radius R, d(s) = R - r
## (r = distance from centre), and E[r] = 2R/3 for a uniform disk, so
## depth_bar(disk) = R/3.
##
## WEIGHTED REFERENCE - reuses radial_concentration_index()'s own machinery
## exactly, via an algebraic identity, not a new derivation: d(s) = R-r is
## DECREASING in r, so (rearrangement inequality) the arrangement that
## MAXIMISES the weighted mean depth is concentric rings densest at the
## centre - the SAME arrangement .annulus_reference_D1() already computes
## (there, because r itself is increasing, so putting heavy weight at
## small r MINIMISES mean distance-to-centre - same sort order, opposite
## reason). Since sum(rho)=1 after normalising and R is the full disk's
## own outer radius (identical to .annulus_reference_D1()'s own outermost
## ring boundary):
##
##   depth_bar_ref = sum_rings rho_k * (R - mean_r_k) = R - D1_ref
##
## so .annulus_reference_depth() below is just R - .annulus_reference_D1(),
## no new closed-form integral needed. weight = NULL collapses this
## exactly to R - (2/3)R = R/3, matching the disk closed form above.
##
## ALGORITHMIC CHOICES: deterministic = TRUE reuses
## radial_concentration_index()'s own .radial_point_cloud()/
## .radial_point_cloud_array() point cloud unchanged (area-adaptive medial
## subdivision), but the bias it corrects runs the OPPOSITE direction from
## radial_concentration_index()'s own: |s-c| is convex everywhere, so a bare
## centroid systematically UNDERestimates a triangle's contribution there.
## d(s) = min_b|s-b| is a MINIMUM of convex functions, which is not convex
## in general - away from the medial axis (where a single boundary edge is
## unambiguously nearest) it's exactly LINEAR, so a bare centroid has ZERO
## bias; only triangles straddling the medial axis itself (where d(s) is a
## concave "tent", min of two linear pieces) collapse badly, and there a
## bare centroid OVERestimates. Verified directly on a deeply-notched star
## (rich medial-axis structure): unsubdivided centroids overshot the true
## mean depth by ~83%, subdivision landed within a fraction of a percent of
## a 200,000-point Monte Carlo reference (see
## vignette("h-understanding-depth-index")'s own "Algorithmic choices"
## section). deterministic = FALSE reuses
## .sample_weighted_points()/.sample_weighted_points_array() unchanged.
## What's genuinely new here (not needed by any of this package's other
## five mesh indices) is the field itself: st_distance() from the point
## cloud to st_boundary(poly_u) - the actual polygon boundary, not
## anything derivable from triangle centroids/areas alone, so poly_u (not
## just the mesh) is threaded through every engine below. Confirmed
## empirically (a real ~50k-vertex-boundary polygon) that this scales
## sub-linearly in boundary vertex count and roughly linearly in sample
## point count - no chunking/memory-safety mechanism needed the way
## convexity_index()'s pairwise cross-product matrix requires.
##
## BOUND HANDLING: like every Monte Carlo/finite-point-cloud estimate in
## this package, a shape whose TRUE index is exactly 1 (a disk) can show a
## small index > 1 from sampling/subdivision noise - documented as an
## expected residual, not clamped (this package's convention throughout:
## fix the estimator, never clamp the output to hide it).

#' @noRd
.disk_reference_depth <- function(area) sqrt(area / pi) / 3

#' Weighted reference: R - D1_ref, an algebraic identity with
#' radial_concentration_index()'s own annulus reference - see file header.
#' @param tri_area numeric vector, each triangle's own area
#' @param weight numeric vector, one entry per triangle, need not sum to 1
#' @return depth_bar_ref
#' @noRd
.annulus_reference_depth <- function(tri_area, weight) {
    R <- sqrt(sum(tri_area) / pi)
    R - .annulus_reference_D1(tri_area, weight)
}

#' Nx2 coordinate matrix -> sfc of POINT geometries, vectorised (no
#' per-row st_point() loop) via st_multipoint()/st_cast() - the point
#' clouds here (deterministic mode especially, via medial subdivision) can
#' run into the tens of thousands of rows per triangle-heavy mesh.
#' @param mat Nx2 coordinate matrix
#' @param crs the CRS to attach
#' @return sfc of N POINT geometries, in row order
#' @noRd
.coords_to_points <- function(mat, crs) {
    sf::st_cast(sf::st_sfc(sf::st_multipoint(mat), crs = crs), "POINT")
}

#' Deterministic mean-depth-to-boundary over a fixed mesh - shared engine
#' behind depth_index().
#' @param tri sf data frame of CDT triangles, with an `area` column
#' @param poly_u the (multi)polygon `tri` triangulates - depth needs the
#'   actual boundary, not just the mesh (see file header)
#' @param weight optional numeric vector, length nrow(tri)
#' @return list(index, mean_depth, ref_depth, area, total_weight, triangles)
#' @noRd
.mesh_depth_index <- function(tri, poly_u, weight = NULL) {
    n <- if (is.null(tri)) 0 else nrow(tri)
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }
    if (n == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, mean_depth = NA_real_, ref_depth = NA_real_,
                    area = NA_real_, total_weight = NA_real_, triangles = tri))
    }

    area      <- sum(tri$area)
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri$area else .normalize_weight(weight)

    cloud <- .radial_point_cloud(tri, w)
    bnd   <- st_boundary(poly_u)
    d     <- as.numeric(st_distance(.coords_to_points(cloud$p, st_crs(poly_u)), bnd))
    mean_depth <- sum(cloud$w * d) / sum(cloud$w)

    ref_depth <- if (is.null(weight)) .disk_reference_depth(area) else .annulus_reference_depth(tri$area, weight)
    index <- mean_depth / ref_depth

    list(index = index, mean_depth = mean_depth, ref_depth = ref_depth, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total, triangles = tri)
}

## Stochastic (deterministic = FALSE): draws `n_lines` points directly from
## the weighted density (.sample_weighted_points(), same as
## radial_concentration_index()'s own Monte Carlo mode) instead of the
## depth^4-per-triangle subdivision cloud .mesh_depth_index() always
## builds - decouples cost from mesh size. Named `n_lines`, not `n_points`,
## for the same shared-argument reason radial_concentration_index() is.

#' Monte Carlo (deterministic = FALSE) mean-depth-to-boundary - internal
#' engine for depth_index()'s deterministic = FALSE. See that function's
#' own roxygen for the user-facing parameter reference; parameters here
#' are the same, without defaults.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param n_lines number of random points to sample, or a function(n_tri)
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon()
#' @param seed optional RNG seed
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`; NULL samples proportional to triangle area
#' @param points optional pre-drawn Nx2 coordinate matrix, from
#'   shape_indices()'s point-sharing mechanism
#' @return list(index, mean_depth, ref_depth, area, total_weight, triangles)
#' @noRd
.random_point_depth_index <- function(poly, n_lines, prep, seed, weight = NULL, points = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly)
    tri    <- prep$tri
    poly_u <- prep$poly
    n_tri  <- if (is.null(tri)) 0 else nrow(tri)

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
        return(list(index = NA_real_, mean_depth = NA_real_, ref_depth = NA_real_,
                    area = NA_real_, total_weight = NA_real_, triangles = tri))
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

    bnd <- st_boundary(poly_u)
    d   <- as.numeric(st_distance(.coords_to_points(coords, st_crs(poly_u)), bnd))
    mean_depth <- mean(d)

    ref_depth <- if (is.null(weight)) .disk_reference_depth(area) else .annulus_reference_depth(tri$area, weight)
    index <- mean_depth / ref_depth

    list(index = index, mean_depth = mean_depth, ref_depth = ref_depth, area = area,
         total_weight = if (is.null(raw_total)) area else raw_total, triangles = tri)
}

#' Depth index: mean distance to the boundary vs an equal-area circle
#'
#' depth_bar(rho) / depth_bar_ref, in (0, 1], where depth_bar(rho) is the
#' rho-weighted mean, over the shape's own interior, of each point's
#' distance to the nearest point on the shape's OWN boundary, and
#' depth_bar_ref is the same quantity for the reference shape (a circle,
#' unweighted; a concentric annulus densest at the centre, weighted) -
#' both provable MAXIMISERS of mean depth, so this ratio is the opposite
#' way up from every other index in this package (see
#' [radial_concentration_index()] for the distance-to-a-single-point
#' analogue this is often confused with). See this function's own source
#' comments for the Brunn-Minkowski proof.
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param deterministic if TRUE (default), compute over a fixed
#'   depth-4-subdivision point cloud (same construction
#'   [radial_concentration_index()] uses) - O(n) in triangle count, with a
#'   large constant factor from that subdivision. If FALSE, Monte Carlo
#'   over `n_lines` points sampled directly from the weighted density
#'   instead.
#' @param n_lines number of random points to sample when
#'   deterministic = FALSE (default 3000), or a function(n_tri) - same
#'   shared argument name as this package's other mesh indices, so one
#'   value passed through shape_indices() sets every one of their sample
#'   counts at once.
#' @param seed optional RNG seed, only used when deterministic = FALSE.
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating when also calling other indices
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, substituting for that triangle's own area/mass throughout.
#'   NULL (default) compares against a circle; otherwise the reference
#'   becomes a concentric annulus matching this weight's own histogram,
#'   densest at the centre (see file header for why that's the
#'   depth-maximising arrangement, not the depth-minimising one
#'   [moment_of_inertia_index()]'s/[radial_concentration_index()]'s own
#'   weighted references use the same sort order for).
#' @param points optional pre-drawn Nx2 coordinate matrix, only meaningful
#'   when deterministic = FALSE - primarily for shape_indices()'s internal
#'   point-sharing mechanism.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly.
#' @return list(index, mean_depth, ref_depth, area, total_weight,
#'   triangles). A shape whose true index is 1 (a disk) can show a small
#'   `index > 1` from finite-sample/subdivision noise - an expected
#'   residual, not clamped (consistent with every other estimator in this
#'   package).
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' depth_index(wake)$index
#' # Camden: a long, narrow county - shallow almost everywhere
#' depth_index(nc[nc$NAME == "Camden", ])$index
#'
#' # deterministic = FALSE: Monte Carlo estimate instead of the full
#' # subdivision cloud - a seed makes it reproducible
#' depth_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#'
#' # weight: substitutes for each triangle's own area - weighting by the
#' # triangle's own area exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' depth_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
depth_index <- function(poly, deterministic = TRUE, n_lines = 3000, seed = NULL,
                         prep = NULL, weight = NULL, points = NULL,
                         simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)

    if (!deterministic) {
        return(.random_point_depth_index(poly, n_lines = n_lines, prep = prep, seed = seed,
                                          weight = weight, points = points))
    }
    if (!is.null(points)) {
        stop("`points` (a pre-drawn sample) has no meaning for deterministic = TRUE, which ",
             "computes over the full subdivision point cloud, not a random sample. Drop ",
             "`points`, or set deterministic = FALSE to use it.")
    }

    .mesh_depth_index(prep$tri, prep$poly, weight = weight)
}

#' depth_index() for every row of an sf data frame
#'
#' Each row indexed independently and unweighted.
#' @param x an sf data frame
#' @param ... passed to depth_index() for every row
#' @return `x` with one new column, depth_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- depth_index_sf(nc[1:5, ])
#' res$depth_index
#' @export
depth_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) depth_index(geoms[i], ...))
    x$depth_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
