## =========================================================================
## 9. Moment isotropy index (shares moment_of_inertia_index()'s
## triangulation and mass-moment machinery)
## =========================================================================
##
## index = lambda_min / lambda_max, the ratio of the two PRINCIPAL moments
## of the mass inertia tensor M = [[Ixx, Ixy], [Ixy, Iyy]] (about the mass
## centroid, same Ixx/Iyy/Ixy moment_of_inertia_index() already computes).
## M is provably positive semi-definite for any 2D mass distribution -
## it's literally an integral of rank-1 PSD matrices [y,x]^T[y,x]
## (mass-weighted) over the shape - so lambda_min >= 0 and lambda_max > 0
## for any polygon with positive area, giving index in (0, 1] with no
## extra derivation: lambda_min = 0 only in the degenerate limit where all
## mass lies on a single line, never a genuine 2D polygon.
##
## NAMING: an earlier version of this called itself "eccentricity_index"
## - a real mistake, not just a style choice. Classical eccentricity of an
## ellipse, e = sqrt(1 - (b/a)^2), INCREASES from 0 (circle) to 1
## (degenerate) as a shape elongates; this index runs the opposite way -
## 1 (isotropic) down to 0 (elongated) - and isn't even the same formula
## (this is b^2/a^2 = 1 - e^2 for an ellipse, not e itself). "Isotropy" is
## the correct word precisely because it already points the right
## direction: high isotropy genuinely does mean round/uniform, matching
## this index's own 1-is-good convention, the same direction every other
## index in this package uses.
##
## index = 1 iff the two principal moments are equal (the mass
## distribution is rotationally isotropic about its own centroid) - a
## disk achieves this, but so does a square, a regular hexagon, or any
## shape with 3-fold or higher rotational symmetry. Unlike
## convexity_index()/moment_of_inertia_index()/span_index()/
## radial_concentration_index(), this is NOT "= 1 iff a disk" - it
## measures pure anisotropy/elongation of the mass distribution, a
## genuinely different question from "how disk-like is this shape". A
## long thin rectangle is perfectly convex but scores low here; an
## isotropic but very non-convex pinwheel shape could score high here
## while scoring low on convexity_index().
##
## NO EXTERNAL REFERENCE SHAPE: the other four compare the actual shape's
## own integral against a provably-optimal disk or concentric-rings
## arrangement (a rearrangement-inequality/bathtub-principle argument).
## This index has no such reference - it compares the shape's own two
## principal moments to each other directly, an elementary PSD-matrix
## fact, not a rearrangement inequality. Weighting is still "free" in the
## same practical sense: Ixx/Iyy/Ixy already come from
## moment_of_inertia_index()'s own rho-weighted core, so a weighted
## isotropy is just the same eigenvalue ratio computed from the weighted
## tensor - there's no separate concentric-rings derivation needed,
## because there's no reference shape to derive one for.

#' Moment isotropy index of a (multi)polygon's mass distribution
#'
#' Ratio of the smaller to larger principal moment of the mass inertia
#' tensor, in `(0, 1]`. `1` means the mass distribution is rotationally
#' isotropic about its own centroid (a disk qualifies, but so does any
#' shape with 3-fold or higher rotational symmetry - this is not the same
#' claim as convexity_index()'s or moment_of_inertia_index()'s "= 1 iff a
#' disk"). Lower values mean the mass is more anisotropic - elongated
#' along one direction - regardless of whether the shape itself is convex
#' or dispersed. See `vignette("g-understanding-moment-isotropy-index")`
#' for the derivation, the bound proof, and how this relates to (but
#' differs from) classical eccentricity.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, used as that triangle's mass instead of its own area.
#'   NULL (default) reproduces the unweighted (area-based) index exactly.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return list(index, lambda_min, lambda_max, Ixx, Iyy, Ixy, area,
#'   total_weight, centroid, triangles). `centroid`/`Ixx`/`Iyy`/`Ixy` are
#'   exactly moment_of_inertia_index()'s own fields, computed about the
#'   same mass centroid.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' moment_isotropy_index(wake)$index
#' # Camden: long and narrow - low isotropy despite being convex
#' moment_isotropy_index(nc[nc$NAME == "Camden", ])$index
#'
#' # weight: substitutes for each triangle's own mass, same as
#' # moment_of_inertia_index() - weighting by the triangle's own area
#' # exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' moment_isotropy_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
moment_isotropy_index <- function(poly, prep = NULL, weight = NULL, simplify_tolerance = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)
    poly <- prep$poly
    tri  <- prep$tri
    n <- if (is.null(tri)) 0 else nrow(tri)
    if (n == 0) {
        warning("Triangulation produced no triangles; index is not defined.")
        return(list(index = NA_real_, lambda_min = NA_real_, lambda_max = NA_real_,
                    Ixx = NA_real_, Iyy = NA_real_, Ixy = NA_real_,
                    area = NA_real_, total_weight = NA_real_,
                    centroid = NULL, triangles = tri))
    }
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }

    tri_area  <- tri$area
    raw_total <- if (is.null(weight)) NULL else sum(weight)
    w         <- if (is.null(weight)) tri_area else .normalize_weight(weight)
    rho       <- w / tri_area

    core <- .moment_of_inertia_core(st_geometry(tri), tri_area, rho, poly)
    eig  <- .principal_moments(core$Ixx, core$Iyy, core$Ixy)

    list(index = eig$lambda_min / eig$lambda_max,
         lambda_min = eig$lambda_min, lambda_max = eig$lambda_max,
         Ixx = core$Ixx, Iyy = core$Iyy, Ixy = core$Ixy,
         area = core$area,
         total_weight = if (is.null(raw_total)) core$total_weight else raw_total,
         centroid = core$centroid, triangles = tri)
}

#' Principal moments (eigenvalues) of the 2x2 mass inertia tensor
#' `[[Ixx, Ixy], [Ixy, Iyy]]` - closed form via the standard 2x2 symmetric
#' eigenvalue formula, same result as eigen() but without the function-call
#' overhead for something this elementary. Mathematically always
#' non-negative (the tensor is PSD - see file header); `max(.., 0)` guards
#' only against machine-epsilon-level floating point noise around an exact
#' zero (e.g. a perfectly isotropic tensor computed via subtraction), not
#' against any real algorithmic bound violation.
#' @param Ixx,Iyy,Ixy scalars, from .moment_of_inertia_core()
#' @return list(lambda_min, lambda_max)
#' @noRd
.principal_moments <- function(Ixx, Iyy, Ixy) {
    mean_i    <- (Ixx + Iyy) / 2
    diff_term <- sqrt(((Ixx - Iyy) / 2)^2 + Ixy^2)
    list(lambda_min = max(mean_i - diff_term, 0), lambda_max = mean_i + diff_term)
}

#' Vectorised wrapper for an sf data frame - each row indexed independently.
#' @param x an sf data frame
#' @param ... passed to moment_isotropy_index() for every row
#' @return `x` with one new column, moment_isotropy_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- moment_isotropy_index_sf(nc[1:5, ])
#' res$moment_isotropy_index
#' @export
moment_isotropy_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) moment_isotropy_index(geoms[i], ...))
    x$moment_isotropy_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
