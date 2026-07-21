#' shapeindices: Shape Indices for (Multi)Polygons
#'
#' Convexity, moment-of-inertia compactness, span (mean-pairwise-distance),
#' and radial concentration (mean-distance-to-geometric-median) indices
#' for (multi)polygons via constrained Delaunay triangulation (CDT), plus
#' classic boundary/hull/bounding-box compactness scores from the
#' redistricting literature (convex-hull ratio, Polsby-Popper,
#' width-length ratio, Reock) that need no triangulation at all. See
#' `vignette("b-understanding-convexity-index", package = "shapeindices")` for the full derivation
#' and worked examples.
#'
#' @keywords internal
#' @import sf
#' @importFrom stats dist runif setNames
#' @importFrom utils combn
"_PACKAGE"
