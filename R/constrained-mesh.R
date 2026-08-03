## =========================================================================
## Public entry point onto the constrained-triangulation mesh
## =========================================================================
##
## .constrained_weighted_mesh() (constrained-triangulation.R) has always
## been internal - every caller inside this package (.shape_indices_sf_grouped())
## only ever needed the index values it feeds into, never the mesh itself.
## Downstream packages building their own per-polygon computations on top of
## the SAME mesh (rather than recomputing it, or reimplementing constrained
## Delaunay triangulation independently) need a stable, documented way in -
## this is that entry point, a thin wrapper with no logic of its own beyond
## projecting `x` first (every internal caller already does this once, up
## front, before touching the mesh) and turning `ok = FALSE` into a real
## error instead of a silent partial result.

#' Build the exact constrained-triangulation mesh for a collection of
#' polygons
#'
#' Triangulates the union of `x`'s rows via constrained Delaunay
#' triangulation, using every kept row's own boundary as a constraint
#' segment - so no triangle ever straddles a row boundary, and each
#' triangle's weight/density exactly reproduces its own row's, with no
#' smearing across boundaries the way a coarse triangulate-then-overlay
#' approach would (see `shape_indices_sf(byrow = FALSE)`'s own weighted
#' mesh for that alternative, and this function's internal counterpart,
#' `.constrained_weighted_mesh()`, for the full derivation).
#'
#' This is the same mesh `shape_indices_sf(byrow = FALSE, weights = ...)`
#' builds internally to evaluate its own indices - exposed here for
#' downstream code that wants to build its own computations directly on top
#' of the triangulation (per-triangle areas/weights/owning-polygon), rather
#' than recomputing an equivalent mesh independently or re-deriving index
#' values already available via `shape_indices_sf()`.
#'
#' @param x an sf data frame of (multi)polygons
#' @param weights NULL (each row weighted by its own area - uniform
#'   density), a column-name string, or a numeric vector, one entry per row
#'   of `x` - see `shape_indices_sf()`'s own `weights` doc. A row with
#'   weight exactly 0 or NA is treated as a hole and excluded from the mesh
#'   entirely - it contributes no triangles, and its row index is skipped
#'   throughout `poly_id`.
#' @param simplify_tolerance optional `st_simplify()` tolerance applied to
#'   the union's own outer boundary only - see `.resolve_union()`'s own doc.
#'   Never simplifies row boundaries themselves, which would reintroduce
#'   the smearing this mesh exists to avoid.
#' @return list(P, T, tri_area, tri_weight, poly_id, poly_u, raw_total, crs):
#'   \describe{
#'     \item{P}{N x 2 point matrix, `RTriangle::triangulate()`'s own output}
#'     \item{T}{M x 3 triangle-vertex-index matrix into `P`, one row per
#'       kept triangle (a triangle whose centroid fell inside some kept row;
#'       triangles covering only a hole's interior are dropped)}
#'     \item{tri_area}{length-M numeric, each triangle's own physical area}
#'     \item{tri_weight}{length-M numeric, each triangle's mass (its row's
#'       density times its own area) - already resolved from `weights`, not
#'       yet normalised to sum to 1}
#'     \item{poly_id}{length-M integer, the 1-based index into `x`'s KEPT
#'       rows (rows with a zero/NA weight are dropped before indexing, same
#'       row-keeping `.resolve_union()` always does) that each triangle
#'       belongs to - never straddles a row boundary, by construction}
#'     \item{poly_u}{the (multi)polygon union of every kept row}
#'     \item{raw_total}{`sum(weights)` before normalisation, NA/hole rows
#'       excluded}
#'     \item{crs}{`x`'s coordinate reference system}
#'   }
#' @examples
#' grid9 <- sf::st_sf(
#'     id = 1:9,
#'     geometry = sf::st_make_grid(sf::st_sfc(
#'         sf::st_polygon(
#'             list(rbind(c(0, 0), c(3, 0), c(3, 3), c(0, 3), c(0, 0)))
#'         ),
#'         crs = 3857
#'     ), n = c(3, 3))
#' )
#' grid9$pop <- c(10, 20, 10, 20, 50, 20, 10, 20, 10)
#' mesh <- constrained_mesh(grid9, weights = "pop")
#' # every triangle's weight, summed per owning polygon, reproduces that
#' # polygon's own `pop` value exactly
#' tapply(mesh$tri_weight, mesh$poly_id, sum)
#' @export
constrained_mesh <- function(x, weights = NULL, simplify_tolerance = NULL) {
    x <- .ensure_projected(x)
    cwm <- .constrained_weighted_mesh(x, weights, simplify_tolerance)
    if (!isTRUE(cwm$ok)) {
        stop("Could not build a constrained-triangulation mesh for `x`: ", cwm$reason)
    }
    cwm[c("P", "T", "tri_area", "tri_weight", "poly_id", "poly_u", "raw_total", "crs")]
}
