# Build the exact constrained-triangulation mesh for a collection of polygons

Triangulates the union of `x`'s rows via constrained Delaunay
triangulation, using every kept row's own boundary as a constraint
segment - so no triangle ever straddles a row boundary, and each
triangle's weight/density exactly reproduces its own row's, with no
smearing across boundaries the way a coarse triangulate-then-overlay
approach would (see `shape_indices_sf(byrow = FALSE)`'s own weighted
mesh for that alternative, and this function's internal counterpart,
`.constrained_weighted_mesh()`, for the full derivation).

## Usage

``` r
constrained_mesh(x, weights = NULL, simplify_tolerance = NULL)
```

## Arguments

- x:

  an sf data frame of (multi)polygons

- weights:

  NULL (each row weighted by its own area - uniform density), a
  column-name string, or a numeric vector, one entry per row of `x` -
  see
  [`shape_indices_sf()`](https://nkaza.github.io/shapeindices/reference/shape_indices_sf.md)'s
  own `weights` doc. A row with weight exactly 0 or NA is treated as a
  hole and excluded from the mesh entirely - it contributes no
  triangles, and its row index is skipped throughout `poly_id`.

- simplify_tolerance:

  optional `st_simplify()` tolerance applied to the union's own outer
  boundary only - see `.resolve_union()`'s own doc. Never simplifies row
  boundaries themselves, which would reintroduce the smearing this mesh
  exists to avoid.

## Value

list(P, T, tri_area, tri_weight, poly_id, poly_u, raw_total, crs):

- P:

  N x 2 point matrix,
  [`RTriangle::triangulate()`](https://rdrr.io/pkg/RTriangle/man/triangulate.html)'s
  own output

- T:

  M x 3 triangle-vertex-index matrix into `P`, one row per kept triangle
  (a triangle whose centroid fell inside some kept row; triangles
  covering only a hole's interior are dropped)

- tri_area:

  length-M numeric, each triangle's own physical area

- tri_weight:

  length-M numeric, each triangle's mass (its row's density times its
  own area) - already resolved from `weights`, not yet normalised to sum
  to 1

- poly_id:

  length-M integer, the 1-based index into `x`'s KEPT rows (rows with a
  zero/NA weight are dropped before indexing, same row-keeping
  `.resolve_union()` always does) that each triangle belongs to - never
  straddles a row boundary, by construction

- poly_u:

  the (multi)polygon union of every kept row

- raw_total:

  `sum(weights)` before normalisation, NA/hole rows excluded

- crs:

  `x`'s coordinate reference system

## Details

This is the same mesh `shape_indices_sf(byrow = FALSE, weights = ...)`
builds internally to evaluate its own indices - exposed here for
downstream code that wants to build its own computations directly on top
of the triangulation (per-triangle areas/weights/owning-polygon), rather
than recomputing an equivalent mesh independently or re-deriving index
values already available via
[`shape_indices_sf()`](https://nkaza.github.io/shapeindices/reference/shape_indices_sf.md).

**`weights` convention, worth stating explicitly for anyone using this
package alongside its raster sibling `gridmorph`**: `weights` here is a
raw, EXTENSIVE per-row total (e.g. a population COUNT column like
`"pop"`) - this function divides by each row's own area internally to
get a density, exactly reproducing the row's raw weight value when
triangle weights are summed back up (see the worked example below).
`gridmorph`'s own raster indices (`weighted = TRUE`) instead require an
already-computed per-cell DENSITY value - a raster cell has no natural
"this cell's own share of a polygon total" the way a row does, so
passing a raw count raster there produces a different, resolution-
dependent result. If you're translating a
`weights = "some_count_column"` workflow from this package (or its
regionalization sibling `reseam`) into a raster one (`gridmorph`, or its
own regionalization sibling `restitch`), divide by area first.

## Examples

``` r
grid9 <- sf::st_sf(
    id = 1:9,
    geometry = sf::st_make_grid(sf::st_sfc(
        sf::st_polygon(
            list(rbind(c(0, 0), c(3, 0), c(3, 3), c(0, 3), c(0, 0)))
        ),
        crs = 3857
    ), n = c(3, 3))
)
grid9$pop <- c(10, 20, 10, 20, 50, 20, 10, 20, 10)
mesh <- constrained_mesh(grid9, weights = "pop")
# every triangle's weight, summed per owning polygon, reproduces that
# polygon's own `pop` value exactly
tapply(mesh$tri_weight, mesh$poly_id, sum)
#>  1  2  3  4  5  6  7  8  9 
#> 10 20 10 20 50 20 10 20 10 
```
