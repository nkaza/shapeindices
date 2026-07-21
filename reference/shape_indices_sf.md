# All indices (or a chosen subset) for every feature in an sf data frame.

All indices (or a chosen subset) for every feature in an sf data frame.

## Usage

``` r
shape_indices_sf(
  x,
  byrow = TRUE,
  which = "all",
  weights = NULL,
  id = NULL,
  parallel_rows = TRUE,
  ...
)
```

## Arguments

- x:

  an sf data frame

- byrow:

  if TRUE (default), index each row separately (`weights` and `id` must
  be NULL). If FALSE, treat every row as a weighted sub-polygon of one
  overall shape (`st_union(x)`) and return a single-row sf.

- which:

  "all" (default), or a character vector naming a subset of these
  thirteen values - each listed here with the function it actually
  calls, since the `which` string and the function name aren't always
  identical:

  - `"convexity"` -
    [`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.md)

  - `"moment_of_inertia"` -
    [`moment_of_inertia_index()`](https://nkaza.github.io/shapeindices/reference/moment_of_inertia_index.md)

  - `"moment_isotropy"` -
    [`moment_isotropy_index()`](https://nkaza.github.io/shapeindices/reference/moment_isotropy_index.md)

  - `"directional_balance"` -
    [`directional_balance_index()`](https://nkaza.github.io/shapeindices/reference/directional_balance_index.md)

  - `"span"` -
    [`span_index()`](https://nkaza.github.io/shapeindices/reference/span_index.md)

  - `"radial_concentration"` -
    [`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)

  - `"depth"` -
    [`depth_index()`](https://nkaza.github.io/shapeindices/reference/depth_index.md)

  - `"hull_ratio"` -
    [`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)

  - `"polsby_popper"` -
    [`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)

  - `"width_length_ratio"` -
    [`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)

  - `"reock"` -
    [`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)

  - `"detour"` -
    [`detour_index()`](https://nkaza.github.io/shapeindices/reference/detour_index.md)

  - `"exchange"` -
    [`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)

  The first seven need a CDT mesh; the last six are classic
  boundary/hull/bounding-box metrics that don't - requesting only those
  six skips triangulation entirely, in both `byrow` modes. An
  unrecognised name in `which` errors immediately, listing all thirteen
  valid values.

- weights:

  only used when `byrow = FALSE`. NULL (default) weights each row by its
  own area; otherwise a numeric vector or column name. Weight 0/NA
  excludes that row as a HOLE (with a warning) rather than
  zero-weighting it - since this changes `poly_u` itself, even
  `hull_ratio_index` can differ between calls that only differ in
  `weights`. When rows genuinely differ in weight, the combined mesh is
  built by an exact constrained triangulation (every kept row's own
  boundary supplied as a triangulation constraint, so no triangle can
  ever straddle a row boundary and density allocation is exact) rather
  than triangulating the union once and area-weight-averaging row
  density onto whatever it overlaps. Falls back to that coarser overlay
  automatically, with a warning, if kept rows genuinely overlap or
  triangulation otherwise fails; `weights = NULL` never needs either
  approach, since every row then weights by its own area and density is
  uniform everywhere.

- id:

  required when `byrow = FALSE`: a single scalar value identifying this
  group in the collapsed single-row result.

- parallel_rows:

  only meaningful when `byrow = TRUE` (errors if TRUE with
  `byrow = FALSE`). If TRUE (default) and furrr is installed, rows run
  across the active
  [`future::plan()`](https://future.futureverse.org/reference/plan.html) -
  pass `deterministic_max_tri` too on real-world polygons first, or
  unbounded per-row cost plus parallelism turns "slow" into "OOM crash".

- ...:

  passed to shape_indices() when `byrow = TRUE`. When `byrow = FALSE`:
  `n_lines`, `seed` are read the same way and set every mesh-based Monte
  Carlo index's sample count together; `n_quad` likewise, but only
  affects convexity/span (radial_concentration's, directional_balance's,
  and depth's deterministic modes have no quadrature refinement of that
  kind to control); `deterministic` picks the method directly (default
  TRUE), shared by convexity, span, radial_concentration,
  directional_balance, and depth alike; `deterministic_max_tri` can push
  `deterministic` from TRUE to FALSE based on the combined triangle
  count, but never overrides an explicit `deterministic = FALSE`; on the
  exact constrained mesh above, convexity's/span's deterministic mode is
  additionally capped by a safety ceiling derived from memory actually
  available on this machine right now (that mesh can be far larger than
  the union-only one, since it's built at the rows' own boundary
  resolution) - left at its default, a mesh above that ceiling falls
  back to the Monte Carlo estimator silently; an explicit
  `deterministic_max_tri` that would still exceed the safe ceiling
  errors instead, rather than either attempting the computation or
  silently overriding what was asked for; `simplify_tolerance`
  simplifies `st_union(x)` itself before triangulating - do not simplify
  the rows yourself beforehand instead: adjacent rows sharing an edge
  get that edge simplified differently on each side, fragmenting the
  union into slivers (see
  [`vignette("a-basic-usage")`](https://nkaza.github.io/shapeindices/articles/a-basic-usage.md)'s
  simplify_tolerance section).

## Value

if `byrow = TRUE`: x with one new `<name>_index` column per requested
index. If `byrow = FALSE`: a single-row sf with id, those same columns,
and total_weight.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)

# byrow = TRUE (default): index each county independently
res <- shape_indices_sf(nc[1:5, ])
#> Warning: parallel_rows = TRUE but no parallel future::plan() is active (still on the default sequential plan) - running in order. Call future::plan(future::multisession, workers = ...) first for actual multi-core speedup.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
res$convexity_index
#> [1] 0.9971652 0.9862805 0.9986526 0.8047873 0.9766696

# a subset, requested the same way as shape_indices()
shape_indices_sf(nc[1:5, ], which = c("hull_ratio", "reock"))
#> Warning: parallel_rows = TRUE but no parallel future::plan() is active (still on the default sequential plan) - running in order. Call future::plan(future::multisession, workers = ...) first for actual multi-core speedup.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
#> Simple feature collection with 5 features and 16 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -81.74107 ymin: 36.07282 xmax: -75.77316 ymax: 36.58965
#> Geodetic CRS:  NAD27
#>    AREA PERIMETER CNTY_ CNTY_ID        NAME  FIPS FIPSNO CRESS_ID BIR74 SID74
#> 1 0.114     1.442  1825    1825        Ashe 37009  37009        5  1091     1
#> 2 0.061     1.231  1827    1827   Alleghany 37005  37005        3   487     0
#> 3 0.143     1.630  1828    1828       Surry 37171  37171       86  3188     5
#> 4 0.070     2.968  1831    1831   Currituck 37053  37053       27   508     1
#> 5 0.153     2.206  1832    1832 Northampton 37131  37131       66  1421     9
#>   NWBIR74 BIR79 SID79 NWBIR79                       geometry hull_ratio_index
#> 1      10  1364     0      19 MULTIPOLYGON (((-81.47276 3...        0.9154096
#> 2      10   542     3      12 MULTIPOLYGON (((-81.23989 3...        0.8256371
#> 3     208  3616     6     260 MULTIPOLYGON (((-80.45634 3...        0.9176972
#> 4     123   830     2     145 MULTIPOLYGON (((-76.00897 3...        0.5467063
#> 5    1066  1606     3    1197 MULTIPOLYGON (((-77.21767 3...        0.8088957
#>   reock_index
#> 1   0.6553073
#> 2   0.4581459
#> 3   0.6269183
#> 4   0.1716233
#> 5   0.3414909

# byrow = TRUE, with deterministic_max_tri forcing the Monte Carlo
# estimator once a row's own mesh exceeds it - `...` passes n_lines/seed
# through to shape_indices() for every row
res_rli <- shape_indices_sf(nc[1:5, ], byrow = TRUE, deterministic_max_tri = 5,
                             n_lines = 2000, seed = 1)
#> Warning: parallel_rows = TRUE but no parallel future::plan() is active (still on the default sequential plan) - running in order. Call future::plan(future::multisession, workers = ...) first for actual multi-core speedup.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 253 triangle-pairs that deterministic = TRUE (23 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: n_lines (2000) is not substantially lower than the 253 triangle-pairs that deterministic = TRUE (23 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 300 triangle-pairs that deterministic = TRUE (25 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: n_lines (2000) is not substantially lower than the 300 triangle-pairs that deterministic = TRUE (25 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 406 triangle-pairs that deterministic = TRUE (29 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: n_lines (2000) is not substantially lower than the 406 triangle-pairs that deterministic = TRUE (29 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 465 triangle-pairs that deterministic = TRUE (31 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: n_lines (2000) is not substantially lower than the 465 triangle-pairs that deterministic = TRUE (31 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: UNRELIABLE VALUE: Future (<unnamed-3>) unexpectedly generated random numbers without specifying argument 'seed'. There is a risk that those random numbers are not statistically sound and the overall results might be invalid. To fix this, specify 'seed=TRUE'. This ensures that proper, parallel-safe random numbers are produced. To disable this check, use 'seed=NULL', or set option 'future.rng.onMisuse' to "ignore". [future <unnamed-3> (ef97d64e929852c0691af3b2fbdef7bd-3); on ef97d64e929852c0691af3b2fbdef7bd@runnervm3jd5f<7785>]
res_rli$convexity_index
#> [1] 0.9949545 0.9844751 0.9984327 0.8028660 0.9753344

# byrow = FALSE: treat several adjacent counties as one weighted shape -
# Wake, Durham, Orange and Chatham (the Research Triangle) are contiguous
triangle <- nc[nc$NAME %in% c("Wake", "Durham", "Orange", "Chatham"), ]
shape_indices_sf(triangle, byrow = FALSE, id = "triangle_by_area")
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.8377, lon_0 = -78.9545) before computing - pass already-projected data instead if you need a specific CRS.
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -54482.01 ymin: -36248.79 xmax: 63283.74 ymax: 44215.23
#> Projected CRS: +proj=laea +lat_0=35.8376681673 +lon_0=-78.9545321393 +datum=WGS84 +units=m +no_defs
#>                 id convexity_index moment_of_inertia_index
#> 1 triangle_by_area        0.983053               0.8154609
#>   moment_isotropy_index directional_balance_index span_index
#> 1             0.5073647                 0.9836302  0.9167346
#>   radial_concentration_index depth_index hull_ratio_index polsby_popper_index
#> 1                  0.9243383   0.6853179        0.8306064           0.5140176
#>   width_length_ratio_index reock_index detour_index exchange_index total_weight
#> 1                0.6856106   0.4828814    0.8297133      0.8073051   5812126377
#>                         geometry
#> 1 POLYGON ((-25588.04 -0.0077...
shape_indices_sf(triangle, byrow = FALSE, weights = "BIR74",
                  id = "triangle_by_births")
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.8377, lon_0 = -78.9545) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: hull_ratio/polsby_popper/width_length_ratio/reock/detour/exchange don't use `weights` in their own formula at all - only whether a row's weight is exactly 0/NA (excluded as a hole, which changes poly_u itself) has any effect. The magnitude of a nonzero weight (5 vs 5000) is completely invisible to any of the six classic metrics; if you meant them to reflect the weighting itself, they won't.
#> Warning: Large mesh: 69 pieces, 2346 pairs, 21114 candidate lines to evaluate at n_quad = 3 (deterministic = TRUE is O(n_quad^2 * n^2)) - this can be slow. If it is, try: n_quad = 1 (drops the 9x quadrature multiplier), deterministic = FALSE (random-line estimate instead of exhaustive), or shape_indices()/shape_indices_sf()'s deterministic_max_tri to switch to deterministic = FALSE automatically above a size threshold.
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -54482.01 ymin: -36248.79 xmax: 63283.74 ymax: 44215.23
#> Projected CRS: +proj=laea +lat_0=35.8376681673 +lon_0=-78.9545321393 +datum=WGS84 +units=m +no_defs
#>                   id convexity_index moment_of_inertia_index
#> 1 triangle_by_births       0.9895831               0.6970804
#>   moment_isotropy_index directional_balance_index span_index
#> 1             0.6430302                 0.9598688  0.8492941
#>   radial_concentration_index depth_index hull_ratio_index polsby_popper_index
#> 1                  0.8519084   0.5074114        0.8306064           0.5140176
#>   width_length_ratio_index reock_index detour_index exchange_index total_weight
#> 1                0.6856106   0.4828814    0.8297133      0.8073051        27264
#>                         geometry
#> 1 POLYGON ((-25588.04 -0.0077...

# weight 0/NA treats a row as a hole - excluded from the shape itself,
# not just zero-weighted - e.g. dropping Chatham (the least urban/most
# rural of the four) out of the footprint entirely
w <- ifelse(triangle$NAME == "Chatham", 0, triangle$BIR74)
shape_indices_sf(triangle, byrow = FALSE, weights = w, id = "triangle_minus_chatham")
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.8377, lon_0 = -78.9545) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: 1 of 4 rows have zero or NA weight; treating them as holes - excluded from the triangulated union entirely (not just zero-weighted area). Pass weights without zeros/NA for those rows if you want that area included with zero weight instead.
#> Warning: hull_ratio/polsby_popper/width_length_ratio/reock/detour/exchange don't use `weights` in their own formula at all - only whether a row's weight is exactly 0/NA (excluded as a hole, which changes poly_u itself) has any effect. The magnitude of a nonzero weight (5 vs 5000) is completely invisible to any of the six classic metrics; if you meant them to reflect the weighting itself, they won't.
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -28530.05 ymin: -35213.33 xmax: 63283.74 ymax: 44215.23
#> Projected CRS: +proj=laea +lat_0=35.8376681673 +lon_0=-78.9545321393 +datum=WGS84 +units=m +no_defs
#>                       id convexity_index moment_of_inertia_index
#> 1 triangle_minus_chatham       0.9771155               0.7165332
#>   moment_isotropy_index directional_balance_index span_index
#> 1             0.4360721                 0.9882301  0.8597266
#>   radial_concentration_index depth_index hull_ratio_index polsby_popper_index
#> 1                  0.8571135   0.5421913        0.8359174           0.5477416
#>   width_length_ratio_index reock_index detour_index exchange_index total_weight
#> 1                  0.73927   0.4845069    0.8362334      0.7674118        25618
#>                         geometry
#> 1 POLYGON ((-3989.542 -26212....

# excluding an INTERIOR row (rather than one on the edge, as above) can
# punch an actual hole through the union, not just trim the boundary. A
# 3x3 grid of unit squares is a solid convex block when every cell is
# kept; excluding only the centre cell turns it into a square ring (CI
# 1 -> 0.86, hull_ratio_index 1 -> 0.89 = 8/9) - two calls differing only in
# `weights` can disagree on hull_ratio_index despite hull_ratio_index having no
# weighted form of its own.
cell <- function(cx, cy) sf::st_polygon(list(rbind(
  c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
  c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
ctr <- expand.grid(x = -1:1, y = -1:1)
grid9 <- sf::st_sf(pop = ifelse(ctr$x == 0 & ctr$y == 0, 0, 10),
                    geometry = sf::st_sfc(mapply(cell, ctr$x, ctr$y, SIMPLIFY = FALSE),
                                           crs = 3857))
shape_indices_sf(grid9, byrow = FALSE, id = "solid_block")            # weights = NULL: keeps all 9
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -1.5 ymin: -1.5 xmax: 1.5 ymax: 1.5
#> Projected CRS: WGS 84 / Pseudo-Mercator
#>            id convexity_index moment_of_inertia_index moment_isotropy_index
#> 1 solid_block               1               0.9549297                     1
#>   directional_balance_index span_index radial_concentration_index depth_index
#> 1                         1  0.9798324                  0.9834955   0.8873007
#>   hull_ratio_index polsby_popper_index width_length_ratio_index reock_index
#> 1                1           0.7853982                        1   0.6366198
#>   detour_index exchange_index total_weight                       geometry
#> 1    0.8862269      0.9094344            9 POLYGON ((1.5 -1.5, 0.5 -1....
shape_indices_sf(grid9, byrow = FALSE, weights = "pop", id = "ring")  # centre cell excluded
#> Warning: 1 of 9 rows have zero or NA weight; treating them as holes - excluded from the triangulated union entirely (not just zero-weighted area). Pass weights without zeros/NA for those rows if you want that area included with zero weight instead.
#> Warning: hull_ratio/polsby_popper/width_length_ratio/reock/detour/exchange don't use `weights` in their own formula at all - only whether a row's weight is exactly 0/NA (excluded as a hole, which changes poly_u itself) has any effect. The magnitude of a nonzero weight (5 vs 5000) is completely invisible to any of the six classic metrics; if you meant them to reflect the weighting itself, they won't.
#> Simple feature collection with 1 feature and 15 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: -1.5 ymin: -1.5 xmax: 1.5 ymax: 1.5
#> Projected CRS: WGS 84 / Pseudo-Mercator
#>     id convexity_index moment_of_inertia_index moment_isotropy_index
#> 1 ring       0.8624193               0.7639437                     1
#>   directional_balance_index span_index radial_concentration_index depth_index
#> 1                         1  0.8750462                  0.8556306   0.4799515
#>   hull_ratio_index polsby_popper_index width_length_ratio_index reock_index
#> 1        0.8888889           0.3926991                        1   0.5658842
#>   detour_index exchange_index total_weight                       geometry
#> 1    0.8355428      0.8399935           80 POLYGON ((1.5 -1.5, 0.5 -1....
```
