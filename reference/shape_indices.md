# All indices (or a chosen subset) for a single (multi)polygon

Shares one triangulation when any CDT-based index is requested.

## Usage

``` r
shape_indices(
  poly,
  which = "all",
  deterministic_max_tri = NULL,
  simplify_tolerance = NULL,
  ...
)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

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
  six skips triangulation entirely. An unrecognised name in `which`
  errors immediately, listing all thirteen valid values.

- deterministic_max_tri:

  if set, meshes with more triangles than this use
  `deterministic = FALSE` automatically - real-world polygons can
  triangulate into thousands of CDT triangles, where the
  exhaustive/full-subdivision methods are minutes-to-hours per polygon.
  Can only push `deterministic` from TRUE/default to FALSE, never
  override an explicit `deterministic = FALSE` back to TRUE. Shared by
  convexity_index()/span_index()/radial_concentration_index()/
  directional_balance_index()/depth_index() alike; has no effect on the
  six classic metrics, which have no mesh to switch at all.

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md) -
  see its own doc. Has no effect when `which` requests only the six
  classic metrics, since no mesh is built for those to begin with.

- ...:

  passed to convexity_index(), span_index(),
  radial_concentration_index(), directional_balance_index(), and
  depth_index() (e.g. `weight` is shared by all five, and also forwarded
  to moment_of_inertia_index()/moment_isotropy_index(), so every
  weighted index stays consistently weighted; `deterministic`/`n_lines`
  are shared by all five mesh-based Monte Carlo indices - `n_lines` sets
  every one of their sample counts with one argument, even though
  span_index()/radial_concentration_index()/directional_balance_index()/
  depth_index() never actually build line geometry; `n_quad` is
  convexity_index()/span_index() only, since the other three mesh
  indices' deterministic modes have no quadrature refinement of that
  kind to control).

## Value

named numeric vector, one entry per requested index, in canonical order

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
shape_indices(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#>            convexity    moment_of_inertia      moment_isotropy 
#>            0.9948142            0.8556213            0.4230179 
#>  directional_balance                 span radial_concentration 
#>            0.9965274            0.9385481            0.9407923 
#>                depth           hull_ratio        polsby_popper 
#>            0.7922441            0.9046430            0.6285667 
#>   width_length_ratio                reock               detour 
#>            0.9595009            0.5379165            0.8643451 
#>             exchange 
#>            0.8233315 

# a subset - skips triangulation entirely, since none of these need a mesh
shape_indices(wake, which = c("hull_ratio", "polsby_popper", "reock"))
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#>    hull_ratio polsby_popper         reock 
#>     0.9046430     0.6285667     0.5379165 

# deterministic_max_tri forces the cheaper Monte Carlo estimator once
# the mesh exceeds it, regardless of how simple the polygon actually is -
# here set deliberately low (5) just to demonstrate the switch on a
# small county
shape_indices(wake, deterministic_max_tri = 5, n_lines = 2000, seed = 1)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#>            convexity    moment_of_inertia      moment_isotropy 
#>            0.9923343            0.8556213            0.4230179 
#>  directional_balance                 span radial_concentration 
#>            0.9909765            0.9447978            0.9331626 
#>                depth           hull_ratio        polsby_popper 
#>            0.7898759            0.9046430            0.6285667 
#>   width_length_ratio                reock               detour 
#>            0.9595009            0.5379165            0.8643451 
#>             exchange 
#>            0.8233315 

# simplify_tolerance trades a small amount of boundary detail (here 50m)
# for a smaller mesh - passed straight through to prepare_polygon()
shape_indices(wake, simplify_tolerance = 50)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#>            convexity    moment_of_inertia      moment_isotropy 
#>            0.9953327            0.8556330            0.4230274 
#>  directional_balance                 span radial_concentration 
#>            0.9965106            0.9385326            0.9408442 
#>                depth           hull_ratio        polsby_popper 
#>            0.7922886            0.9047569            0.6286460 
#>   width_length_ratio                reock               detour 
#>            0.9595009            0.5379842            0.8643995 
#>             exchange 
#>            0.8233065 
```
