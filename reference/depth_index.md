# Depth index: mean distance to the boundary vs an equal-area circle

depth_bar(rho) / depth_bar_ref, in (0, 1\], where depth_bar(rho) is the
rho-weighted mean, over the shape's own interior, of each point's
distance to the nearest point on the shape's OWN boundary, and
depth_bar_ref is the same quantity for the reference shape (a circle,
unweighted; a concentric annulus densest at the centre, weighted) - both
provable MAXIMISERS of mean depth, so this ratio is the opposite way up
from every other index in this package (see
[`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)
for the distance-to-a-single-point analogue this is often confused
with). See this function's own source comments for the Brunn-Minkowski
proof.

## Usage

``` r
depth_index(
  poly,
  deterministic = TRUE,
  n_lines = 3000,
  seed = NULL,
  prep = NULL,
  weight = NULL,
  points = NULL,
  simplify_tolerance = NULL
)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

- deterministic:

  if TRUE (default), compute over a fixed depth-4-subdivision point
  cloud (same construction
  [`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)
  uses) - O(n) in triangle count, with a large constant factor from that
  subdivision. If FALSE, Monte Carlo over `n_lines` points sampled
  directly from the weighted density instead.

- n_lines:

  number of random points to sample when deterministic = FALSE (default
  3000), or a function(n_tri) - same shared argument name as this
  package's other mesh indices, so one value passed through
  shape_indices() sets every one of their sample counts at once.

- seed:

  optional RNG seed, only used when deterministic = FALSE.

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating when also calling other indices

- weight:

  optional numeric vector, one entry per triangle in `prep$tri`,
  substituting for that triangle's own area/mass throughout. NULL
  (default) compares against a circle; otherwise the reference becomes a
  concentric annulus matching this weight's own histogram, densest at
  the centre (see file header for why that's the depth-maximising
  arrangement, not the depth-minimising one
  [`moment_of_inertia_index()`](https://nkaza.github.io/shapeindices/reference/moment_of_inertia_index.md)'s/[`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)'s
  own weighted references use the same sort order for).

- points:

  optional pre-drawn Nx2 coordinate matrix, only meaningful when
  deterministic = FALSE - primarily for shape_indices()'s internal
  point-sharing mechanism.

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly.

## Value

list(index, mean_depth, ref_depth, area, total_weight, triangles). A
shape whose true index is 1 (a disk) can show a small `index > 1` from
finite-sample/subdivision noise - an expected residual, not clamped
(consistent with every other estimator in this package).

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
depth_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.7922441
# Camden: a long, narrow county - shallow almost everywhere
depth_index(nc[nc$NAME == "Camden", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4012, lon_0 = -76.234) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.5075611

# deterministic = FALSE: Monte Carlo estimate instead of the full
# subdivision cloud - a seed makes it reproducible
depth_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.7784485

# weight: substitutes for each triangle's own area - weighting by the
# triangle's own area exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
depth_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.7922441
```
