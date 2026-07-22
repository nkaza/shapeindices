# Radial concentration index: mean distance to the geometric median vs an equal-area circle

D1_ref/D1, in (0, 1\], where D1 is the mean distance from random
interior points to the shape's own geometric median (the point
minimising that mean distance - not the centroid) and D1_ref is the same
quantity for the reference shape (a circle, unweighted; a concentric
annulus, weighted) - both provable minimisers of D1. See
[`vignette("f-understanding-radial-concentration-index")`](https://nkaza.github.io/shapeindices/articles/f-understanding-radial-concentration-index.md)
for the derivation and the proofs.

## Usage

``` r
radial_concentration_index(
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
  cloud (4^4 = 256 points per CDT triangle) - O(n) in triangle count,
  but with a large constant factor from that 256x expansion. If FALSE,
  Monte Carlo over `n_lines` points sampled directly from the weighted
  density instead - decouples cost from mesh size, useful when the mesh
  is large enough that the 256x subdivision cloud itself becomes the
  bottleneck.

- n_lines:

  number of random points to sample when deterministic = FALSE (default
  3000), or a function(n_tri) - same argument name as
  convexity_index()'s/span_index()'s Monte Carlo mode, so one value
  passed through shape_indices() sets every mesh index's sample count at
  once, even though this draws single points, not lines or point-pairs.

- seed:

  optional RNG seed, only used when deterministic = FALSE.

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating when also calling other indices

- weight:

  optional numeric vector, one entry per triangle in `prep$tri`,
  substituting for that triangle's own area/mass throughout. NULL
  (default) compares against a circle; otherwise the reference becomes a
  concentric annulus matching this weight's own histogram.

- points:

  optional pre-drawn Nx2 coordinate matrix, only meaningful when
  deterministic = FALSE - primarily for shape_indices()'s internal
  point-sharing mechanism (draw one sample and reuse it across
  convexity/span/radial_concentration instead of each drawing
  independently), not typically supplied directly.

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly, since simplification only happens while
  building `prep` in the first place.

## Value

list(index, D1, D1_ref, area, total_weight, center, triangles). `center`
is the geometric median found by Weiszfeld's algorithm. For symmetric
shapes the geometric median can be non-unique (the minimising "point" is
a whole segment - e.g. between two identical separated blobs), so
`center` may land anywhere along it, including inside a hole or the gap
between multi-part pieces; the index value itself is unaffected.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
radial_concentration_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9407923
# Camden: a long, narrow county - far from a circle
radial_concentration_index(nc[nc$NAME == "Camden", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4012, lon_0 = -76.234) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.606436

# deterministic = FALSE: Monte Carlo estimate instead of the full
# subdivision cloud - a seed makes it reproducible
radial_concentration_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9294861

# weight: substitutes for each triangle's own area - weighting by the
# triangle's own area exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
radial_concentration_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.9407923
```
