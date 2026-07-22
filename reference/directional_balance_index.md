# Directional balance index of a (multi)polygon's mass distribution

`1 - R`, where `R` is the mean resultant length of the (weighted)
BEARING of the shape's own interior mass, viewed from its own mass
centroid, in `[0, 1]`. `1` means the mass is directionally balanced - no
bearing pulls harder than the rest; lower values mean the mass leans
toward one direction (a one-sided appendage, an off-centre lobe). One
caveat to know before using it: balanced does not mean round - a
symmetric dumbbell (equal lobes on opposite sides) scores 1 here, same
as a disk, because opposite pulls cancel exactly. See
[`vignette("h-understanding-directional-balance-index")`](https://nkaza.github.io/shapeindices/articles/h-understanding-directional-balance-index.md)
for the derivation, the bound proof, and that blind spot worked through.

## Usage

``` r
directional_balance_index(
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

  if TRUE (default), compute over the same area-adaptive subdivision
  point cloud
  [`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)
  uses. If FALSE, Monte Carlo over `n_lines` points sampled directly
  from the weighted density instead - see the section below on the small
  downward bias this carries for near-balanced shapes, and how to
  compute an unbiased estimate of the squared-scale variant from the
  returned `R`.

- n_lines:

  number of random points to sample when deterministic = FALSE (default
  3000), or a function(n_tri) - same argument name as the package's
  other Monte Carlo mesh indices, so one value passed through
  shape_indices() sets every sample count at once.

- seed:

  optional RNG seed, only used when deterministic = FALSE.

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating when also calling other indices

- weight:

  optional numeric vector, one entry per triangle in `prep$tri`,
  substituting for that triangle's own area/mass throughout. NULL
  (default) reproduces the unweighted (area-based) index exactly.

- points:

  optional pre-drawn Nx2 coordinate matrix, only meaningful when
  deterministic = FALSE - primarily for shape_indices()'s internal
  point-sharing mechanism, not typically supplied directly.

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly, since simplification only happens while
  building `prep` in the first place.

## Value

list(index, R, mean_angle, area, total_weight, centroid, triangles). `R`
is the mean resultant length the index is built from (`index = 1 - R`) -
kept in the output so the unbiased squared-scale estimate described
above can be computed from it. `mean_angle` (radians,
[`atan2()`](https://rdrr.io/r/base/Trig.html) convention) is the bearing
of the net directional pull - which way the mass leans, if it leans at
all; not meaningful when `R` is near 0 (any direction is roughly as good
as any other). `centroid` is the exact mass centroid, always found in
closed form regardless of `deterministic`.

## Monte Carlo bias, and computing an unbiased estimate yourself

With `deterministic = FALSE`, the returned `index` is a biased estimate
of `1 - R`: a finite sample's resultant length is never exactly zero, so
it systematically overestimates a small true `R`, and near-balanced
shapes score slightly LOW - by up to roughly `1/sqrt(n_lines)` for an
exactly balanced shape (about 0.018 at the default `n_lines = 3000`),
shrinking as the true `R` grows or as `n_lines` increases. This cannot
be corrected away on the `1 - R` scale: no unbiased estimator of `1 - R`
exists for any finite sample, for the same structural reason a sample
standard deviation cannot be unbiased even though a sample variance can.

What does admit an exactly unbiased estimate is the closely related
`1 - R^2` (the same quantity on a squared scale: same `[0, 1]` range,
same direction, same shapes scoring exactly 1). Compute it yourself from
the returned `R`: `1 - (n_lines * R^2 - 1) / (n_lines - 1)`. This is
exactly unbiased for `1 - R^2` - any weight distribution, any
`n_lines >= 2` - and is the right quantity to average or compare across
many near-balanced shapes in a batch. It can land slightly above 1 for
near-balanced shapes; that is a necessary feature, not an error (no
estimator confined to `[0, 1]` can be unbiased at the boundary), so do
not truncate it back to 1 if unbiasedness is the point.

`deterministic = TRUE` carries no statistical bias at all - only
deterministic quadrature error, which shrinks with mesh resolution.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
directional_balance_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9965274

# deterministic = FALSE: Monte Carlo estimate - a seed makes it reproducible
directional_balance_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.987398

# weight: substitutes for each triangle's own area - weighting by the
# triangle's own area exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
directional_balance_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.9965274
```
