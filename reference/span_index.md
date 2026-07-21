# Span index: mean pairwise interior distance vs. an equal-area circle

D_ref/D, in (0, 1\], where D is the mean distance between two random
interior points and D_ref is the same quantity for the reference shape
(a circle, unweighted; a concentric annulus, weighted) - both are
provable minimisers of D. Distinct from moment_of_inertia_index(), which
a squared-distance version of this index would just collapse to. See
[`vignette("d-understanding-span-index")`](https://nkaza.github.io/shapeindices/articles/d-understanding-span-index.md)
for the derivation and both proofs.

## Usage

``` r
span_index(
  poly,
  deterministic = TRUE,
  n_quad = 3,
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

  if TRUE (default), compute over the polygon's own fixed CDT quadrature
  mesh - O(n^2). If FALSE, Monte Carlo over `n_lines` random interior
  point pairs instead - O(n_lines), useful when the mesh is too large.
  `n_quad` is unused (and an error if passed explicitly) when
  deterministic = FALSE.

- n_quad:

  quadrature points per triangle when deterministic = TRUE (1 or 3),
  used for the cross-triangle term only (the same-triangle term is
  always computed separately, by subdividing each triangle). 3 (default)
  is more accurate; drop to 1 on meshes of more than a few hundred
  triangles if that matters.

- n_lines:

  number of random point pairs to sample when deterministic = FALSE
  (default 3000), or a function(n_tri) - same argument name as
  convexity_index()'s Monte Carlo mode, so one value passed through
  shape_indices() sets every mesh index's sample count at once, though
  no line geometry is actually built here. Warns if not substantially
  lower than what deterministic = TRUE would need for the same polygon.

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

  optional pre-drawn (2\*n_lines)x2 coordinate matrix, only meaningful
  when deterministic = FALSE - primarily for shape_indices()'s internal
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

list(index, D, D_ref, area, total_weight, triangles).

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
span_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9385481
# Camden: a long, narrow county - far from a circle
span_index(nc[nc$NAME == "Camden", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4012, lon_0 = -76.234) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.6125629

# deterministic = FALSE: Monte Carlo estimate instead of the exhaustive mesh
span_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> [1] 0.9447978

# weight: substitutes for each triangle's own area - weighting by the
# triangle's own area exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
span_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.9385479
```
