# Convexity/dispersal index of a (multi)polygon

1 minus the expected fraction of a random interior line lying outside
the polygon. 1 = convex; lower means more concave and/or more spatially
dispersed.

## Usage

``` r
convexity_index(
  poly,
  deterministic = TRUE,
  n_quad = 3,
  n_lines = 3000,
  seed = NULL,
  plot = FALSE,
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

  if TRUE (default), compute over a fixed quadrature grid on all pairs
  of the polygon's own CDT triangles - O(n^2). If FALSE, use a Monte
  Carlo estimate over `n_lines` random interior lines instead -
  O(n_lines), useful when the mesh is too large, and also the more
  accurate option for shapes with small concavities relative to triangle
  size, since a fixed quadrature grid can miss a concavity that falls
  between its sample points (see
  [`vignette("b-understanding-convexity-index")`](https://nkaza.github.io/shapeindices/articles/b-understanding-convexity-index.md)).
  `n_quad` is unused (and an error if passed explicitly) when
  deterministic = FALSE.

- n_quad:

  quadrature points per triangle when deterministic = TRUE (1 or 3). 3
  (default) uses a 3-point Hammer-Stroud rule, closing most of the gap
  to the random-line index at ~9x the cost of n_quad = 1; drop to 1 on
  meshes of more than a few hundred triangles if that matters.

- n_lines:

  number of random lines to sample when deterministic = FALSE (default
  3000), or a function(n_tri) for callers whose polygons vary widely in
  complexity. Warns if not substantially lower than what deterministic =
  TRUE would need for the same polygon.

- seed:

  optional RNG seed, only used when deterministic = FALSE.

- plot:

  draw a diagnostic plot (needs a graphics device)

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating when also calling moment_of_inertia_index()

- weight:

  optional numeric vector, one entry per triangle in `prep$tri`,
  substituting for that triangle's own area/mass throughout. NULL
  (default) reproduces the unweighted index exactly.

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

list(index, triangles, edges). index is in `[0, 1]`, 1 = fully convex.
`triangles` holds the CDT mesh. `edges` holds one row per evaluated
line.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]

# Wake: a simple, nearly-convex Piedmont county
convexity_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9948142
# Dare: Outer Banks mainland piece plus a separated barrier-island strip
convexity_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.7375487

# deterministic = FALSE: Monte Carlo (random-line) estimate instead of
# the exhaustive all-pairs method - a seed makes it reproducible
convexity_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> [1] 0.9923343

# weight: substitutes for each triangle's own area - weighting by the
# triangle's own area exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
convexity_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.9948142
```
