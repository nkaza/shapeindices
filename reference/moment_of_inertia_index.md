# Polar moment-of-inertia compactness/dispersal index of a (multi)polygon.

Polar moment-of-inertia compactness/dispersal index of a (multi)polygon.

## Usage

``` r
moment_of_inertia_index(
  poly,
  prep = NULL,
  weight = NULL,
  simplify_tolerance = NULL
)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating

- weight:

  optional numeric vector, one entry per triangle in `prep$tri`, used as
  that triangle's mass instead of its own area. NULL (default)
  reproduces the unweighted index exactly. See the concentric-rings note
  above for why the reference changes once density is non-uniform.

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly, since simplification only happens while
  building `prep` in the first place.

## Value

list(index, J, Ixx, Iyy, Ixy, J_ref, area, total_weight, centroid,
triangles). `centroid` is the plain geometric centroid when weight =
NULL, or the true mass centroid when weighted - J/Ixx/Iyy/Ixy are always
computed about this point (`Ixy`, the product of inertia, isn't used by
this index at all - it's here for
[`moment_isotropy_index()`](https://nkaza.github.io/shapeindices/reference/moment_isotropy_index.md),
which shares this same triangulation and mass centroid). `total_weight`
is sum(weight) as supplied (before normalisation), or area when weight =
NULL.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
moment_of_inertia_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.8556213
# Camden: a long, narrow county - convex, but far from a disk
moment_of_inertia_index(nc[nc$NAME == "Camden", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4012, lon_0 = -76.234) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.3196296

# weight: substitutes for each triangle's own mass - weighting by the
# triangle's own area exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
moment_of_inertia_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.8556213
```
