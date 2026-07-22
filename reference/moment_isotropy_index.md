# Moment isotropy index of a (multi)polygon's mass distribution

Ratio of the smaller to larger principal moment of the mass inertia
tensor, in `(0, 1]`. `1` means the mass distribution is rotationally
isotropic about its own centroid (a disk qualifies, but so does any
shape with 3-fold or higher rotational symmetry - this is not the same
claim as convexity_index()'s or moment_of_inertia_index()'s "= 1 iff a
disk"). Lower values mean the mass is more anisotropic - elongated along
one direction - regardless of whether the shape itself is convex or
dispersed. See
[`vignette("g-understanding-moment-isotropy-index")`](https://nkaza.github.io/shapeindices/articles/g-understanding-moment-isotropy-index.md)
for the derivation, the bound proof, and how this relates to (but
differs from) classical eccentricity.

## Usage

``` r
moment_isotropy_index(
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
  reproduces the unweighted (area-based) index exactly.

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly, since simplification only happens while
  building `prep` in the first place.

## Value

list(index, lambda_min, lambda_max, Ixx, Iyy, Ixy, area, total_weight,
centroid, triangles). `centroid`/`Ixx`/`Iyy`/`Ixy` are exactly
moment_of_inertia_index()'s own fields, computed about the same mass
centroid.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
wake <- nc[nc$NAME == "Wake", ]
moment_isotropy_index(wake)$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.4230179
# Camden: long and narrow - low isotropy despite being convex
moment_isotropy_index(nc[nc$NAME == "Camden", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4012, lon_0 = -76.234) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.03639716

# weight: substitutes for each triangle's own mass, same as
# moment_of_inertia_index() - weighting by the triangle's own area
# exactly reproduces the unweighted index
prep <- prepare_polygon(wake)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
moment_isotropy_index(wake, prep = prep, weight = prep$tri$area)$index
#> [1] 0.4230179
```
