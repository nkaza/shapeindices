# Adaptive subdivision of a CDT mesh into a finer triangle mesh.

Standalone utility, not currently called by any of this package's own
index functions -
[`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)
does the same area-adaptive subdivision internally, but collapses each
sub-triangle straight to a weighted centroid rather than keeping
triangle geometries around, so it doesn't go through this function. See
[`convex_decompose()`](https://nkaza.github.io/shapeindices/reference/convex_decompose.md)
for the opposite direction (fewer, larger convex pieces rather than
more, smaller triangles).

## Usage

``` r
subdivide_mesh(poly, prep = NULL, max_depth = 4, simplify_tolerance = NULL)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating

- max_depth:

  subdivision depth ceiling, applied to the mesh's own largest triangle
  only (4^max_depth children for that triangle; fewer for smaller ones,
  adaptively - see `.adaptive_tri_depth()`)

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly, since simplification only happens while
  building `prep` in the first place.

## Value

an sf data frame shaped like
[`cdt_triangles()`](https://nkaza.github.io/shapeindices/reference/cdt_triangles.md)'s
output (tri_id, area, geometry), or NULL if the polygon triangulates to
no pieces

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
tri  <- cdt_triangles(nc[nc$NAME == "Dare", ])
fine <- subdivide_mesh(nc[nc$NAME == "Dare", ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
nrow(tri)    # coarse CDT triangles
#> [1] 15
nrow(fine)   # more, smaller triangles
#> [1] 1296
```
