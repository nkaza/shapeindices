# Convex decomposition of a (multi)polygon (Hertel-Mehlhorn).

Not currently called by any of this package's own index functions -
every index works directly off
[`cdt_triangles()`](https://nkaza.github.io/shapeindices/reference/cdt_triangles.md)'s
triangles (optionally refined finer by
[`subdivide_mesh()`](https://nkaza.github.io/shapeindices/reference/subdivide_mesh.md)),
never off merged convex pieces. This is a standalone utility for a
caller who wants an approximate convex decomposition for their own
purposes; see
[`subdivide_mesh()`](https://nkaza.github.io/shapeindices/reference/subdivide_mesh.md)
for the opposite direction (more, smaller triangles rather than fewer,
larger convex pieces).

## Usage

``` r
convex_decompose(poly, prep = NULL, simplify_tolerance = NULL)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

- prep:

  optional pre-computed list(poly, tri) from prepare_polygon(), to skip
  re-triangulating

- simplify_tolerance:

  passed to
  [`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md)
  when `prep` is NULL - see its own doc. Ignored (with no error) if
  `prep` is supplied directly, since simplification only happens while
  building `prep` in the first place.

## Value

an sf data frame shaped like cdt_triangles()'s output (piece_id, area,
geometry), or NULL if the polygon triangulates to no pieces

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
tri    <- cdt_triangles(nc[nc$NAME == "Dare", ])
pieces <- convex_decompose(nc[nc$NAME == "Dare", ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
nrow(tri)     # many small triangles
#> [1] 15
nrow(pieces)  # fewer, larger convex pieces
#> [1] 10
```
