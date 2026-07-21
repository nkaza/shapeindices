# Constrained Delaunay triangulation of a (multi)polygon

Triangulates each part separately so disjoint parts never bridge, and
hole boundaries are always respected as constraints.

## Usage

``` r
cdt_triangles(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON (length 1); does not run
  st_make_valid() itself - callers are expected to have cleaned it up

## Value

an sf data frame with one row per triangle (tri_id, area, geometry), or
NULL if the polygon has no triangulatable area

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
tri <- cdt_triangles(nc[nc$NAME == "Wake", ])
nrow(tri)
#> [1] 24
plot(sf::st_geometry(tri), border = "grey40")
```
