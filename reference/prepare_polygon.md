# Validate and triangulate a polygon once

Builds the `prep` object every mesh-based index accepts, so one
triangulation can be shared across several index calls on the same
polygon (see shape_indices(), which does this automatically).

## Usage

``` r
prepare_polygon(poly, simplify_tolerance = NULL)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

- simplify_tolerance:

  optional st_simplify() tolerance. NULL (default) skips simplification
  entirely. Applied after cleaning, before triangulating - trades a
  small, bounded amount of boundary detail for a smaller mesh (fewer CDT
  triangles and fewer total boundary edges, both of which reduce every
  downstream index's cost and memory use). Always in the geometry's OWN
  CRS unit at the point simplification runs, NOT necessarily metres:
  geographic (lon/lat) input gets auto-projected to a metric CRS first
  (always metres), but input already in a projected CRS keeps whatever
  linear unit that CRS uses - many US State Plane CRSs are US survey
  feet, for example. Warns if that unit isn't metres, since a
  silently-wrong unit makes the same tolerance value mean a very
  different amount of simplification than intended.

## Value

list(poly, tri) - cleaned, planar geometry and its triangle mesh

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
prep <- prepare_polygon(nc[nc$NAME == "Wake", ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
nrow(prep$tri)
#> [1] 24
```
