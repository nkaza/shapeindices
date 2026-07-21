# Convex hull area ratio of a (multi)polygon

Convex hull area ratio of a (multi)polygon

## Usage

``` r
hull_ratio_index(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

## Value

list(index, area, hull_area, hull). index is in (0, 1\], 1 = the polygon
is convex (equals its own convex hull).

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
hull_ratio_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.904643
hull_ratio_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.2647488
```
