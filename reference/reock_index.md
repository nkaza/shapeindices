# Reock compactness score of a (multi)polygon

area/area(minimum bounding circle), in (0, 1\], 1 = a circle.

## Usage

``` r
reock_index(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

## Value

list(index, area, mbc_area, mbc). `mbc` is the minimum bounding circle
itself, as an sfc POLYGON, for plotting.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
reock_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.5379165
reock_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.09032822
```
