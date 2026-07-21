# Polsby-Popper compactness score of a (multi)polygon

4*pi*area/perimeter^2, in (0, 1\], 1 = a circle. Perimeter is the full
boundary length (outer ring plus any holes).

## Usage

``` r
polsby_popper_index(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

## Value

list(index, area, perimeter)

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
polsby_popper_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.6285667
polsby_popper_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.08229314
```
