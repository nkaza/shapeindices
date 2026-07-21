# polsby_popper_index() for every row of an sf data frame

polsby_popper_index() for every row of an sf data frame

## Usage

``` r
polsby_popper_index_sf(x)
```

## Arguments

- x:

  an sf data frame

## Value

`x` with one new column, polsby_popper_index

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
res <- polsby_popper_index_sf(nc[1:5, ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
res$polsby_popper_index
#> [1] 0.71217996 0.53389426 0.69442497 0.09600441 0.42538411
```
