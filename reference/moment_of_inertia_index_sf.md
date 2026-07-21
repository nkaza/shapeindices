# Vectorised wrapper for an sf data frame - each row indexed independently.

Vectorised wrapper for an sf data frame - each row indexed
independently.

## Usage

``` r
moment_of_inertia_index_sf(x, ...)
```

## Arguments

- x:

  an sf data frame

- ...:

  passed to moment_of_inertia_index() for every row

## Value

`x` with one new column, moi_index

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
res <- moment_of_inertia_index_sf(nc[1:5, ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
res$moi_index
#> [1] 0.9543717 0.7390698 0.9188184 0.3589253 0.6646385
```
