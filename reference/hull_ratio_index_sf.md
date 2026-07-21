# hull_ratio_index() for every row of an sf data frame

Each row indexed independently (see shape_indices_sf() for a version
that also computes the other indices).

## Usage

``` r
hull_ratio_index_sf(x)
```

## Arguments

- x:

  an sf data frame

## Value

`x` with one new column, hull_ratio_index

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
res <- hull_ratio_index_sf(nc[1:5, ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
res$hull_ratio_index
#> [1] 0.9154096 0.8256371 0.9176972 0.5467063 0.8088957
```
