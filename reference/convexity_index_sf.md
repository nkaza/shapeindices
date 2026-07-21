# convexity_index() for every row of an sf data frame

Each row indexed independently and unweighted (see shape_indices_sf()
for weighting a collection of rows as one shape).

## Usage

``` r
convexity_index_sf(x, ...)
```

## Arguments

- x:

  an sf data frame

- ...:

  passed to convexity_index() for every row (e.g. deterministic, n_quad,
  n_lines)

## Value

`x` with one new column, convexity_index

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
res <- convexity_index_sf(nc[1:5, ])
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
res$convexity_index
#> [1] 0.9971652 0.9862805 0.9986526 0.8047873 0.9766696

# `...` passes through to convexity_index() for every row - e.g. force
# the Monte Carlo estimator (deterministic = FALSE) for all five counties
res_rli <- convexity_index_sf(nc[1:5, ], deterministic = FALSE, n_lines = 2000, seed = 1)
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4314, lon_0 = -81.4982) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 276 triangle-pairs that deterministic = TRUE (24 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4911, lon_0 = -81.1251) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 253 triangle-pairs that deterministic = TRUE (23 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4125, lon_0 = -80.6857) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 300 triangle-pairs that deterministic = TRUE (25 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4071, lon_0 = -76.0272) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 406 triangle-pairs that deterministic = TRUE (29 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 36.4224, lon_0 = -77.4105) before computing - pass already-projected data instead if you need a specific CRS.
#> Warning: n_lines (2000) is not substantially lower than the 465 triangle-pairs that deterministic = TRUE (31 triangles) would evaluate for this same polygon; deterministic = FALSE is meant as a cheaper approximation for meshes too large to enumerate exhaustively - consider deterministic = TRUE instead, or a smaller n_lines.
res_rli$convexity_index
#> [1] 0.9949545 0.9844751 0.9984327 0.8028660 0.9753344
```
