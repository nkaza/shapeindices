# Detour compactness score of a (multi)polygon

The ratio of the equal-area circle's perimeter to the perimeter of the
polygon's own convex hull, in (0, 1\], 1 = a circle (whose hull is
itself). Angel, Parent & Civco (2010) introduce this to measure how hard
a shape is to circumnavigate as an obstacle - a smooth, rounded hull
scores high regardless of how convoluted the actual boundary inside that
hull is, since only the hull's own perimeter enters the ratio.

## Usage

``` r
detour_index(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

## Value

list(index, area, hull_perimeter, hull). `hull` is the convex hull
itself, as an sfc POLYGON, for plotting.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
detour_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.8643451
detour_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.4126559
```
