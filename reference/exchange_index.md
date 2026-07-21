# Exchange compactness score of a (multi)polygon

The share of the polygon's own area that falls inside the equal-area
circle centred at its centroid, in `[0, 1]`, 1 = a circle. Angel, Parent
& Civco (2010) introduce this as a natural metric for gerrymandering: a
district that "reaches out" to grab distant voters, or excludes nearby
ones, has most of its area outside that circle.

## Usage

``` r
exchange_index(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

## Value

list(index, area, circle_area, circle). `circle` is the equal-area
reference circle itself, as an sfc POLYGON, for plotting.

## Details

KNOWN LIMITATION: for a multi-part polygon, the reference circle is
centred at the *overall* centroid, which for well-separated parts sits
in the empty space between them. Once the circle's radius is smaller
than the distance to the nearest part, the intersection - and so the
index - is exactly 0, not just low; two individually-compact pieces a
few units apart can already score under 0.05. Angel et al. themselves
note the same failure mode on real districts split by water. Included
despite this, the same way
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
ships a real, useful, but documented-blind-spot metric rather than
omitting it.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
exchange_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.8233315
exchange_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.6477935
```
