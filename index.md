# shapeindices

The goal of shapeindices is to create indices for sf polygons that
describe its shape. Only indices that are (0,1\] are considered. It does
so by triangulating the polygon as well as creating some bounding or
reference shapes. It explicitly accounts for holes and multi-polygons

## Installation

You can install the development version of shapeindices

``` r

devtools::install_github(nkaza/shapeindices)
```

## Example

Here is a basic example of usage.

``` r

library(shapeindices)

square <- sf::st_sfc(sf::st_polygon(list(rbind(c(0,0), c(10,0), c(10,10), c(0,10), c(0,0)))), crs = "+proj=cartesian")
moment_of_inertia_index(square)

#> $index
#> [1] 0.9549297
#> 
#> $J
#> [1] 1666.667
#> 
#> $Ixx
#> [1] 833.3333
#> 
#> $Iyy
#> [1] 833.3333
#> 
#> $Ixy
#> [1] 2.220446e-14
#> 
#> $J_ref
#> [1] 1591.549
#> 
#> $area
#> [1] 100
#> 
#> $total_weight
#> [1] 100
#> 
#> $centroid
#> Geometry set for 1 feature 
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 5 ymin: 5 xmax: 5 ymax: 5
#> CRS:           NA
#> POINT (5 5)
#> 
#> $triangles
#> Simple feature collection with 2 features and 2 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 0 ymin: 0 xmax: 10 ymax: 10
#> CRS:           NA
#>   tri_id area                       geometry
#> 1      1   50 POLYGON ((0 10, 0 0, 10 0, ...
#> 2      2   50 POLYGON ((10 0, 10 10, 0 10...
```
