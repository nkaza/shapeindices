# 3-point interior (Hammer-Stroud) quadrature rule for a triangle

3-point interior (Hammer-Stroud) quadrature rule for a triangle

## Usage

``` r
tri_quad_points(coords3x2)
```

## Arguments

- coords3x2:

  a 3x2 matrix of the triangle's corner coordinates

## Value

a 3x2 matrix of quadrature point coordinates

## Examples

``` r
triangle <- rbind(c(0, 0), c(4, 0), c(0, 3))
tri_quad_points(triangle)
#>           [,1] [,2]
#> [1,] 0.6666667  0.5
#> [2,] 2.6666667  0.5
#> [3,] 0.6666667  2.0
```
