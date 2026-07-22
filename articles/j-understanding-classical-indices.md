# 10. Understanding the Classical Compactness Metrics

Code

``` r

library(shapeindices)
library(sf)
library(ggplot2)

theme_set(theme_minimal(base_size = 11))
theme_gallery <- theme_void(base_size = 10) +
  theme(strip.text = element_text(size = 9, face = "bold"))
```

## 1 Introduction

Six of this package’s thirteen indices need no triangulation at all -
just a polygon’s own area, perimeter, convex hull, or bounding shapes,
all cheap closed-form calculations:
[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md),
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md),
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md),
[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md),
[`detour_index()`](https://nkaza.github.io/shapeindices/reference/detour_index.md),
and
[`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md).
This vignette covers all six together, the way
[`vignette("k-nc-counties-comparison")`](https://nkaza.github.io/shapeindices/articles/k-nc-counties-comparison.md)
covers the other ten together on real data - not deriving each from
scratch (their math is well-established, cited below), but showing
concretely how they behave, and where each one breaks.

Each of these has a real history, reviewed comprehensively in Frolov
(1975),[^1] MacEachren (1985),[^2] and, most recently, Murray
(2025).[^3] MacEachren’s own categorisation - perimeter/area,
circle-related, comparison-to-a-reference-shape,
dispersion-around-a-central-point - still organises the field:
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
is a perimeter/area measure,
[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)
and
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
are circle- and box-related respectively,
[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)
and
[`detour_index()`](https://nkaza.github.io/shapeindices/reference/detour_index.md)
compare the shape to its own convex hull, and
[`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)
compares it to a reference circle placed at its centroid.

- **[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)**:
  area of the polygon over the area of its convex hull. First formalised
  by Dori & Ben-Bassat (1983);[^4] MacEachren’s “comparison to shape”
  category.
- **[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)**:
  the isoperimetric quotient, $`4\pi A/p^2`$. The underlying fact - a
  circle has the least perimeter for a given area - goes back to
  Zenodorus and, per some accounts, Pythagoras;[^5] the specific
  redistricting application is Polsby & Popper (1991),[^6] itself
  building on Ritter’s 1822 proposal reviewed in Frolov.
- **[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)**:
  the shorter over the longer side of the polygon’s minimum-area
  bounding rectangle, at any rotation
  ([`sf::st_minimum_rotated_rectangle()`](https://r-spatial.github.io/sf/reference/geos_unary.html)) -
  not the axis-aligned bounding box, which would make the score depend
  on orientation rather than shape alone. A common,
  informally-attributed redistricting metric with no single definitive
  origin in the reviews consulted here - included because it’s simple,
  cheap, and widely used, not because of a strong literature pedigree.
- **[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)**:
  area of the polygon over the area of its minimum bounding circle.
  Commonly attributed to Reock (1961),[^7] with origins traced to
  Ehrenburg (1892) by Frolov.
- **[`detour_index()`](https://nkaza.github.io/shapeindices/reference/detour_index.md)**
  and
  **[`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)**:
  both Angel, Parent & Civco (2010),[^8] covered together with the
  package’s other ten indices in
  [`vignette("e-understanding-span-index")`](https://nkaza.github.io/shapeindices/articles/e-understanding-span-index.md)’s
  Introduction (Cohesion) and in this vignette’s own sections below.

## 2 Mathematical formulae

Each of the six is a ratio comparing the actual polygon to some
reference shape or bound. All are already provably in $`(0, 1]`$
([`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)’s
ceiling too), 1 = a circle, except where noted - the proofs are
established results, not derived here.

| Index | Formula | 1 means |
|:---|:---|:---|
| \`hull_ratio_index()\` | \$A(P) / A(\text{hull}(P))\$ | the polygon equals its own convex hull |
| \`polsby_popper_index()\` | \$4\pi A(P) / p(P)^2\$ | the polygon is a circle |
| \`width_length_ratio_index()\` | \$\min(w, h) / \max(w, h)\$ | the minimum bounding rectangle is a square |
| \`reock_index()\` | \$A(P) / A(\text{MBC})\$ | the polygon equals its own minimum bounding circle |
| \`detour_index()\` | \$2\sqrt{\pi A(P)} / p(\text{hull}(P))\$ | the convex hull is a circle |
| \`exchange_index()\` | \$A(P \cap C) / A(P)\$ | the polygon equals its own equal-area circle |

$`p`$ is perimeter, $`w`$/$`h`$ the bounding-box width/height, MBC the
minimum bounding circle (radius found by the classic Welzl-style
incremental algorithm), and $`C`$ the circle of area $`A(P)`$ centred at
$`P`$’s own centroid
([`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)’s
reference). None of the six take a `weight` argument - see “Holes,
multi-part shapes, and what weighting actually does” below for why.

## 3 Illustrations: canonical shapes

Code

``` r

make_square <- function(half = 5) st_polygon(list(rbind(
  c(-half, -half), c(half, -half), c(half, half), c(-half, half), c(-half, -half))))
make_rectangle <- function(w, h) st_polygon(list(rbind(
  c(0, 0), c(w, 0), c(w, h), c(0, h), c(0, 0))))
make_disk <- function(r = 5, n = 60) st_buffer(st_sfc(st_point(c(0, 0))), dist = r, nQuadSegs = n)[[1]]
make_star <- function(n_points, r_outer, r_inner) {
  n <- n_points * 2
  angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n + 1)[1:n]
  radii  <- rep(c(r_outer, r_inner), n_points)
  x <- radii * cos(angles); y <- radii * sin(angles)
  st_polygon(list(rbind(cbind(x, y), c(x[1], y[1]))))
}

canonical <- list(
  square         = make_square(5),
  disk           = make_disk(5.64),   # same area as the square, for a fair side-by-side
  `rect 2:1`     = make_rectangle(14.14, 7.07),
  `rect 10:1`    = make_rectangle(31.6, 3.16),
  `hexagon`      = st_polygon(list(rbind(t(sapply(seq(0, 300, 60) * pi / 180,
                     function(a) c(5.5 * cos(a), 5.5 * sin(a)))), c(5.5, 0)))),
  `star (mild)`  = make_star(6, 5.64, 3.5),
  `star (sharp)` = make_star(6, 5.64, 0.6)
)
```

| shape | name | hull_ratio_index | polsby_popper_index | width_length_ratio_index | reock_index | detour_index | exchange_index |
|:---|:---|---:|---:|---:|---:|---:|---:|
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAGFBMVEUAAABMTExJSUlNTU1LS0tPT09cXFy/v78JXKQ0AAAACHRSTlMAGxwhIv///1AOH2YAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAA0SURBVDiNY2AYBQMFmFgwABM2dYysbBiAlRGLQmY2dgzAxjyqcFThcFFIdFYgOnONAroAAF6zC6n/x/coAAAAAElFTkSuQmCC) | square | 1.000 | 0.785 | 1.000 | 0.637 | 0.886 | 0.909 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAQlBMVEUAAABAQEBJSUlOTk5OTk5PT09MTExMTExOTk5NTU1OTk5RUVFVVVVaWlpubm58fHx9fX2lpaWzs7O9vb2+vr6/v7+vRdNrAAAAFnRSTlMABAcNGh0ef4C9z////////////////88EeAAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAJNJREFUOI3tk9kKgzAQRWPcJqtZnP//1YqUoiFpL4j0xfN8mHAnd4R4+DdymuRvayBtnDOahq9aRzbEzJxjsNS1vV75xG+SV31znlpW/rAuqjVz9gdvM/3cyGETn0i2nogCFwSqijqWYtQ1T5pcitnUNj+60mN24xURfhoOg68HXjj8hXgp4JrhxRXwKexgx/VwMy+E+hiBDo0ywQAAAABJRU5ErkJggg==) | disk | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 | 1.000 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAIVBMVEUAAAB/f39QUFBNTU1MTExOTk5OTk5SUlJcXFxzc3O/v7/BGbJLAAAAC3RSTlMAAiAhSk5V/////9mht9oAAAAJcEhZcwAADsMAAA7DAcdvqGQAAAA2SURBVDiNY2AYBaOAcsDIyoYBWBixKGRi58QA7ExYFDJzcGEADuZRhaQqJDrAiY7CUTAKSAUAu8gIpx9CeEEAAAAASUVORK5CYII=) | rect 2:1 | 1.000 | 0.698 | 0.500 | 0.509 | 0.836 | 0.742 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAFVBMVEUAAABmZmZNTU1NTU1NTU1TU1OCgoLbzv4MAAAAB3RSTlMABSFxd///xRPWjAAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAClJREFUOI1jYBgFo2BkAEZmFgzAzIhFIRMrGwZgZaJEIdFWj4JRMBwBABTDAZVN9DDcAAAAAElFTkSuQmCC) | rect 10:1 | 1.000 | 0.260 | 0.100 | 0.126 | 0.510 | 0.352 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAilBMVEUAAAAAAABmZmZVVVVRUVFKSkpSUlJOTk5OTk5NTU1PT09NTU1MTExOTk5NTU1MTExOTk5MTExNTU1NTU1NTU1MTExOTk5NTU1NTU1NTU1NTU1NTU1YWFhjY2NkZGRvb29wcHCAgICBgYGCgoKRkZGSkpKTk5OYmJijo6Ourq6wsLC5ubm+vr6/v79Qw+KJAAAALnRSTlMAAQUGFhgZMTQ1RFZXWHR1dpqbnMHF2drd4Orr////////////////////////eQVSMAAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAJpJREFUOI3t0l8TQkAUh+ETlZIShZQifyLs9/96bkx22Z8x01WN9/qZ2TlnD9Hcf6W7lzZXH3Oql5dtuaeOwGPIPkUmdsvru4PlbQ3hKWFcsY2cFlQ8rIItgOcXE8qchdTt77UI64chc4pfsF6Fr3wDyRg+vQPDZCJMwTDD9WzkjsiOefi0kOt/4QpCMqMOhgfspp/Z5MOd+8EaTFIrx5fzJT0AAAAASUVORK5CYII=) | hexagon | 1.000 | 0.907 | 0.866 | 0.827 | 0.952 | 0.963 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAh1BMVEUAAABJSUlSUlJMTExJSUlPT09NTU1MTExPT09OTk5MTExNTU1NTU1MTExMTExOTk5NTU1MTExNTU1NTU1NTU1MTExNTU1NTU1NTU1NTU1dXV1eXl5gYGBhYWFjY2N/f3+GhoaHh4eOjo6WlpaXl5ebm5uoqKipqamqqqq0tLS9vb2+vr6/v79lMrknAAAALXRSTlMABxkbHDc/QERFW2BxgpOur7GztcHFxujp6v////////////////////////8uyDJeAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAAnklEQVQ4je2TORKCQBBFAcWFTRBFAWdYRlCGvv/5tIrASdr6CUXCi1/13pa1sji2DYpxjHmeEB7iOff+nW8B8VQT1ee/iusHYZRITaRlEoWB7zJiKqumVQN9GVTbVDJlxH2hyWAsd1zuy9MU1ZUt8vgYjYDiwLeTdT+xy3gPFuHUaDPweOCBwyucAI5iwrn1r3wDiPDh4q+AP9fKjHwAvlYegAH4WesAAAAASUVORK5CYII=) | star (mild) | 0.717 | 0.524 | 0.866 | 0.593 | 0.806 | 0.880 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAbFBMVEUAAABVVVVJSUlMTExLS0tJSUlOTk5MTExNTU1MTExNTU1MTExNTU1MTExMTExNTU1MTExNTU1NTU1MTExMTExOTk5NTU1NTU1NTU1NTU1MTExNTU1MTExNTU1NTU1NTU1NTU1NTU1NTU1wcHBud25aAAAAJHRSTlMABgcKERUuTU9eZnJ3fJCRlpmfoKqur73M29zl5ujp6uv7///EPhhmAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAAl0lEQVQ4je2TuRKDMAxEARNiQiCASbgv8f//iIBxabxNhoYtrOaNZK0kx7l1uTwPBOMYBD8JCKocBH9fCxCII5blEUVgANM+2tGq2rGoT00pn+34ZrSp2aJXO4Tm4n5GjE4jY1Q8Tv8pO+rmmR95im1JFS0LKd/GwSBYWjfjWppB7YENh0eoZV0KLXjN4MWFTwE+rlt/1Aoz7AoMxqWCzwAAAABJRU5ErkJggg==) | star (sharp) | 0.123 | 0.034 | 0.866 | 0.102 | 0.334 | 0.487 |

The square and disk anchor the two ends familiar from every other
vignette in this package: the disk reads at or near 1 on all six (its
minimum bounding rectangle is the one exception to “1 = circle” - a
circle’s minimum bounding rectangle is a square regardless of how it’s
drawn, so
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
reads 1 for both), while the square is already informative on its own -
[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)/[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
read exactly 1 (a square is its own hull and its own minimum bounding
rectangle), but
[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)
is capped at exactly $`2/\pi \approx 0.637`$ regardless of size, since a
square’s minimum bounding circle always wastes its four corners.
Elongating the rectangle moves every index down together but at
different rates, and the two stars - same six points, shrinking only the
inner radius - separate the six most clearly; see “Key pathological
properties” below for why.

## 4 Holes, multi-part shapes, and what weighting actually does

``` r

make_square_with_hole <- function(outer_half = 5, hole_frac = 0.3) {
  outer <- rbind(c(-outer_half, -outer_half), c(outer_half, -outer_half),
                 c(outer_half, outer_half), c(-outer_half, outer_half), c(-outer_half, -outer_half))
  hh <- outer_half * sqrt(hole_frac)
  hole <- rbind(c(-hh, -hh), c(-hh, hh), c(hh, hh), c(hh, -hh), c(-hh, -hh))
  st_polygon(list(outer, hole))
}
make_dumbbell_gap <- function(gap) {
  sq1 <- st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
  sq2 <- st_polygon(list(rbind(c(2 + gap, 0), c(4 + gap, 0), c(4 + gap, 2), c(2 + gap, 2), c(2 + gap, 0))))
  st_union(st_sfc(sq1, sq2))
}
```

| shape | hole % | hull_ratio_index | polsby_popper_index | width_length_ratio_index | reock_index | detour_index | exchange_index |
|:---|:---|---:|---:|---:|---:|---:|---:|
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAG1BMVEUAAABMTExJSUlNTU1LS0tPT09cXFyqqqq/v78eyHlVAAAACXRSTlMAGxwhIv////9bgcpyAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAAPElEQVQ4jWNgGAUDBZhYMAATNnWMrGwYgJURi0JmNg4MwMY8qpB8hezs1FY4FHw9eBUSnRWIzlyjgC4AAF/vDOm9LfriAAAAAElFTkSuQmCC) | 0% | 1.0 | 0.785 | 1 | 0.637 | 0.886 | 0.909 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAKlBMVEUAAABMTExJSUlNTU1LS0tNTU1NTU1NTU1PT09cXFyBgYGCgoKkpKS/v792rOg+AAAADnRSTlMAGxwhInR3tf///////3UCcl0AAAAJcEhZcwAADsMAAA7DAcdvqGQAAABQSURBVDiN7dM3DgAgCEBRO2K5/3UtKyYwWRLfRvIXJSj1nWIcYVadjkhEvQgtVgLt22HJU2HDBKGDxId+DP7uUPoY8fe8sGs+FJ+C+Li+LRpYrRKFW3E5twAAAABJRU5ErkJggg==) | 10% | 0.9 | 0.408 | 1 | 0.573 | 0.841 | 0.849 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAKlBMVEUAAABMTExJSUlNTU1LS0tLS0tPT09MTExPT09cXFxkZGRmZmaCgoK/v79zT+I1AAAADnRSTlMAGxwhIjM6Yf///////yLl7PkAAAAJcEhZcwAADsMAAA7DAcdvqGQAAABQSURBVDiNY2AYBQMFmFgwABM2dYwcnBiAgxGLQmZOXgzAyUw/hTxccMCDVyE3OysUsHPjV8gG47KNKqREIbEBTnQUDkwyIzorEJ25RgFdAACI1w919PY3QwAAAABJRU5ErkJggg==) | 30% | 0.7 | 0.230 | 1 | 0.446 | 0.741 | 0.571 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAALVBMVEUAAABMTExJSUlNTU1LS0tMTExNTU1NTU1PT09cXFyZmZmampqcnJy7u7u/v7+TQiysAAAAD3RSTlMAGxwhIqqw5f////////8i9w2GAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAAUklEQVQ4je3TOQrAMAxEUS+xZXm7/3HTCiTDVDEEvVL8chSCuyU9SrK62Fhp0Qgzb4WzHa7RhbGO4aQi0DyHVR6qh38O0VHAM8OGC78C/FzuEy/RZQstynbDgQAAAABJRU5ErkJggg==) | 50% | 0.5 | 0.135 | 1 | 0.318 | 0.627 | 0.091 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAJ1BMVEUAAABMTExJSUlNTU1LS0tNTU1NTU1NTU1PT09cXFyAgICCgoKkpKR3id1VAAAADXRSTlMAGxwhInF3s///////+QJcYgAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAElJREFUOI3t0zsOACAIREH/CHr/81ormGylDVOSVy4huF9SUZLVxc5Kj0aYeY7D5GyGQm1DcgvrfqkeeoiH6MzA4cKvAD+Xe2IBDh8H0Q619v4AAAAASUVORK5CYII=) | 70% | 0.3 | 0.070 | 1 | 0.191 | 0.485 | 0.000 |

| shape | gap | hull_ratio_index | polsby_popper_index | width_length_ratio_index | reock_index | detour_index | exchange_index |
|:---|---:|---:|---:|---:|---:|---:|---:|
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAJ1BMVEUAAABQUFBNTU1OTk5MTExNTU1NTU1YWFhcXFx0dHR+fn6oqKi/v7/xYwIwAAAADXRSTlMAECFFbsTM////////Tzqj4AAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAEJJREFUOI1jYBgFo4COgJGVDQqYmWEsVkYsCpnYuSGAk4WFE8pkZ8KmkIMHArhYWLigTI7hr5DY4CE6wEfBKKAZAAAN0AeBCG6eVgAAAABJRU5ErkJggg==) | 0.5 | 0.889 | 0.393 | 0.444 | 0.420 | 0.771 | 0.617 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAKlBMVEUAAABOTk5NTU1MTExMTExNTU1NTU1MTExYWFhcXFxlZWVubm6hoaG/v78WOSP+AAAADnRSTlMADSElSrO5u////////3o6S28AAAAJcEhZcwAADsMAAA7DAcdvqGQAAABGSURBVDiN7dA7DgAgCERBFfwgev/rmsXWwk4LpnvJNhCCcz+KuUIhBBWLHA/D1AZ0RnC3aOk0lAm6h2ohvw1vj7l+j3PvLdgJBPHkfgXLAAAAAElFTkSuQmCC) | 2.0 | 0.667 | 0.393 | 0.333 | 0.255 | 0.627 | 0.242 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAKlBMVEUAAABGRkZOTk5NTU1KSkpOTk5NTU1MTExWVlZYWFhcXFxeXl6ZmZm/v79sNda3AAAADnRSTlMACw0hJqGiqv///////96oRUcAAAAJcEhZcwAADsMAAA7DAcdvqGQAAABISURBVDiNY2AYBaNgSANGVnZ2djYmGJeJDchlZcSikJmDh4eHkwXGZeEEcjmYsSnk4uXl5UYo5AZyuShTSKzVRHtmFIyCIQQAIoACicKt1rwAAAAASUVORK5CYII=) | 5.0 | 0.444 | 0.393 | 0.222 | 0.120 | 0.456 | 0.000 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAG1BMVEUAAAAAAABJSUlKSkpOTk5OTk5MTExSUlJgYGBp0GzRAAAACXRSTlMAARUfS1V4//9PnCXHAAAACXBIWXMAAA7DAAAOwwHHb6hkAAAAMElEQVQ4jWNgGAWjYGQARhZWJlQRJlYWRiwKmdk52FBF2DjYmSlRSLTVo2AUDEcAANO3AJH2v+FsAAAAAElFTkSuQmCC) | 20.0 | 0.167 | 0.393 | 0.083 | 0.018 | 0.193 | 0.000 |

A hole never touches the convex hull or the minimum bounding rectangle,
so
[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)/[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
respond only through the numerator:
[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)
tracks `1 - hole_frac` exactly (the table above confirms this
precisely), and
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
doesn’t move at all.
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
takes a double hit instead - the hole’s own boundary adds directly to
the total perimeter on top of the area loss - which is why it falls
fastest of the six here.
[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)
degrades at the same rate as
[`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)
(same mechanism, fixed denominator), and
[`detour_index()`](https://nkaza.github.io/shapeindices/reference/detour_index.md)
more gently still (its denominator is hull *perimeter*, unaffected by
the hole, while its numerator only shrinks with $`\sqrt{\text{area}}`$).
[`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)
falls hardest of all and reaches exactly 0 by a 70% hole here - a large,
roughly-centred hole empties out exactly the region the reference circle
would occupy, so the circle it draws around the centroid finds nothing
left to overlap.

For multi-part shapes,
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
shows something worth flagging on its own: it barely moves at all as the
gap widens (0.393 at gap 0.5 vs. 0.393 at gap 20) - both area and
perimeter are simple sums over disjoint pieces, indifferent to how far
apart those pieces actually are, so
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
is structurally blind to dispersal as long as the pieces don’t touch.
[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)
moves the most: an ever-more-distant pair of squares needs an
ever-larger minimum bounding circle to enclose both, so its area
balloons quadratically with the gap even while the true combined area
stays fixed.
[`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)
shows its own already-documented pathology here - a hard collapse to
exactly 0, not merely a low score, once the reference circle can no
longer reach either piece.

**What weighting actually does to these six.** None of them take a
`weight` argument - there’s nothing for a per-triangle or per-row weight
to attach to beyond a binary in/out decision (see
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)’s
and
[`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)’s
own roxygen for why a magnitude-based weighted form isn’t mathematically
defensible for either). The only way `weights` can move any of the six,
via `shape_indices_sf(byrow = FALSE)`, is that a row weighted exactly 0
or `NA` is excluded as a hole *before* the union is built - changing
which polygon is measured, not how it’s measured:

``` r

sq1 <- st_polygon(list(rbind(c(0, 0), c(2, 0), c(2, 2), c(0, 2), c(0, 0))))
sq2 <- st_polygon(list(rbind(c(2, 0), c(4, 0), c(4, 2), c(2, 2), c(2, 0))))
x <- st_sf(name = c("a", "b"), geometry = st_sfc(sq1, sq2, crs = 3857))

equal_weights   <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, weights = c(1, 1),
                                                       which = c("hull_ratio", "detour", "exchange"), id = "z"))
unequal_weights <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, weights = c(1, 1000),
                                                       which = c("hull_ratio", "detour", "exchange"), id = "z"))
zeroed_weight   <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, weights = c(1, 0),
                                                       which = c("hull_ratio", "detour", "exchange"), id = "z"))

data.frame(
  case = c("weights = c(1, 1)", "weights = c(1, 1000)", "weights = c(1, 0)"),
  rbind(st_drop_geometry(equal_weights)[, c("hull_ratio_index", "detour_index", "exchange_index")],
        st_drop_geometry(unequal_weights)[, c("hull_ratio_index", "detour_index", "exchange_index")],
        st_drop_geometry(zeroed_weight)[, c("hull_ratio_index", "detour_index", "exchange_index")])
) |> knitr::kable(format = "html", digits = 4, row.names = FALSE)
```

| case                 | hull_ratio_index | detour_index | exchange_index |
|:---------------------|-----------------:|-------------:|---------------:|
| weights = c(1, 1)    |                1 |       0.8355 |         0.7420 |
| weights = c(1, 1000) |                1 |       0.8355 |         0.7420 |
| weights = c(1, 0)    |                1 |       0.8862 |         0.9094 |

`weights = c(1, 1)` and `weights = c(1, 1000)` give *identical*
results - a thousandfold difference in magnitude has zero effect, since
both squares stay in the union either way. `weights = c(1, 0)` gives
different numbers entirely, because the second square is dropped as a
hole and only the first remains - binarisation, not weighting in any
graduated sense.

## 5 Key pathological properties

### 5.1 Spikes: `width_length_ratio_index()`’s blind spot

| shape | r_inner / r_outer | hull_ratio_index | polsby_popper_index | width_length_ratio_index | reock_index | detour_index | exchange_index |
|:---|:---|---:|---:|---:|---:|---:|---:|
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAjVBMVEUAAAB/f39mZmZGRkZVVVVVVVVQUFBLS0tHR0dMTExQUFBNTU1OTk5OTk5MTExOTk5NTU1MTExNTU1NTU1MTExNTU1NTU1NTU1NTU1NTU1PT09TU1NUVFRVVVVcXFxfX19vb292dnZ3d3eAgICenp6kpKSmpqanp6eoqKirq6usrKy3t7e7u7u+vr6/v7/IxYJXAAAAL3RSTlMAAgULDA8QERIUI0lSb3KAmarb5ebo7e7v+////////////////////////////8rltfgAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACmSURBVDiN7dPJDoIwEIDhARQUEIVWccEiAi4sff/HUxNj0qatkxDjhf/8XTqdARj7e7aNYk5EaeR8ZZMFOd3u5/VyanYhLRr+rCloaGDuitX83TWP5zoXkLLjn7qSBBqYVFyoSjRw24qwTdXO2vci7HeWEs6OXCrzlNBnMmT+MOhlMjy4wx4DqTyejdrhB47+QmEpLoalAPSaAX5xAX0Kr5DHNfbTHuaNLdjNZy1RAAAAAElFTkSuQmCC) | 0.90 | 1.000 | 0.938 | 0.900 | 0.859 | 0.969 | 0.972 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAjVBMVEUAAABJSUlMTExVVVVKSkpLS0tOTk5NTU1PT09NTU1NTU1OTk5MTExNTU1OTk5MTExNTU1MTExNTU1OTk5NTU1NTU1NTU1NTU1NTU1WVlZXV1dZWVlaWlpdXV1eXl53d3d6enqBgYGCgoKFhYWKioqLi4uhoaGjo6OkpKS1tbW3t7e4uLi5ubm7u7u/v7/dvvJYAAAAL3RSTlMABwoMJik0NTc4TFVkZnZ/i4ybq6yts97g/////////////////////////////0qjfLAAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACjSURBVDiN7ZPJDoJAEERBXJBdQHRUYFgFlf7/z9OEA+OhSV2MF975JZWu7ta0hb+j66AYRZhnZZmFeKtTez8bgHjIifJgVjFtx/V80RP1wvdcxzYZMUmlLMqOPnRlIWWaMOL2OpDCcNtz2XGtig3f0vrymLyn2PDjhNUkViHvweJX9Gsm+tiow9Qx5+3QeuDC4RWOAEcxAp8ZfLj4K+DPtfBD3rZyHx+jMhtRAAAAAElFTkSuQmCC) | 0.60 | 0.693 | 0.490 | 0.866 | 0.573 | 0.793 | 0.871 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAk1BMVEUAAAAAAABVVVVJSUlHR0dMTExQUFBOTk5MTExMTExPT09OTk5NTU1MTExNTU1MTExMTExNTU1NTU1MTExOTk5MTExNTU1NTU1MTExNTU1NTU1NTU1NTU1NTU1NTU1PT09TU1NVVVVXV1dYWFhbW1tcXFxpaWlzc3OKioqZmZmfn5+hoaGkpKSmpqaoqKi4uLi/v785uh9yAAAAMXRSTlMAAQMHEhQjJCUvRFhdYWNkZXqFjJegt7m7vL3H6e7//////////////////////////kraQwAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAKxJREFUOI3t07kSgjAQBuDlEgFBuUVujYRDwfd/OkeGwobM3zg2/E1SfDOb3SREW/4eRQGh54HwkmBOynMJgvs708XC0ObF4tyaN5qxAsMydlQiv219ItWJy3AFykHRVJGdDkNqR1VTBPJqcZeNHa+nqebdyFzRMQ/Xx2vO83YUN2Rm/cf1mSl2MPwqfRI5sBl4PPDA4StcojO2E4sl8DOj5Iw5/CvAn2vLD/MG0gQTSTCasxQAAAAASUVORK5CYII=) | 0.30 | 0.346 | 0.138 | 0.866 | 0.286 | 0.560 | 0.705 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAhFBMVEUAAABmZmZVVVVJSUlAQEBVVVVMTExJSUlQUFBMTExOTk5MTExMTExMTExNTU1NTU1NTU1MTExNTU1MTExMTExOTk5NTU1MTExNTU1MTExOTk5NTU1NTU1NTU1NTU1NTU1NTU1NTU1OTk5NTU1NTU1NTU1NTU1NTU1NTU1NTU1OTk5oaGhwQ4CfAAAALHRSTlMABQYHCAkKDhAULkpNXmZxd3yPkJaXmZ2eqqusucfV1uHi4+Tl6Onq+P///8IDhUsAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACjSURBVDiNY2AYBQMOmJmJVCgoSKRCYSEiFUqIE6lQWoqAAi42CC0jA6HZuHAoFFHlZwXRsrIgkpVfVQSXkTzK6gJApfJyDAxMfCqqvLgt5xDTVBNgUVZi4lPSlOTE605uRU0FDQ0FTUVuvMqAgF1US1tbS5SdkDqiFRJpNcwzjAQ8gxI8yriDh+gAJzoKYYBgooABopMZ0QmX6KxAdOYaBTQEAE2HDJPr9XTpAAAAAElFTkSuQmCC) | 0.10 | 0.115 | 0.031 | 0.866 | 0.095 | 0.324 | 0.476 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAk1BMVEUAAAAAAABAQEBVVVVJSUlVVVVRUVFOTk5KSkpPT09LS0tOTk5NTU1NTU1LS0tMTExOTk5NTU1OTk5NTU1OTk5MTExOTk5NTU1NTU1NTU1NTU1OTk5OTk5NTU1MTExNTU1NTU1MTExNTU1MTExNTU1NTU1OTk5NTU1NTU1NTU1MTExOTk5NTU1MTExNTU1NTU1RUVFr/mCgAAAAMXRSTlMAAQQGBwkTFxgdIi41OD1ASFlveoCCg4SFh4iKjY+QkZKTn6Cio6SzxMvS2eHw8///rFV85wAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAKNJREFUOI1jYBgFAw5YWIhUKCxMpEIxUSIVSogTqVBKkoACblYILS0NoVm5cShU1BFkBNEyMiCSUVBHEZeR/FqafEBKThZI8GhoC+C2nElIV52LQVGBgVNFT4QZrzvZlQwkVVUlDZTY8SoDAV5NfX1NXoLKiFbIJm8gr6YsYqDCgVcZsZ4hNniIDnCioxAGCCYKGCA6mRGdcInOCkRnrlFAQwAAAuUMfHwG17YAAAAASUVORK5CYII=) | 0.02 | 0.023 | 0.005 | 0.866 | 0.019 | 0.145 | 0.244 |

A regular 6-point star with fixed outer radius and shrinking inner
radius goes from nearly solid to six needle-thin spikes - every index
here falls toward 0 as the notches deepen, except one:
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
barely moves (0.900 to 0.866). The six outer points never move, so the
minimum bounding rectangle never changes shape regardless of how thin
the arms between them get -
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
only ever sees the envelope, never how much of it is actually filled.
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
falls hardest of the rest, since deepening notches add a great deal of
perimeter for comparatively little area loss.

### 5.2 Boundary/fractal noise: `polsby_popper_index()` doesn’t converge

| shape | resolution | hull_ratio_index | polsby_popper_index | width_length_ratio_index | reock_index | detour_index | exchange_index |
|:---|:---|---:|---:|---:|---:|---:|---:|
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAQlBMVEUAAABAQEBJSUlOTk5SUlJOTk5PT09MTExMTExNTU1NTU1RUVFVVVVaWlptbW18fHx9fX2mpqazs7O9vb2+vr6/v78ectFoAAAAFnRSTlMABAcNGRodHn+9zv//////////////Q/srlgAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAJNJREFUOI3tk8kKwzAMRB1nk9d4if7/V1tKKa2R24EQesk7P0ZYHil18W/0uurf1kTWheAsTV+1gXzKlbnm5Gnoe6OJhZ+UaMZuntl2frFvppdJ8c27m5Fkb/aFPyh+lgMTNyQ50uZWzFbytKutWJ20+SW0HnNYjojwaPgx+HrghcNfiJcCrhleXAWfwgPsuC5O5ga4OhiRoFjq5wAAAABJRU5ErkJggg==) | true (unpixelated) | 1.000 | 1.000 | 1 | 1.000 | 1.000 | 1.000 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAM1BMVEUAAAAAAAB/f39QUFBNTU1LS0tMTExOTk5MTExOTk5NTU1SUlJcXFxycnJzc3OSkpK/v790uRlhAAAAEXRSTlMAAQIgISJKTlFVi////////yQEbmUAAAAJcEhZcwAADsMAAA7DAcdvqGQAAABiSURBVDiN7dO7DoAgDEBRXgKlIPz/14oEdZCQOggLZ2typ0IZW2bi2413Q42+Qt0NbUhVsONCYQDAPaHLoxGNUOKexSuM54SyESqfXrxa4deQvHDyExZzvllBPgXycS1/OwAAOROp4vK61wAAAABJRU5ErkJggg==) | cell = 4.00 | 0.857 | 0.589 | 1 | 0.764 | 0.899 | 0.910 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAS1BMVEUAAABRUVFOTk5JSUlNTU1LS0tMTExOTk5OTk5NTU1NTU1NTU1NTU1MTExPT09cXFxeXl5gYGCFhYWGhoaJiYm2tra3t7e7u7u/v7+P1xhoAAAAGXRSTlMAFhocISJrbHbe7O7v8P//////////////A9ZuRAAAAAlwSFlzAAAOwwAADsMBx2+oZAAAAIBJREFUOI3t08kKgCAUQFFTc2jQMlP//0sLgyjSepsIwrs+qzcgVPo2UsfIk6uUjqnqlrG2H0Js6FuWd0Ia5zfonZEi57i0c9ibreQZ2EwuHHJTk3a4O7lVdjjp1OjP0I8qJakOlzT9OwSPBzxw+ArhRwE+M/jhIvArIPhzlV5uAclcHPVbc9fLAAAAAElFTkSuQmCC) | cell = 1.00 | 0.945 | 0.638 | 1 | 0.895 | 0.961 | 0.977 |
| ![](data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAACgAAAAoCAMAAAC7IEhfAAAAYFBMVEUAAAAAAABVVVVJSUlKSkpQUFBNTU1LS0tNTU1PT09NTU1NTU1OTk5NTU1NTU1NTU1NTU1RUVFTU1NcXFxhYWFwcHBxcXF5eXmfn5+ioqKoqKipqam2tra3t7e+vr6/v78j7VQ0AAAAIHRSTlMAAQYHHyAhIistkZKXmNbb4P///////////////////6UQ/zwAAAAJcEhZcwAADsMAAA7DAcdvqGQAAACeSURBVDiN7dPZDsIgEEBR6EL3lk3ainT+/y+riTFKBp3EGF96n09CGAbGjv4cF1Ul+EdWD6M2Ro9D/ZZlrXJLALgsJ9VmaVdK6+He2coy5XI5b/Bom2WegJ19cldpO9w1ysNLXjUo7B1EuR5zfFpjuE7YPAsdYhh0gUBhYgdgxDeQfDT5MuTx0AdOfkL6UpDXjL64jPwVbhE/19GP2wGj/yRnkBKNXQAAAABJRU5ErkJggg==) | cell = 0.25 | 0.976 | 0.619 | 1 | 0.964 | 0.985 | 0.993 |

Rasterizing a smooth circle at finer and finer grid resolutions is a
stand-in for real vectorized raster/classification boundaries. Five of
the six indices here converge cleanly back toward 1 as the pixel grid
gets finer - they’re all area-based, and the stairstep noise this
introduces is a vanishing fraction of total area at fine resolution.
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
does not: 0.589 at the coarsest resolution shown, 0.619 at the finest -
stuck well below 1, not steadily improving. This is the classical
staircase phenomenon: a stairstep boundary’s total perimeter does not
converge to the smooth curve’s true length as the steps shrink, because
each step’s own edges never point in the curve’s actual tangent
direction, however small the steps get. Any perimeter-based measure
inherits this -
[`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
is the one index in this package’s whole collection genuinely vulnerable
to it.
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
stays at exactly 1 throughout, for an unrelated reason: a circle’s
minimum bounding rectangle is a square at any pixelation coarse enough
to still resemble the original extent.

### 5.3 Holes and multi-part shapes: `exchange_index()`’s two distinct failure modes

[`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)
already showed the worse of its two pathologies above: an exact collapse
to 0 once a multi-part shape’s pieces separate far enough that the
reference circle can’t reach either one
([`vignette("j-understanding-classical-indices")`](https://nkaza.github.io/shapeindices/articles/j-understanding-classical-indices.md)’s
own hole/multi-part tables above; Angel et al. note the identical
failure mode for real districts split by water). The second, shown in
the hole table above, is distinct in mechanism though similar in effect:
a large, roughly-centred hole can empty out exactly the region the
shrinking reference circle would occupy, driving the index to 0 there
too - not because the shape has separated into pieces, but because the
one piece it has has a void precisely where the index looks. Both are
documented, both are kept rather than omitted, matching this package’s
existing precedent with
[`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)’s
own multi-part/hole blind spots.

## 6 Key takeaways

- All six classical metrics are cheap and closed-form - area, perimeter,
  convex hull, or a bounding circle/box - with no CDT triangulation and
  no weighted form; a `weight` can only exclude a row as a hole before
  union, never scale an index’s own magnitude.
- [`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)
  and
  [`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
  respond to holes only through their numerator
  ([`hull_ratio_index()`](https://nkaza.github.io/shapeindices/reference/hull_ratio_index.md)
  tracks `1 - hole_frac` exactly;
  [`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
  doesn’t move at all, since neither the hull nor the minimum bounding
  rectangle ever reaches into a hole).
- [`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
  is the most exposed to boundary/fractal noise of any index in this
  package - a stairstepped (rasterized) boundary inflates its perimeter
  and keeps it inflated even at very fine pixel resolution, a real
  mathematical fact about approximating curves with axis-aligned steps,
  not a fixable artifact of coarse data.
- [`polsby_popper_index()`](https://nkaza.github.io/shapeindices/reference/polsby_popper_index.md)
  is also structurally blind to multi-part dispersal specifically: area
  and perimeter are simple sums over disjoint pieces, unaffected by how
  far apart those pieces sit, as long as they don’t touch.
- [`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
  is blind to notch/spike depth in the opposite sense - a shape can go
  from solid to needle-thin spiked without moving its own minimum
  bounding rectangle at all.
- All six are orientation-invariant -
  [`width_length_ratio_index()`](https://nkaza.github.io/shapeindices/reference/width_length_ratio_index.md)
  uses the minimum-area bounding rectangle at any rotation, not the
  axis-aligned bounding box, specifically so that rotating a shape in
  place, without changing it at all otherwise, can’t move the index.
- [`exchange_index()`](https://nkaza.github.io/shapeindices/reference/exchange_index.md)
  has two distinct failure modes, both ending in an exact 0 rather than
  merely a low score: multi-part separation past the reference circle’s
  reach, and a large hole occupying the same ground the reference circle
  would.
- [`reock_index()`](https://nkaza.github.io/shapeindices/reference/reock_index.md)
  degrades the most severely of the six under multi-part dispersal,
  since its own reference (the minimum bounding circle) must grow to
  enclose ever more distant pieces, ballooning in area quadratically
  with the separation.

[^1]: Frolov, Y. S. (1975). Measuring the shape of geographical
    phenomena: a history of the issue. *Soviet Geography*, 16(10),
    676-687.

[^2]: MacEachren, A. M. (1985). Compactness of geographic shape:
    comparison and evaluation of measures. *Geografiska Annaler: Series
    B, Human Geography*, 67(1), 53-67.

[^3]: Murray, A. T. (2025). Geographical Compactness in Shape
    Assessment. *Geographical Analysis*, 57, 88-113.

[^4]: Dori, D., & Ben-Bassat, M. (1983). Efficient nonpolynomial
    computation of the ratio of area to inertia of a shape. *IEEE
    Transactions on Pattern Analysis and Machine Intelligence*, 5(4),
    471-475.

[^5]: Murray, A. T. (2025). Geographical Compactness in Shape
    Assessment. *Geographical Analysis*, 57, 88-113.

[^6]: Polsby, D. D., & Popper, R. D. (1991). The third criterion:
    compactness as a procedural safeguard against partisan
    gerrymandering. *Yale Law and Policy Review*, 9(2), 301-353.

[^7]: Reock, E. C. (1961). A note: measuring compactness as a
    requirement of legislative apportionment. *Midwest Journal of
    Political Science*, 5(1), 70-74.

[^8]: Angel, S., Parent, J., & Civco, D. L. (2010). Ten compactness
    properties of circles: measuring shape in geography. *Canadian
    Geographer*, 54(4), 441-461.
