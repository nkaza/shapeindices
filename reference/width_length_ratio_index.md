# Width-length ratio of a (multi)polygon's minimum-area bounding rectangle

The shorter of the bounding rectangle's two side lengths over the longer
(`length` is the longer one, matching everyday usage, not the shorter
one), in (0, 1\], 1 = a square bounding rectangle. Uses the MINIMUM-AREA
rectangle at any rotation
([`sf::st_minimum_rotated_rectangle()`](https://r-spatial.github.io/sf/reference/geos_unary.html),
GEOS-backed - no new package dependency), not the axis-aligned bounding
box: an earlier version of this function used the axis-aligned box,
which meant the same shape could score anywhere from its true ratio up
to a spurious 1 depending purely on which way it happened to be drawn
relative to the coordinate axes - orientation is not a property of the
shape itself, so an index built on it shouldn't depend on it either.
Verified directly: a fixed 2:1 rectangle rotated from 0 to 90 degrees
now returns the same 0.5 throughout, rather than swinging up to 1.0 at
45 degrees the way the axis-aligned version did.

## Usage

``` r
width_length_ratio_index(poly)
```

## Arguments

- poly:

  a single sfg/sfc (MULTI)POLYGON

## Value

list(index, length, width)

## Details

KNOWN LIMITATION: this looks only at the bounding rectangle, so it's
blind to both holes and multi-part dispersal in ways
`hull_ratio_index`/`reock_index` are not. A solid square and the same
square with a large hole punched through it score identically (a hole
never extends past the outer boundary, so it never moves the bounding
rectangle). Likewise, several small disconnected pieces sitting near the
corners of a roughly square overall extent can score close to 1 despite
being scattered, not compact - the bounding rectangle is taken over
every part combined, independent of how much of it is actually filled.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
width_length_ratio_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9595009
width_length_ratio_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.4047457
```
