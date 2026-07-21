# Width-length ratio of a (multi)polygon's bounding box

The shorter of the bounding box's x/y extents over the longer (`length`
is the longer one, matching everyday usage, not the shorter one), in (0,
1\], 1 = a square bounding box. Axis-aligned, not the minimum bounding
rectangle at any rotation - a diagonally-oriented elongated shape can
score deceptively high, the classic limitation of this score.

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

KNOWN LIMITATION: this looks only at the bounding box, so it's blind to
both holes and multi-part dispersal in ways
`hull_ratio_index`/`reock_index` are not. A solid square and the same
square with a large hole punched through it score identically (a hole
never extends past the outer boundary, so it never moves the bbox).
Likewise, several small disconnected pieces sitting near the corners of
a roughly square overall extent can score close to 1 despite being
scattered, not compact - the bbox is taken over every part combined,
independent of how much of it is actually filled.

## Examples

``` r
nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
width_length_ratio_index(nc[nc$NAME == "Wake", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7844, lon_0 = -78.653) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.9006867
width_length_ratio_index(nc[nc$NAME == "Dare", ])$index
#> Input is in geographic (lon/lat) coordinates; auto-projecting to a local azimuthal-equal-area CRS centred on the data (lat_0 = 35.7353, lon_0 = -75.8095) before computing - pass already-projected data instead if you need a specific CRS.
#> [1] 0.4431556
```
