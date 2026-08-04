# 11. Comparing Shape Indices: North Carolina Counties

Code

``` r

library(shapeindices)
library(sf)
library(ggplot2)
library(ggalign)
library(dplyr)

theme_set(theme_minimal(base_size = 11))
```

## 1 Introduction

The other vignettes in this package derive each index in isolation, on
synthetic shapes built to isolate one property at a time. This one runs
some of the indices on a real dataset: the North Carolina county
boundaries bundled with `sf`) [^1] The vignette is primarily to show how
the indices perform on a real dataset with varying shapes and weights.

## 2 Indices point to stories about geography and history

County lines are surveyed straight lines and coastlines, not organic
growth boundaries, so the shapes here test something specific: how much
does an index notice a barrier island, a diagonal sliver, or a county
that’s convex but stretched?

[TABLE]

**Caswell** is the baseline: a simple, nearly-rectangular Piedmont
county with no natural features at the boundary. The county is created
from bisecting Orange county and then again to form Caswell adn Pearson.
Convexity, hull ratio, and width-length ratio are all at or near 1. Even
here, Reock only reaches 0.63 (a perfect square’s Reock score is exactly
$`2/\pi \approx 0.637`$). Truly, disk shaped counties are hard to come
by and Reock penalises any deviation from the circle more heavily than
other indices.

**Dare and Currituck** are Outer Banks counties where a compact mainland
piece connects to a long barrier-island strip across an open sound.
Almost index penalises this geographical arrangement: convexity
0.74/0.80, hull ratio 0.26/0.55, Polsby-Popper as low as 0.08 for Dare.

Width-length ratio agrees too, now: 0.40 for Dare and 0.43 for
Currituck, both low like every other elongation-sensitive measure here.
That agreement is itself worth a note: this index used to score
Currituck at 0.96 - nearly a perfect square - because its old
axis-aligned bounding box happened to line up close to the coordinate
axes, an accident of how the county is drawn on the map, not a property
of its actual shape. Switching to the minimum-area bounding rectangle at
any rotation
([`sf::st_minimum_rotated_rectangle()`](https://r-spatial.github.io/sf/reference/geos_unary.html))
fixes this: Currituck’s true elongation, independent of orientation, is
close to Dare’s own, not the 0.96 the old orientation-sensitive version
reported.

**Camden** is the classic split: convexity 0.99, essentially convex, no
real notches for a line to cross. Span/radial concentration agree
(0.61/0.61): Camden is a long, narrow diagonal sliver, and
convexity/span/radial-concentration agree broadly. On the other hand,
moment isotropy as well as moment of inertia isolates that elongation
directly and puts a sharp number on it: 0.04/0.32. Camden’s own
principal moments differ by roughly 27-to-1.

**Lincoln** sharpens that same point with no boundary complexity at all
involved: convexity is a perfect 1.00 and hull ratio 0.99. Yet Reock is
0.34, width-length ratio 0.30, and moment isotropy 0.11; all three
elongation-sensitive measures still penalise it for being stretched,
each via a different route (bounding circle, bounding box, mass tensor).
**Convex and compact are different properties**, and a shape can score
perfectly on one while scoring poorly on the other.

## 3 Rescaling convexity for one corpus

[`vignette("c-understanding-convexity-index")`](https://nkaza.github.io/shapeindices/articles/c-understanding-convexity-index.md)’s
“Why the index compresses toward convex” section shows that CI stays
near 1 for any shape that keeps a substantial convex core, which is most
real polygons. This dataset is a clean demonstration of exactly that:

``` r

summary(nc_df$convexity_index)
```

       Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
     0.7375  0.9918  0.9970  0.9854  0.9992  1.0000 

Raw convexity never drops below 0.74 (Dare, the barrier-island county
above) across all 100 counties, and 96% of them score above 0.9. Median
convexity is 0.9970, and even the 90th percentile only reaches 0.99993 -
practically a perfect square. That’s correct: these are surveyed county
boundaries, not hand-drawn stars, and almost none of them lose their
convex core the way it takes to pull CI down (see the linked section
above). But it makes raw CI a poor scale for a task like a choropleth
map or a ranked list, where 96 counties reading as “about 0.99” is not a
useful map even though it’s an accurate one.

### 3.1 Percentile rescaling

The fix isn’t to alter `convexity_index` itself - a value of 0.9992
versus 0.9999 is a real, if tiny, geometric difference, and rescaling
the definition to manufacture spread would misrepresent it (the same
reasoning
[`vignette("c-understanding-convexity-index")`](https://nkaza.github.io/shapeindices/articles/c-understanding-convexity-index.md)
gives for not clamping or rescaling CI’s compressed region directly).
Instead, compute a second, corpus-relative statistic on top of it: each
county’s percentile rank within this specific set of 100.

``` r

nc_df$convexity_percentile <- rank(nc_df$convexity_index, ties.method = "average") / nrow(nc_df)
```

That spreads the same 100 counties uniformly across `[0, 1]` by
construction - one unit of the rescaled value always means “10 more
counties,” regardless of how bunched the raw values are:

Code

``` r

rescale_counties <- c("Dare", "Currituck", "Camden", "Lincoln", "Caswell")
knitr::kable(
    nc_df[match(rescale_counties, nc_df$NAME), c("NAME", "convexity_index", "convexity_percentile")],
    digits = c(NA, 6, 3), col.names = c("County", "Raw convexity", "Percentile (this corpus)"),
    row.names = FALSE
)
```

| County    | Raw convexity | Percentile (this corpus) |
|:----------|--------------:|-------------------------:|
| Dare      |      0.737548 |                     0.01 |
| Currituck |      0.804788 |                     0.02 |
| Camden    |      0.987415 |                     0.21 |
| Lincoln   |      0.999996 |                     0.91 |
| Caswell   |      1.000000 |                     0.96 |

Dare and Currituck, at raw 0.738 and 0.805, land at the 1th and 2th
percentiles - clearly the least convex counties in the state, which the
raw scores already showed well enough on their own. The more useful case
is the middle of the distribution: Camden, at a raw 0.987 that reads as
“basically convex” next to Lincoln’s and Caswell’s near-perfect scores,
actually sits at only the 21th percentile - noticeably less convex than
most of the state, a distinction raw CI compresses almost out of view.
Lincoln and Caswell themselves make the same point at the opposite
extreme: both round to 1.000 at three decimal places, yet sit five
percentile points apart, because so much of the state is packed into
that last ten-thousandth of raw score.

### 3.2 When to use which

Keep raw `convexity_index` for anything that needs to be comparable
outside this dataset: a fixed threshold, a cross-study comparison, or
any claim about the shape itself (“this county is 99.9% as convex as
physically possible”). Reach for the percentile version only for
within-corpus tasks where the audience is these 100 counties
specifically - a choropleth that needs to use its full colour ramp, or a
ranked table of “most/least convex NC counties.” Percentile rescaling is
corpus-relative by definition: add a state with more fragmented coastal
counties to the dataset and every percentile shifts, even though no
county’s own geometry changed.

## 4 How do the indices relate to each other?

![](k-nc-counties-comparison_files/figure-html/nc-corr-1.png)

The above correlation matrix uses unweighted indices and hierarchical
clustering. Across these 100 counties, the indices condense into three
primary behaviors dominated by whether they respond to elongation or
boundary deformities. The core **Elongation Cluster**—comprising moment
of inertia, span, radial concentration, and Reock—behaves almost
interchangeably, as county-level elongation swamps the theoretical
differences between their formulas. Other measures fold into this group
via distinct mechanics: `exchange_index` reads elongation via area
overlap (0.96 with MOI), `detour_index` captures hull-level stretching
(0.94 with span), `depth_index` aligns with boundary complexity (0.93
with Polsby-Popper), and `width_length_ratio_index` - since being
corrected to use the minimum-area bounding rectangle at any rotation
rather than an axis-aligned one - now tracks the same elongation cleanly
(0.80 with moment isotropy, its own strongest tie in the whole matrix).
Moment isotropy belongs to this family but sits slightly apart (0.74
with MOI) because it isolates pure axis elongation, remaining
structurally blind to dispersal or boundary noise.

**Convexity and hull ratio** form a second, distinct cluster (0.84) that
isolates concave boundary departures while ignoring elongation entirely.
The near-zero correlation between moment isotropy and convexity (0.24)
provides clean empirical proof that elongation and concavity are
independent structural properties on real geography rather than
redundant descriptions. Similarly, directional balance tracks
first-harmonic bearing preference rather than axis shape; its weak tie
to moment isotropy (0.46) confirms that bearing preference and
mass-distribution symmetry operate as independent axes despite both
belonging to the elongation-adjacent family.

Connecting these groups, Polsby-Popper acts as a hybrid bridge between
the convexity and elongation clusters (0.90 with hull ratio; 0.84 with
MOI) because perimeter length is inherently inflated by both boundary
roughness and stretched geometries.

`width_length_ratio_index` is no longer the outlier it once was. An
earlier, axis-aligned version of this index correlated weakly with
almost everything in this dataset, because its own orientation noise -
which way a county happens to be drawn relative to the coordinate axes,
not a property of its shape - swamped the real elongation signal
underneath. With that noise removed, its strongest ties sit squarely
inside the elongation cluster: 0.80 with moment isotropy, 0.78 with
exchange, 0.76 with moment of inertia. It still correlates more loosely
with the boundary-complexity side of things (0.31 with convexity, 0.31
with hull ratio) - a bounding rectangle, however it’s computed, still
can’t see concavity the way a convex-hull-based measure can - but it is
a member of a real cluster now, not an island of its own.

## 5 Weighting a collection of polygons

`byrow = FALSE` treats every row instead as a weighted sub-polygon of
**one** overall shape (`st_union(x)`), returning a single-row result.
`weights = NULL` (the default) weights each row by its own physical
area, which reproduces the plain unweighted indices on the union; a
column name or numeric vector weights rows by something else instead -
here, births (`BIR74`) as a population proxy for four contiguous
Research Triangle counties:

|                            | weighted_by_area | weighted_by_births |
|:---------------------------|-----------------:|-------------------:|
| convexity_index            |             0.98 |               0.99 |
| moment_of_inertia_index    |             0.82 |               0.70 |
| moment_isotropy_index      |             0.51 |               0.64 |
| directional_balance_index  |             0.98 |               0.96 |
| span_index                 |             0.92 |               0.85 |
| radial_concentration_index |             0.92 |               0.85 |
| depth_index                |             0.69 |               0.51 |
| hull_ratio_index           |             0.83 |               0.83 |
| polsby_popper_index        |             0.51 |               0.51 |
| width_length_ratio_index   |             0.69 |               0.69 |
| reock_index                |             0.48 |               0.48 |
| detour_index               |             0.83 |               0.83 |
| exchange_index             |             0.81 |               0.81 |
| total_weight               |    5811314465.44 |           27264.00 |

Same union geometry results in the same `hull_ratio_index` (it has no
weighted form - it’s a property of the convex hull’s boundary, not a sum
over pieces) - but `convexity_index`, `moment_of_inertia_index`, and
`span_index` all shift once births, rather than raw land area, decide
how much each county’s shape “counts”.

Weighted by each county’s own area (`weights = NULL`, the default), the
four-county union scores convexity 0.983 and MOI 0.815. Switching the
weight to 1974 births, so that Wake dominates in one part of the union,
moves MOI, span, and radial concentration all *down* together (0.697,
0.849, 0.852). Convexity moves the other way, slightly *up* (0.990). It
does not depend on where the centre of mass is at all and therefore
behaves as expected.

Reweighting the mesh reveals that moment isotropy and directional
balance capture fundamentally different spatial mechanisms when mass
shifts toward a dominant county like Wake. **Moment isotropy moves up**
(0.507 to 0.643), shifting opposite to MOI, span, and radial
concentration. While center-relative metrics penalize how far mass
drifts from the union’s geometric middle, moment isotropy evaluates only
the *shape* of the distribution around its own center—concentrating
weight into Wake (a relatively round county within a long four-county
strip) makes the core mass footprint more isotropic even as its center
drifts. Conversely, **directional balance moves down** (0.984 to 0.960),
tracking alongside MOI and span because shifting weight toward Wake
creates an asymmetric angular pull from the new centroid. Moment
isotropy acts as the lone outlier here precisely because it ignores
center-relative displacement, proving it operates as an independent axis
rather than a redundant metric. *(Note: `total_weight` reflects raw area
in* $`m^2`$ versus total births, so row totals represent different
scales).

## 6 Conclusion

No single shape index here is redundant: real-world geographies split
metrics into distinct structural behaviors rather than a simple
continuum. While Moment of Inertia (MOI), span, and radial concentration
cluster tightly around overall elongation, measures like convexity and
hull ratio track a fundamentally different property—concavity—explaining
why a shape can be nearly convex (e.g., Lincoln) while scoring poorly on
circle-referenced compactness. Moment isotropy sits loosely in the
elongation group but acts as an independent axis: it is virtually immune
to boundary concavity and can shift in the opposite direction from its
cluster under population reweighting because it ignores centroid offset
entirely, evaluating only the symmetry of the mass footprint itself.

Directional balance is about angular bearing preference (whether mass
leans toward one direction) rather than axis preference (whether mass
aligns along a line), yielding a remarkably weak correlation with moment
isotropy (0.46). Under mesh reweighting, directional balance aligns with
MOI and span because all three respond to centroid-relative mass
drift—making it uniquely sensitive to split geometries like Dare and
Currituck’s mainland-and-barrier-island shapes. Ultimately,
`width_length_ratio_index` - once corrected to use the minimum-area
bounding rectangle at any rotation rather than an axis-aligned one -
joins the elongation cluster rather than sitting apart from it, and
population-weighting doesn’t just tweak scores across the mesh-based
indices—it fundamentally alters which spatial property each one reveals.

[^1]: It is probably useful to think about the history of how the shapes
    came to be and how they changed over time: see Kelly, S. R. 2015.
    “The Boundary Hunters Uncovering North Carolina’s Lost Borders.”
    Technology. The Atlantic, October 6.
    https://www.theatlantic.com/technology/archive/2015/10/north-carolina-lost-county-lines/409090/.
