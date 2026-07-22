# 1. Basic Usage

Code

``` r

library(shapeindices)
library(sf)
library(ggplot2)
library(dplyr)

theme_set(theme_minimal(base_size = 11))
```

This vignette is a practical, code-first tour of the package
`shapeindices`. The goal of the package is to provide indices for `sf`,
`sfc` objects that describe the shape/mass configurations. In this
package, I demonstrate some key indices and approaches. For deep dive
into specific index, see the corresponding vignettes.

## 1 Three basic shapes

Three small synthetic polygons, built with base `sf` calls, cover a
useful range of boundary behaviour: a **star** (deep, evenly-spaced
notches), a **spiral** (a single winding corridor, no notches at all),
and a **blob with a hole** (an irregular ring with an interior void).

Code

``` r

make_star <- function(n_points, r_outer = 1, r_inner = 0.5, center = c(0, 0)) {
  n <- n_points * 2
  angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n + 1)[1:n]
  radii  <- rep(c(r_outer, r_inner), n_points)
  x <- center[1] + radii * cos(angles)
  y <- center[2] + radii * sin(angles)
  coords <- rbind(cbind(x, y), c(x[1], y[1]))
  st_polygon(list(coords))
}

make_turtle_path <- function(n_steps, angle_deg, step0 = 1, step_growth = 0) {
  angle <- 0
  pos   <- c(0, 0)
  pts   <- matrix(pos, ncol = 2)
  step  <- step0
  for (i in seq_len(n_steps)) {
    angle <- angle + angle_deg * pi / 180
    pos   <- pos + step * c(cos(angle), sin(angle))
    pts   <- rbind(pts, pos)
    step  <- step + step_growth
  }
  st_linestring(pts)
}

make_spiral <- function(n_steps = 48, angle_deg = 24, step0 = 0.35,
                         step_growth = 0.09, width = 0.35) {
  st_buffer(make_turtle_path(n_steps, angle_deg, step0, step_growth),
            dist = width, endCapStyle = "FLAT", joinStyle = "MITRE", mitreLimit = 3)
}

make_blob_hole <- function(n = 14, seed = 1, roughness = 0.55, hole_frac = 0.3) {
  set.seed(seed)
  angles <- sort(runif(n, 0, 2 * pi))
  radii  <- 1 + roughness * (runif(n) - 0.5) * 2
  outer  <- cbind(radii * cos(angles), radii * sin(angles))
  outer  <- rbind(outer, outer[1, ])

  # the hole is the SAME irregular outline scaled down by sqrt(hole_frac)
  # (so hole area / outer area ~= hole_frac) and traced in reverse, so sf
  # treats it as an interior ring rather than a second shell. Scaling every
  # radius by the same constant keeps the hole strictly inside the outer
  # boundary at every angle, so the result is always a valid polygon.
  hole <- outer[nrow(outer):1, ] * sqrt(hole_frac)

  st_polygon(list(outer, hole))
}

star      <- make_star(6, r_outer = 1, r_inner = 0.4)
spiral    <- make_spiral()
blob_hole <- make_blob_hole()

basic_shapes <- list(star = star, spiral = spiral, "blob with hole" = blob_hole)

# geom_sf() + facet_wrap() can't use free scales, so normalise each shape
# onto a common bounding box (center it, scale to unit span) before faceting
normalize_geom <- function(geom, center, scale) (geom - center) * scale
shapes_norm <- lapply(names(basic_shapes), function(nm) {
  g      <- st_sfc(basic_shapes[[nm]])
  bb     <- st_bbox(g)
  center <- unname(c((bb["xmin"] + bb["xmax"]) / 2, (bb["ymin"] + bb["ymax"]) / 2))
  span   <- max(bb["xmax"] - bb["xmin"], bb["ymax"] - bb["ymin"])
  st_sf(shape = nm, geometry = normalize_geom(g, center, 1 / span))
})
shapes_sf <- do.call(rbind, shapes_norm)

ggplot(shapes_sf) +
  geom_sf(fill = "steelblue", alpha = 0.6, color = "grey20") +
  facet_wrap(~ shape) +
  theme_void(base_size = 11) +
  theme(strip.text = element_text(face = "bold"))
```

![](a-basic-usage_files/figure-html/shapes-1.png)

## 2 The thirteen indices

Each index has its own function `*_index()` that takes, sf or sfc
objects as input. The following table briefly describes the conceptual
framing of the index and shows, as a thumbnail, a shape/mass
distribution that achieves close to the maximum value. These are just
some examples and the emphasis is on showing some unexpected profiles
that score high on these indices. While it is tempting to think that
there is an unique shape and a mass profile configuration
(e.g. center-weighted disk) that maximises the index, the following
table shows the folly of such naivete.

[TABLE]

The three thumbnails on `width_length_ratio_index`’s own row aren’t a
*good* example, but a caveat. A hole punched clean through a shape,
needle-thin spikes, guitar picks, dendritic shapes with square convex
hull and a shape split into several pieces all still score close to `1`
here, because the index only ever looks at the minimum-area bounding
rectangle. See
[`vignette("j-understanding-classical-indices", package = "shapeindices")`](https://nkaza.github.io/shapeindices/articles/j-understanding-classical-indices.md)
for the other blind spots this family of indices shares.

Called directly on a shape, for example,:

``` r

convexity_index(star)$index
```

    [1] 0.9444644

``` r

moment_of_inertia_index(spiral)$index
```

    [1] 0.220872

``` r

hull_ratio_index(blob_hole)$index
```

    [1] 0.5117867

Running the seven mesh-based indices separately re-triangulates the
polygon up to seven times.
[`shape_indices()`](https://nkaza.github.io/shapeindices/reference/shape_indices.md)
triangulates once and returns all thirteen as a named vector, which
might have implications for speed.

``` r

round(sapply(basic_shapes, shape_indices), digits = 3)
```

                          star spiral blob with hole
    convexity            0.944  0.300          0.695
    moment_of_inertia    0.761  0.221          0.492
    moment_isotropy      1.000  0.828          0.770
    directional_balance  1.000  0.997          0.947
    span                 0.882  0.470          0.705
    radial_concentration 0.900  0.467          0.684
    depth                0.512  0.102          0.318
    hull_ratio           0.462  0.270          0.512
    polsby_popper        0.224  0.018          0.129
    width_length_ratio   0.866  0.922          0.901
    reock                0.382  0.231          0.343
    detour               0.647  0.511          0.657
    exchange             0.772  0.182          0.540

The spiral scores low on every index (a long, narrow, winding corridor
is far from convex and far from disk-like); the blob’s hole drags its
indices down relative to a solid blob of the same outer boundary; the
star sits in between, close to its `hull_ratio_index` twin but lower on
`convexity_index`, which - unlike `hull_ratio_index` - is sensitive to
the *evenly-spaced* notches cut into it.

`polsby_popper_index` in particular penalises the blob’s hole twice
over, not once: `area` is already net of the hole (every index here
measures net area), but `perimeter` grows too, since
[`st_boundary()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
returns *every* ring - the hole’s own boundary, not just the outer one -
and
[`st_length()`](https://r-spatial.github.io/sf/reference/geos_measures.html)
sums across all of them. A shape with a hole cut out of it genuinely is
less compact than the solid version, so both terms moving the same
direction is the right behaviour here - but it’s not a universal
convention. Some Polsby-Popper implementations in the redistricting
literature measure only the *outer* ring’s length, treating interior
holes (a small enclave, a data-digitisation artifact) as
compactness-irrelevant. This package prizes consistency of definition
and leaves this data cleaning exercise to the user.

## 3 Deterministic vs. random estimation

- [`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.md),
  [`span_index()`](https://nkaza.github.io/shapeindices/reference/span_index.md),
  [`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md),
  [`directional_balance_index()`](https://nkaza.github.io/shapeindices/reference/directional_balance_index.md),
  and
  [`depth_index()`](https://nkaza.github.io/shapeindices/reference/depth_index.md)
  all have a `deterministic`
  argument.[`moment_of_inertia_index()`](https://nkaza.github.io/shapeindices/reference/moment_of_inertia_index.md),
- [`moment_isotropy_index()`](https://nkaza.github.io/shapeindices/reference/moment_isotropy_index.md),
  and the six classic metrics don’t have this argument. They are always
  computed directly off the CDT mesh / convex hull / bounding circle.

`deterministic = TRUE` (the default for all five) computes the index
over a **fixed grid** on the CDT mesh - `O(n^2)` in piece count for
convexity/span, `O(n)` for
radial_concentration/directional_balance/depth but with each triangle’s
own subdivision depth (up to 256 points) scaled to its area, so only the
mesh’s largest triangles pay the full cost. It’s called `deterministic`,
not `exact`, deliberately: for shapes with small concavities relative to
triangle size, or (for radial_concentration/directional_balance/depth)
an integral with no closed form at all, it’s still only an approximation
of the true value, just a non-random one (see
[`vignette("c-understanding-convexity-index")`](https://nkaza.github.io/shapeindices/articles/c-understanding-convexity-index.md)’s
convexity section for why).

`deterministic = FALSE` instead computes a **Monte Carlo estimate**: all
five functions draw `n_lines` independent random samples and average
some property of them: -
[`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.md)
the fraction of a connecting line’s length that falls outside the
polygon, -
[`span_index()`](https://nkaza.github.io/shapeindices/reference/span_index.md)
the line’s length itself -
[`radial_concentration_index()`](https://nkaza.github.io/shapeindices/reference/radial_concentration_index.md)
the distance from each sampled point to the shape’s own geometric
median -
[`directional_balance_index()`](https://nkaza.github.io/shapeindices/reference/directional_balance_index.md)
the bearing of each sampled point from the shape’s own centroid, -
[`depth_index()`](https://nkaza.github.io/shapeindices/reference/depth_index.md)
the distance from each sampled point to the shape’s own boundary.

Only
[`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.md)
ever actually builds line geometry. Sharing the argument name means one
`n_lines` the sample count at once for all the five functions.

Code

``` r

library(purrr)
library(dplyr)
library(knitr)

# Non-deterministic parameters helper
rand_opts <- list(deterministic = FALSE, n_lines = 3000, seed = 1)

# Helper function to compute metrics for a single shape 
compute_shape_metrics <- function(nm, g) {
  ci_det   <- convexity_index(g, deterministic = TRUE)$index
  ci_rand  <- do.call(convexity_index, c(list(g), rand_opts))$index
  
  span_det  <- span_index(g, deterministic = TRUE)$index
  span_rand <- do.call(span_index, c(list(g), rand_opts))$index
  
  rci_det  <- radial_concentration_index(g, deterministic = TRUE)$index
  rci_rand <- do.call(radial_concentration_index, c(list(g), rand_opts))$index
  
  db_det   <- directional_balance_index(g, deterministic = TRUE)$index
  db_rand  <- do.call(directional_balance_index, c(list(g), rand_opts))$index
  
  depth_det  <- depth_index(g, deterministic = TRUE)$index
  depth_rand <- do.call(depth_index, c(list(g), rand_opts))$index

  tibble::tibble(
    shape               = nm,
    ci_deterministic    = ci_det,
    ci_random_line      = ci_rand,
    span_deterministic  = span_det,
    span_random_pair    = span_rand,
    rci_deterministic   = rci_det,
    rci_random_point    = rci_rand,
    db_deterministic    = db_det,
    db_random_point     = db_rand,
    depth_deterministic = depth_det,
    depth_random_point  = depth_rand,
    moi                 = moment_of_inertia_index(g)$index,
    moment_isotropy     = moment_isotropy_index(g)$index,
    hull_ratio          = hull_ratio_index(g)$index
  )
}

# 1. Loop over shape names to avoid type confusion with `basic_shapes[[nm]]`
det_compare <- map_dfr(names(basic_shapes), function(nm) {
  compute_shape_metrics(nm, basic_shapes[[nm]])
})

# 2. Extract matrix cleanly
det_mat <- as.matrix(det_compare[, -1])
rownames(det_mat) <- det_compare$shape

# 3. Create thumbnail headers safely
headers <- map_chr(names(basic_shapes), function(nm) {
  paste0(shape_thumb(basic_shapes[[nm]]), nm)
})

# 4. Render Table
kable(
  t(det_mat),
  format = "html", 
  digits = 2, 
  row.names = TRUE,
  col.names = headers, 
  escape = FALSE
)
```

[TABLE]

The Monte Carlo estimate is noisy, but close to the deterministic one,
and that noise shrinks as `n_lines` grows.

[`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.md)
always works from the polygon’s own CDT triangle mesh - `n_quad`
controls the quadrature refinement on that mesh (see above), but there’s
no alternative tessellation to choose. Two standalone mesh utilities,
[`convex_decompose()`](https://nkaza.github.io/shapeindices/reference/convex_decompose.md)
and
[`subdivide_mesh()`](https://nkaza.github.io/shapeindices/reference/subdivide_mesh.md),
are also available for downstream geometry work (neither is wired into
any of the thirteen indices) - see
[`vignette("b-understanding-triangulation")`](https://nkaza.github.io/shapeindices/articles/b-understanding-triangulation.md)
for both, alongside the CDT mesh they both start from.

## 4 Applying it to an `sf` object, row by row

Real polygons usually arrive as rows of an `sf` object.
[`shape_indices_sf()`](https://nkaza.github.io/shapeindices/reference/shape_indices_sf.md)
with `byrow = TRUE` (the default) runs
[`shape_indices()`](https://nkaza.github.io/shapeindices/reference/shape_indices.md)
on every row independently and appends four columns. Two North Carolina
counties from the `sf` package - Wake (a simple, nearly-convex Piedmont
county) and Dare (a mainland piece plus a separated barrier-island
strip) - show the range:

``` r

nc <- st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE) %>%
  st_transform(32119)   # NC state plane (meters) - nc.shp ships in NAD27 lon/lat

pair     <- nc[nc$NAME %in% c("Wake", "Dare"), ]
pair_res <- shape_indices_sf(pair)

pair_res %>%
  st_drop_geometry() %>%
  select(NAME, convexity_index, moment_of_inertia_index, span_index, hull_ratio_index) %>%
  knitr::kable(format = "html", digits = 3, row.names = FALSE)
```

| NAME | convexity_index | moment_of_inertia_index | span_index | hull_ratio_index |
|:-----|----------------:|------------------------:|-----------:|-----------------:|
| Wake |           0.995 |                   0.856 |      0.939 |            0.905 |
| Dare |           0.738 |                   0.228 |      0.545 |            0.265 |

``` r

ggplot(pair_res) +
  geom_sf(aes(fill = convexity_index), color = "grey20") +
  scale_fill_viridis_c(limits = c(0, 1)) +
  labs(fill = "CI")
```

![](a-basic-usage_files/figure-html/nc-byrow-plot-1.png)

[`shape_indices_sf()`](https://nkaza.github.io/shapeindices/reference/shape_indices_sf.md)
also supports `byrow = FALSE`, treating a whole collection of rows as
sub-polygons of one weighted overall shape rather than scoring each row
independently - and weighting them changes the CDT mesh itself, not just
the resulting index values, down to a design decision that treats a row
weighted `0` or `NA` as an excluded hole rather than a zero-weighted
piece. See
[`vignette("b-understanding-triangulation")`](https://nkaza.github.io/shapeindices/articles/b-understanding-triangulation.md)
for both mechanics, and
[`vignette("k-nc-counties-comparison")`](https://nkaza.github.io/shapeindices/articles/k-nc-counties-comparison.md)’s
“Weighting a collection of polygons” for how the weighting choice moves
the index values themselves.

## 5 Parallelisation and efficiency

This package has one parallelism knob, at the row level:
[`shape_indices_sf()`](https://nkaza.github.io/shapeindices/reference/shape_indices_sf.md)’s
`parallel_rows` argument. There is no parallelism *within* a single
index call - a single very complex polygon can’t be split across cores,
only a collection of rows can.

With `byrow = TRUE`, `parallel_rows = TRUE` distributes whole rows
across the active
[`future::plan()`](https://future.futureverse.org/reference/plan.html)
via `furrr` - this dispatch is unconditional (no fast-path
short-circuit), so it delivers real multi-core speedups. Row-level
parallelism is only worth turning on once each row’s own cost is
bounded - pass `deterministic_max_tri` (see
[`shape_indices()`](https://nkaza.github.io/shapeindices/reference/shape_indices.md))
so no single county can blow up into an expensive deterministic-mode
computation (`O(n^2)` for convexity/span, `O(n)` with a large constant
for radial_concentration) inside a worker. Raising
`deterministic_max_tri` and `n_lines` below (relative to the row-by-row
example earlier) gives every county enough work to make the per-worker
overhead worth paying:

``` r

library(furrr)

future::plan(future::sequential)
t_seq_rows <- system.time(
  res_seq <- shape_indices_sf(nc, deterministic_max_tri = 60, n_lines = 6000, seed = 1,
                               parallel_rows = FALSE)
)

future::plan(future::multisession, workers = 4)
t_par_rows <- system.time(
  res_par_rows <- shape_indices_sf(nc, deterministic_max_tri = 60, n_lines = 6000, seed = 1,
                                    parallel_rows = TRUE)
)
future::plan(future::sequential)

timing <- data.frame(
  mode    = c("parallel_rows = FALSE", "parallel_rows = TRUE (4 workers)"),
  elapsed = c(t_seq_rows[["elapsed"]], t_par_rows[["elapsed"]])
)
timing$speedup <- timing$elapsed[1] / timing$elapsed
knitr::kable(format = "html", timing, digits = 2, row.names = FALSE)
```

| mode                             | elapsed | speedup |
|:---------------------------------|--------:|--------:|
| parallel_rows = FALSE            |   24.81 |    1.00 |
| parallel_rows = TRUE (4 workers) |   12.49 |    1.99 |

Both modes agree on every county’s indices (not shown - parallelising
changes *how* the 100 rows get computed, not the values themselves); the
gain is in `elapsed`. Overhead from spinning up workers is fixed per
call, so it matters less as the batch gets bigger or each row gets more
expensive - not worth turning on for a handful of small, cheap polygons,
but pays off exactly on the kind of workload this package targets: many
real-world rows, each potentially complex enough to need
`deterministic_max_tri`’s random-line fallback. For a real workload,
size `workers` to
[`future::availableCores()`](https://parallelly.futureverse.org/reference/availableCores.html)
(or fewer, to leave headroom) and prefer
[`future::multicore`](https://future.futureverse.org/reference/multicore.html)
over
[`future::multisession`](https://future.futureverse.org/reference/multisession.html)
on Mac/Linux, which forks rather than starting fresh R sessions per
worker.

## 6 `simplify_tolerance`: reducing boundary detail, not just mesh size

[`prepare_polygon()`](https://nkaza.github.io/shapeindices/reference/prepare_polygon.md),
[`shape_indices()`](https://nkaza.github.io/shapeindices/reference/shape_indices.md),
and
[`shape_indices_sf()`](https://nkaza.github.io/shapeindices/reference/shape_indices_sf.md)
all accept a `simplify_tolerance` argument: an
[`st_simplify()`](https://r-spatial.github.io/sf/reference/geos_unary.html)
pass, applied once the polygon is already cleaned and (for
`byrow = FALSE`) already unioned, before triangulating. Its units are
whatever the geometry’s own CRS uses at that point - metres for
geographic (lon/lat) input, since these functions auto-project that to a
metric CRS themselves, but whatever unit an *already*-projected input’s
CRS uses otherwise (many US State Plane CRSs are US survey feet, not
metres) - a warning fires if that unit isn’t metres, since the same
numeric tolerance means a very different amount of simplification
depending on which.

It’s worth reaching for on genuinely detailed real-world boundaries
(block/parcel data with hundreds of thousands of vertices), for a reason
that isn’t obvious from mesh size alone:
[`convexity_index()`](https://nkaza.github.io/shapeindices/reference/convexity_index.md)’s
`deterministic = FALSE` mode tests candidate lines against *every*
boundary edge of the polygon, and that edge count doesn’t shrink just
because the triangulated mesh does - simplifying the boundary itself is
what actually bounds it.

**With `byrow = FALSE`, simplify the union, never the rows beforehand.**
`shape_indices_sf(byrow = FALSE)` applies `simplify_tolerance` to
`st_union(x)` itself, after every row’s own boundary has already merged
with its neighbours - internal edges shared between adjacent rows are
already gone by that point, so there’s nothing there for simplification
to misalign. Simplifying each row independently *first* - a
natural-looking alternative, and the wrong one - simplifies the same
shared edge differently on either side of it, turning what used to be a
clean tessellation into thousands of sliver gaps between rows that no
longer quite touch. As a concrete example, on one real 50,424-row
Census-block dataset: unioning the raw rows gives one clean,
~50,000-vertex polygon, while independently simplifying every row first
and *then* unioning gives 1,851 disjoint sliver fragments instead -
*more* total boundary than not simplifying at all, which at that scale
is the difference between a computation that fits in memory and one that
doesn’t. Passing `simplify_tolerance` here instead keeps the union at
one part and cuts its boundary edge count by an order of magnitude.
Weight accuracy is unaffected by design: the row-level weight overlay
always uses the *unsimplified* original rows - only the mesh comes from
the simplified union - so a tolerance modest relative to the shape’s own
scale still captures upward of 99% of the true row-level weight total.
