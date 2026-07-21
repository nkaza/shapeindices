test_that(".ellip_E matches numerical integration of its own definition", {
  for (k in c(0, 0.3, 0.5, 0.7, 1 / sqrt(2), 0.9, 0.999, 1)) {
    ref <- stats::integrate(function(phi) sqrt(1 - k^2 * sin(phi)^2), 0, pi / 2)$value
    expect_equal(shapeindices:::.ellip_E(k), ref, tolerance = 1e-6)
  }
})

test_that(".gauss_legendre is exact for polynomials up to degree 2m - 1", {
  for (m in 2:8) {
    gl <- shapeindices:::.gauss_legendre(m)
    got   <- sum(gl$w * gl$x^(2 * m - 2))
    exact <- 2 / (2 * m - 1)
    expect_equal(got, exact, tolerance = 1e-10)
  }
})

test_that(".mean_chord matches known closed forms", {
  # same circle (r1 = r2 = r): known "circle line picking" mean chord, 4r/pi
  expect_equal(shapeindices:::.mean_chord(2.5, 2.5), 4 * 2.5 / pi, tolerance = 1e-8)
  # point at the centre (r2 = 0): distance to any point on the r1-circle is r1
  expect_equal(shapeindices:::.mean_chord(3.7, 0), 3.7, tolerance = 1e-8)
})

test_that("span_index is ~1 for a disk", {
  disk <- make_disk(5, n_seg = 60)
  expect_equal(span_index(wc(disk))$index, 1, tolerance = 0.01)
})

test_that("span_index falls as a shape elongates", {
  square_span <- span_index(wc(make_square()))$index
  rect_span   <- span_index(wc(make_rectangle(20, 1)))$index
  expect_lt(square_span, 1)
  expect_lt(rect_span, square_span)
})

test_that("span_index is scale invariant", {
  # tolerance is looser than a pure-math scale-invariance guarantee would
  # need, because sfdct's own triangulation isn't perfectly scale-invariant
  # (the 10x-scaled star splits into slightly differently-shaped triangles,
  # not just a uniformly-scaled copy of the original mesh) - the self-term's
  # area-adaptive depth then quantizes each mesh's own triangles slightly
  # differently, on top of that pre-existing triangulation difference
  star <- make_star(7, 1, 0.4)
  base <- span_index(wc(star))$index
  expect_equal(span_index(wc(star * 10))$index, base, tolerance = 1e-4)
})

test_that("weight = triangle area reproduces the unweighted index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  plain <- span_index(star, prep = prep)$index
  wt    <- span_index(star, prep = prep, weight = prep$tri$area)$index
  expect_equal(wt, plain, tolerance = 1e-6)
})

test_that("span_index total_weight reports raw (pre-normalisation) sum", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area * 3
  res <- span_index(star, prep = prep, weight = w)
  expect_equal(res$total_weight, sum(w), tolerance = 1e-8)
})

test_that("uniform density (weight proportional to area) collapses the annulus reference to the disk reference", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  res <- span_index(star, prep = prep, weight = prep$tri$area)
  expect_equal(res$D_ref, shapeindices:::.disk_reference_D(res$area), tolerance = 1e-3)
})

test_that("deterministic = TRUE and deterministic = FALSE roughly agree", {
  square <- wc(make_square(5))
  ex <- span_index(square)
  mc <- suppressWarnings(span_index(square, deterministic = FALSE, n_lines = 20000, seed = 1))
  expect_equal(mc$D, ex$D, tolerance = 0.02)
})

test_that("deterministic = FALSE errors if n_quad is passed explicitly", {
  expect_error(span_index(wc(make_square()), deterministic = FALSE, n_quad = 3), "n_quad")
})

test_that("n_quad must be 1 or 3", {
  expect_error(span_index(wc(make_square()), n_quad = 2), "n_quad")
})

test_that("weight length is validated", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  expect_error(span_index(star, prep = prep, weight = 1:3), "one entry per triangle")
})

test_that("concentrating weight near the shape's own centroid scores closer to 1 than concentrating it at the edge", {
  # regression/sanity check for the Riesz-rearrangement claim behind the
  # annulus reference: since D_ref is the same for both (same weight
  # histogram, just relabelled which physical triangle gets which value),
  # only D_actual should differ, and it should differ in the direction
  # that rewards mass placed near the shape's own centre
  # a square CDTs into just 2 (symmetric) triangles, too coarse to tell
  # "centre" from "edge" apart - a star has enough triangles of enough
  # different distances from the centroid for the two weightings to differ
  star <- wc(make_star(8, 1, 0.4))
  prep <- prepare_polygon(star)
  cen <- sf::st_coordinates(sf::st_centroid(sf::st_union(prep$poly)))[1, 1:2]
  tc  <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(prep$tri)))
  d   <- sqrt((tc[, 1] - cen[1])^2 + (tc[, 2] - cen[2])^2)

  center_res <- span_index(star, prep = prep, weight = exp(-d / mean(d)))
  edge_res   <- span_index(star, prep = prep, weight = exp(d / mean(d)))

  expect_gt(center_res$index, edge_res$index)
})

test_that("sorting by density descending minimises the annulus reference (Riesz rearrangement check)", {
  # .annulus_reference_D() always sorts internally; this directly checks
  # that sorted order beats a shuffled order using the same ring-pair sum,
  # empirically confirming the rearrangement-inequality claim in the file
  # header rather than just trusting it
  set.seed(1)
  n <- 8
  tri_area <- runif(n, 0.5, 2)
  weight   <- runif(n, 0.1, 10)

  ring_pair_sum <- function(order_idx) {
    r_hi <- sqrt(cumsum(tri_area[order_idx]) / pi)
    r_lo <- c(0, r_hi[-length(r_hi)])
    W <- (weight / sum(weight))[order_idx]
    gl <- shapeindices:::.gauss_legendre(8)
    nodes <- Map(shapeindices:::.radial_nodes, r_lo, r_hi, MoreArgs = list(gl = gl))
    r_all <- unlist(lapply(nodes, `[[`, "r"))
    p_all <- unlist(Map(function(nd, w_i) nd$p * w_i, nodes, W))
    sum(outer(p_all, p_all) * outer(r_all, r_all, shapeindices:::.mean_chord))
  }

  sorted_order   <- order(weight / tri_area, decreasing = TRUE)
  shuffled_order <- sample(seq_len(n))

  expect_lte(ring_pair_sum(sorted_order), ring_pair_sum(shuffled_order))
})

test_that("span_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- span_index_sf(x))
  expect_equal(res$span_index[1], span_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$span_index[2], span_index(x[2, ])$index, tolerance = 1e-8)
})

## -- .annulus_reference_D() ring-count regression (was O(n_triangles^2)) --
##
## Was one Gauss-Legendre ring per triangle, an O((8*n)^2) pairwise sum
## computed unconditionally whenever `weight` is supplied - regardless of
## `deterministic`, since this is the reference value, not the actual
## index. Crashed with an out-of-memory error on a real 2169-triangle
## mesh, computing nothing else. Fixed by capping the ring count and
## building rings via a density-ratio-bounded greedy pass rather than
## one ring per triangle.

.grid_mesh_tri <- function(n_side) {
  sq <- function(cx, cy) sf::st_polygon(list(rbind(
    c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
    c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
  cells <- expand.grid(x = seq_len(n_side), y = seq_len(n_side))
  polys <- mapply(sq, cells$x, cells$y, SIMPLIFY = FALSE)
  u <- wc(sf::st_union(sf::st_sfc(polys)))
  list(u = u, prep = prepare_polygon(u))
}

test_that(".annulus_reference_D completes quickly on a mesh well past the old O(n^2) breaking point", {
  m <- .grid_mesh_tri(50)
  tri <- m$prep$tri
  skip_if(nrow(tri) < 120, "grid mesh too small for this test")
  set.seed(1)
  w <- runif(nrow(tri), 1, 100)
  t <- system.time(D_ref <- shapeindices:::.annulus_reference_D(tri$area, w))
  expect_true(is.finite(D_ref))
  expect_lt(t[["elapsed"]], 5)  # was minutes-to-crash territory before the fix
})

test_that("span_index stays <= 1 for an adversarial weight spike concentrated on a single (small-area) triangle", {
  # a single triangle with a wildly disproportionate weight relative to
  # its own area used to dilute across many merged neighbours under
  # naive rank/area binning, inflating the reference past the actual
  # value (verified during development: index = 1.22 without the fix)
  m <- .grid_mesh_tri(50)
  tri <- m$prep$tri
  skip_if(nrow(tri) < 120, "grid mesh too small for this test")
  w <- rep(1, nrow(tri))
  w[which.min(tri$area)] <- 1e6
  res <- span_index(m$u, prep = m$prep, weight = w, deterministic = TRUE)
  expect_lte(res$index, 1 + 1e-9)
  expect_true(is.finite(res$index))
})

test_that("span_index stays <= 1 across many random adversarial weight-spike trials", {
  m <- .grid_mesh_tri(50)
  tri <- m$prep$tri
  skip_if(nrow(tri) < 120, "grid mesh too small for this test")
  set.seed(42)
  for (trial in 1:20) {
    w <- rep(1, nrow(tri))
    n_spikes <- sample(1:5, 1)
    w[sample(nrow(tri), n_spikes)] <- 10^runif(n_spikes, 2, 8)
    res <- span_index(m$u, prep = m$prep, weight = w, deterministic = TRUE)
    expect_lte(res$index, 1 + 1e-9)
  }
})

test_that(".tri_self_mean_distance depth 0 collapses to a single centroid (distance 0)", {
  v <- rbind(c(0, 0), c(4, 0), c(0, 3))
  expect_equal(shapeindices:::.tri_self_mean_distance(v, depth = 0), 0)
})

test_that(".tri_self_mean_distance converges toward the same value as depth increases", {
  v <- rbind(c(0, 0), c(4, 0), c(0, 3))
  d1 <- shapeindices:::.tri_self_mean_distance(v, depth = 1)
  d4 <- shapeindices:::.tri_self_mean_distance(v, depth = 4)
  d6 <- shapeindices:::.tri_self_mean_distance(v, depth = 6)
  expect_true(abs(d6 - d4) < abs(d4 - d1))
})

test_that("span_index's self-term gives smaller triangles less adaptive depth than the mesh's largest", {
  poly <- make_spiky_square()
  prep <- prepare_polygon(poly)
  tri  <- prep$tri
  w    <- tri$area
  depth_i <- shapeindices:::.adaptive_tri_depth(w^2 * sqrt(tri$area), 4L)
  expect_equal(depth_i[which.max(tri$area)], 4L)
  expect_lt(depth_i[which.min(tri$area)], 4L)
})

test_that("span_index's adaptive self-term point count is substantially lower than fixed depth 4 on a mesh with real size variation", {
  poly <- make_spiky_square()
  prep <- prepare_polygon(poly)
  tri  <- prep$tri
  w    <- tri$area
  depth_i <- shapeindices:::.adaptive_tri_depth(w^2 * sqrt(tri$area), 4L)
  expect_lt(sum(4^depth_i), nrow(tri) * 4^4)
})

test_that("span_index gives a result close to the old fixed-depth-4 self-term on a mesh with real size variation", {
  poly <- make_spiky_square()
  prep <- prepare_polygon(poly)
  tri  <- prep$tri
  w    <- tri$area

  old_self_sum <- sum(vapply(seq_len(nrow(tri)), function(i) {
    v <- st_coordinates(st_geometry(tri)[i])[1:3, 1:2, drop = FALSE]
    w[i]^2 * shapeindices:::.tri_self_mean_distance(v, depth = 4)
  }, numeric(1)))
  depth_i <- shapeindices:::.adaptive_tri_depth(w^2 * sqrt(tri$area), 4L)
  new_self_sum <- sum(vapply(seq_len(nrow(tri)), function(i) {
    v <- st_coordinates(st_geometry(tri)[i])[1:3, 1:2, drop = FALSE]
    w[i]^2 * shapeindices:::.tri_self_mean_distance(v, depth = depth_i[i])
  }, numeric(1)))
  expect_equal(new_self_sum, old_self_sum, tolerance = 0.02)
})

test_that(".annulus_reference_D matches the exact per-triangle result below the ring cap (no behaviour change for typical meshes)", {
  star <- wc(make_star(8, 1, 0.4))
  prep <- prepare_polygon(star)
  skip_if(nrow(prep$tri) >= 100, "mesh too large - would trigger ring binning, not what this test checks")
  set.seed(1)
  w <- runif(nrow(prep$tri), 1, 50)
  # exact reference: one ring per triangle, no binning - the pre-fix formula
  exact_D_ref <- local({
    W <- shapeindices:::.normalize_weight(w)
    rho <- W / prep$tri$area
    ord <- order(rho, decreasing = TRUE)
    r_hi <- sqrt(cumsum(prep$tri$area[ord]) / pi)
    r_lo <- c(0, r_hi[-length(r_hi)])
    W <- W[ord]
    gl <- shapeindices:::.gauss_legendre(8)
    nodes <- Map(shapeindices:::.radial_nodes, r_lo, r_hi, MoreArgs = list(gl = gl))
    r_all <- unlist(lapply(nodes, `[[`, "r"))
    p_all <- unlist(Map(function(nd, w_i) nd$p * w_i, nodes, W))
    sum(outer(p_all, p_all) * outer(r_all, r_all, shapeindices:::.mean_chord))
  })
  expect_equal(shapeindices:::.annulus_reference_D(prep$tri$area, w), exact_D_ref, tolerance = 1e-10)
})
