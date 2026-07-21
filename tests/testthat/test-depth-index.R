# a thin annulus, for the annulus-vs-thin-rectangle indistinguishability
# test below - not in helper-shapes.R since no other test file needs it
make_annulus <- function(r_outer, r_inner, center = c(0, 0), n_seg = 90) {
  outer <- sf::st_buffer(sf::st_sfc(sf::st_point(center)), dist = r_outer, nQuadSegs = n_seg)
  inner <- sf::st_buffer(sf::st_sfc(sf::st_point(center)), dist = r_inner, nQuadSegs = n_seg)
  sf::st_difference(outer, inner)[[1]]
}

test_that("depth_index is ~1 for a disk", {
  disk <- make_disk(5, n_seg = 60)
  expect_equal(depth_index(wc(disk))$index, 1, tolerance = 0.01)
})

test_that("depth_index falls as a shape elongates", {
  square_di <- depth_index(wc(make_square()))$index
  rect_di   <- depth_index(wc(make_rectangle(20, 1)))$index
  expect_lt(square_di, 1)
  expect_lt(rect_di, square_di)
})

test_that("depth_index is scale invariant", {
  star <- make_star(7, 1, 0.4)
  base <- depth_index(wc(star))$index
  expect_equal(depth_index(wc(star * 10))$index, base, tolerance = 1e-3)
})

test_that("depth_index against a direct Monte Carlo cross-check", {
  # for a square [-h,h]^2, d(x,y) = h - max(|x|,|y|) - ground truth via
  # direct simulation, not through any package machinery
  h <- 5
  sq <- make_square(h)
  set.seed(1)
  n <- 300000
  x <- runif(n, -h, h); y <- runif(n, -h, h)
  mean_depth_true <- mean(h - pmax(abs(x), abs(y)))
  res <- depth_index(wc(sq))
  expect_equal(res$mean_depth, mean_depth_true, tolerance = 0.01)
})

test_that("depth_index handles a hole (shallower than the same square without one)", {
  outer <- rbind(c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5))
  hole  <- rbind(c(-2, -2), c(-2, 2), c(2, 2), c(2, -2), c(-2, -2))
  square_hole <- wc(sf::st_polygon(list(outer, hole)))
  res <- depth_index(square_hole)
  expect_true(is.finite(res$index))
  expect_lt(res$index, depth_index(wc(sf::st_polygon(list(outer))))$index)
})

test_that("depth_index handles a multi-part shape without hanging", {
  sq1 <- make_square(1, center = c(0, 0))
  sq2 <- make_square(1, center = c(10, 0))
  dispersed <- wc(sf::st_union(sf::st_sfc(sq1, sq2)))
  res <- depth_index(dispersed)
  expect_true(is.finite(res$index))
})

test_that("depth_index is unaffected by how far apart two identical parts are, holding their combined area fixed", {
  # each point's own nearest edge is local, so mean_depth doesn't react to
  # separation - UNLIKE span_index()/radial_concentration_index()/
  # convexity_index(), which all connect distant parts via lines/pairs and
  # so do react. This is a different claim from the previous test: total
  # area (and hence the disk reference) must stay fixed for this to hold -
  # comparing against a single square of HALF the area would conflate
  # fragmentation with the reference simply scaling up.
  two_parts <- function(gap) {
    sq1 <- make_square(1, center = c(0, 0))
    sq2 <- make_square(1, center = c(2 + gap, 0))
    wc(sf::st_union(sf::st_sfc(sq1, sq2)))
  }
  near <- depth_index(two_parts(0.01))$index
  far  <- depth_index(two_parts(100))$index
  expect_equal(far, near, tolerance = 1e-6)
})

test_that("an annulus and a thin rectangle of the same width and area are indistinguishable", {
  # documented characteristic, not a bug: depth only sees LOCAL distance to
  # the nearest edge, so global topology (a ring vs. a strip) is invisible
  # to it once the width matches
  width <- 1
  r_out <- 5; r_in <- r_out - width
  ann <- wc(make_annulus(r_out, r_in))
  ann_area <- as.numeric(sf::st_area(sf::st_sfc(ann)))
  rect <- wc(make_rectangle(ann_area / width, width))
  expect_equal(depth_index(ann)$index, depth_index(rect)$index, tolerance = 0.01)
})

test_that("weight = triangle area reproduces the unweighted index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  plain <- depth_index(star, prep = prep)$index
  wt    <- depth_index(star, prep = prep, weight = prep$tri$area)$index
  expect_equal(wt, plain, tolerance = 1e-6)
})

test_that("depth_index total_weight reports raw (pre-normalisation) sum", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area * 3
  res <- depth_index(star, prep = prep, weight = w)
  expect_equal(res$total_weight, sum(w), tolerance = 1e-8)
})

test_that("uniform density (weight proportional to area) collapses the annulus reference to the disk reference", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  res <- depth_index(star, prep = prep, weight = prep$tri$area)
  expect_equal(res$ref_depth, shapeindices:::.disk_reference_depth(res$area), tolerance = 1e-8)
})

test_that(".disk_reference_depth matches the closed form R/3", {
  area <- 78.53981634  # pi * 5^2
  R <- sqrt(area / pi)
  expect_equal(shapeindices:::.disk_reference_depth(area), R / 3, tolerance = 1e-8)
})

test_that("weight concentrated near the shape's own centre scores higher than weight concentrated at the edge", {
  # per file header: densest-at-centre MAXIMISES weighted mean depth, the
  # opposite reason radial_concentration_index()'s own reference sorts the
  # same way (there, to minimise mean distance-to-centre)
  disk <- wc(make_disk(5, n_seg = 60))
  prep <- prepare_polygon(disk)
  tc <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(prep$tri)))
  d  <- sqrt(tc[, 1]^2 + tc[, 2]^2)

  center_res <- depth_index(disk, prep = prep, weight = prep$tri$area * exp(-d / mean(d)))
  edge_res   <- depth_index(disk, prep = prep, weight = prep$tri$area * exp(d / mean(d)))
  expect_gt(center_res$index, edge_res$index)
})

## -- 0/NA weight consistency with the other mesh indices -------------------

test_that("depth_index's own `weight` errors on NA, same as the other indices", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area
  w[1] <- NA
  expect_error(depth_index(star, prep = prep, weight = w), "contains NA")
})

test_that("depth_index's own `weight` allows an individual 0 (zero density, not a hole)", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area
  w[1] <- 0
  expect_no_error(res <- depth_index(star, prep = prep, weight = w))
  expect_true(is.finite(res$index))
})

test_that("weight length is validated", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  expect_error(depth_index(star, prep = prep, weight = 1:3), "one entry per triangle")
})

test_that("depth_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- depth_index_sf(x))
  expect_equal(res$depth_index[1], depth_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$depth_index[2], depth_index(x[2, ])$index, tolerance = 1e-8)
})

## -- deterministic = FALSE (Monte Carlo) mode ------------------------------

test_that("deterministic = TRUE and deterministic = FALSE roughly agree", {
  star <- wc(make_star(8, 1, 0.4))
  det_val <- depth_index(star, deterministic = TRUE)$index
  reps <- vapply(1:5, function(s) {
    depth_index(star, deterministic = FALSE, n_lines = 8000, seed = s)$index
  }, numeric(1))
  expect_equal(mean(reps), det_val, tolerance = 0.02)
})

test_that("deterministic = TRUE and deterministic = FALSE roughly agree, weighted", {
  star <- wc(make_star(8, 1, 0.4))
  prep <- prepare_polygon(star)
  set.seed(1)
  w <- runif(nrow(prep$tri), 1, 10)
  det_val <- depth_index(star, prep = prep, weight = w, deterministic = TRUE)$index
  reps <- vapply(1:5, function(s) {
    depth_index(star, prep = prep, weight = w, deterministic = FALSE, n_lines = 8000, seed = s)$index
  }, numeric(1))
  expect_equal(mean(reps), det_val, tolerance = 0.03)
})

test_that("deterministic = FALSE is reproducible with a seed", {
  star <- wc(make_star(8, 1, 0.4))
  a <- depth_index(star, deterministic = FALSE, n_lines = 500, seed = 42)$index
  b <- depth_index(star, deterministic = FALSE, n_lines = 500, seed = 42)$index
  expect_identical(a, b)
})

test_that("`points` (a pre-drawn sample) with deterministic = TRUE errors", {
  star <- wc(make_star(6, 1, 0.4))
  expect_error(
    depth_index(star, deterministic = TRUE, points = matrix(1:4, 2, 2)),
    "no meaning for deterministic = TRUE"
  )
})

test_that("`points` supplied to deterministic = FALSE skips internal sampling and is used directly", {
  star <- wc(make_star(6, 1, 0.4))
  # a single point at the centre, repeated - mean_depth must equal that
  # point's own distance to the boundary exactly, regardless of n_lines/seed
  prep <- prepare_polygon(star)
  ctr  <- sf::st_coordinates(sf::st_centroid(sf::st_union(prep$poly)))[1, 1:2]
  pts  <- matrix(rep(ctr, 4), ncol = 2, byrow = TRUE)
  d_true <- as.numeric(sf::st_distance(sf::st_sfc(sf::st_point(ctr), crs = sf::st_crs(prep$poly)),
                                        sf::st_boundary(prep$poly)))
  res <- depth_index(star, prep = prep, deterministic = FALSE, n_lines = 999, seed = 999, points = pts)
  expect_equal(res$mean_depth, d_true, tolerance = 1e-6)
})

test_that("shape_indices() with multiple Monte Carlo mesh indices (including depth) shares one point draw", {
  star <- wc(make_star(8, 1, 0.3))
  r <- suppressWarnings(shape_indices(star, which = c("radial_concentration", "depth"),
                                       deterministic = FALSE, n_lines = 3000, seed = 1))
  expect_true(all(is.finite(r)))
  expect_true(all(r > 0 & r <= 1.05))
})

test_that("shape_indices_sf(byrow = FALSE) with multiple Monte Carlo mesh indices (including depth) shares one point draw", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  res <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, id = "z",
                           which = c("radial_concentration", "depth"),
                           deterministic = FALSE, n_lines = 3000, seed = 1))
  expect_true(is.finite(res$radial_concentration_index))
  expect_true(is.finite(res$depth_index))
})

test_that("shape_indices_sf(byrow = FALSE, weights = ...) computes depth on the exact constrained mesh", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  res <- shape_indices_sf(x, byrow = FALSE, id = "z", which = "depth", weights = c(5, 50))
  expect_true(is.finite(res$depth_index))
  expect_gt(res$depth_index, 0)
})
