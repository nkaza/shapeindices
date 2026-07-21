test_that(".geometric_median matches the centroid for a symmetric point cloud", {
  # a square grid of equal-weight points: by symmetry the geometric median
  # is the centroid exactly
  p <- as.matrix(expand.grid(x = c(-1, 1), y = c(-1, 1)))
  gm <- shapeindices:::.geometric_median(p, rep(1, 4))
  expect_equal(unname(gm$center), c(0, 0), tolerance = 1e-6)
})

test_that(".geometric_median converges to the connecting midpoint for two equal weighted points (degenerate case)", {
  # classical degeneracy: for two equal point masses, every point on the
  # connecting segment minimises the objective equally - Weiszfeld should
  # still converge (to the segment's midpoint, since it starts there) and
  # not hang or error
  p <- rbind(c(0, 0), c(10, 0))
  gm <- shapeindices:::.geometric_median(p, c(1, 1))
  expect_equal(unname(gm$center), c(5, 0), tolerance = 1e-6)
  expect_equal(gm$D1, 5, tolerance = 1e-6)
})

test_that("radial_concentration_index is ~1 for a disk", {
  disk <- make_disk(5, n_seg = 60)
  expect_equal(radial_concentration_index(wc(disk))$index, 1, tolerance = 0.01)
})

test_that("radial_concentration_index falls as a shape elongates", {
  square_rci <- radial_concentration_index(wc(make_square()))$index
  rect_rci   <- radial_concentration_index(wc(make_rectangle(20, 1)))$index
  expect_lt(square_rci, 1)
  expect_lt(rect_rci, square_rci)
})

test_that("radial_concentration_index is scale invariant", {
  # looser tolerance than an exact closed form would need: CDT of the
  # scaled and unscaled star can differ by tiny floating-point tie-breaks
  # (a different, equally valid triangulation), same order-of-magnitude
  # effect as span_index()'s cross-triangulation tolerance
  star <- make_star(7, 1, 0.4)
  base <- radial_concentration_index(wc(star))$index
  expect_equal(radial_concentration_index(wc(star * 10))$index, base, tolerance = 1e-4)
})

test_that("radial_concentration_index against a direct Monte Carlo cross-check", {
  # for a square, the geometric median is exactly the centre (4-fold
  # symmetry), so ground truth is just mean distance from a uniform
  # random point to the centre
  sq <- make_square(5) # half = 5 -> a 10x10 square centred at the origin
  set.seed(1)
  n <- 300000
  x <- runif(n, -5, 5); y <- runif(n, -5, 5)
  D1_true <- mean(sqrt(x^2 + y^2))
  res <- radial_concentration_index(wc(sq))
  expect_equal(res$D1, D1_true, tolerance = 0.01)
})

test_that("radial_concentration_index handles a hole (centre falls inside it, by symmetry)", {
  outer <- rbind(c(-5, -5), c(5, -5), c(5, 5), c(-5, 5), c(-5, -5))
  hole  <- rbind(c(-2, -2), c(-2, 2), c(2, 2), c(2, -2), c(-2, -2))
  square_hole <- wc(sf::st_polygon(list(outer, hole)))
  res <- radial_concentration_index(square_hole)
  expect_true(is.finite(res$index))
  expect_lt(res$index, 1)
  expect_equal(unname(sf::st_coordinates(res$center)[1, ]), c(0, 0), tolerance = 1e-6)
})

test_that("radial_concentration_index handles a multi-part (dispersed) shape without hanging", {
  # degenerate geometric median case (two identical, separated parts) -
  # the centre should land exactly on the connecting midpoint and the
  # computation should still terminate promptly
  sq1 <- make_square(1, center = c(0, 0))
  sq2 <- make_square(1, center = c(10, 0))
  dispersed <- wc(sf::st_union(sf::st_sfc(sq1, sq2)))
  res <- radial_concentration_index(dispersed)
  expect_true(is.finite(res$index))
  expect_equal(unname(sf::st_coordinates(res$center)[1, ]), c(5, 0), tolerance = 1e-4)
  # dispersing should score much lower than the solid union would
  solid <- wc(make_square(1, center = c(0, 0)))
  expect_lt(res$index, radial_concentration_index(solid)$index)
})

test_that("weight = triangle area reproduces the unweighted index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  plain <- radial_concentration_index(star, prep = prep)$index
  wt    <- radial_concentration_index(star, prep = prep, weight = prep$tri$area)$index
  expect_equal(wt, plain, tolerance = 1e-6)
})

test_that("radial_concentration_index total_weight reports raw (pre-normalisation) sum", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area * 3
  res <- radial_concentration_index(star, prep = prep, weight = w)
  expect_equal(res$total_weight, sum(w), tolerance = 1e-8)
})

test_that("uniform density (weight proportional to area) collapses the annulus reference to the disk reference", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  res <- radial_concentration_index(star, prep = prep, weight = prep$tri$area)
  expect_equal(res$D1_ref, shapeindices:::.disk_reference_D1(res$area), tolerance = 1e-8)
})

test_that("weight concentrated near the shape's own centre scores higher than weight concentrated at the edge", {
  star <- wc(make_star(8, 1, 0.4))
  prep <- prepare_polygon(star)
  cen <- sf::st_coordinates(sf::st_centroid(sf::st_union(prep$poly)))[1, 1:2]
  tc  <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(prep$tri)))
  d   <- sqrt((tc[, 1] - cen[1])^2 + (tc[, 2] - cen[2])^2)

  center_res <- radial_concentration_index(star, prep = prep, weight = exp(-d / mean(d)))
  edge_res   <- radial_concentration_index(star, prep = prep, weight = exp(d / mean(d)))
  expect_gt(center_res$index, edge_res$index)
})

## -- 0/NA weight consistency with the other three indices -----------------

test_that("radial_concentration_index's own `weight` errors on NA, same as the other indices", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area
  w[1] <- NA
  expect_error(radial_concentration_index(star, prep = prep, weight = w), "contains NA")
  # same error, from the same .normalize_weight() call, as convexity/MOI/span
  expect_error(convexity_index(star, prep = prep, weight = w), "contains NA")
  expect_error(moment_of_inertia_index(star, prep = prep, weight = w), "contains NA")
  expect_error(span_index(star, prep = prep, weight = w), "contains NA")
})

test_that("radial_concentration_index's own `weight` allows an individual 0 (zero density, not a hole)", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area
  w[1] <- 0
  expect_no_error(res <- radial_concentration_index(star, prep = prep, weight = w))
  expect_true(is.finite(res$index))
})

test_that("weight length is validated", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  expect_error(radial_concentration_index(star, prep = prep, weight = 1:3), "one entry per triangle")
})

test_that("radial_concentration_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- radial_concentration_index_sf(x))
  expect_equal(res$radial_concentration_index[1], radial_concentration_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$radial_concentration_index[2], radial_concentration_index(x[2, ])$index, tolerance = 1e-8)
})

## -- deterministic = FALSE (Monte Carlo) mode ------------------------------

test_that(".radial_point_cloud preserves total weight and point count regardless of vectorisation", {
  star <- wc(make_star(8, 1, 0.4))
  prep <- prepare_polygon(star)
  tri  <- prep$tri
  w    <- tri$area
  cloud <- shapeindices:::.radial_point_cloud(tri, w, max_depth = 4)
  # the star's own triangles vary somewhat in size, so this is no longer
  # necessarily exactly n_tri * 4^4 (see the adaptive-depth tests below) -
  # what must still hold is the upper bound and weight conservation
  expect_lte(nrow(cloud$p), nrow(tri) * 4^4)
  expect_gte(nrow(cloud$p), nrow(tri))
  expect_equal(sum(cloud$w), sum(w), tolerance = 1e-8)
  # every sub-triangle centroid must lie within its own parent triangle's
  # bounding box - a cheap sanity check that vertices weren't scrambled
  # across triangles during the batched rbind() stacking
  bb <- sf::st_bbox(star)
  expect_true(all(cloud$p[, 1] >= bb["xmin"] - 1e-6 & cloud$p[, 1] <= bb["xmax"] + 1e-6))
  expect_true(all(cloud$p[, 2] >= bb["ymin"] - 1e-6 & cloud$p[, 2] <= bb["ymax"] + 1e-6))
})

## -- area-adaptive subdivision depth ----------------------------------------
##
## A fixed depth (4^4 = 256 points) for every triangle regardless of its
## own size was wasteful on real meshes: a complex boundary produces many
## small triangles alongside a handful of large ones, and only the large
## ones actually need deep subdivision to control the Jensen's-inequality
## centroid-collapse bias (see file header). Verified during development:
## 62-65x fewer points on real Census-block-derived meshes (area ratios of
## tens of millions to one between the largest and smallest triangle),
## with the resulting index changing by <0.05% relative to a much
## higher-resolution reference - invisible at the 3 decimal places this
## package reports indices to.

test_that("smaller triangles get less subdivision depth than larger ones", {
  poly <- make_spiky_square()
  prep <- prepare_polygon(poly)
  tri  <- prep$tri
  skip_if(nrow(tri) < 3, "mesh too small for this test")
  depth_i <- shapeindices:::.adaptive_tri_depth(tri$area, 4L)
  expect_equal(depth_i[which.max(tri$area)], 4L)
  expect_lt(depth_i[which.min(tri$area)], 4L)
})

test_that("adaptive depth substantially reduces point count on a mesh with real size variation", {
  poly <- make_spiky_square()
  prep <- prepare_polygon(poly)
  tri  <- prep$tri
  skip_if(nrow(tri) < 3, "mesh too small for this test")
  w <- tri$area
  cloud <- shapeindices:::.radial_point_cloud(tri, w)
  old_style_points <- nrow(tri) * 4^4
  expect_lt(nrow(cloud$p), old_style_points)
  expect_equal(sum(cloud$w), sum(w), tolerance = 1e-8)
})

test_that("adaptive depth gives an index close to a much higher-resolution reference", {
  poly <- make_spiky_square()
  prep <- prepare_polygon(poly)
  tri  <- prep$tri
  skip_if(nrow(tri) < 3, "mesh too small for this test")
  w <- tri$area

  cloud_adaptive <- shapeindices:::.radial_point_cloud(tri, w)
  gm_adaptive <- shapeindices:::.geometric_median(cloud_adaptive$p, cloud_adaptive$w)

  # depth=6 for every triangle, far finer than the adaptive default ever
  # uses even for the largest triangle - a "ground truth" for comparison
  hi_res <- local({
    corner_mat <- vapply(sf::st_geometry(tri), function(g) {
      v <- sf::st_coordinates(g)[1:3, 1:2, drop = FALSE]
      c(v[1, ], v[2, ], v[3, ])
    }, numeric(6))
    A <- t(corner_mat[1:2, , drop = FALSE]); B <- t(corner_mat[3:4, , drop = FALSE]); C <- t(corner_mat[5:6, , drop = FALSE])
    orig_idx <- seq_len(nrow(tri))
    for (d in seq_len(6)) {
      sub <- shapeindices:::.subdivide_tri_batch(A, B, C)
      A <- sub$A; B <- sub$B; C <- sub$C
      orig_idx <- rep(orig_idx, times = 4)
    }
    list(p = (A + B + C) / 3, w = (w / 4^6)[orig_idx])
  })
  gm_hi <- shapeindices:::.geometric_median(hi_res$p, hi_res$w)

  expect_equal(gm_adaptive$D1, gm_hi$D1, tolerance = 0.01)  # within 1%, generous vs. the <0.05% seen empirically
})

test_that("deterministic = TRUE and deterministic = FALSE roughly agree", {
  star <- wc(make_star(8, 1, 0.4))
  det_val <- radial_concentration_index(star, deterministic = TRUE)$index
  reps <- vapply(1:5, function(s) {
    radial_concentration_index(star, deterministic = FALSE, n_lines = 8000, seed = s)$index
  }, numeric(1))
  expect_equal(mean(reps), det_val, tolerance = 0.02)
})

test_that("deterministic = TRUE and deterministic = FALSE roughly agree, weighted", {
  star <- wc(make_star(8, 1, 0.4))
  prep <- prepare_polygon(star)
  set.seed(1)
  w <- runif(nrow(prep$tri), 1, 10)
  det_val <- radial_concentration_index(star, prep = prep, weight = w, deterministic = TRUE)$index
  reps <- vapply(1:5, function(s) {
    radial_concentration_index(star, prep = prep, weight = w, deterministic = FALSE, n_lines = 8000, seed = s)$index
  }, numeric(1))
  expect_equal(mean(reps), det_val, tolerance = 0.03)
})

test_that("deterministic = FALSE is reproducible with a seed", {
  star <- wc(make_star(8, 1, 0.4))
  a <- radial_concentration_index(star, deterministic = FALSE, n_lines = 500, seed = 42)$index
  b <- radial_concentration_index(star, deterministic = FALSE, n_lines = 500, seed = 42)$index
  expect_identical(a, b)
})

test_that("`points` (a pre-drawn sample) with deterministic = TRUE errors", {
  star <- wc(make_star(6, 1, 0.4))
  expect_error(
    radial_concentration_index(star, deterministic = TRUE, points = matrix(1:4, 2, 2)),
    "no meaning for deterministic = TRUE"
  )
})

test_that("`points` supplied to deterministic = FALSE skips internal sampling and is used directly", {
  star <- wc(make_star(6, 1, 0.4))
  # 4 points at the same spot -> geometric median must land there exactly,
  # regardless of n_lines/seed (neither should matter once points is given)
  pts <- matrix(c(1, 1, 1, 1, 2, 2, 2, 2), ncol = 2)
  res <- radial_concentration_index(star, deterministic = FALSE, n_lines = 999, seed = 999, points = pts)
  expect_equal(unname(sf::st_coordinates(res$center)[1, ]), c(1, 2), tolerance = 1e-8)
})

test_that("shape_indices() with multiple Monte Carlo mesh indices shares one point draw and still gives finite, sensible results", {
  star <- wc(make_star(8, 1, 0.3))
  # suppressWarnings(): this test's tiny mesh legitimately trips the
  # unrelated "n_lines not substantially lower than exhaustive" advisory
  # warning at n_lines = 3000 - not what's being tested here
  r <- suppressWarnings(shape_indices(star, which = c("convexity", "span", "radial_concentration"),
                                       deterministic = FALSE, n_lines = 3000, seed = 1))
  expect_true(all(is.finite(r)))
  expect_true(all(r > 0 & r <= 1))
})

test_that("shape_indices_sf(byrow = FALSE) with multiple Monte Carlo mesh indices shares one point draw", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  res <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, id = "z",
                           which = c("convexity", "span", "radial_concentration"),
                           deterministic = FALSE, n_lines = 3000, seed = 1))
  expect_true(is.finite(res$convexity_index))
  expect_true(is.finite(res$span_index))
  expect_true(is.finite(res$radial_concentration_index))
})
