test_that("convexity_index is 1 for convex shapes", {
  expect_equal(convexity_index(wc(make_square()))$index, 1, tolerance = 1e-8)
  expect_equal(convexity_index(wc(make_rectangle(10, 2)))$index, 1, tolerance = 1e-8)
  expect_equal(convexity_index(wc(make_star(6, 1, 1)))$index, 1, tolerance = 1e-8) # ratio 1 -> regular hexagon
  expect_equal(convexity_index(wc(make_disk(1, n_seg = 12)))$index, 1, tolerance = 1e-3)
})

test_that("convexity_index drops below 1 for a notched star, deeper notches score lower", {
  # ratio = 0.9 (very shallow notch) is inside CI's own documented blind
  # spot on symmetric few-point stars (see the "Understanding Convexity
  # Index" vignette) - the mesh isn't fine enough to detect a notch that
  # shallow, so CI == 1 there is expected, not a bug. Use a moderate ratio
  # that's reliably detected instead.
  shallow <- convexity_index(wc(make_star(6, 1, 0.5)))$index
  deep    <- convexity_index(wc(make_star(6, 1, 0.1)))$index
  expect_lt(shallow, 1)
  expect_lt(deep, shallow)
})

test_that("convexity_index is approximately scale and translation invariant", {
  # "approximately", not exactly: scaling/shifting changes the exact
  # floating-point coordinates CDT triangulates, which can shift the mesh
  # topology slightly (a different, but equally valid, triangulation) -
  # the index is invariant in the underlying geometry, not bit-for-bit in
  # a specific mesh realisation of it
  base <- make_star(7, 1, 0.4)
  ci_base <- convexity_index(wc(base))$index

  scaled <- base * 5
  expect_equal(convexity_index(wc(scaled))$index, ci_base, tolerance = 0.01)

  shifted <- base + c(100, -50)
  expect_equal(convexity_index(wc(shifted))$index, ci_base, tolerance = 0.01)
})

test_that("deterministic = TRUE and deterministic = FALSE (random-line) roughly agree", {
  star <- wc(make_star(8, 1, 0.4))
  det_val    <- convexity_index(star, deterministic = TRUE)$index
  approx_val <- convexity_index(star, deterministic = FALSE, n_lines = 4000, seed = 1)$index
  expect_equal(approx_val, det_val, tolerance = 0.05)
})

test_that("n_quad passed explicitly together with deterministic = FALSE errors", {
  star <- wc(make_star(6, 1, 0.4))
  expect_error(convexity_index(star, deterministic = FALSE, n_quad = 1), "deterministic = FALSE")
})

test_that("convexity_index no longer accepts mesh/parallel arguments", {
  # regression guard: convexity_index() dropped mesh (convex decomposition
  # is no longer wired into it, though convex_decompose() itself remains a
  # standalone utility - see test-convex-decompose.R) and parallel/
  # min_pairs_for_parallel (pair-level parallelism removed entirely - see
  # test-utils.R for the row-level parallelism that replaced it)
  star <- wc(make_star(6, 1, 0.3))
  expect_error(convexity_index(star, mesh = "decompose"), "unused argument")
  expect_error(convexity_index(star, parallel = TRUE), "unused argument")
  expect_error(convexity_index(star, min_pairs_for_parallel = 10), "unused argument")
})

test_that("weight = triangle area reproduces the unweighted index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  ci_plain  <- convexity_index(star, prep = prep)$index
  ci_wtarea <- convexity_index(star, prep = prep, weight = prep$tri$area)$index
  expect_equal(ci_wtarea, ci_plain, tolerance = 1e-8)
})

test_that("deterministic = TRUE rejects a wrong-length weight instead of silently recycling it", {
  # regression: .mesh_convexity_index() used to have no length check at all
  # (unlike moment_of_inertia_index()/the random-line path, which both
  # already validated this) - a wrong-length weight silently recycled via
  # ordinary R vector arithmetic instead of erroring, observed to come out
  # as index = NA with no warning for one concrete case
  sq <- wc(make_square())
  prep <- prepare_polygon(sq)
  n_tri <- nrow(prep$tri)
  expect_error(convexity_index(sq, prep = prep, weight = 5), "one entry per triangle")
  expect_error(convexity_index(sq, prep = prep, weight = rep(1, n_tri + 1)), "one entry per triangle")
  expect_no_error(convexity_index(sq, prep = prep, weight = rep(1, n_tri)))
})

test_that("a single-triangle polygon is vacuously convex", {
  tri_poly <- wc(sf::st_polygon(list(rbind(c(0, 0), c(4, 0), c(0, 3), c(0, 0)))))
  expect_equal(convexity_index(tri_poly)$index, 1)
})

test_that("convexity_index_sf preserves CRS and matches row-by-row convexity_index", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- convexity_index_sf(x))
  expect_equal(res$convexity_index[1], convexity_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$convexity_index[2], convexity_index(x[2, ])$index, tolerance = 1e-8)
})

test_that(".frac_outside_geos computes fraction-outside correctly (sequential GEOS fallback)", {
  # a 10x10 square centred at the origin - one line straddling it exactly
  # half in/half out, one line fully interior
  square <- make_square(5)
  seg <- sf::st_sfc(
    sf::st_linestring(rbind(c(0, 0), c(10, 0))),   # half inside (0-5), half outside (5-10)
    sf::st_linestring(rbind(c(-2, -2), c(2, 2)))   # fully interior
  )
  len_total <- as.numeric(sf::st_length(seg))

  result <- shapeindices:::.frac_outside_geos(seg, len_total, idx_out = 1:2, poly_u = square)
  expect_equal(result, c(0.5, 0), tolerance = 1e-6)
})
