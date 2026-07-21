# helper: an nx*ny grid of unit squares, one row per cell, with a weight
# function so tests can control how much density varies cell-to-cell -
# real per-row density variance is what actually exercises the
# constrained-triangulation path (see shape-indices.R's dispatch: it's
# skipped entirely when weights = NULL, since density is then uniform)
make_weighted_grid <- function(nx, ny, weight_fn = function(i, j) 1 + i + j * 3) {
  polys <- vector("list", nx * ny)
  w <- numeric(nx * ny)
  k <- 1
  for (i in 0:(nx - 1)) for (j in 0:(ny - 1)) {
    polys[[k]] <- make_square(0.5, center = c(i, j))
    w[k] <- weight_fn(i, j)
    k <- k + 1
  }
  list(x = sf::st_sf(pop = w, geometry = sf::st_sfc(polys, crs = TEST_CRS)), w = w)
}

test_that(".extract_pslg() handles holes and multi-part features without dropping rings", {
  square_hole <- make_square_with_hole(outer_half = 5, hole_frac = 0.25) # 2 rings: outer + hole
  second_part <- make_square(1, center = c(20, 20))                     # 1 ring, separate part
  mp <- sf::st_multipolygon(list(unclass(square_hole), unclass(second_part)))
  geoms <- sf::st_sfc(mp, crs = TEST_CRS)

  res <- shapeindices:::.extract_pslg(geoms)
  expect_equal(res$n_rings, 3) # outer + hole + second part's own ring
  expect_equal(nrow(res$S), nrow(res$P)) # every ring is closed: #segments == #vertices
})

test_that(".extract_pslg() shares vertices between touching rows via coordinate rounding", {
  sq1 <- make_square(1, center = c(-1, 0))
  sq2 <- make_square(1, center = c(1, 0)) # shares the x=0 edge with sq1
  geoms <- sf::st_sfc(sq1, sq2, crs = TEST_CRS)

  res <- shapeindices:::.extract_pslg(geoms)
  # the shared edge's 2 endpoints must collapse to 2 (not 4) unique vertices:
  # 8 boundary points total, minus 2 shared -> 6 unique
  expect_equal(nrow(res$P), 6)
})

test_that(".constrained_weighted_mesh() recovers the true total weight exactly", {
  g <- make_weighted_grid(3, 3)
  res <- shapeindices:::.constrained_weighted_mesh(g$x, weights = "pop")
  expect_true(res$ok)
  expect_equal(sum(res$tri_weight), sum(g$w), tolerance = 1e-8)
  expect_equal(res$raw_total, sum(g$w))
})

test_that(".constrained_weighted_mesh() excludes hole area from weight and triangle area", {
  outer <- make_square(5)
  hole_poly <- make_square(1) # a real hole cut from the outer row via st_difference
  row_geom <- sf::st_difference(sf::st_sfc(outer, crs = TEST_CRS), sf::st_sfc(hole_poly, crs = TEST_CRS))
  x <- sf::st_sf(pop = 100, geometry = row_geom)

  res <- shapeindices:::.constrained_weighted_mesh(x, weights = "pop")
  expect_true(res$ok)
  expect_equal(sum(res$tri_area), as.numeric(sf::st_area(row_geom)), tolerance = 1e-6)
  # the hole must NOT contribute area even though it's inside the outer square's bbox
  expect_lt(sum(res$tri_area), as.numeric(sf::st_area(sf::st_sfc(outer, crs = TEST_CRS))))
})

test_that(".constrained_weighted_mesh() reports ok = FALSE for genuinely overlapping rows", {
  sq1 <- make_square(2, center = c(0, 0))
  sq2 <- make_square(2, center = c(1, 0)) # overlaps sq1
  x <- sf::st_sf(pop = c(10, 20), geometry = wc(sf::st_sfc(sq1, sq2)))

  res <- shapeindices:::.constrained_weighted_mesh(x, weights = "pop")
  expect_false(res$ok)
  expect_match(res$reason, "overlap")
})

test_that(".constrained_weighted_mesh() does not flag touching (non-overlapping) rows as unsafe", {
  sq1 <- make_square(1, center = c(-1, 0))
  sq2 <- make_square(1, center = c(1, 0))
  x <- sf::st_sf(pop = c(10, 20), geometry = wc(sf::st_sfc(sq1, sq2)))

  res <- shapeindices:::.constrained_weighted_mesh(x, weights = "pop")
  expect_true(res$ok)
})

test_that(".constrained_weighted_mesh() density is scale-invariant under a constant weight multiplier", {
  g <- make_weighted_grid(3, 3)
  res1 <- shapeindices:::.constrained_weighted_mesh(g$x, weights = "pop")
  x2 <- g$x
  x2$pop <- x2$pop * 1e6
  res2 <- shapeindices:::.constrained_weighted_mesh(x2, weights = "pop")

  # per-triangle DENSITY (weight / area) must be identical up to the same
  # constant multiplier - not just the total
  expect_equal(res2$tri_weight / res1$tri_weight, rep(1e6, length(res1$tri_weight)), tolerance = 1e-8)
})

test_that(".safe_deterministic_tri_ceiling() ceiling grows with available memory and shrinks with n_quad", {
  testthat::local_mocked_bindings(.available_memory_mb = function() 1000, .package = "shapeindices")
  small_mem_ceiling <- shapeindices:::.safe_deterministic_tri_ceiling(3, "convexity")

  testthat::local_mocked_bindings(.available_memory_mb = function() 100000, .package = "shapeindices")
  large_mem_ceiling <- shapeindices:::.safe_deterministic_tri_ceiling(3, "convexity")

  expect_gt(large_mem_ceiling, small_mem_ceiling)

  ceiling_q1 <- shapeindices:::.safe_deterministic_tri_ceiling(1, "convexity")
  ceiling_q3 <- shapeindices:::.safe_deterministic_tri_ceiling(3, "convexity")
  expect_gt(ceiling_q1, ceiling_q3)
})

test_that(".safe_deterministic_tri_ceiling() inverts each index's own cost formula correctly", {
  testthat::local_mocked_bindings(.available_memory_mb = function() 500, .package = "shapeindices")
  n_quad <- 3

  n_conv <- shapeindices:::.safe_deterministic_tri_ceiling(n_quad, "convexity")
  budget_units <- 500 * 1024^2 * 0.2 / 200
  expect_lte(choose(n_conv, 2) * n_quad^2, budget_units)
  expect_gt(choose(n_conv + 1, 2) * n_quad^2, budget_units)

  n_span <- shapeindices:::.safe_deterministic_tri_ceiling(n_quad, "span")
  expect_lte((n_span * n_quad)^2, budget_units)
  expect_gt(((n_span + 1) * n_quad)^2, budget_units)
})

test_that("array-native moment-of-inertia matches the sf-based core exactly on a real constrained mesh", {
  g <- make_weighted_grid(3, 3)
  cwm <- shapeindices:::.constrained_weighted_mesh(g$x, weights = "pop")
  expect_true(cwm$ok)

  arr <- shapeindices:::.moment_of_inertia_core_array(
    cwm$P, cwm$T, cwm$tri_area, cwm$tri_weight / cwm$tri_area, cwm$crs)

  pieces <- shapeindices:::.array_mesh_to_sf_pieces(cwm$P, cwm$T, cwm$tri_area, cwm$crs)
  sfbased <- shapeindices:::.moment_of_inertia_core(
    sf::st_geometry(pieces), pieces$area, cwm$tri_weight / cwm$tri_area, cwm$poly_u)

  expect_equal(arr$index, sfbased$index, tolerance = 1e-10)
  expect_equal(arr$Ixx, sfbased$Ixx, tolerance = 1e-10)
  expect_equal(arr$Iyy, sfbased$Iyy, tolerance = 1e-10)
  expect_equal(arr$Ixy, sfbased$Ixy, tolerance = 1e-10)
})

test_that("array-native directional-balance/radial-concentration match their sf-based mirrors", {
  g <- make_weighted_grid(3, 3)
  cwm <- shapeindices:::.constrained_weighted_mesh(g$x, weights = "pop")
  expect_true(cwm$ok)

  pieces <- shapeindices:::.array_mesh_to_sf_pieces(cwm$P, cwm$T, cwm$tri_area, cwm$crs)

  arr_db <- shapeindices:::.mesh_directional_balance_index_array(
    cwm$P, cwm$T, cwm$tri_area, cwm$tri_weight, cwm$crs)
  sf_db  <- shapeindices:::.mesh_directional_balance_index(pieces, weight = cwm$tri_weight)
  expect_equal(arr_db$index, sf_db$index, tolerance = 1e-8)

  arr_rc <- shapeindices:::.mesh_radial_concentration_index_array(
    cwm$P, cwm$T, cwm$tri_area, cwm$tri_weight, cwm$crs)
  sf_rc  <- shapeindices:::.mesh_radial_concentration_index(pieces, weight = cwm$tri_weight)
  expect_equal(arr_rc$index, sf_rc$index, tolerance = 1e-8)
})
