test_that(".normalize_weight sums to 1 and preserves relative proportions", {
  w <- shapeindices:::.normalize_weight(c(1, 2, 3, 4))
  expect_equal(sum(w), 1)
  expect_equal(w, c(1, 2, 3, 4) / 10)
})

test_that(".normalize_weight validates input", {
  expect_error(shapeindices:::.normalize_weight(factor(c("a", "b"))), "factor")
  expect_error(shapeindices:::.normalize_weight("a"), "numeric")
  expect_error(shapeindices:::.normalize_weight(c(1, NA)), "NA")
  expect_error(shapeindices:::.normalize_weight(c(1, Inf)), "finite")
  expect_error(shapeindices:::.normalize_weight(c(0, 0)), "positive")
  expect_error(shapeindices:::.normalize_weight(c(-1, -2)), "non-negative")
})

test_that(".normalize_weight rejects a negative entry even when the total sums positive", {
  # regression: c(-5, 10) sums to a positive total (5), so the old
  # sum(weight) <= 0 check alone let it straight through despite one
  # triangle having a nonsensical negative weight
  expect_error(shapeindices:::.normalize_weight(c(-5, 10)), "non-negative")
  expect_no_error(shapeindices:::.normalize_weight(c(0, 10)))
})

test_that("convexity_index rejects a weight vector with a negative entry", {
  sq <- wc(make_square())
  prep <- prepare_polygon(sq)
  n_tri <- nrow(prep$tri)
  w <- rep(10, n_tri); w[1] <- -5
  expect_error(convexity_index(sq, prep = prep, weight = w), "non-negative")
})

test_that(".ensure_projected passes already-planar data through unchanged", {
  sq <- sf::st_sfc(make_square(), crs = TEST_CRS)
  out <- shapeindices:::.ensure_projected(sq)
  expect_equal(sf::st_crs(out), sf::st_crs(sq))
})

test_that(".ensure_projected warns (not errors) when CRS is missing", {
  sq <- sf::st_sfc(make_square()) # no CRS
  expect_warning(out <- shapeindices:::.ensure_projected(sq), "no CRS")
  expect_equal(out, sq)
})

test_that(".ensure_projected auto-projects geographic input with a message", {
  ll <- sf::st_sfc(make_square(0.01, center = c(-79, 35)), crs = 4326)
  expect_message(out <- shapeindices:::.ensure_projected(ll), "auto-projecting")
  expect_false(sf::st_is_longlat(out))
})

test_that(".ensure_projected errors on empty geometry", {
  empty <- sf::st_sfc(sf::st_polygon())
  expect_error(shapeindices:::.ensure_projected(empty), "empty")
})

test_that(".make_valid_warn is silent for already-valid geometry", {
  ok <- sf::st_sfc(make_square(), crs = TEST_CRS)
  expect_no_warning(shapeindices:::.make_valid_warn(ok))
})

test_that(".make_valid_warn warns for genuinely invalid (self-intersecting) geometry", {
  bt <- sf::st_sfc(make_bowtie(), crs = TEST_CRS)
  expect_warning(shapeindices:::.make_valid_warn(bt), "not a valid, simple polygon")
})

test_that(".make_valid_warn pluralises its message correctly for multi-row input", {
  x <- sf::st_sfc(make_square(), make_bowtie(), make_bowtie(), crs = TEST_CRS)
  expect_warning(shapeindices:::.make_valid_warn(x), "2 of 3")
})

test_that("%||% picks the right-hand side only when the left is NULL", {
  `%||%` <- shapeindices:::`%||%`
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
  expect_equal(FALSE %||% 5, FALSE) # only NULL triggers the default, not falsy values
})
