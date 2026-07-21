test_that("shape_indices_sf byrow = TRUE appends ten columns matching shape_indices row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  # not testing parallel dispatch here (covered below) - parallel_rows = FALSE
  # avoids the "no active future::plan()" advisory warning
  res <- shape_indices_sf(x, parallel_rows = FALSE)
  expect_true(all(c("convexity_index", "moment_of_inertia_index", "moment_isotropy_index",
                     "directional_balance_index", "span_index",
                     "radial_concentration_index", "depth_index", "hull_ratio_index", "polsby_popper_index",
                     "width_length_ratio_index", "reock_index") %in% names(res)))
  expect_equal(res$convexity_index[1], unname(shape_indices(x[1, ])["convexity"]), tolerance = 1e-8)
  expect_equal(res$convexity_index[2], unname(shape_indices(x[2, ])["convexity"]), tolerance = 1e-8)
  expect_equal(res$span_index[1], unname(shape_indices(x[1, ])["span"]), tolerance = 1e-8)
  expect_equal(res$span_index[2], unname(shape_indices(x[2, ])["span"]), tolerance = 1e-8)
  expect_equal(res$radial_concentration_index[1],
               unname(shape_indices(x[1, ])["radial_concentration"]), tolerance = 1e-8)
  expect_equal(res$radial_concentration_index[2],
               unname(shape_indices(x[2, ])["radial_concentration"]), tolerance = 1e-8)
  expect_equal(res$depth_index[1], unname(shape_indices(x[1, ])["depth"]), tolerance = 1e-8)
  expect_equal(res$depth_index[2], unname(shape_indices(x[2, ])["depth"]), tolerance = 1e-8)
  expect_equal(res$reock_index[1], unname(shape_indices(x[1, ])["reock"]), tolerance = 1e-8)
  expect_equal(res$reock_index[2], unname(shape_indices(x[2, ])["reock"]), tolerance = 1e-8)
})

test_that("shape_indices_sf byrow = TRUE `which` appends only the requested columns", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  res <- shape_indices_sf(x, which = c("hull_ratio", "reock"), parallel_rows = FALSE)
  expect_true(all(c("hull_ratio_index", "reock_index") %in% names(res)))
  expect_false(any(c("convexity_index", "span_index", "moment_of_inertia_index",
                      "radial_concentration_index", "polsby_popper_index",
                      "width_length_ratio_index") %in% names(res)))
})

test_that("shape_indices_sf byrow = TRUE `which` naming only classic metrics skips triangulation", {
  testthat::local_mocked_bindings(
    prepare_polygon = function(...) stop("prepare_polygon() should not be called"),
    .package = "shapeindices"
  )
  x <- sf::st_sf(name = "sq", geometry = sf::st_sfc(make_square(), crs = TEST_CRS))
  expect_no_error(shape_indices_sf(x, which = c("hull_ratio", "polsby_popper"), parallel_rows = FALSE))
})

test_that("shape_indices_sf byrow = TRUE preserves CRS (no spurious 'no CRS' warning)", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  # parallel_rows = FALSE to isolate the CRS check from the (unrelated,
  # separately expected) "no active future::plan()" warning
  expect_no_warning(shape_indices_sf(x, parallel_rows = FALSE))
})

test_that("shape_indices_sf byrow = TRUE rejects weights/id", {
  x <- sf::st_sf(name = "sq", geometry = sf::st_sfc(make_square(), crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = TRUE, weights = 1), "byrow = FALSE")
  expect_error(shape_indices_sf(x, byrow = TRUE, id = "a"), "byrow = FALSE")
})

test_that("shape_indices_sf byrow = TRUE with deterministic_max_tri passes through to every row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(10, 1, 0.3), crs = TEST_CRS))
  # parallel_rows = FALSE: this test is about deterministic_max_tri's per-row
  # switch, not parallel dispatch (covered separately below) - row-level
  # futures sampling RNG (the random-line/random-pair paths) without an
  # explicit future.seed otherwise trips future's own "unreliable value"
  # safety warning, unrelated to what's being tested here. n_lines
  # are small since deterministic_max_tri = 1 forces even these tiny meshes onto
  # the random estimators - large values relative to the (also tiny)
  # deterministic pair count would trip convexity_index()'s/span_index()'s
  # own "not substantially lower than" warnings.
  res <- shape_indices_sf(x, deterministic_max_tri = 1, n_lines = 10, seed = 1,
                           parallel_rows = FALSE)
  expect_true(all(is.finite(res$convexity_index)))
  expect_true(all(is.finite(res$span_index)))
  expect_true(all(is.finite(res$radial_concentration_index)))
})

test_that("shape_indices_sf byrow = TRUE dispatches through furrr and matches sequential", {
  skip_if_not_installed("furrr")
  skip_if_not_installed("future")

  # future::sequential (the default plan) still exercises furrr::future_map()'s
  # dispatch machinery, just in-process - unlike multisession (separate R
  # sessions), it doesn't require shapeindices to be a properly installed
  # package findable by library() in a worker, so this works the same way
  # whether run via devtools::test()/load_all() or against an installed copy
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))

  res_seq <- shape_indices_sf(x, parallel_rows = FALSE)
  res_par <- suppressWarnings(shape_indices_sf(x, parallel_rows = TRUE))

  expect_equal(res_par$convexity_index, res_seq$convexity_index, tolerance = 1e-8)
})
