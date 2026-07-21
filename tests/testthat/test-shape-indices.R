test_that("shape_indices returns all thirteen named indices, sharing one triangulation", {
  sq <- wc(make_square())
  r <- shape_indices(sq)
  expect_named(r, c("convexity", "moment_of_inertia", "moment_isotropy", "directional_balance",
                     "span", "radial_concentration", "depth",
                     "hull_ratio", "polsby_popper", "width_length_ratio", "reock", "detour", "exchange"))
  expect_equal(unname(r["convexity"]), 1, tolerance = 1e-8)
  expect_equal(unname(r["hull_ratio"]), 1, tolerance = 1e-8)
  expect_equal(unname(r["moment_isotropy"]), 1, tolerance = 1e-8)
  expect_equal(unname(r["directional_balance"]), 1, tolerance = 1e-8)
  expect_equal(unname(r["detour"]), sqrt(pi) / 2, tolerance = 1e-8)
  expect_equal(unname(r["exchange"]), exchange_index(sq)$index, tolerance = 1e-8)
})

test_that("shape_indices matches calling the thirteen index functions directly", {
  star <- wc(make_star(6, 1, 0.4))
  r <- shape_indices(star)
  expect_equal(unname(r["convexity"]), convexity_index(star)$index, tolerance = 1e-8)
  expect_equal(unname(r["hull_ratio"]), hull_ratio_index(star)$index, tolerance = 1e-8)
  expect_equal(unname(r["polsby_popper"]), polsby_popper_index(star)$index, tolerance = 1e-8)
  expect_equal(unname(r["width_length_ratio"]), width_length_ratio_index(star)$index, tolerance = 1e-8)
  expect_equal(unname(r["reock"]), reock_index(star)$index, tolerance = 1e-8)
  expect_equal(unname(r["detour"]), detour_index(star)$index, tolerance = 1e-8)
  expect_equal(unname(r["exchange"]), exchange_index(star)$index, tolerance = 1e-8)
  # moment_of_inertia_index()/moment_isotropy_index()/directional_balance_index()/
  # span_index()/radial_concentration_index()/depth_index() re-triangulate
  # independently but on the same (deterministic) geometry, so this should
  # still match closely
  expect_equal(unname(r["moment_of_inertia"]), moment_of_inertia_index(star)$index, tolerance = 1e-6)
  expect_equal(unname(r["moment_isotropy"]), moment_isotropy_index(star)$index, tolerance = 1e-6)
  expect_equal(unname(r["directional_balance"]), directional_balance_index(star)$index, tolerance = 1e-6)
  expect_equal(unname(r["span"]), span_index(star)$index, tolerance = 1e-6)
  expect_equal(unname(r["radial_concentration"]), radial_concentration_index(star)$index, tolerance = 1e-6)
  expect_equal(unname(r["depth"]), depth_index(star)$index, tolerance = 1e-6)
})

test_that("shape_indices deterministic_max_tri switches to the random estimators above the threshold", {
  star <- wc(make_star(10, 1, 0.3))
  prep <- prepare_polygon(star)
  n_tri <- nrow(prep$tri)
  skip_if(n_tri < 2, "star mesh too small for this test")

  ci_det    <- convexity_index(star)$index
  span_det  <- span_index(star)$index
  rci_det   <- radial_concentration_index(star)$index
  depth_det <- depth_index(star)$index
  # small n_lines: forcing deterministic_max_tri just below n_tri switches
  # this small mesh onto the random estimators for all four mesh indices
  # that have one (convexity, span, radial_concentration, depth alike), and
  # a large n_lines relative to its (also small) deterministic pair count
  # would trip the respective "not substantially lower than" warnings.
  # One shared n_lines sets the Monte Carlo sample count for all four.
  r_forced <- shape_indices(star, deterministic_max_tri = n_tri - 1, n_lines = 50, seed = 1)
  expect_equal(unname(r_forced["convexity"]), ci_det, tolerance = 0.1)
  expect_equal(unname(r_forced["span"]), span_det, tolerance = 0.1)
  expect_equal(unname(r_forced["radial_concentration"]), rci_det, tolerance = 0.1)
  expect_equal(unname(r_forced["depth"]), depth_det, tolerance = 0.1)
})

test_that("shape_indices: an explicit deterministic = FALSE overrides deterministic_max_tri, not the other way round", {
  # deterministic_max_tri only ever pushes deterministic TRUE -> FALSE (its
  # guardrail purpose); an explicit deterministic = FALSE must win outright
  # even when deterministic_max_tri is generous enough that the size check
  # alone would have picked deterministic = TRUE. Detected via the
  # random-only "not substantially lower than" warning, same trick as
  # elsewhere in this file.
  sq <- wc(make_square())
  expect_warning(
    shape_indices(sq, deterministic = FALSE, deterministic_max_tri = 1000, n_lines = 50, seed = 1),
    "not substantially lower"
  )
  # sanity: deterministic_max_tri can still force TRUE -> FALSE when
  # deterministic isn't explicitly FALSE (the guardrail use case, unaffected
  # by this change)
  expect_warning(
    shape_indices(sq, deterministic_max_tri = 1, n_lines = 50, seed = 1),
    "not substantially lower"
  )
  # and plain deterministic = TRUE/default with a generous deterministic_max_tri stays on
  # the deterministic path (no random-estimator warning)
  expect_no_warning(shape_indices(sq, deterministic_max_tri = 1000))
})

test_that("shape_indices `weight` forwards consistently to convexity, IMI, moment isotropy, directional balance, span, radial_concentration, and depth", {
  # shape_indices() builds its own prep internally (that's the whole point
  # of it sharing one triangulation) and doesn't accept `prep` itself -
  # `weight` must line up with ITS internal triangulation, so derive it
  # from a separate prepare_polygon() call on the same geometry instead of
  # passing prep through
  star <- wc(make_star(6, 1, 0.4))
  w <- prepare_polygon(star)$tri$area
  r_wt    <- shape_indices(star, weight = w)
  r_plain <- shape_indices(star)
  expect_equal(unname(r_wt["convexity"]), unname(r_plain["convexity"]), tolerance = 1e-6)
  expect_equal(unname(r_wt["moment_of_inertia"]), unname(r_plain["moment_of_inertia"]), tolerance = 1e-6)
  expect_equal(unname(r_wt["moment_isotropy"]), unname(r_plain["moment_isotropy"]), tolerance = 1e-6)
  expect_equal(unname(r_wt["directional_balance"]), unname(r_plain["directional_balance"]), tolerance = 1e-6)
  expect_equal(unname(r_wt["span"]), unname(r_plain["span"]), tolerance = 1e-6)
  expect_equal(unname(r_wt["radial_concentration"]), unname(r_plain["radial_concentration"]), tolerance = 1e-6)
  expect_equal(unname(r_wt["depth"]), unname(r_plain["depth"]), tolerance = 1e-6)
})

test_that("shape_indices `n_lines` sets both convexity's and span's Monte Carlo sample count", {
  # convexity_index() and span_index() share the `n_lines` argument name
  # (even though span_index() never builds line geometry - see its own
  # file header for why) precisely so one shared value controls both,
  # rather than needing two separately-named arguments
  sq <- wc(make_square())
  r <- suppressWarnings(
    shape_indices(sq, deterministic = FALSE, n_lines = 500, seed = 1, which = c("convexity", "span"))
  )
  expect_equal(unname(r["convexity"]),
               suppressWarnings(convexity_index(sq, deterministic = FALSE, n_lines = 500, seed = 1)$index),
               tolerance = 1e-8)
  expect_equal(unname(r["span"]),
               suppressWarnings(span_index(sq, deterministic = FALSE, n_lines = 500, seed = 1)$index),
               tolerance = 1e-8)
})

## -- which -------------------------------------------------------------

test_that("shape_indices `which` returns exactly the requested subset, in canonical order", {
  sq <- wc(make_square())
  r <- shape_indices(sq, which = c("reock", "convexity", "hull_ratio"))
  # canonical order (convexity, ..., hull_ratio, ..., reock), not caller's order
  expect_named(r, c("convexity", "hull_ratio", "reock"))
})

test_that("shape_indices `which` matches the full result for the same indices", {
  star <- wc(make_star(6, 1, 0.4))
  full <- shape_indices(star)
  subset <- shape_indices(star, which = c("span", "polsby_popper"))
  expect_equal(unname(subset["span"]), unname(full["span"]), tolerance = 1e-6)
  expect_equal(unname(subset["polsby_popper"]), unname(full["polsby_popper"]), tolerance = 1e-8)
})

test_that("shape_indices `which` rejects unknown names with a clear error listing valid choices", {
  sq <- wc(make_square())
  expect_error(shape_indices(sq, which = "not_a_real_index"), "Unknown index")
  expect_error(shape_indices(sq, which = "not_a_real_index"), "reock")
})

test_that("shape_indices `which` naming only classic metrics skips triangulation entirely", {
  testthat::local_mocked_bindings(
    prepare_polygon = function(...) stop("prepare_polygon() should not be called"),
    .package = "shapeindices"
  )
  sq <- wc(make_square())
  expect_no_error(
    r <- shape_indices(sq, which = c("hull_ratio", "polsby_popper", "width_length_ratio", "reock", "detour", "exchange"))
  )
  expect_named(r, c("hull_ratio", "polsby_popper", "width_length_ratio", "reock", "detour", "exchange"))
})
