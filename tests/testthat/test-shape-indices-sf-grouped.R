test_that("byrow = FALSE requires id", {
  x <- sf::st_sf(name = "a", geometry = sf::st_sfc(make_square(), crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = FALSE), "id")
})

test_that("byrow = FALSE requires a single-length id", {
  x <- sf::st_sf(name = "a", geometry = sf::st_sfc(make_square(), crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = FALSE, id = c("a", "b")), "single value")
})

test_that("byrow = FALSE with weights = NULL reproduces plain indices on the union", {
  sq1 <- make_square(1, center = c(-1, 0))
  sq2 <- make_square(1, center = c(1, 0)) # touching -> union is a solid 4x2 rectangle
  x <- sf::st_sf(name = c("left", "right"), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))

  res <- shape_indices_sf(x, byrow = FALSE, id = "grp")
  union_geom <- sf::st_union(sf::st_geometry(x))
  plain <- shape_indices(union_geom)

  expect_equal(res$convexity_index, unname(plain["convexity"]), tolerance = 1e-6)
  expect_equal(res$moment_of_inertia_index, unname(plain["moment_of_inertia"]), tolerance = 1e-6)
  # looser tolerance than convexity/MOI: unlike those two, span_index()'s
  # self-term and radial_concentration_index()'s point cloud are both
  # subdivision-based approximations (see their own file headers), so two
  # independently-built triangulations of the same union polygon can
  # differ at this order even though they're the same shape
  expect_equal(res$span_index, unname(plain["span"]), tolerance = 1e-4)
  expect_equal(res$radial_concentration_index, unname(plain["radial_concentration"]), tolerance = 1e-4)
  expect_equal(res$depth_index, unname(plain["depth"]), tolerance = 1e-4)
  expect_equal(res$hull_ratio_index, unname(plain["hull_ratio"]), tolerance = 1e-8)
  expect_equal(res$polsby_popper_index, unname(plain["polsby_popper"]), tolerance = 1e-8)
  expect_equal(res$width_length_ratio_index, unname(plain["width_length_ratio"]), tolerance = 1e-8)
  expect_equal(res$reock_index, unname(plain["reock"]), tolerance = 1e-6)
})

test_that("byrow = FALSE `which` returns only the requested columns", {
  sq1 <- make_square(1, center = c(-1, 0))
  sq2 <- make_square(1, center = c(1, 0))
  x <- sf::st_sf(name = c("left", "right"), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))

  res <- shape_indices_sf(x, byrow = FALSE, which = c("hull_ratio", "reock"), id = "grp")
  expect_true(all(c("hull_ratio_index", "reock_index", "total_weight") %in% names(res)))
  expect_false(any(c("convexity_index", "moment_of_inertia_index", "span_index",
                      "radial_concentration_index", "polsby_popper_index",
                      "width_length_ratio_index") %in% names(res)))
})

test_that("byrow = FALSE `which` naming only classic metrics skips triangulation entirely", {
  testthat::local_mocked_bindings(
    cdt_triangles = function(...) stop("cdt_triangles() should not be called"),
    .package = "shapeindices"
  )
  sq1 <- make_square(1, center = c(-1, 0))
  sq2 <- make_square(1, center = c(1, 0))
  x <- sf::st_sf(name = c("left", "right"), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))
  expect_no_error(
    shape_indices_sf(x, byrow = FALSE, which = c("hull_ratio", "polsby_popper", "reock"), id = "grp")
  )
})

test_that("byrow = FALSE moment_of_inertia_index uses the true mass centroid when weighted", {
  # regression (same underlying fix as moment_of_inertia_index() directly):
  # a heavily right-weighted row must pull the moment-of-inertia reference
  # point toward it, not leave it at the union's plain area centroid.
  # Cross-validated against a direct moment_of_inertia_index() call using
  # the SAME per-triangle weight (area-weighted overlay of the row weights,
  # matching what .weighted_mesh() does internally) - if shape_indices_sf()
  # silently kept using the plain area centroid, this would disagree with
  # the direct call, which is already regression-tested against the fix.
  sq1 <- make_square(1, center = c(-1, 0)) # left square, weight 1
  sq2 <- make_square(1, center = c(1, 0))  # right square, weight 100
  x <- sf::st_sf(name = c("left", "right"), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))

  # which narrowed to just what's being tested - avoids the unrelated
  # "hull_ratio/width_length_ratio don't use weights" warning `which = "all"`
  # would otherwise trigger here, and skips computing the other 7 indices
  res_area <- shape_indices_sf(x, byrow = FALSE, which = "moment_of_inertia", id = "by_area")
  res_wt   <- shape_indices_sf(x, byrow = FALSE, which = "moment_of_inertia",
                                weights = c(1, 100), id = "by_weight")

  # heavily right-weighted MOI must differ from the unweighted case - if the
  # reference point silently stayed at the area centroid (the bug), these
  # would be numerically identical instead
  expect_false(isTRUE(all.equal(res_area$moment_of_inertia_index, res_wt$moment_of_inertia_index)))

  # cross-check: the combined mesh's per-triangle weight is just each row's
  # weight spread uniformly across its own triangles here (no straddling,
  # since the mesh is built from the union and each triangle falls cleanly
  # within one square) - so a direct call using that same per-triangle
  # weight should match shape_indices_sf()'s result exactly
  union_geom <- sf::st_union(sf::st_geometry(x))
  prep <- prepare_polygon(union_geom)
  tri_cx <- vapply(sf::st_geometry(prep$tri), function(t) mean(sf::st_coordinates(t)[1:3, 1]), numeric(1))
  tri_w <- ifelse(tri_cx < 0, 1, 100)
  direct <- moment_of_inertia_index(union_geom, prep = prep, weight = tri_w)
  expect_equal(res_wt$moment_of_inertia_index, direct$index, tolerance = 1e-6)
})

test_that("byrow = FALSE treats weight 0 and weight NA identically, as holes", {
  left   <- make_square(1, center = c(-2, 0))
  middle <- make_square(1, center = c(0, 0))
  right  <- make_square(1, center = c(2, 0))
  x <- sf::st_sf(name = c("left", "middle", "right"),
                  geometry = sf::st_sfc(left, middle, right, crs = TEST_CRS))

  # each call also emits the (harmless, expected) "hull_ratio/width_length_ratio
  # don't use weights" warning alongside the hole one, since `which = "all"`
  # includes hull_ratio and a nonzero weight is present here too -
  # expect_warning() only requires a match, not an exclusive one
  expect_warning(r_zero <- shape_indices_sf(x, byrow = FALSE, weights = c(10, 0, 10), id = "z"),
                 "zero or NA weight")
  expect_warning(r_na <- shape_indices_sf(x, byrow = FALSE, weights = c(10, NA, 10), id = "n"),
                 "zero or NA weight")

  # weight 0 and weight NA must give IDENTICAL results across every index,
  # not just convexity/MOI - both are resolved to the same "hole" before
  # any index ever sees the mesh, so nothing downstream can tell them apart
  expect_equal(r_zero$convexity_index, r_na$convexity_index, tolerance = 1e-8)
  expect_equal(r_zero$moment_of_inertia_index, r_na$moment_of_inertia_index, tolerance = 1e-8)
  expect_equal(r_zero$span_index, r_na$span_index, tolerance = 1e-8)
  expect_equal(r_zero$radial_concentration_index, r_na$radial_concentration_index, tolerance = 1e-8)
  expect_equal(r_zero$depth_index, r_na$depth_index, tolerance = 1e-8)
  expect_equal(r_zero$hull_ratio_index, r_na$hull_ratio_index, tolerance = 1e-8)
  expect_equal(r_zero$total_weight, r_na$total_weight, tolerance = 1e-8)

  # excluding the middle square should leave a real gap, i.e. a dispersed
  # (non-convex) union of just left + right, strictly less convex than the
  # solid, gapless union of all three
  r_all <- shape_indices_sf(x, byrow = FALSE, id = "all")
  expect_lt(r_zero$convexity_index, r_all$convexity_index)

  expected_area <- as.numeric(sf::st_area(sf::st_union(sf::st_geometry(x)[c(1, 3)])))
  expect_equal(as.numeric(sf::st_area(r_zero)), expected_area, tolerance = 1e-6)
})

test_that("byrow = FALSE errors cleanly when every row is a hole", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-2, 0)), make_square(center = c(2, 0)),
                                         crs = TEST_CRS))
  expect_error(suppressWarnings(shape_indices_sf(x, byrow = FALSE, weights = c(0, 0), id = "z")),
               "nothing left to triangulate")
})

test_that("byrow = FALSE total_weight sums resolved weights, excluding holes", {
  x <- sf::st_sf(name = c("a", "b", "c"),
                  geometry = sf::st_sfc(make_square(center = c(-2, 0)), make_square(center = c(0, 0)),
                                         make_square(center = c(2, 0)), crs = TEST_CRS))
  res <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, weights = c(5, NA, 7), id = "z"))
  expect_equal(res$total_weight, 12)
})

test_that("byrow = FALSE weights column-name lookup works and rejects unknown columns", {
  x <- sf::st_sf(name = c("a", "b"), pop = c(5, 15),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  # which narrowed to avoid the unrelated hull_ratio/width_length_ratio
  # weight-magnitude warning - this test is about total_weight, not those
  res <- shape_indices_sf(x, byrow = FALSE, which = "convexity", weights = "pop", id = "z")
  expect_equal(res$total_weight, 20)
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = "nope", id = "z"), "not a column")
})

test_that("byrow = FALSE weights naming the geometry column itself gives a specific error", {
  # regression: "geometry" IS technically in names(x), so the naive check
  # let it through, then st_drop_geometry(x)[["geometry"]] silently
  # returned NULL and failed later with a confusing "must be numeric, got
  # NULL" instead of a clear, specific error
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = "geometry", id = "z"),
               "geometry column")
})

test_that("byrow = FALSE weights rejects non-numeric input with a graceful, specific error", {
  x <- sf::st_sf(name = c("a", "b", "c"),
                  geometry = sf::st_sfc(make_square(center = c(-2, 0)), make_square(center = c(0, 0)),
                                         make_square(center = c(2, 0)), crs = TEST_CRS))

  # direct vectors
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = c("a", "b", "c"), id = "z"),
               "must be numeric")
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = c(TRUE, FALSE, TRUE), id = "z"),
               "must be numeric")
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = factor(c("a", "b", "c")), id = "z"),
               "is a factor")

  # column-name lookups
  x$char_col <- c("a", "b", "c")
  x$fac_col  <- factor(c("a", "b", "c"))
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = "char_col", id = "z"), "must be numeric")

  # regression: the factor column-lookup error used to suggest invalid R
  # code (as.numeric(as.character(column "fac_col"))) - the fix-it snippet
  # must now be a real, syntactically valid accessor into x, not the
  # human-readable "column \"fac_col\"" description
  err <- tryCatch(shape_indices_sf(x, byrow = FALSE, weights = "fac_col", id = "z"),
                   error = function(e) conditionMessage(e))
  expect_match(err, "is a factor")
  expect_match(err, 'x\\[\\["fac_col"\\]\\]')
})

test_that("byrow = FALSE hull_ratio is unaffected by weights, and warns to say so", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  r_area <- shape_indices_sf(x, byrow = FALSE, id = "a")
  expect_warning(
    r_weight <- shape_indices_sf(x, byrow = FALSE, weights = c(1, 100), id = "b"),
    "don't use `weights`"
  )
  expect_equal(r_area$hull_ratio_index, r_weight$hull_ratio_index, tolerance = 1e-8)
})

test_that("byrow = FALSE weight-magnitude warning fires for any of the six classic metrics, not the CDT-based ones", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  for (idx in c("hull_ratio", "polsby_popper", "width_length_ratio", "reock", "detour", "exchange")) {
    expect_warning(
      shape_indices_sf(x, byrow = FALSE, which = idx, weights = c(1, 100), id = "z"),
      "don't use `weights`"
    )
  }
  # narrowing `which` to only the CDT-based indices keeps it quiet
  expect_no_warning(
    shape_indices_sf(x, byrow = FALSE,
                      which = c("convexity", "moment_of_inertia", "span", "radial_concentration"),
                      weights = c(1, 100), id = "z")
  )
})

test_that("byrow = FALSE does not warn about overlap for touching, non-overlapping rows", {
  sq1 <- make_square(1, center = c(-1, 0))
  sq2 <- make_square(1, center = c(1, 0))
  x <- sf::st_sf(name = c("a", "b"), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))
  expect_no_warning(shape_indices_sf(x, byrow = FALSE, id = "z"))
})

test_that("byrow = FALSE stays well-defined when rows overlap (no warning, correct union geometry)", {
  # triangulating the union directly (rather than per-row triangles pooled
  # together) means an overlapping pair of rows no longer risks producing
  # literally overlapping triangles - st_union() already resolves the
  # overlap correctly, so this is silent and geometrically sound, just with
  # weight double-counted in the overlap region (see .weighted_mesh())
  sq1 <- make_square(2, center = c(0, 0))
  sq2 <- make_square(2, center = c(1, 0)) # overlaps sq1
  x <- sf::st_sf(name = c("a", "b"), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))

  expect_no_warning(res <- shape_indices_sf(x, byrow = FALSE, id = "z"))
  expect_true(is.finite(res$convexity_index))
  expect_equal(as.numeric(sf::st_area(res)),
               as.numeric(sf::st_area(sf::st_union(sf::st_geometry(x)))),
               tolerance = 1e-6)
})

test_that("byrow = FALSE, parallel_rows = TRUE is rejected", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = FALSE, id = "z", parallel_rows = TRUE),
               "no meaning when byrow = FALSE")
  # leaving it at the default (unset) is fine - only an explicit TRUE errors
  expect_no_error(shape_indices_sf(x, byrow = FALSE, id = "z"))
  expect_no_error(shape_indices_sf(x, byrow = FALSE, id = "z", parallel_rows = FALSE))
})

test_that("byrow = FALSE weights vector rejects negative and non-finite values", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = c(-1, 2), id = "z"), "non-negative")
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = c(Inf, 2), id = "z"), "finite")
})

test_that("byrow = FALSE weights vector must have one entry per row", {
  x <- sf::st_sf(name = c("a", "b", "c"),
                  geometry = sf::st_sfc(make_square(center = c(-2, 0)), make_square(center = c(0, 0)),
                                         make_square(center = c(2, 0)), crs = TEST_CRS))
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = c(1, 2), id = "z"),
               "one entry per row.*\\(3\\).*got 2")
  expect_error(shape_indices_sf(x, byrow = FALSE, weights = c(1, 2, 3, 4), id = "z"),
               "one entry per row.*\\(3\\).*got 4")
  expect_no_error(suppressWarnings(shape_indices_sf(x, byrow = FALSE, weights = c(1, 2, 3), id = "z")))
})

test_that("byrow = FALSE deterministic_max_tri switches the combined mesh to the random-line estimator", {
  # enough rows/complexity to get a mesh worth switching away from
  squares <- lapply(0:5, function(i) make_square(1, center = c(i * 2, 0)))
  x <- sf::st_sf(name = paste0("s", seq_along(squares)),
                  geometry = sf::st_sfc(squares, crs = TEST_CRS))
  # small n_lines: deterministic_max_tri = 1 forces even this small mesh
  # onto the random estimators, and large values relative to its (also
  # small) deterministic pair count would trip convexity_index()'s/span_index()'s
  # own "not substantially lower than" warnings
  res <- shape_indices_sf(x, byrow = FALSE, id = "z", deterministic_max_tri = 1,
                           n_lines = 10, seed = 1)
  expect_true(is.finite(res$convexity_index))
  expect_true(is.finite(res$span_index))
  # unaffected by deterministic_max_tri (no deterministic/Monte Carlo split),
  # but still computed correctly alongside the switched indices
  expect_true(is.finite(res$radial_concentration_index))
  expect_true(is.finite(res$depth_index))
})

test_that("byrow = FALSE respects an explicit deterministic = FALSE (regression: used to be silently ignored)", {
  # .shape_indices_sf_grouped() used to only ever consult deterministic_max_tri -
  # a directly-passed `deterministic` fell into `...` and was never read, so
  # deterministic = FALSE (and the n_lines that goes with it) silently had no
  # effect and the exhaustive method always ran regardless. Detect the
  # random-line path actually running via its own "n_lines vs deterministic
  # pairs" warning, which only fires from that code path.
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  expect_warning(
    shape_indices_sf(x, byrow = FALSE, id = "z", deterministic = FALSE, n_lines = 100, seed = 1),
    "not substantially lower"
  )
  # deterministic = TRUE (or unset) must NOT trip that random-line-only warning
  expect_no_warning(shape_indices_sf(x, byrow = FALSE, id = "z", deterministic = TRUE))
  expect_no_warning(shape_indices_sf(x, byrow = FALSE, id = "z"))
})

test_that("byrow = FALSE deterministic_max_tri can force deterministic = TRUE down to FALSE (the guardrail case)", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  # deterministic = TRUE, but deterministic_max_tri = 1 forces the small mesh onto the
  # random-line path anyway - same detection trick as above
  expect_warning(
    shape_indices_sf(x, byrow = FALSE, id = "z", deterministic = TRUE, deterministic_max_tri = 1,
                      n_lines = 100, seed = 1),
    "not substantially lower"
  )
})

test_that("byrow = FALSE: an explicit deterministic = FALSE overrides deterministic_max_tri, not the other way round", {
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  # deterministic = FALSE with a generous deterministic_max_tri (would pick deterministic = TRUE on
  # size alone) must still take the random-line path - deterministic_max_tri only
  # ever pushes TRUE -> FALSE, never FALSE -> TRUE
  expect_warning(
    shape_indices_sf(x, byrow = FALSE, id = "z", deterministic = FALSE, deterministic_max_tri = 1000,
                      n_lines = 100, seed = 1),
    "not substantially lower"
  )
})

## -- exact constrained-triangulation mesh dispatch -----------------------
##
## weights = NULL never triggers the exact mesh: .resolve_row_weights()
## makes every row's weight equal to its own area, so density is uniform
## everywhere and .weighted_mesh()'s coarse overlay has nothing to smear -
## see constrained-triangulation.R's own file header.

test_that("byrow = FALSE with weights = NULL never attempts the exact constrained mesh", {
  testthat::local_mocked_bindings(
    .constrained_weighted_mesh = function(...) stop(".constrained_weighted_mesh() should not be called"),
    .package = "shapeindices"
  )
  x <- sf::st_sf(name = c("a", "b"),
                  geometry = sf::st_sfc(make_square(center = c(-1, 0)), make_square(center = c(1, 0)),
                                         crs = TEST_CRS))
  expect_no_error(shape_indices_sf(x, byrow = FALSE, id = "z"))
})

test_that("byrow = FALSE with real weights uses the exact constrained mesh and recovers total weight exactly", {
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000),
                  geometry = sf::st_sfc(squares, crs = TEST_CRS))
  res <- shape_indices_sf(x, byrow = FALSE, id = "z", which = "moment_of_inertia", weights = "pop")
  expect_equal(res$total_weight, sum(x$pop))
})

test_that("byrow = FALSE falls back to the coarse mesh (with a warning) when weighted rows overlap", {
  sq1 <- make_square(2, center = c(0, 0))
  sq2 <- make_square(2, center = c(1, 0)) # overlaps sq1
  x <- sf::st_sf(pop = c(10, 1000), geometry = sf::st_sfc(sq1, sq2, crs = TEST_CRS))

  expect_warning(
    res <- shape_indices_sf(x, byrow = FALSE, id = "z", which = "moment_of_inertia", weights = "pop"),
    "Falling back to the coarse"
  )
  expect_true(is.finite(res$moment_of_inertia_index))
})

test_that("byrow = FALSE convexity/span silently fall back to Monte Carlo above the memory-safe ceiling at default settings", {
  testthat::local_mocked_bindings(
    .safe_deterministic_tri_ceiling = function(...) 1L,
    .package = "shapeindices"
  )
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000), geometry = sf::st_sfc(squares, crs = TEST_CRS))

  # deterministic_max_tri left at its default (NULL) - must not error, and
  # must not go anywhere near .mesh_convexity_index()/.mesh_span_index()
  # (the O(n^2) path), while moment_of_inertia stays on its own
  # (unceilinged, O(n)) deterministic path
  testthat::local_mocked_bindings(
    .mesh_convexity_index = function(...) stop(".mesh_convexity_index() should not be called above the ceiling"),
    .mesh_span_index = function(...) stop(".mesh_span_index() should not be called above the ceiling"),
    .package = "shapeindices"
  )
  res <- expect_no_error(
    shape_indices_sf(x, byrow = FALSE, id = "z", weights = "pop",
                      which = c("convexity", "span", "moment_of_inertia"), n_lines = 20, seed = 1)
  )
  expect_true(is.finite(res$convexity_index))
  expect_true(is.finite(res$span_index))
  expect_true(is.finite(res$moment_of_inertia_index))
})

test_that("byrow = FALSE convexity/span fail gracefully when an explicit deterministic_max_tri exceeds the memory-safe ceiling", {
  testthat::local_mocked_bindings(
    .safe_deterministic_tri_ceiling = function(...) 1L,
    .package = "shapeindices"
  )
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000), geometry = sf::st_sfc(squares, crs = TEST_CRS))

  expect_error(
    shape_indices_sf(x, byrow = FALSE, id = "z", weights = "pop",
                      which = "convexity", deterministic_max_tri = 1000),
    "memory-safe ceiling"
  )
})

test_that("byrow = FALSE convexity/span skip the memory-safe ceiling entirely under an explicit deterministic = FALSE", {
  testthat::local_mocked_bindings(
    .safe_deterministic_tri_ceiling = function(...) 1L,
    .package = "shapeindices"
  )
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000), geometry = sf::st_sfc(squares, crs = TEST_CRS))

  res <- expect_no_error(
    shape_indices_sf(x, byrow = FALSE, id = "z", weights = "pop",
                      which = "convexity", deterministic = FALSE, n_lines = 20, seed = 1)
  )
  expect_true(is.finite(res$convexity_index))
})

test_that("byrow = FALSE convexity/span deterministic_max_tri below n_tri gives ordinary Monte Carlo, no error, even with a generous ceiling", {
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000), geometry = sf::st_sfc(squares, crs = TEST_CRS))

  res <- expect_no_error(
    shape_indices_sf(x, byrow = FALSE, id = "z", weights = "pop",
                      which = "convexity", deterministic_max_tri = 1, n_lines = 20, seed = 1)
  )
  expect_true(is.finite(res$convexity_index))
})

test_that("byrow = FALSE exact-mesh indices are scale-invariant under a constant weight multiplier", {
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000), geometry = sf::st_sfc(squares, crs = TEST_CRS))
  x_scaled <- x
  x_scaled$pop <- x_scaled$pop * 1e6

  which_cdt <- c("convexity", "moment_of_inertia", "moment_isotropy", "directional_balance",
                 "span", "radial_concentration")
  r1 <- shape_indices_sf(x, byrow = FALSE, id = "z", which = which_cdt, weights = "pop", seed = 1)
  r2 <- shape_indices_sf(x_scaled, byrow = FALSE, id = "z", which = which_cdt, weights = "pop", seed = 1)

  expect_equal(r1$convexity_index, r2$convexity_index, tolerance = 1e-8)
  expect_equal(r1$moment_of_inertia_index, r2$moment_of_inertia_index, tolerance = 1e-8)
  expect_equal(r1$moment_isotropy_index, r2$moment_isotropy_index, tolerance = 1e-8)
  expect_equal(r1$directional_balance_index, r2$directional_balance_index, tolerance = 1e-8)
  expect_equal(r1$span_index, r2$span_index, tolerance = 1e-8)
  expect_equal(r1$radial_concentration_index, r2$radial_concentration_index, tolerance = 1e-8)
})

test_that("byrow = FALSE exact-mesh moment_of_inertia_index matches a direct call built from the same array mesh", {
  squares <- lapply(0:3, function(i) make_square(0.5, center = c(i, 0)))
  x <- sf::st_sf(pop = c(5, 50, 500, 5000), geometry = sf::st_sfc(squares, crs = TEST_CRS))

  res <- shape_indices_sf(x, byrow = FALSE, id = "z", which = "moment_of_inertia", weights = "pop")

  cwm <- shapeindices:::.constrained_weighted_mesh(x, weights = "pop")
  expect_true(cwm$ok)
  pieces <- shapeindices:::.array_mesh_to_sf_pieces(cwm$P, cwm$T, cwm$tri_area, cwm$crs)
  direct <- moment_of_inertia_index(cwm$poly_u, prep = list(poly = cwm$poly_u, tri = pieces), weight = cwm$tri_weight)

  expect_equal(res$moment_of_inertia_index, direct$index, tolerance = 1e-8)
})
