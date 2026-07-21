test_that("convex_decompose produces no more pieces than triangles, same total area", {
  star   <- wc(make_star(10, 1, 0.3))
  tri    <- cdt_triangles(star)
  pieces <- convex_decompose(star)
  expect_true(nrow(pieces) <= nrow(tri))
  expect_equal(sum(pieces$area), sum(tri$area), tolerance = 1e-6)
})

test_that("convex_decompose collapses an already-convex shape to a single piece", {
  sq <- wc(make_square())
  pieces <- convex_decompose(sq)
  expect_equal(nrow(pieces), 1)
})

test_that("convex_decompose accepts a precomputed prep to skip re-triangulating", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  pieces <- convex_decompose(star, prep = prep)
  expect_true(nrow(pieces) >= 1)
})
