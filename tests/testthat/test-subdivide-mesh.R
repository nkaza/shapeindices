test_that("subdivide_mesh produces more pieces than the original CDT, same total area", {
  star  <- wc(make_star(10, 1, 0.3))
  tri   <- cdt_triangles(star)
  fine  <- subdivide_mesh(star)
  expect_true(nrow(fine) >= nrow(tri))
  expect_equal(sum(fine$area), sum(tri$area), tolerance = 1e-6)
})

test_that("subdivide_mesh accepts a precomputed prep to skip re-triangulating", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  fine <- subdivide_mesh(star, prep = prep)
  expect_true(nrow(fine) >= nrow(prep$tri))
})

test_that("subdivide_mesh with max_depth = 0 returns the CDT mesh unchanged", {
  star <- wc(make_star(8, 1, 0.5))
  tri  <- cdt_triangles(star)
  fine <- subdivide_mesh(star, max_depth = 0)
  expect_equal(nrow(fine), nrow(tri))
  expect_equal(sum(fine$area), sum(tri$area), tolerance = 1e-6)
})

test_that("subdivide_mesh gives every triangle at least 4^max_depth children in a uniform mesh", {
  # a single square triangulates to 2 congruent triangles - no size variation,
  # so adaptive depth degenerates to the fixed-depth case for both
  sq   <- wc(make_square())
  tri  <- cdt_triangles(sq)
  fine <- subdivide_mesh(sq, max_depth = 2)
  expect_equal(nrow(fine), nrow(tri) * 4^2)
})

test_that("subdivide_mesh gives smaller triangles proportionally fewer children than larger ones", {
  # one big square plus one tiny sliver triangle far smaller in area
  big   <- rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))
  spike <- rbind(c(5, 10), c(5.05, 10), c(5.025, 10.1), c(5, 10))
  poly  <- wc(sf::st_union(sf::st_sfc(sf::st_polygon(list(big)), sf::st_polygon(list(spike)))))
  tri   <- cdt_triangles(poly)
  fine  <- subdivide_mesh(poly, max_depth = 4)

  # every source triangle contributes at least 1 child, the largest gets 4^4
  expect_true(nrow(fine) >= nrow(tri))
  expect_true(nrow(fine) < nrow(tri) * 4^4)
  expect_equal(sum(fine$area), sum(tri$area), tolerance = 1e-6)
})

test_that("subdivide_mesh returns NULL when the polygon has no triangulatable area", {
  degenerate <- wc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(2, 0), c(0, 0)))))
  expect_warning(result <- subdivide_mesh(degenerate), "not a valid, simple polygon")
  expect_null(result)
})
