test_that("moment_isotropy_index is ~1 for a disk", {
  disk <- make_disk(1, n_seg = 60)
  expect_equal(moment_isotropy_index(wc(disk))$index, 1, tolerance = 0.01)
})

test_that("moment_isotropy_index is exactly 1 for a square (isotropic, not a disk)", {
  expect_equal(moment_isotropy_index(wc(make_square()))$index, 1, tolerance = 1e-8)
})

test_that("moment_isotropy_index is ~1 for a regular hexagon (3-fold+ rotational symmetry)", {
  ang <- seq(0, 2 * pi, length.out = 7)[1:6]
  hex <- sf::st_polygon(list(rbind(cbind(cos(ang), sin(ang)), c(1, 0))))
  expect_equal(moment_isotropy_index(wc(hex))$index, 1, tolerance = 1e-8)
})

test_that("moment_isotropy_index falls as a shape elongates, matching the closed form for a rectangle", {
  square_ecc <- moment_isotropy_index(wc(make_square()))$index
  rect_ecc   <- moment_isotropy_index(wc(make_rectangle(20, 1)))$index
  expect_equal(square_ecc, 1, tolerance = 1e-8)
  expect_lt(rect_ecc, square_ecc)
  # closed form: for a WxH rectangle, lambda_min/lambda_max = min(H^2,W^2)/max(H^2,W^2)
  expect_equal(rect_ecc, 1^2 / 20^2, tolerance = 1e-8)
})

test_that("moment_isotropy_index is rotation invariant", {
  rect <- make_rectangle(20, 2)
  base <- moment_isotropy_index(wc(rect))$index
  theta <- 0.37
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  coords_rot <- sf::st_coordinates(rect)[, 1:2] %*% t(rot)
  rect_rot <- sf::st_polygon(list(coords_rot))
  expect_equal(moment_isotropy_index(wc(rect_rot))$index, base, tolerance = 1e-8)
})

test_that("moment_isotropy_index is scale invariant", {
  star <- make_star(7, 1, 0.4)
  base <- moment_isotropy_index(wc(star))$index
  expect_equal(moment_isotropy_index(wc(star * 10))$index, base, tolerance = 1e-6)
})

test_that("weight = triangle area reproduces the unweighted index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  ecc_plain <- moment_isotropy_index(star, prep = prep)$index
  ecc_wt    <- moment_isotropy_index(star, prep = prep, weight = prep$tri$area)$index
  expect_equal(ecc_wt, ecc_plain, tolerance = 1e-8)
})

test_that("weighting mass toward the far ends lowers moment isotropy; toward the centre raises it", {
  # a symmetric grid so any symmetric reweighting leaves the mass centroid
  # fixed (an asymmetric mesh would confound the comparison - the centroid
  # itself would move under reweighting, verified during development)
  cell <- function(cx, cy) sf::st_polygon(list(rbind(
    c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
    c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
  cells <- expand.grid(x = seq(-19.5, 19.5, by = 1), y = seq(-1.5, 1.5, by = 1))
  polys <- mapply(cell, cells$x, cells$y, SIMPLIFY = FALSE)
  rect  <- wc(sf::st_union(sf::st_sfc(polys)))
  prep  <- prepare_polygon(rect)
  tri   <- prep$tri
  cen   <- sf::st_coordinates(sf::st_centroid(sf::st_geometry(tri)))

  idx_plain <- moment_isotropy_index(rect, prep = prep)$index
  idx_ends  <- moment_isotropy_index(rect, prep = prep, weight = cen[, 1]^2 + 1)$index
  idx_mid   <- moment_isotropy_index(rect, prep = prep, weight = 1 / (cen[, 1]^2 + 1))$index

  expect_lt(idx_ends, idx_plain)
  expect_gt(idx_mid, idx_plain)
})

test_that("moment_isotropy_index shares Ixx/Iyy/Ixy/centroid with moment_of_inertia_index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  moi <- moment_of_inertia_index(star, prep = prep)
  ecc <- moment_isotropy_index(star, prep = prep)
  expect_equal(ecc$Ixx, moi$Ixx)
  expect_equal(ecc$Iyy, moi$Iyy)
  expect_equal(sf::st_coordinates(ecc$centroid), sf::st_coordinates(moi$centroid))
  expect_equal(ecc$lambda_min + ecc$lambda_max, moi$J, tolerance = 1e-8)  # trace = J
})

test_that("moment_isotropy_index warns and returns NA on an empty triangulation", {
  degenerate <- wc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(2, 0), c(0, 0)))))
  # two warnings fire here: prepare_polygon()'s own validity repair, then
  # moment_isotropy_index()'s own "no triangles" check - suppress the first
  # since only the second is what this test is actually checking
  suppressWarnings(expect_warning(res <- moment_isotropy_index(degenerate), "not defined"))
  expect_true(is.na(res$index))
})

test_that("moment_isotropy_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- moment_isotropy_index_sf(x))
  expect_equal(res$moment_isotropy_index[1], moment_isotropy_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$moment_isotropy_index[2], moment_isotropy_index(x[2, ])$index, tolerance = 1e-8)
})

test_that("shape_indices_sf(byrow = FALSE) computes moment isotropy consistently with the direct call", {
  cell <- function(cx, cy) sf::st_polygon(list(rbind(
    c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
    c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
  cells <- expand.grid(x = 1:6, y = 1:2)
  geoms <- mapply(cell, cells$x, cells$y, SIMPLIFY = FALSE)
  x <- sf::st_sf(pop = 1, geometry = sf::st_sfc(geoms, crs = TEST_CRS))

  res <- shape_indices_sf(x, byrow = FALSE, id = "g", which = c("moment_of_inertia", "moment_isotropy"))
  poly_u <- sf::st_union(sf::st_geometry(x))
  direct <- moment_isotropy_index(wc(poly_u[[1]]))$index
  expect_equal(res$moment_isotropy_index, direct, tolerance = 1e-8)
})
