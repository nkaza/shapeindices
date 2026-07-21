test_that("hull_ratio_index is 1 for convex shapes", {
  expect_equal(hull_ratio_index(wc(make_square()))$index, 1, tolerance = 1e-8)
  expect_equal(hull_ratio_index(wc(make_disk(1, n_seg = 12)))$index, 1, tolerance = 1e-8)
})

test_that("hull_ratio_index is < 1 for a non-convex star", {
  expect_lt(hull_ratio_index(wc(make_star(6, 1, 0.3)))$index, 1)
})

test_that("hull_ratio_index is scale and translation invariant", {
  star <- make_star(7, 1, 0.4)
  hr_base <- hull_ratio_index(wc(star))$index
  expect_equal(hull_ratio_index(wc(star * 5))$index, hr_base, tolerance = 1e-8)
  expect_equal(hull_ratio_index(wc(star + c(50, -20)))$index, hr_base, tolerance = 1e-8)
})

test_that("hull_ratio_index warns on a non-simple polygon, same as prepare_polygon", {
  bt <- sf::st_sfc(make_bowtie(), crs = TEST_CRS)
  expect_warning(hull_ratio_index(bt), "not a valid, simple polygon")
})

test_that("hull_ratio_index_sf preserves CRS and matches row-by-row hull_ratio_index", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- hull_ratio_index_sf(x))
  expect_equal(res$hull_ratio_index[1], hull_ratio_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$hull_ratio_index[2], hull_ratio_index(x[2, ])$index, tolerance = 1e-8)
})

## -- polsby_popper_index ------------------------------------------------------

test_that("polsby_popper_index is 1 for a circle and falls for non-circular shapes", {
  expect_equal(polsby_popper_index(wc(make_disk(1, n_seg = 90)))$index, 1, tolerance = 1e-4)
  expect_lt(polsby_popper_index(wc(make_square()))$index, 1)
  expect_lt(polsby_popper_index(wc(make_rectangle(20, 1)))$index, polsby_popper_index(wc(make_square()))$index)
})

test_that("polsby_popper_index matches the exact known value for a square", {
  # 4*pi*s^2 / (4s)^2 = pi/4
  expect_equal(polsby_popper_index(wc(make_square(5)))$index, pi / 4, tolerance = 1e-8)
})

test_that("polsby_popper_index is scale and translation invariant", {
  star <- make_star(7, 1, 0.4)
  base <- polsby_popper_index(wc(star))$index
  expect_equal(polsby_popper_index(wc(star * 5))$index, base, tolerance = 1e-8)
  expect_equal(polsby_popper_index(wc(star + c(50, -20)))$index, base, tolerance = 1e-8)
})

test_that("polsby_popper_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- polsby_popper_index_sf(x))
  expect_equal(res$polsby_popper_index[1], polsby_popper_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$polsby_popper_index[2], polsby_popper_index(x[2, ])$index, tolerance = 1e-8)
})

## -- width_length_ratio_index --------------------------------------------------

test_that("width_length_ratio_index is 1 for a square and matches the known ratio for a rectangle", {
  expect_equal(width_length_ratio_index(wc(make_square()))$index, 1, tolerance = 1e-8)
  expect_equal(width_length_ratio_index(wc(make_rectangle(20, 1)))$index, 1 / 20, tolerance = 1e-8)
})

test_that("width_length_ratio_index is invariant to which axis is longer", {
  expect_equal(width_length_ratio_index(wc(make_rectangle(20, 1)))$index,
               width_length_ratio_index(wc(make_rectangle(1, 20)))$index, tolerance = 1e-8)
})

test_that("width_length_ratio_index is invariant to rotation (minimum rotated bounding rectangle, not axis-aligned)", {
  rect <- make_rectangle(2, 1)
  rotate <- function(g, deg) {
    th <- deg * pi / 180
    rot <- matrix(c(cos(th), sin(th), -sin(th), cos(th)), 2, 2)
    (g - sf::st_centroid(g)) * rot + sf::st_centroid(g)
  }
  base_index <- width_length_ratio_index(wc(rect))$index
  for (deg in c(15, 30, 45, 60, 75, 90)) {
    rotated <- wc(rotate(sf::st_sfc(rect), deg)[[1]])
    expect_equal(width_length_ratio_index(rotated)$index, base_index, tolerance = 1e-6,
                 info = sprintf("rotation = %d degrees", deg))
  }
})

test_that("width_length_ratio_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "rect"),
                  geometry = sf::st_sfc(make_square(), make_rectangle(20, 1), crs = TEST_CRS))
  expect_no_warning(res <- width_length_ratio_index_sf(x))
  expect_equal(res$width_length_ratio_index[1], width_length_ratio_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$width_length_ratio_index[2], width_length_ratio_index(x[2, ])$index, tolerance = 1e-8)
})

## -- reock_index and the minimum enclosing circle -----------------------------

test_that(".min_enclosing_circle matches the known exact value for a square", {
  # circumcircle through all 4 corners: radius = side * sqrt(2) / 2
  pts <- rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10))
  mec <- shapeindices:::.min_enclosing_circle(pts)
  expect_equal(unname(mec$center), c(5, 5), tolerance = 1e-8)
  expect_equal(mec$r, 10 * sqrt(2) / 2, tolerance = 1e-8)
})

test_that(".min_enclosing_circle matches the known exact value for an equilateral triangle", {
  s <- 10
  pts <- rbind(c(0, 0), c(s, 0), c(s / 2, s * sqrt(3) / 2))
  mec <- shapeindices:::.min_enclosing_circle(pts)
  expect_equal(mec$r, s / sqrt(3), tolerance = 1e-6)
})

test_that(".min_enclosing_circle handles an obtuse triangle (diameter-circle case)", {
  # a very flat, obtuse triangle: MEC is the circle over its two farthest
  # points, NOT the (larger) circumcircle - regression guard for the
  # "try each pair as diameter first" logic in .circle_from_3()
  pts <- rbind(c(0, 0), c(10, 0), c(5, 0.5))
  mec <- shapeindices:::.min_enclosing_circle(pts)
  expect_equal(unname(mec$center), c(5, 0), tolerance = 1e-6)
  expect_equal(mec$r, 5, tolerance = 1e-6)
  # every point must actually be inside (or on) the returned circle
  d <- sqrt((pts[, 1] - mec$center[1])^2 + (pts[, 2] - mec$center[2])^2)
  expect_true(all(d <= mec$r + 1e-6))
})

test_that("reock_index matches the known exact value for a square and an equilateral triangle", {
  expect_equal(reock_index(wc(make_square(5)))$index, 2 / pi, tolerance = 1e-6)
  s <- 10
  tri <- sf::st_polygon(list(rbind(c(0, 0), c(s, 0), c(s / 2, s * sqrt(3) / 2), c(0, 0))))
  expect_equal(reock_index(wc(tri))$index, 3 * sqrt(3) / (4 * pi), tolerance = 1e-6)
})

test_that("reock_index is ~1 for a circle and scale/translation invariant", {
  expect_equal(reock_index(wc(make_disk(1, n_seg = 90)))$index, 1, tolerance = 1e-3)
  star <- make_star(7, 1, 0.4)
  base <- reock_index(wc(star))$index
  expect_equal(reock_index(wc(star * 5))$index, base, tolerance = 1e-6)
  expect_equal(reock_index(wc(star + c(50, -20)))$index, base, tolerance = 1e-6)
})

test_that("reock_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- reock_index_sf(x))
  expect_equal(res$reock_index[1], reock_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$reock_index[2], reock_index(x[2, ])$index, tolerance = 1e-8)
})

## -- detour_index ---------------------------------------------------------------

test_that("detour_index is ~1 for a circle and matches the known exact value for a square", {
  expect_equal(detour_index(wc(make_disk(1, n_seg = 90)))$index, 1, tolerance = 1e-3)
  # square: hull = itself, perimeter 4s; circle perimeter 2*sqrt(pi*s^2) = 2s*sqrt(pi)
  # ratio = 2s*sqrt(pi) / 4s = sqrt(pi)/2, independent of s (scale-invariant)
  expect_equal(detour_index(wc(make_square()))$index, sqrt(pi) / 2, tolerance = 1e-8)
})

test_that("detour_index matches the known exact value for an equilateral triangle", {
  # convex, so hull = itself: perimeter 3*s, area sqrt(3)/4*s^2
  s <- 10
  tri <- sf::st_polygon(list(rbind(c(0, 0), c(s, 0), c(s / 2, s * sqrt(3) / 2), c(0, 0))))
  area <- sqrt(3) / 4 * s^2
  expect_equal(detour_index(wc(tri))$index, 2 * sqrt(pi * area) / (3 * s), tolerance = 1e-6)
})

test_that("detour_index is <= 1 and falls for elongated shapes", {
  expect_lt(detour_index(wc(make_rectangle(20, 1)))$index, detour_index(wc(make_square()))$index)
  expect_lte(detour_index(wc(make_star(6, 1, 0.3)))$index, 1)
})

test_that("detour_index's hull is untouched by a hole - only the shrunken area moves the index", {
  solid <- wc(make_square(5))
  holed <- wc(make_square_with_hole(5, 0.3))
  # the hole is fully interior, so the convex hull (and its perimeter) is
  # identical to the solid square's - only `area` (and so circle_perimeter)
  # shrinks by the hole's own area fraction
  expect_equal(detour_index(holed)$hull_perimeter, detour_index(solid)$hull_perimeter, tolerance = 1e-8)
  expect_equal(detour_index(holed)$index, detour_index(solid)$index * sqrt(1 - 0.3), tolerance = 1e-6)
  expect_lt(detour_index(holed)$index, detour_index(solid)$index)
})

test_that("detour_index degrades smoothly (no discontinuity) as multi-part pieces separate", {
  make_dumbbell_gap <- function(gap) {
    sq1 <- make_square(1, center = c(0, 0))
    sq2 <- make_square(1, center = c(2 + gap, 0))
    wc(sf::st_union(sf::st_sfc(sq1, sq2)))
  }
  gaps <- c(0.5, 2, 5, 10, 20)
  vals <- vapply(gaps, function(g) detour_index(make_dumbbell_gap(g))$index, numeric(1))
  expect_true(all(is.finite(vals)) && all(vals > 0))
  expect_true(all(diff(vals) < 0))  # strictly decreasing as the gap widens
})

test_that("detour_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- detour_index_sf(x))
  expect_equal(res$detour_index[1], detour_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$detour_index[2], detour_index(x[2, ])$index, tolerance = 1e-8)
})

## -- exchange_index --------------------------------------------------------------

test_that("exchange_index is ~1 for a circle and less than 1 for an elongated shape", {
  expect_equal(exchange_index(wc(make_disk(1, n_seg = 90)))$index, 1, tolerance = 1e-2)
  expect_lt(exchange_index(wc(make_rectangle(20, 1)))$index, exchange_index(wc(make_square()))$index)
})

test_that("exchange_index is unaffected by a hole only insofar as the reference circle shrinks", {
  # the hole is centred and small relative to the square, so it doesn't
  # change WHICH points are near the centroid - just the circle's radius
  solid <- wc(make_square(5))
  holed <- wc(make_square_with_hole(5, 0.1))
  expect_lt(exchange_index(holed)$index, exchange_index(solid)$index)
})

test_that("exchange_index collapses toward, and then exactly to, 0 as multi-part pieces separate - the documented pathology", {
  make_dumbbell_gap <- function(gap) {
    sq1 <- make_square(1, center = c(0, 0))
    sq2 <- make_square(1, center = c(2 + gap, 0))
    wc(sf::st_union(sf::st_sfc(sq1, sq2)))
  }
  vals <- vapply(c(0.5, 2, 5, 20), function(g) exchange_index(make_dumbbell_gap(g))$index, numeric(1))
  expect_true(vals[1] > vals[2])  # still falling while the circle partially reaches both squares
  expect_equal(vals[3], 0, tolerance = 1e-8)  # gap = 5: reference circle no longer reaches either square
  expect_equal(vals[4], 0, tolerance = 1e-8)  # gap = 20: same, further still
})

test_that("exchange_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- exchange_index_sf(x))
  expect_equal(res$exchange_index[1], exchange_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$exchange_index[2], exchange_index(x[2, ])$index, tolerance = 1e-8)
})
