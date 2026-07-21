test_that("moment_of_inertia_index is ~1 for a disk", {
  disk <- make_disk(1, n_seg = 60)
  expect_equal(moment_of_inertia_index(wc(disk))$index, 1, tolerance = 0.01)
})

test_that("moment_of_inertia_index falls as a shape elongates", {
  square_moi <- moment_of_inertia_index(wc(make_square()))$index
  rect_moi   <- moment_of_inertia_index(wc(make_rectangle(20, 1)))$index
  expect_lt(square_moi, 1)
  expect_lt(rect_moi, square_moi)
})

test_that("moment_of_inertia_index is scale invariant", {
  star <- make_star(7, 1, 0.4)
  moi_base <- moment_of_inertia_index(wc(star))$index
  expect_equal(moment_of_inertia_index(wc(star * 10))$index, moi_base, tolerance = 1e-6)
})

test_that("weight = triangle area reproduces the unweighted index", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  moi_plain <- moment_of_inertia_index(star, prep = prep)$index
  moi_wt    <- moment_of_inertia_index(star, prep = prep, weight = prep$tri$area)$index
  expect_equal(moi_wt, moi_plain, tolerance = 1e-8)
})

test_that("moment_of_inertia_index total_weight reports raw (pre-normalisation) sum", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  w <- prep$tri$area * 3 # arbitrary non-unit scale
  res <- moment_of_inertia_index(star, prep = prep, weight = w)
  expect_equal(res$total_weight, sum(w), tolerance = 1e-8)
})

test_that("weight = uniform density reproduces the plain area centroid", {
  # rho uniform (weight proportional to area) is the case where the mass
  # centroid and the plain geometric centroid are mathematically identical
  # (verified separately to floating-point precision) - regression guard
  # that the internally-computed mass centroid still agrees with GEOS's
  # own st_centroid() in this case
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  res <- moment_of_inertia_index(star, prep = prep, weight = prep$tri$area)
  geos_centroid <- sf::st_coordinates(sf::st_centroid(sf::st_union(prep$poly)))[1, 1:2]
  mass_centroid <- sf::st_coordinates(res$centroid)[1, ]
  expect_equal(unname(mass_centroid), unname(geos_centroid), tolerance = 1e-8)
})

test_that("moment_of_inertia_index uses the true mass centroid, not the plain area centroid, when weighted", {
  # regression: an earlier version always used st_centroid(st_union(poly))
  # (the plain, unweighted geometric centroid) as the reference point, even
  # when `weight` made density non-uniform - by the parallel axis theorem,
  # computing J about any point other than the true mass centroid strictly
  # inflates J (by mass * distance^2), silently understating the index.
  # Two adjacent unit squares, heavily weighted toward the right one: the
  # true mass centroid must sit well right of the rectangle's plain area
  # centroid (which stays fixed at the midpoint, x = 1).
  left  <- sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
  right <- sf::st_polygon(list(rbind(c(1, 0), c(2, 0), c(2, 1), c(1, 1), c(1, 0))))
  rect2 <- sf::st_union(sf::st_sfc(left, right))
  prep <- prepare_polygon(wc(rect2))
  tri <- prep$tri
  cx <- vapply(sf::st_geometry(tri), function(t) mean(sf::st_coordinates(t)[1:3, 1]), numeric(1))
  w <- ifelse(cx < 1, 1, 100) # heavy weight on the right square

  res <- moment_of_inertia_index(wc(rect2), prep = prep, weight = w)
  mass_centroid <- sf::st_coordinates(res$centroid)[1, ]

  # the mass centroid must have moved meaningfully right of the plain area
  # centroid (x = 1, the rectangle's midpoint) - not just be near it
  expect_gt(mass_centroid[1], 1.3)

  # independent cross-check via the parallel axis theorem: J computed about
  # the WRONG (plain area) centroid must equal J (about the correct mass
  # centroid, what the function now returns) plus total_mass * distance^2 -
  # verifies the returned centroid is not just "moved", but the exact
  # right point, not merely a step in the right direction
  area_centroid <- sf::st_coordinates(sf::st_centroid(rect2))[1, 1:2]
  d2 <- sum((area_centroid - mass_centroid)^2)
  rho_norm <- (w / sum(w)) / tri$area
  moms <- vapply(sf::st_geometry(tri), function(t) {
    v <- sf::st_coordinates(t)[1:3, 1:2, drop = FALSE]
    v[, 1] <- v[, 1] - area_centroid[1]; v[, 2] <- v[, 2] - area_centroid[2]
    signed2A <- sum(v[, 1] * v[c(2, 3, 1), 2] - v[c(2, 3, 1), 1] * v[, 2])
    if (signed2A < 0) v <- v[c(1, 3, 2), ]
    x <- v[, 1]; y <- v[, 2]; xn <- x[c(2, 3, 1)]; yn <- y[c(2, 3, 1)]; cross <- x * yn - xn * y
    c(sum((y^2 + y * yn + yn^2) * cross) / 12, sum((x^2 + x * xn + xn^2) * cross) / 12)
  }, numeric(2))
  J_about_area_centroid <- sum(rho_norm * moms[1, ]) + sum(rho_norm * moms[2, ])
  expect_equal(res$J + d2, J_about_area_centroid, tolerance = 1e-8) # total_mass (normalised) = 1
})

test_that("moment_of_inertia_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- moment_of_inertia_index_sf(x))
  expect_equal(res$moi_index[1], moment_of_inertia_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$moi_index[2], moment_of_inertia_index(x[2, ])$index, tolerance = 1e-8)
})
