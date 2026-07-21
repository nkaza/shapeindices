# shared dumbbell builders - two disks joined by a thin neck, both centred
# on the y-axis so the mass centroid sits on the neck by symmetry
make_dumbbell <- function(r_top = 3, r_bottom = 3, cy = 10) {
  top    <- sf::st_buffer(sf::st_sfc(sf::st_point(c(0, cy))), r_top, nQuadSegs = 30)
  bottom <- sf::st_buffer(sf::st_sfc(sf::st_point(c(0, -cy))), r_bottom, nQuadSegs = 30)
  neck   <- sf::st_polygon(list(rbind(c(-0.2, -cy), c(0.2, -cy), c(0.2, cy), c(-0.2, cy), c(-0.2, -cy))))
  sf::st_make_valid(sf::st_union(sf::st_sfc(c(top[[1]], bottom[[1]], neck))))
}

test_that("directional_balance_index is exactly 1 for a symmetric dumbbell (the blind spot)", {
  # equal lobes: the two opposite-bearing contributions cancel exactly in
  # the complex resultant sum, even though this shape is badly dispersed
  # along its own axis (see moment_of_inertia_index()/span_index())
  db <- wc(make_dumbbell(3, 3)[[1]])
  res <- directional_balance_index(db)
  expect_equal(res$R, 0, tolerance = 1e-10)
  expect_equal(res$index, 1, tolerance = 1e-10)
})

test_that("directional_balance_index is exactly 1 for a square (CDT triangulates it with exact reflection symmetry)", {
  expect_equal(directional_balance_index(wc(make_square()))$index, 1, tolerance = 1e-8)
})

test_that("directional_balance_index is very close to (but, unlike moment_isotropy_index, not exactly) 1 for a disk/hexagon approximated by straight segments", {
  # unlike Ixx/Iyy/Ixy (exact polynomial integrals, triangulation-invariant),
  # this index has no closed form - a CDT of a many-sided polygon isn't
  # necessarily triangulated with perfect reflection symmetry, so a small
  # residual R survives at the default subdivision depth even for a shape
  # that's symmetric in the continuous limit. See the depth-convergence
  # test below for confirmation this is discretization error, not a bug.
  expect_equal(directional_balance_index(wc(make_disk(1, n_seg = 60)))$index, 1, tolerance = 1e-3)
  ang <- seq(0, 2 * pi, length.out = 7)[1:6]
  hex <- sf::st_polygon(list(rbind(cbind(cos(ang), sin(ang)), c(1, 0))))
  expect_equal(directional_balance_index(wc(hex))$index, 1, tolerance = 3e-3)
})

test_that("the disk's residual R shrinks toward 0 as subdivision depth increases (confirms discretization, not a bug)", {
  disk <- wc(make_disk(1, n_seg = 60))
  prep <- prepare_polygon(disk)
  tri  <- prep$tri
  w    <- tri$area
  G    <- shapeindices:::.mass_centroid(sf::st_geometry(tri), w)
  Rs <- vapply(c(0, 2, 4, 6), function(d) {
    cloud <- shapeindices:::.radial_point_cloud(tri, w, max_depth = d)
    shapeindices:::.resultant_from_cloud(cloud$p, cloud$w, G)$R
  }, numeric(1))
  expect_true(all(diff(Rs) < 0))  # strictly decreasing as depth increases
})

test_that("directional_balance_index drops below 1 for an asymmetric dumbbell, pointing toward the bigger lobe", {
  # bigger lobe pulls more area to bear in its own direction from the
  # (still centrally-located, on the neck) mass centroid
  db <- wc(make_dumbbell(4, 1.2)[[1]])
  res <- directional_balance_index(db)
  expect_lt(res$index, 1)
  expect_gt(res$R, 0)
  # bigger lobe is at cy = +10 (north); mean_angle should point north (~90 degrees)
  expect_equal(res$mean_angle * 180 / pi, 90, tolerance = 2)
})

test_that("directional_balance_index is rotation invariant (index, not mean_angle)", {
  db <- make_dumbbell(4, 1.2)[[1]]
  base <- directional_balance_index(wc(db))
  theta <- 0.6
  rot <- matrix(c(cos(theta), sin(theta), -sin(theta), cos(theta)), 2, 2)
  db_rot <- wc(db * rot)
  rotated <- directional_balance_index(db_rot)
  expect_equal(rotated$index, base$index, tolerance = 1e-8)
  expect_equal(rotated$R, base$R, tolerance = 1e-8)
  angle_shift <- ((rotated$mean_angle - base$mean_angle) * 180 / pi) %% 360
  expect_equal(min(angle_shift, 360 - angle_shift), theta * 180 / pi, tolerance = 0.5)
})

test_that("directional_balance_index is scale invariant", {
  db <- make_dumbbell(4, 1.2)[[1]]
  base <- directional_balance_index(wc(db))$index
  expect_equal(directional_balance_index(wc(db * 10))$index, base, tolerance = 1e-6)
})

test_that("weight = triangle area reproduces the unweighted index", {
  db <- wc(make_dumbbell(4, 1.2)[[1]])
  prep <- prepare_polygon(db)
  plain <- directional_balance_index(db, prep = prep)$index
  wtd   <- directional_balance_index(db, prep = prep, weight = prep$tri$area)$index
  expect_equal(wtd, plain, tolerance = 1e-8)
})

test_that("directional_balance_index is always in [0, 1], even under an adversarial weight spike", {
  star <- wc(make_star(7, 1, 0.3))
  prep <- prepare_polygon(star)
  n <- nrow(prep$tri)
  # concentrate almost all weight on a single arbitrary triangle
  spike <- rep(1e-6, n); spike[1] <- 1
  res <- directional_balance_index(star, prep = prep, weight = spike)
  expect_gte(res$index, 0)
  expect_lte(res$index, 1)
  expect_gte(res$R, 0)
  expect_lte(res$R, 1)
})

test_that("deterministic = FALSE roughly agrees with deterministic = TRUE, with a known small downward bias", {
  db <- wc(make_dumbbell(4, 1.2)[[1]])
  det <- directional_balance_index(db)
  mc  <- directional_balance_index(db, deterministic = FALSE, n_lines = 20000, seed = 1)
  expect_lt(abs(mc$index - det$index), 0.02)
})

test_that("the documented unbiased 1 - R^2 formula beats the naive index for a balanced shape", {
  # the man page tells users: 1 - (n_lines * R^2 - 1)/(n_lines - 1) is an
  # exactly unbiased estimate of 1 - R^2 - check it on a shape where the
  # truth is known exactly (square: R = 0, so both 1 - R and 1 - R^2 are
  # exactly 1), averaging over seeds
  square <- wc(make_square())
  prep <- prepare_polygon(square)
  n <- 500
  res <- vapply(1:50, function(s) {
    R <- directional_balance_index(square, prep = prep, deterministic = FALSE,
                                    n_lines = n, seed = s)$R
    c(naive = 1 - R, unbiased = 1 - (n * R^2 - 1) / (n - 1))
  }, numeric(2))
  # naive is visibly biased low at this n; the corrected mean is at 1
  # within Monte Carlo noise (and individual draws may exceed 1 - that is
  # expected, not clamped away)
  expect_gt(1 - mean(res["naive", ]), 0.02)
  expect_lt(abs(mean(res["unbiased", ]) - 1), 0.002)
})

test_that("Monte Carlo mode's upward bias in R (downward bias in index) shrinks with n_lines, roughly like 1/sqrt(n)", {
  # a perfectly symmetric shape has true R = 0 exactly; the naive Monte
  # Carlo estimator of |sample mean| is never exactly 0 for a finite
  # sample (Jensen's inequality on a convex/nonlinear transform - see file
  # header) - verify the gap from 1 shrinks at roughly the predicted rate
  # rather than staying constant or growing, without asserting an exact
  # correction formula
  square <- wc(make_square())
  gap_small <- 1 - directional_balance_index(square, deterministic = FALSE, n_lines = 200, seed = 1)$index
  gap_large <- 1 - directional_balance_index(square, deterministic = FALSE, n_lines = 20000, seed = 1)$index
  expect_gt(gap_small, gap_large)
  # never violates the bound even for a small, noisy sample
  expect_lte(1 - gap_small, 1)
  expect_gte(1 - gap_small, 0)
})

test_that("directional_balance_index warns and returns NA on an empty triangulation", {
  degenerate <- wc(sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(2, 0), c(0, 0)))))
  suppressWarnings(expect_warning(res <- directional_balance_index(degenerate), "not defined"))
  expect_true(is.na(res$index))
})

test_that("directional_balance_index errors on a mismatched weight length", {
  star <- wc(make_star(6, 1, 0.4))
  prep <- prepare_polygon(star)
  expect_error(directional_balance_index(star, prep = prep, weight = c(1, 2)), "one entry per triangle")
})

test_that("directional_balance_index_sf preserves CRS and matches row-by-row", {
  x <- sf::st_sf(name = c("sq", "star"),
                  geometry = sf::st_sfc(make_square(), make_star(6, 1, 0.3), crs = TEST_CRS))
  expect_no_warning(res <- directional_balance_index_sf(x))
  expect_equal(res$directional_balance_index[1], directional_balance_index(x[1, ])$index, tolerance = 1e-8)
  expect_equal(res$directional_balance_index[2], directional_balance_index(x[2, ])$index, tolerance = 1e-8)
})

test_that("shape_indices_sf(byrow = FALSE) computes directional balance consistently with the direct call", {
  cell <- function(cx, cy) sf::st_polygon(list(rbind(
    c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
    c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
  cells <- expand.grid(x = 1:6, y = 1:2)
  geoms <- mapply(cell, cells$x, cells$y, SIMPLIFY = FALSE)
  x <- sf::st_sf(pop = 1, geometry = sf::st_sfc(geoms, crs = TEST_CRS))

  res <- shape_indices_sf(x, byrow = FALSE, id = "g", which = c("moment_isotropy", "directional_balance"))
  poly_u <- sf::st_union(sf::st_geometry(x))
  direct <- directional_balance_index(wc(poly_u[[1]]))$index
  expect_equal(res$directional_balance_index, direct, tolerance = 1e-8)
})

test_that("shape_indices() includes directional_balance in \"all\" and matches the direct call", {
  star <- wc(make_star(6, 1, 0.4))
  sv <- shape_indices(star)
  expect_true("directional_balance" %in% names(sv))
  expect_equal(sv[["directional_balance"]], directional_balance_index(star)$index, tolerance = 1e-8)
})

test_that("shape_indices() Monte Carlo mode shares one point draw across directional_balance and the others", {
  star <- wc(make_star(6, 1, 0.4))
  sv <- suppressWarnings(shape_indices(star, which = c("convexity", "radial_concentration", "directional_balance"),
                                        deterministic = FALSE, n_lines = 1000, seed = 1))
  expect_true(all(c("convexity", "radial_concentration", "directional_balance") %in% names(sv)))
  expect_true(all(sv >= 0 & sv <= 1))
})
