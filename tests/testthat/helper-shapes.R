# Shared shape constructors + fixtures for tests. Files named "helper-*"
# are sourced by testthat before any test file runs, without being treated
# as tests themselves.

make_square <- function(half = 1, center = c(0, 0)) {
  h <- half
  sf::st_polygon(list(rbind(
    c(center[1] - h, center[2] - h), c(center[1] + h, center[2] - h),
    c(center[1] + h, center[2] + h), c(center[1] - h, center[2] + h),
    c(center[1] - h, center[2] - h)
  )))
}

make_rectangle <- function(width, height, center = c(0, 0)) {
  hw <- width / 2; hh <- height / 2
  sf::st_polygon(list(rbind(
    c(center[1] - hw, center[2] - hh), c(center[1] + hw, center[2] - hh),
    c(center[1] + hw, center[2] + hh), c(center[1] - hw, center[2] + hh),
    c(center[1] - hw, center[2] - hh)
  )))
}

make_star <- function(n_points, r_outer = 1, r_inner = 0.5, center = c(0, 0)) {
  n <- n_points * 2
  angles <- seq(pi / 2, pi / 2 + 2 * pi, length.out = n + 1)[1:n]
  radii  <- rep(c(r_outer, r_inner), n_points)
  x <- center[1] + radii * cos(angles)
  y <- center[2] + radii * sin(angles)
  coords <- rbind(cbind(x, y), c(x[1], y[1]))
  sf::st_polygon(list(coords))
}

make_disk <- function(radius = 1, center = c(0, 0), n_seg = 60) {
  sf::st_buffer(sf::st_sfc(sf::st_point(center)), dist = radius, nQuadSegs = n_seg)[[1]]
}

make_square_with_hole <- function(outer_half = 5, hole_frac = 0.3) {
  outer <- rbind(c(-outer_half, -outer_half), c(outer_half, -outer_half),
                 c(outer_half,  outer_half),  c(-outer_half,  outer_half),
                 c(-outer_half, -outer_half))
  hh <- outer_half * sqrt(hole_frac)
  hole <- rbind(c(-hh, -hh), c(-hh, hh), c(hh, hh), c(hh, -hh), c(-hh, -hh))
  sf::st_polygon(list(outer, hole))
}

# self-intersecting ring - not a valid simple polygon
make_bowtie <- function() {
  sf::st_polygon(list(rbind(c(0, 0), c(10, 10), c(10, 0), c(0, 10), c(0, 0))))
}

# a big square with one tiny spike - CDT gives a handful of large triangles
# for the square body and several tiny ones near the spike, a huge
# (self-contained, no external file dependency) area ratio - used by
# radial_concentration_index() and span_index() adaptive-depth tests alike
make_spiky_square <- function() {
  big   <- rbind(c(0, 0), c(100, 0), c(100, 100), c(0, 100), c(0, 0))
  spike <- rbind(c(50, 100), c(50.3, 100), c(50.15, 100.5), c(50, 100))
  poly  <- sf::st_union(sf::st_sfc(sf::st_polygon(list(big)), sf::st_polygon(list(spike))))
  wc(poly)
}

# a fixed projected (metres) CRS for tests that need one - anything planar
# works here, it just has to be non-geographic so .ensure_projected()
# treats it as already-planar rather than triggering auto-projection
TEST_CRS <- 3857

# wrap a bare sfg in an sfc carrying TEST_CRS - for tests that aren't
# themselves about CRS handling, this keeps output free of the (correct,
# expected, but irrelevant-to-what's-being-tested) "no CRS set" warning
# .ensure_projected() raises on bare geometry. Tests that DO care about
# missing/geographic CRS behaviour construct their own sfc directly
# instead of using this.
wc <- function(g, crs = TEST_CRS) sf::st_sfc(g, crs = crs)
