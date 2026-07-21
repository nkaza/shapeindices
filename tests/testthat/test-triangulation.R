test_that("cdt_triangles triangulates a simple convex polygon exactly (area preserved)", {
  sq <- wc(make_square())
  tri <- cdt_triangles(sq)
  expect_s3_class(tri, "sf")
  expect_true(nrow(tri) >= 1)
  expect_true(all(c("tri_id", "area", "geometry") %in% names(tri)))
  expect_equal(sum(tri$area), as.numeric(sf::st_area(sq)), tolerance = 1e-6)
})

test_that("cdt_triangles handles a polygon with a hole (area preserved)", {
  h <- wc(make_square_with_hole())
  tri <- cdt_triangles(h)
  expect_true(nrow(tri) > 0)
  expect_equal(sum(tri$area), as.numeric(sf::st_area(sf::st_make_valid(h))), tolerance = 1e-6)
})

test_that("cdt_triangles keeps disjoint multipolygon parts from bridging", {
  left  <- make_square(1, center = c(-10, 0))
  right <- make_square(1, center = c(10, 0))
  # unclass()ing an sfg POLYGON gives its "list of rings" structure, which
  # is exactly the shape st_multipolygon() expects per part
  mp <- sf::st_multipolygon(list(unclass(left), unclass(right)))
  tri <- cdt_triangles(sf::st_sfc(mp))
  # no triangle should span the 18-unit gap between the two squares
  bb <- lapply(sf::st_geometry(tri), sf::st_bbox)
  widths <- vapply(bb, function(b) b["xmax"] - b["xmin"], numeric(1))
  expect_true(all(widths < 3))
})

test_that("tri_quad_points returns 3 interior points for a triangle", {
  triangle <- rbind(c(0, 0), c(4, 0), c(0, 3))
  pts <- tri_quad_points(triangle)
  expect_equal(dim(pts), c(3, 2))
  expect_true(all(is.finite(pts)))
})

test_that("prepare_polygon returns cleaned geometry + triangle mesh", {
  sq <- wc(make_square())
  prep <- prepare_polygon(sq)
  expect_named(prep, c("poly", "tri"))
  expect_true(nrow(prep$tri) >= 1)
})

test_that("prepare_polygon errors on empty geometry", {
  empty <- sf::st_polygon()
  expect_error(prepare_polygon(empty), "empty")
})

test_that("prepare_polygon messages and auto-projects geographic (lon/lat) input", {
  ll <- sf::st_sfc(make_square(0.01, center = c(-79, 35)), crs = 4326)
  expect_message(prep <- prepare_polygon(ll), "auto-projecting")
  expect_false(sf::st_is_longlat(prep$poly))
})

test_that("prepare_polygon warns on a non-simple (self-intersecting) polygon", {
  bt <- sf::st_sfc(make_bowtie(), crs = TEST_CRS)
  expect_warning(prepare_polygon(bt), "not a valid, simple polygon")
})

test_that("prepare_polygon is silent (no warning) on already-valid, already-projected input", {
  sq <- sf::st_sfc(make_square(), crs = TEST_CRS)
  expect_no_warning(prepare_polygon(sq))
})
