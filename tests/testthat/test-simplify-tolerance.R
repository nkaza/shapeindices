# a jagged (non-rectangular) arrangement of unit grid cells, one sf row per
# cell - the union's OUTER boundary is a staircase with plenty of detail
# for simplification to smooth away, while the cells themselves tessellate
# with zero gaps (like the real Census-block scenario this exists for)
.jagged_grid_sf <- function(n_side = 12, seed = 1) {
  cell <- function(cx, cy) sf::st_polygon(list(rbind(
    c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
    c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
  set.seed(seed)
  cells <- expand.grid(x = seq_len(n_side), y = seq_len(n_side))
  # trim a random jagged edge instead of a clean rectangle
  keep <- with(cells, y <= n_side - floor(n_side / 3 * (x / n_side)) + sample(-1:1, nrow(cells), replace = TRUE))
  cells <- cells[keep, ]
  geoms <- mapply(cell, cells$x, cells$y, SIMPLIFY = FALSE)
  sf::st_sf(pop = seq_len(nrow(cells)), geometry = sf::st_sfc(geoms, crs = TEST_CRS))
}

test_that("simplify_tolerance (post-union) does not fragment a clean tessellation", {
  x <- .jagged_grid_sf()
  ru_plain <- shapeindices:::.resolve_union(x, weights = NULL)
  ru_simp  <- shapeindices:::.resolve_union(x, weights = NULL, simplify_tolerance = 0.3)

  parts_plain <- shapeindices:::.split_polygon_parts(ru_plain$poly_u)
  parts_simp  <- shapeindices:::.split_polygon_parts(ru_simp$poly_u)
  expect_equal(length(parts_simp), length(parts_plain))
})

test_that("simplify_tolerance (post-union) reduces total boundary edge count", {
  x <- .jagged_grid_sf()
  ru_plain <- shapeindices:::.resolve_union(x, weights = NULL)
  ru_simp  <- shapeindices:::.resolve_union(x, weights = NULL, simplify_tolerance = 0.3)

  E_plain <- nrow(shapeindices:::.boundary_edges(ru_plain$poly_u)$p1)
  E_simp  <- nrow(shapeindices:::.boundary_edges(ru_simp$poly_u)$p1)
  expect_lt(E_simp, E_plain)
})

test_that("simplify_tolerance preserves at least 95% of true weight for a modest tolerance", {
  x <- .jagged_grid_sf()
  ru <- shapeindices:::.resolve_union(x, weights = "pop", simplify_tolerance = 0.3)
  row_area <- as.numeric(sf::st_area(ru$geoms))
  row_density <- ru$w_row / row_area
  tri <- cdt_triangles(ru$poly_u)
  rows_kept <- sf::st_sf(density = row_density, geometry = ru$geoms)
  tri_sfc <- sf::st_sf(tri_id = seq_len(nrow(tri)), geometry = sf::st_geometry(tri))
  ov <- suppressWarnings(sf::st_intersection(tri_sfc, rows_kept))
  captured <- sum(as.numeric(sf::st_area(ov)) * ov$density)
  expect_gt(captured / ru$raw_total, 0.95)
})

test_that(".weighted_mesh with simplify_tolerance gives a smaller mesh with weight still summing to 1", {
  x <- .jagged_grid_sf()
  wm_plain <- shapeindices:::.weighted_mesh(x, weights = "pop")
  wm_simp  <- shapeindices:::.weighted_mesh(x, weights = "pop", simplify_tolerance = 0.3)
  expect_lt(nrow(wm_simp$pieces), nrow(wm_plain$pieces))
  expect_equal(sum(wm_simp$pieces$weight), 1, tolerance = 1e-8)
})

test_that("shape_indices_sf(byrow = FALSE, simplify_tolerance = ...) runs end-to-end and stays close to the unsimplified result", {
  x <- .jagged_grid_sf()
  res_plain <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, id = "plain", weights = "pop",
                                                  which = c("moment_of_inertia", "hull_ratio")))
  res_simp  <- suppressWarnings(shape_indices_sf(x, byrow = FALSE, id = "simp", weights = "pop", simplify_tolerance = 0.3,
                                                  which = c("moment_of_inertia", "hull_ratio")))
  expect_equal(res_simp$moment_of_inertia_index, res_plain$moment_of_inertia_index, tolerance = 0.1)
  expect_equal(res_simp$hull_ratio_index, res_plain$hull_ratio_index, tolerance = 0.1)
})

test_that("prepare_polygon(simplify_tolerance = ...) reduces triangle count on a detailed single polygon", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))  # a finely-segmented circle - lots of near-collinear detail
  prep_plain <- prepare_polygon(disk)
  prep_simp  <- prepare_polygon(disk, simplify_tolerance = 0.2)
  expect_lt(nrow(prep_simp$tri), nrow(prep_plain$tri))
})

test_that("shape_indices(simplify_tolerance = ...) is passed through to prepare_polygon() when a mesh is needed", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  r <- shape_indices(disk, which = c("moment_of_inertia", "hull_ratio"), simplify_tolerance = 0.2)
  expect_true(is.finite(unname(r["moment_of_inertia"])))
  expect_equal(unname(r["hull_ratio"]), 1, tolerance = 0.02)
})

test_that("simplify_tolerance has no effect when only classic metrics are requested (no mesh built)", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  r_plain <- shape_indices(disk, which = "hull_ratio")
  r_simp  <- shape_indices(disk, which = "hull_ratio", simplify_tolerance = 0.2)
  expect_equal(r_simp, r_plain)
})

test_that("convexity_index(simplify_tolerance = ...) matches building prep separately, in both deterministic modes", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  prep_simp <- prepare_polygon(disk, simplify_tolerance = 0.2)

  direct_det <- convexity_index(disk, simplify_tolerance = 0.2)
  via_prep_det <- convexity_index(disk, prep = prep_simp)
  expect_equal(direct_det$index, via_prep_det$index)
  expect_equal(nrow(direct_det$triangles), nrow(prep_simp$tri))

  direct_mc <- suppressWarnings(convexity_index(disk, deterministic = FALSE, n_lines = 50, seed = 1, simplify_tolerance = 0.2))
  via_prep_mc <- suppressWarnings(convexity_index(disk, deterministic = FALSE, n_lines = 50, seed = 1, prep = prep_simp))
  expect_equal(direct_mc$index, via_prep_mc$index)
})

test_that("span_index(simplify_tolerance = ...) matches building prep separately, in both deterministic modes", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  prep_simp <- prepare_polygon(disk, simplify_tolerance = 0.2)

  direct_det <- span_index(disk, simplify_tolerance = 0.2)
  via_prep_det <- span_index(disk, prep = prep_simp)
  expect_equal(direct_det$index, via_prep_det$index)

  direct_mc <- suppressWarnings(span_index(disk, deterministic = FALSE, n_lines = 50, seed = 1, simplify_tolerance = 0.2))
  via_prep_mc <- suppressWarnings(span_index(disk, deterministic = FALSE, n_lines = 50, seed = 1, prep = prep_simp))
  expect_equal(direct_mc$index, via_prep_mc$index)
})

test_that("radial_concentration_index(simplify_tolerance = ...) matches building prep separately, in both deterministic modes", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  prep_simp <- prepare_polygon(disk, simplify_tolerance = 0.2)

  direct_det <- radial_concentration_index(disk, simplify_tolerance = 0.2)
  via_prep_det <- radial_concentration_index(disk, prep = prep_simp)
  expect_equal(direct_det$index, via_prep_det$index)

  direct_mc <- suppressWarnings(radial_concentration_index(disk, deterministic = FALSE, n_lines = 50, seed = 1, simplify_tolerance = 0.2))
  via_prep_mc <- suppressWarnings(radial_concentration_index(disk, deterministic = FALSE, n_lines = 50, seed = 1, prep = prep_simp))
  expect_equal(direct_mc$index, via_prep_mc$index)
})

test_that("moment_of_inertia_index(simplify_tolerance = ...) matches building prep separately", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  prep_simp <- prepare_polygon(disk, simplify_tolerance = 0.2)
  expect_equal(moment_of_inertia_index(disk, simplify_tolerance = 0.2)$index,
               moment_of_inertia_index(disk, prep = prep_simp)$index)
})

test_that("convex_decompose(simplify_tolerance = ...) matches building prep separately", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  prep_plain <- prepare_polygon(disk)
  prep_simp  <- prepare_polygon(disk, simplify_tolerance = 0.2)
  direct <- convex_decompose(disk, simplify_tolerance = 0.2)
  via_prep <- convex_decompose(disk, prep = prep_simp)
  expect_equal(nrow(direct), nrow(via_prep))
  expect_lt(nrow(direct), nrow(convex_decompose(disk, prep = prep_plain)))
})

test_that("subdivide_mesh(simplify_tolerance = ...) matches building prep separately", {
  disk <- wc(make_disk(radius = 5, n_seg = 200))
  prep_plain <- prepare_polygon(disk)
  prep_simp  <- prepare_polygon(disk, simplify_tolerance = 0.2)
  direct <- subdivide_mesh(disk, simplify_tolerance = 0.2)
  via_prep <- subdivide_mesh(disk, prep = prep_simp)
  expect_equal(nrow(direct), nrow(via_prep))
  expect_lt(nrow(prep_simp$tri), nrow(prep_plain$tri))
})

test_that("prepare_polygon warns when simplify_tolerance is used on a non-metric projected CRS", {
  # EPSG:2264 - NC State Plane, US survey feet
  sq_feet <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 2264)
  expect_warning(prepare_polygon(sq_feet, simplify_tolerance = 50), "US survey foot")
})

test_that("prepare_polygon does not warn when simplify_tolerance is used on a metric CRS", {
  sq_m <- wc(make_square(500))
  expect_no_warning(prepare_polygon(sq_m, simplify_tolerance = 50))
})

test_that("prepare_polygon does not warn when simplify_tolerance is used on geographic input (auto-projects to metres)", {
  disk <- sf::st_sfc(make_disk(radius = 0.05, n_seg = 50), crs = 4326)
  expect_no_warning(prepare_polygon(disk, simplify_tolerance = 50))
})

test_that("prepare_polygon does not warn about units when simplify_tolerance is NULL, even on a non-metric CRS", {
  sq_feet <- sf::st_sfc(sf::st_polygon(list(rbind(
    c(0, 0), c(1000, 0), c(1000, 1000), c(0, 1000), c(0, 0)))), crs = 2264)
  expect_no_warning(prepare_polygon(sq_feet))
})

test_that("shape_indices_sf(byrow = FALSE) warns when simplify_tolerance is used on a non-metric CRS", {
  cell <- function(cx, cy) sf::st_polygon(list(rbind(
    c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
    c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
  cells <- expand.grid(x = 1:3, y = 1:3)
  geoms <- mapply(cell, cells$x, cells$y, SIMPLIFY = FALSE)
  x <- sf::st_sf(pop = 1, geometry = sf::st_sfc(geoms, crs = 2264))
  expect_warning(
    shape_indices_sf(x, byrow = FALSE, id = "g", simplify_tolerance = 0.5,
                      which = "moment_of_inertia"),
    "US survey foot"
  )
})
