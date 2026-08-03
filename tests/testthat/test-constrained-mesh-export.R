test_that("constrained_mesh() is a thin pass-through of .constrained_weighted_mesh()", {
    g <- make_weighted_grid(3, 3)
    direct <- shapeindices:::.constrained_weighted_mesh(g$x, weights = "pop")
    expect_true(direct$ok)

    via_export <- constrained_mesh(g$x, weights = "pop")

    expect_equal(via_export$P, direct$P)
    expect_equal(via_export$T, direct$T)
    expect_equal(via_export$tri_area, direct$tri_area)
    expect_equal(via_export$tri_weight, direct$tri_weight)
    expect_equal(via_export$poly_id, direct$poly_id)
    expect_equal(via_export$raw_total, direct$raw_total)
})

test_that("constrained_mesh()'s poly_id reproduces each input polygon's own area exactly", {
    g <- make_weighted_grid(4, 3)
    mesh <- constrained_mesh(g$x, weights = "pop")

    area_by_poly <- tapply(mesh$tri_area, mesh$poly_id, sum)
    expect_equal(length(area_by_poly), nrow(g$x))
    expect_equal(as.numeric(area_by_poly), as.numeric(sf::st_area(g$x)), tolerance = 1e-8)
})

test_that("constrained_mesh()'s poly_id reproduces each input polygon's own weight exactly", {
    g <- make_weighted_grid(4, 3)
    mesh <- constrained_mesh(g$x, weights = "pop")

    weight_by_poly <- tapply(mesh$tri_weight, mesh$poly_id, sum)
    expect_equal(as.numeric(weight_by_poly), g$w, tolerance = 1e-8)
})

test_that("constrained_mesh() errors (not a partial result) when rows genuinely overlap", {
    sq1 <- make_square(2, center = c(0, 0))
    sq2 <- make_square(2, center = c(1, 0)) # overlaps sq1
    x <- sf::st_sf(pop = c(10, 20), geometry = wc(sf::st_sfc(sq1, sq2)))

    expect_error(constrained_mesh(x, weights = "pop"), "overlap")
})

test_that("constrained_mesh() auto-projects geographic input like every other entry point", {
    nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
    mesh <- constrained_mesh(nc[1:5, ], weights = "BIR74")
    expect_false(sf::st_is_longlat(mesh$crs))
    expect_equal(length(unique(mesh$poly_id)), 5)
})
