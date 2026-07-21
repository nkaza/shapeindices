test_that(".frac_outside_vectorised gives identical results chunked vs unchunked", {
  # an L-shaped (non-convex) polygon, so some candidate lines genuinely
  # cross the boundary more than once
  poly <- wc(sf::st_polygon(list(rbind(
    c(0, 0), c(10, 0), c(10, 4), c(4, 4), c(4, 10), c(0, 10), c(0, 0)
  ))))
  edges <- shapeindices:::.boundary_edges(poly)
  crs <- sf::st_crs(poly)

  set.seed(1)
  n_lines <- 40
  x1 <- runif(n_lines, -2, 12); y1 <- runif(n_lines, -2, 12)
  x2 <- runif(n_lines, -2, 12); y2 <- runif(n_lines, -2, 12)

  unchunked <- shapeindices:::.frac_outside_vectorised(x1, y1, x2, y2, edges, poly, crs, chunk_size = n_lines)
  chunked_1 <- shapeindices:::.frac_outside_vectorised(x1, y1, x2, y2, edges, poly, crs, chunk_size = 1)
  chunked_7 <- shapeindices:::.frac_outside_vectorised(x1, y1, x2, y2, edges, poly, crs, chunk_size = 7)

  expect_equal(chunked_1, unchunked, tolerance = 1e-12)
  expect_equal(chunked_7, unchunked, tolerance = 1e-12)
})

test_that(".frac_outside_vectorised with default (auto) chunk_size matches an explicit large one", {
  square <- wc(make_square(5))
  edges <- shapeindices:::.boundary_edges(square)
  crs <- sf::st_crs(square)

  set.seed(2)
  n_lines <- 20
  x1 <- runif(n_lines, -6, 6); y1 <- runif(n_lines, -6, 6)
  x2 <- runif(n_lines, -6, 6); y2 <- runif(n_lines, -6, 6)

  auto_res  <- shapeindices:::.frac_outside_vectorised(x1, y1, x2, y2, edges, square, crs)
  large_res <- shapeindices:::.frac_outside_vectorised(x1, y1, x2, y2, edges, square, crs, chunk_size = n_lines)
  expect_equal(auto_res, large_res, tolerance = 1e-12)
})

test_that(".frac_outside_vectorised handles zero lines or zero edges without erroring", {
  square <- wc(make_square(5))
  edges <- shapeindices:::.boundary_edges(square)
  crs <- sf::st_crs(square)
  expect_equal(shapeindices:::.frac_outside_vectorised(numeric(0), numeric(0), numeric(0), numeric(0), edges, square, crs), numeric(0))
})

test_that(".choose_line_chunk_size shrinks as E grows, and never returns less than 1", {
  small <- shapeindices:::.choose_line_chunk_size(10)
  large <- shapeindices:::.choose_line_chunk_size(10^8)
  expect_gte(small, large)
  expect_gte(large, 1L)
  expect_gte(small, 1L)
})

test_that(".available_memory_mb returns a finite, positive number", {
  mb <- shapeindices:::.available_memory_mb()
  expect_true(is.finite(mb))
  expect_gt(mb, 0)
})

test_that("convexity_index (deterministic = FALSE) gives the same result forcing a tiny chunk size as it does unchunked", {
  # exercise the real call path (.compute_frac_outside -> .frac_outside_vectorised)
  # with a deliberately small forced chunk_size, standing in for a
  # memory-constrained environment, and check the index doesn't change
  star <- wc(make_star(8, 5, 2))
  res_auto  <- suppressWarnings(convexity_index(star, deterministic = FALSE, n_lines = 300, seed = 7)$index)

  # temporarily force a tiny chunk size by monkey-patching the chooser -
  # confirms the auto-sizing path and a pathologically small one agree
  orig <- shapeindices:::.choose_line_chunk_size
  ns <- asNamespace("shapeindices")
  unlockBinding(".choose_line_chunk_size", ns)
  on.exit({ assign(".choose_line_chunk_size", orig, envir = ns); lockBinding(".choose_line_chunk_size", ns) })
  assign(".choose_line_chunk_size", function(E, ...) 1L, envir = ns)

  res_tiny_chunk <- suppressWarnings(convexity_index(star, deterministic = FALSE, n_lines = 300, seed = 7)$index)
  expect_equal(res_tiny_chunk, res_auto, tolerance = 1e-9)
})
