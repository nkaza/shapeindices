## =========================================================================
## Convexity / dispersal index for (multi)polygons via constrained
## Delaunay triangulation (CDT)
##
## Definition: CI(P) = 1 - E[f(X,Y)], where X, Y are independent points
## drawn uniformly at random from polygon P and f(X,Y) is the fraction of
## segment XY lying outside P. CI = 1 exactly when P is convex.
##
## Method (deterministic = TRUE, the default): CDT the polygon into triangles
## T_1..T_n; connect every pair with q^2 sub-pair lines (q = n_quad, a
## Hammer-Stroud quadrature rule per triangle instead of just the
## centroid); average each pair's fraction-outside, weighted by the
## PRODUCT of the two triangles' areas (not their sum - sum-weighting
## under-penalises pairs touching a tiny triangle), plus an exact
## same-piece diagonal term (sum(area^2)/2) added to the denominator only,
## accounting for the zero-outside probability mass of both points landing
## in the same (convex) triangle. `weight` can replace triangle area
## throughout - exact, not a heuristic, since the derivation only needs
## g_ii = 0 and g_ij = g_ji.
##
## deterministic = FALSE replaces the exhaustive sum with a Monte Carlo
## estimate over `n_lines` random interior lines - useful when the mesh is
## too large to enumerate.
##
## Cost: deterministic = TRUE is O(n_quad * n^2). Per-line cost is reduced via
## vectorised segment-vs-boundary-edge intersection instead of one GEOS
## call per line (~4x faster, falls back to the GEOS loop on unexpected
## geometry). Parallelism in this package is row-level only
## (shape_indices_sf(byrow = TRUE)'s `parallel_rows`) - pair-level dispatch
## inside convexity_index() isn't worth the added failure modes for this
## rarely-exercised fallback path.
##
## KNOWN LIMITATION: on spiky/comb-like shapes, the all-pairs average can
## be misleading. Too few triangles (e.g. a 3-5 point star) can let every
## evaluated line miss the notches entirely, giving a false 1.0. Symmetric
## even-pointed stars have a geometric floor: the diametrically-opposite
## point pair is genuinely mutually visible through the shared centre,
## capping how low the index can go regardless of weighting. Treat the
## index as a lower bound on non-convexity for such shapes, not an exact
## one - which is also why this argument is called `deterministic`, not
## `exact`: it's a FIXED quadrature rule (q^2 sample lines per triangle
## pair), not a literal closed-form answer, and a small concavity relative
## to triangle size can still slip between those fixed sample points even
## with plenty of triangles. deterministic = FALSE with a large `n_lines`
## converges to the true value; deterministic = TRUE at the default n_quad
## does not, in general - verified empirically: refining the quadrature
## (n_quad 1 -> 3 -> a much finer manual grid) moves the deterministic
## value monotonically toward the Monte Carlo one, not the other way
## round.

## -- the convexity index itself --------------------------------------
##
## SAME-PIECE DIAGONAL TERM: the all-pairs sum over DIFFERENT pieces omits
## the case where both random points land in the SAME piece
## (frac_outside = 0 there for any convex piece) - fine for the numerator,
## but omitting it from the denominator too biases the index down, since
## that's real zero-outside probability mass. The exact fix is
## sum(weight^2)/2 added to the denominator (product-of-areas weighting
## extended to i = j, halved for unordered points) - only valid when every
## piece is guaranteed convex, which every mesh this package ever builds
## (CDT triangles, convex_decompose()'s Hertel-Mehlhorn pieces) is.

#' Deterministic (deterministic = TRUE) index over a fixed mesh of convex
#' pieces - the shared engine behind convexity_index() and
#' shape_indices_sf(byrow = FALSE)'s combined mesh.
#' @param pieces sf data frame of CDT triangles, with an `area` column and
#'   POLYGON geometry
#' @param poly the (multi)polygon `pieces` triangulates, used to build the
#'   reference union tested against
#' @param n_quad quadrature points per triangle (1 or 3)
#' @param plot draw a diagnostic plot (needs a graphics device)
#' @param weight optional numeric vector, length nrow(pieces); NULL uses
#'   `pieces$area`
#' @return list(index, triangles, edges) - see convexity_index()'s own
#'   @return for the field meanings
#' @noRd
.mesh_convexity_index <- function(pieces, poly, n_quad, plot, weight = NULL) {
    n <- if (is.null(pieces)) 0 else nrow(pieces)
    # checked before the n == 0/1 short-circuits - a wrong-length weight
    # would otherwise silently recycle into a plausible-looking wrong
    # index, or NA with no warning
    if (!is.null(weight) && length(weight) != n) {
        stop("`weight` must have one entry per triangle (", n, "), got ", length(weight), ".")
    }
    if (n == 0) {
        warning("Triangulation/decomposition produced no pieces; index is not defined.")
        return(list(index = NA_real_, triangles = pieces, edges = NULL))
    }
    if (n == 1) {
        return(list(index = 1, triangles = pieces, edges = NULL))  # vacuously convex
    }
    # cost tracks n_pairs * n_quad^2, not n_pairs alone - a large mesh can
    # OOM the vectorised line-clipping method, which then falls back
    # (gracefully, sequentially) to a per-line GEOS loop
    n_pairs          <- choose(n, 2)
    n_candidate_lines <- n_pairs * n_quad^2
    if (n_candidate_lines > 20000) {
        warning(
            "Large mesh: ", n, " pieces, ", n_pairs, " pairs, ", n_candidate_lines,
            " candidate lines to evaluate at n_quad = ", n_quad, " (deterministic = TRUE ",
            "is O(n_quad^2 * n^2)) - this can be slow. If it is, try: n_quad = 1 (drops ",
            "the ", n_quad^2, "x quadrature multiplier), deterministic = FALSE (random-line ",
            "estimate instead of exhaustive), or shape_indices()/shape_indices_sf()'s ",
            "deterministic_max_tri to switch to deterministic = FALSE automatically above a ",
            "size threshold."
        )
    }

    crs    <- st_crs(poly)
    poly_u <- st_union(poly)
    # weight already normalised to sum to 1 by the caller; NULL keeps
    # physical (unnormalised) area
    piece_weight <- if (is.null(weight)) pieces$area else weight
    q            <- n_quad

    qmat <- if (q == 1) {
        st_coordinates(st_centroid(st_geometry(pieces)))[, 1:2, drop = FALSE]
    } else {
        do.call(rbind, lapply(st_geometry(pieces), function(g) {
            tri_quad_points(st_coordinates(g)[1:3, 1:2, drop = FALSE])
        }))
    }

    res   <- .eval_mesh_pairs(qmat, q, n, poly_u, crs)
    edges <- res$edges
    g_ij  <- res$g_ij

    # product weighting (not sum) - see file header for why
    w <- piece_weight[edges[, 1]] * piece_weight[edges[, 2]]
    diag_term <- sum(piece_weight^2) / 2   # same-piece term - see file header
    index <- 1 - sum(w * g_ij) / (sum(w) + diag_term)   # 1 = fully convex

    edges_df <- setNames(as.data.frame(edges), c("from", "to"))
    # diagnostic geometry always at centroid level, regardless of n_quad
    cen      <- st_centroid(st_geometry(pieces))
    plot_seg <- mapply(function(i, j) st_linestring(rbind(cen[[i]], cen[[j]])),
                        edges_df$from, edges_df$to, SIMPLIFY = FALSE)
    plot_seg <- st_sfc(plot_seg, crs = crs)
    edges_sf <- st_sf(edges_df, frac_outside = g_ij, weight = w, geometry = plot_seg)

    if (plot) {
        plot(poly, col = "grey90", border = "grey40")
        plot(st_geometry(pieces), add = TRUE, border = "grey75")
        ok <- g_ij == 0
        if (any(ok))  plot(st_geometry(edges_sf)[ok], add = TRUE, col = "steelblue")
        if (any(!ok)) plot(st_geometry(edges_sf)[!ok], add = TRUE, col = "red", lwd = 2)
        plot(cen, add = TRUE, pch = 20, cex = 0.6)
    }

    list(index = index, triangles = pieces, edges = edges_sf)
}

## Stochastic (deterministic = FALSE) index: draws `n_lines` independent
## lines, each with two fresh i.i.d. random endpoints inside the polygon,
## and returns 1 minus the mean fraction of each line's length outside -
## reuses deterministic = TRUE's own scoring rule so the two are directly
## comparable. weight = NULL samples uniformly (st_sample()); weight
## supplied samples from that piecewise-constant density instead (see
## .sample_weighted_points()).

#' Draw points from a per-triangle weighted density (pick triangle i with
#' probability W_i / sum(W), then uniform within it) - exact direct
#' sampling, cheaper than importance-weight reweighting when the density
#' varies by orders of magnitude.
#' @param tri the triangle mesh (sf with POLYGON geometry, one row per
#'   triangle)
#' @param weight numeric vector, one entry per triangle, need not sum to 1
#' @param n number of points to draw
#' @return n x 2 coordinate matrix
#' @noRd
.sample_weighted_points <- function(tri, weight, n) {
    n_tri <- nrow(tri)
    # 6 x n_tri: Ax,Ay,Bx,By,Cx,Cy per triangle, built once so sampling is
    # fully vectorised
    corner_mat <- vapply(st_geometry(tri), function(g) {
        v <- st_coordinates(g)[1:3, 1:2, drop = FALSE]
        c(v[1, 1], v[1, 2], v[2, 1], v[2, 2], v[3, 1], v[3, 2])
    }, numeric(6))

    tri_idx <- sample.int(n_tri, size = n, replace = TRUE, prob = weight)
    sel     <- corner_mat[, tri_idx, drop = FALSE]

    # uniform point within a triangle via barycentric coords, folding
    # points outside the (r1+r2 <= 1) triangle back in
    r1 <- runif(n); r2 <- runif(n)
    flip     <- (r1 + r2) > 1
    r1[flip] <- 1 - r1[flip]
    r2[flip] <- 1 - r2[flip]

    Ax <- sel[1, ]; Ay <- sel[2, ]; Bx <- sel[3, ]; By <- sel[4, ]; Cx <- sel[5, ]; Cy <- sel[6, ]
    cbind(x = Ax + r1 * (Bx - Ax) + r2 * (Cx - Ax),
          y = Ay + r1 * (By - Ay) + r2 * (Cy - Ay))
}

#' The random-line index (RLI) - internal engine for convexity_index()'s
#' deterministic = FALSE. See that function's own roxygen for the
#' user-facing parameter reference; parameters here are the same, without
#' defaults.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param n_lines number of random lines to sample, or a function(n_tri)
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon()
#' @param seed optional RNG seed
#' @param plot draw a diagnostic plot (needs a graphics device)
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`; NULL samples uniformly
#' @param points optional pre-drawn (2*n_lines)x2 coordinate matrix, from
#'   shape_indices()'s point-sharing mechanism (see its own comments) -
#'   when supplied, sampling is skipped entirely
#' @return list(index, triangles, edges) - see convexity_index()'s own
#'   @return for the field meanings
#' @noRd
.random_line_index <- function(poly, n_lines, prep, seed, plot, weight = NULL, points = NULL) {
    if (is.null(prep)) prep <- prepare_polygon(poly)
    poly_geom <- prep$poly
    tri       <- prep$tri
    n_tri     <- if (is.null(tri)) 0 else nrow(tri)

    if (!is.null(weight)) {
        if (n_tri == 0) {
            stop("`weight` needs a triangle mesh to sample from, but this polygon ",
                 "triangulated to no triangles.")
        }
        if (length(weight) != n_tri) {
            stop("`weight` must have one entry per triangle (", n_tri, "), got ", length(weight), ".")
        }
    }

    # n_lines can be a function(n_tri), evaluated once triangle count is
    # known, so complexity-varying callers aren't stuck with one fixed size
    if (is.function(n_lines)) n_lines <- max(1L, round(n_lines(n_tri)))

    if (n_tri >= 2) {
        n_pairs_deterministic <- choose(n_tri, 2)
        if (n_lines >= n_pairs_deterministic / 2) {
            warning(sprintf(
                paste("n_lines (%d) is not substantially lower than the %d triangle-pairs",
                      "that deterministic = TRUE (%d triangles) would",
                      "evaluate for this same polygon; deterministic = FALSE is meant as a",
                      "cheaper approximation for meshes too large to enumerate exhaustively -",
                      "consider deterministic = TRUE instead, or a smaller n_lines."),
                n_lines, n_pairs_deterministic, n_tri))
        }
    }

    poly_u <- st_union(poly_geom)
    crs    <- st_crs(poly_geom)

    needed <- 2 * n_lines

    if (!is.null(points)) {
        if (nrow(points) < needed) {
            stop("`points` has ", nrow(points), " rows but ", needed, " are needed for n_lines = ", n_lines, ".")
        }
        coords <- points[seq_len(needed), , drop = FALSE]
    } else if (!is.null(weight)) {
        if (!is.null(seed)) set.seed(seed)
        coords <- .sample_weighted_points(tri, weight, needed)
    } else {
        if (!is.null(seed)) set.seed(seed)
        # i.i.d. uniform interior points, topping up if st_sample() returns
        # slightly fewer than requested
        pts <- st_sample(poly_u, size = needed, type = "random", exact = TRUE)
        tries <- 0
        while (length(pts) < needed && tries < 10) {
            pts <- c(pts, st_sample(poly_u, size = needed, type = "random", exact = TRUE))
            tries <- tries + 1
        }
        if (length(pts) < needed) {
            warning("Could not sample enough interior points; index is not defined.")
            return(list(index = NA_real_, triangles = tri, edges = NULL))
        }
        pts    <- pts[seq_len(needed)]
        coords <- st_coordinates(pts)[, 1:2, drop = FALSE]
    }

    # consecutive pairs -> every line has two points appearing in no other line
    x1 <- coords[seq(1, needed, 2), , drop = FALSE]
    x2 <- coords[seq(2, needed, 2), , drop = FALSE]

    if (requireNamespace("sfheaders", quietly = TRUE)) {
        seg_df <- data.frame(
            id = rep(seq_len(n_lines), each = 2),
            x  = as.vector(rbind(x1[, 1], x2[, 1])),
            y  = as.vector(rbind(x1[, 2], x2[, 2]))
        )
        seg <- sfheaders::sfc_linestring(seg_df, x = "x", y = "y", linestring_id = "id")
        st_crs(seg) <- crs
    } else {
        seg <- mapply(function(i) st_linestring(rbind(x1[i, ], x2[i, ])), seq_len(n_lines), SIMPLIFY = FALSE)
        seg <- st_sfc(seg, crs = crs)
    }

    len_total <- as.numeric(st_length(seg))
    inside  <- lengths(suppressWarnings(st_covered_by(seg, poly_u))) > 0
    idx_out <- which(!inside & len_total > 0)

    frac_outside <- .compute_frac_outside(x1, x2, idx_out, seg, len_total, poly_u, crs)

    index    <- 1 - mean(frac_outside)
    edges_sf <- st_sf(line_id = seq_len(n_lines), frac_outside = frac_outside, geometry = seg)

    if (plot) {
        plot(poly_geom, col = "grey90", border = "grey40")
        ok <- frac_outside == 0
        if (any(ok))  plot(st_geometry(edges_sf)[ok], add = TRUE, col = "steelblue")
        if (any(!ok)) plot(st_geometry(edges_sf)[!ok], add = TRUE, col = "red", lwd = 2)
    }

    list(index = index, triangles = tri, edges = edges_sf)
}

#' Convexity/dispersal index of a (multi)polygon
#'
#' 1 minus the expected fraction of a random interior line lying outside
#' the polygon. 1 = convex; lower means more concave and/or more
#' spatially dispersed.
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param deterministic if TRUE (default), compute over a fixed quadrature
#'   grid on all pairs of the polygon's own CDT triangles - O(n^2). If
#'   FALSE, use a Monte Carlo estimate over `n_lines` random interior lines
#'   instead - O(n_lines), useful when the mesh is too large, and also the
#'   more accurate option for shapes with small concavities relative to
#'   triangle size, since a fixed quadrature grid can miss a concavity that
#'   falls between its sample points (see
#'   `vignette("c-understanding-convexity-index")`). `n_quad` is unused
#'   (and an error if passed explicitly) when deterministic = FALSE.
#' @param n_quad quadrature points per triangle when deterministic = TRUE
#'   (1 or 3). 3 (default) uses a 3-point Hammer-Stroud rule, closing most
#'   of the gap to the random-line index at ~9x the cost of n_quad = 1;
#'   drop to 1 on meshes of more than a few hundred triangles if that
#'   matters.
#' @param n_lines number of random lines to sample when deterministic =
#'   FALSE (default 3000), or a function(n_tri) for callers whose polygons
#'   vary widely in complexity. Warns if not substantially lower than what
#'   deterministic = TRUE would need for the same polygon.
#' @param seed optional RNG seed, only used when deterministic = FALSE.
#' @param plot draw a diagnostic plot (needs a graphics device)
#' @param prep optional pre-computed list(poly, tri) from prepare_polygon(),
#'   to skip re-triangulating when also calling moment_of_inertia_index()
#' @param weight optional numeric vector, one entry per triangle in
#'   `prep$tri`, substituting for that triangle's own area/mass
#'   throughout. NULL (default) reproduces the unweighted index exactly.
#' @param points optional pre-drawn (2*n_lines)x2 coordinate matrix, only
#'   meaningful when deterministic = FALSE - primarily for
#'   shape_indices()'s internal point-sharing mechanism (draw one sample
#'   and reuse it across convexity/span/radial_concentration instead of
#'   each drawing independently), not typically supplied directly.
#' @param simplify_tolerance passed to `prepare_polygon()` when `prep` is
#'   NULL - see its own doc. Ignored (with no error) if `prep` is supplied
#'   directly, since simplification only happens while building `prep` in
#'   the first place.
#' @return list(index, triangles, edges). index is in `[0, 1]`, 1 = fully
#'   convex. `triangles` holds the CDT mesh. `edges` holds one row per
#'   evaluated line.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#'
#' # Wake: a simple, nearly-convex Piedmont county
#' convexity_index(wake)$index
#' # Dare: Outer Banks mainland piece plus a separated barrier-island strip
#' convexity_index(nc[nc$NAME == "Dare", ])$index
#'
#' # deterministic = FALSE: Monte Carlo (random-line) estimate instead of
#' # the exhaustive all-pairs method - a seed makes it reproducible
#' convexity_index(wake, deterministic = FALSE, n_lines = 2000, seed = 1)$index
#'
#' # weight: substitutes for each triangle's own area - weighting by the
#' # triangle's own area exactly reproduces the unweighted index
#' prep <- prepare_polygon(wake)
#' convexity_index(wake, prep = prep, weight = prep$tri$area)$index
#' @export
convexity_index <- function(poly, deterministic = TRUE, n_quad = 3, n_lines = 3000, seed = NULL,
                             plot = FALSE, prep = NULL, weight = NULL, points = NULL,
                             simplify_tolerance = NULL) {
    n_quad_given <- !missing(n_quad)
    if (!is.null(weight)) weight <- .normalize_weight(weight)
    if (is.null(prep)) prep <- prepare_polygon(poly, simplify_tolerance = simplify_tolerance)

    if (!deterministic) {
        if (n_quad_given) {
            stop("`n_quad` selects the quadrature refinement used by deterministic = TRUE; ",
                 "it has no meaning for deterministic = FALSE (the random-line index), which ",
                 "samples fresh random points rather than quadrature points on a fixed ",
                 "mesh. Drop the `n_quad` argument, or set deterministic = TRUE to use it.")
        }
        return(.random_line_index(poly, n_lines = n_lines, prep = prep, seed = seed,
                                   plot = plot, weight = weight, points = points))
    }
    if (!is.null(points)) {
        stop("`points` (a pre-drawn sample) has no meaning for deterministic = TRUE, which ",
             "computes over the full quadrature grid, not a random sample. Drop `points`, ",
             "or set deterministic = FALSE to use it.")
    }

    if (!n_quad %in% c(1, 3)) {
        stop("n_quad must be 1 (centroid only) or 3 (Hammer-Stroud rule); other ",
             "quadrature orders aren't implemented.")
    }
    poly_geom <- prep$poly

    .mesh_convexity_index(prep$tri, poly_geom, n_quad = n_quad, plot = plot, weight = weight)
}

#' convexity_index() for every row of an sf data frame
#'
#' Each row indexed independently and unweighted (see shape_indices_sf()
#' for weighting a collection of rows as one shape).
#' @param x an sf data frame
#' @param ... passed to convexity_index() for every row (e.g. deterministic,
#'   n_quad, n_lines)
#' @return `x` with one new column, convexity_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- convexity_index_sf(nc[1:5, ])
#' res$convexity_index
#'
#' # `...` passes through to convexity_index() for every row - e.g. force
#' # the Monte Carlo estimator (deterministic = FALSE) for all five counties
#' res_rli <- convexity_index_sf(nc[1:5, ], deterministic = FALSE, n_lines = 2000, seed = 1)
#' res_rli$convexity_index
#' @export
convexity_index_sf <- function(x, ...) {
    geoms <- st_geometry(x)
    # geoms[i], not lapply(geoms, .) directly - the latter strips to a bare
    # sfg with no CRS of its own
    res <- lapply(seq_along(geoms), function(i) convexity_index(geoms[i], ...))
    x$convexity_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
