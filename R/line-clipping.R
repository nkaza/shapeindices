## -- vectorised line-clipping engine -------------------------------------
##
## Shared "what fraction of this line's length lies outside the polygon"
## machinery, used by both the exhaustive mesh method (.eval_mesh_pairs())
## and the random-line method (.random_line_index()) in convexity-index.R.
##
## QUADRATURE REFINEMENT (n_quad): a single centroid-to-centroid line is a
## coarse proxy for a triangle-pair's true area-averaged visibility.
## Replacing it with a 3-point Hammer-Stroud rule per triangle and
## averaging the 9 resulting sub-pair lines closes most of the remaining
## gap to the random-line index.

#' 3-point interior (Hammer-Stroud) quadrature rule for a triangle
#'
#' @param coords3x2 a 3x2 matrix of the triangle's corner coordinates
#' @return a 3x2 matrix of quadrature point coordinates
#' @examples
#' triangle <- rbind(c(0, 0), c(4, 0), c(0, 3))
#' tri_quad_points(triangle)
#' @export
tri_quad_points <- function(coords3x2) {
    alpha <- 2/3; beta <- 1/6
    bary <- rbind(c(alpha, beta, beta), c(beta, alpha, beta), c(beta, beta, alpha))
    bary %*% coords3x2
}

## -- vectorised replacement for GEOS's per-segment st_difference() loop --
##
## GEOS has no batched many-lines-vs-one-polygon op, so calling
## st_difference() once per line is dominated by per-call interop overhead.
## Since every polygon here is a simple ring set (no self-intersections),
## "fraction outside" reduces to one vectorised segment-vs-boundary-edge
## intersection test (all lines x all edges at once) - ~4x faster,
## verified to ~5e-14 of the GEOS result. Falls back to the GEOS loop if
## the vectorised path errors on unexpected geometry.
##
## CHUNKED OVER LINES, NOT DONE IN ONE SHOT: the naive version builds
## several M-lines x E-edges intermediate vectors at once. E is the
## polygon's TOTAL boundary edge count across every ring/part of poly_u -
## for a single simply-shaped county this is small (hundreds to low
## thousands), but shape_indices_sf(byrow = FALSE) can hand this a poly_u
## that's the union of many rows, and if those rows were simplified
## independently before unioning (breaking what used to be shared edges
## between adjacent rows into thousands of tiny mismatched slivers), E can
## reach the hundreds of thousands even though the triangulated MESH stays
## small - the two aren't coupled, since triangulation only needs the
## outer shape. At real-world scale (verified: E ~ 215,000, M ~ 5,000 from
## one genuine failing case) M*E reaches past a billion, and the ~15
## per-pair intermediate vectors needed past 100GB - a real crash, not a
## theoretical one. Processing lines in chunks against all E edges at once
## keeps the same per-chunk vectorised speed while bounding peak memory to
## O(chunk_size * E) regardless of M or E individually - chunk_size is
## picked automatically from currently-available memory
## (.choose_line_chunk_size()), not fixed, so a small poly_u still gets one
## unchunked pass (no speed cost) while a pathological one degrades
## gracefully into more, smaller passes instead of allocating past
## whatever RAM exists.

#' Every boundary edge of a (multi)polygon (outer + holes) as coordinate
#' pairs; winding direction doesn't matter for the crossing-count logic.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(p1, p2) - Ex2 matrices of edge start/end coordinates
#' @noRd
.boundary_edges <- function(poly) {
    g     <- st_geometry(poly)[[1]]
    rings <- if (inherits(g, "MULTIPOLYGON")) do.call(c, unclass(g)) else unclass(g)
    p1 <- do.call(rbind, lapply(rings, function(r) r[-nrow(r), 1:2, drop = FALSE]))
    p2 <- do.call(rbind, lapply(rings, function(r) r[-1,        1:2, drop = FALSE]))
    list(p1 = p1, p2 = p2)
}

#' Fraction of each candidate line's length outside the polygon, via a
#' vectorised M-lines x E-edges intersection pass chunked over lines (see
#' file header for why), plus a per-crossing-line sort-and-sum step (pure
#' arithmetic, one batched point-in-polygon call per chunk).
#' @param x1,y1,x2,y2 coordinate vectors, one entry per candidate line
#' @param edges list(p1, p2) boundary edges from .boundary_edges(), in the
#'   SAME (planar) CRS as the candidate lines
#' @param poly_u the (unioned) polygon, used only for the batched
#'   start-point-inside test
#' @param crs CRS to attach to the start-point geometries passed to
#'   st_intersects()
#' @param chunk_size lines processed per vectorised pass; NULL (default)
#'   picks automatically from currently-available memory
#'   (.choose_line_chunk_size()) - exposed mainly so tests can force a small
#'   value without needing a pathologically large E
#' @return numeric vector of length(x1), fraction outside in `[0, 1]`
#' @noRd
.frac_outside_vectorised <- function(x1, y1, x2, y2, edges, poly_u, crs, chunk_size = NULL) {
    M <- length(x1)
    E <- nrow(edges$p1)
    frac_outside <- numeric(M)
    if (M == 0 || E == 0) return(frac_outside)

    if (is.null(chunk_size)) chunk_size <- .choose_line_chunk_size(E)
    chunk_size <- max(1L, min(as.integer(chunk_size), M))

    ex1 <- edges$p1[, 1]; ey1 <- edges$p1[, 2]
    ex2 <- edges$p2[, 1]; ey2 <- edges$p2[, 2]
    sx  <- ex2 - ex1;     sy  <- ey2 - ey1

    for (cs in seq(1L, M, by = chunk_size)) {
        idx <- cs:min(cs + chunk_size - 1L, M)
        m   <- length(idx)

        X1 <- rep(x1[idx], each = E);              Y1 <- rep(y1[idx], each = E)
        RX <- rep(x2[idx] - x1[idx], each = E);     RY <- rep(y2[idx] - y1[idx], each = E)
        EX1 <- rep(ex1, times = m);                 EY1 <- rep(ey1, times = m)
        SX  <- rep(sx,  times = m);                  SY  <- rep(sy,  times = m)

        denom <- RX * SY - RY * SX
        dx <- EX1 - X1; dy <- EY1 - Y1
        t <- (dx * SY - dy * SX) / denom
        u <- (dx * RY - dy * RX) / denom

        valid <- is.finite(denom) & abs(denom) > 1e-12 &
            t > 1e-9 & t < 1 - 1e-9 & u >= -1e-9 & u <= 1 + 1e-9
        local_id <- rep(seq_len(m), each = E)

        cross_local <- unique(local_id[valid])
        if (length(cross_local) == 0) next
        cross_global <- idx[cross_local]

        # start-point in/out for just this chunk's crossing lines - one
        # batched call per chunk, not one for the whole M at once
        starts  <- st_sfc(lapply(cross_global, function(i) st_point(c(x1[i], y1[i]))), crs = crs)
        inside0 <- lengths(suppressWarnings(st_intersects(starts, poly_u))) > 0

        # per crossing line: sort crossing parameters, sum "outside"
        # sub-intervals by alternating in/out status from the known
        # start-point status
        ts_by_line <- split(t[valid], local_id[valid])
        for (k in seq_along(cross_local)) {
            tv <- sort(unique(ts_by_line[[as.character(cross_local[k])]]))
            bounds     <- c(0, tv, 1)
            seg_status <- inside0[k] != (seq_len(length(bounds) - 1) %% 2 == 0)
            frac_outside[cross_global[k]] <- min(max(sum(diff(bounds)[!seg_status]), 0), 1)
        }
    }
    frac_outside
}

#' GEOS-based fallback: one st_difference() call per line, run sequentially
#' - the safety net .compute_frac_outside() drops back to when the
#' vectorised path errors. Not parallelised across future.apply workers:
#' that only ever helps on this rare fallback path, and isn't worth its
#' own failure modes there. Row-level parallelism via shape_indices_sf()'s
#' `parallel_rows` is the one parallelism knob this package offers.
#' @param seg sfc of candidate line segments (all of them, not just idx_out)
#' @param len_total numeric vector, each line's total length (same order as
#'   `seg`)
#' @param idx_out integer indices into `seg`/`len_total` for the lines that
#'   actually need the expensive exact-geometry check (already known not to
#'   be fully inside)
#' @param poly_u the (unioned) polygon to clip against
#' @return numeric vector of length(len_total), fraction outside in `[0, 1]`
#'   (0 for every index not in `idx_out`)
#' @noRd
.frac_outside_geos <- function(seg, len_total, idx_out, poly_u) {
    frac_outside <- numeric(length(len_total))
    compute_frac <- function(k) {
        out <- suppressWarnings(st_length(st_difference(seg[k], poly_u)))
        out <- if (length(out) == 0 || is.na(out)) 0 else as.numeric(out)
        min(max(out / len_total[k], 0), 1)
    }
    frac_outside[idx_out] <- vapply(idx_out, compute_frac, numeric(1))
    frac_outside
}

#' Shared dispatcher used by both .eval_mesh_pairs() and
#' .random_line_index(): fills in `frac_outside[idx_out]` using the fast
#' vectorised path, falling back to the GEOS loop on error.
#' @param x1,x2 full Nx2 coordinate matrices (all candidate lines, not
#'   just idx_out) - matches the layout both callers already build for
#'   segment construction
#' @param idx_out integer indices (into x1/x2/seg/len_total) of lines that
#'   need the expensive check; every other index is left at frac_outside = 0
#' @param seg sfc of candidate line segments (all of them), used only by
#'   the GEOS fallback path
#' @param len_total numeric vector, each line's total length (same order as
#'   `seg`), used only by the GEOS fallback path
#' @param poly_u the (unioned) polygon to test against
#' @param crs CRS for any geometries built internally
#' @return numeric vector of length nrow(x1), fraction outside in `[0, 1]`
#' @noRd
.compute_frac_outside <- function(x1, x2, idx_out, seg, len_total, poly_u, crs) {
    frac_outside <- numeric(nrow(x1))
    if (length(idx_out) == 0) return(frac_outside)
    result <- tryCatch({
        edges <- .boundary_edges(poly_u)
        .frac_outside_vectorised(x1[idx_out, 1], x1[idx_out, 2],
                                  x2[idx_out, 1], x2[idx_out, 2],
                                  edges, poly_u, crs)
    }, error = function(e) NULL)
    if (is.null(result)) {
        frac_outside <- .frac_outside_geos(seg, len_total, idx_out, poly_u)
    } else {
        frac_outside[idx_out] <- result
    }
    frac_outside
}

#' Shared vectorised core: given an (n*q) x 2 matrix of quadrature points
#' (q per piece, triangle-major - row (i-1)*q+a is piece i's a-th point),
#' evaluate every piece-PAIR's mean fraction-outside over its q^2 sub-pair
#' lines. Used by .mesh_convexity_index() below.
#' @param qmat (n*q) x 2 matrix of quadrature point coordinates
#' @param q quadrature points per triangle (1 or 3)
#' @param n number of triangles
#' @param poly_u the (unioned) polygon to test lines against
#' @param crs CRS for internally-built geometries
#' @return list(edges, g_ij) - edges is a Mx2 matrix of piece-index pairs
#'   (M = choose(n, 2)), g_ij is a length-M vector of each pair's mean
#'   fraction outside across its q^2 sub-pair lines
#' @noRd
.eval_mesh_pairs <- function(qmat, q, n, poly_u, crs) {
    edges <- t(combn(n, 2))
    M   <- nrow(edges)
    qsq <- q * q

    ab       <- expand.grid(a = seq_len(q), b = seq_len(q))
    pair_rep <- rep(seq_len(M), each = qsq)
    row1     <- (edges[pair_rep, 1] - 1) * q + rep(ab$a, times = M)
    row2     <- (edges[pair_rep, 2] - 1) * q + rep(ab$b, times = M)
    x1 <- qmat[row1, , drop = FALSE]
    x2 <- qmat[row2, , drop = FALSE]
    m  <- nrow(x1)   # = M * q^2, total segments across all pairs

    # fast, vectorised segment construction
    if (requireNamespace("sfheaders", quietly = TRUE)) {
        seg_df <- data.frame(
            id = rep(seq_len(m), each = 2),
            x  = as.vector(rbind(x1[, 1], x2[, 1])),
            y  = as.vector(rbind(x1[, 2], x2[, 2]))
        )
        seg <- sfheaders::sfc_linestring(seg_df, x = "x", y = "y", linestring_id = "id")
        st_crs(seg) <- crs
    } else {
        seg <- mapply(function(k) st_linestring(rbind(x1[k, ], x2[k, ])), seq_len(m), SIMPLIFY = FALSE)
        seg <- st_sfc(seg, crs = crs)
    }

    len_total <- as.numeric(st_length(seg))

    # cheap short-circuit: flag lines already fully inside (prepared-geometry
    # predicate, safe to vectorise), so the expensive exact check below only
    # runs on the boundary-crossing minority
    inside  <- lengths(suppressWarnings(st_covered_by(seg, poly_u))) > 0
    idx_out <- which(!inside & len_total > 0)

    # results always assigned back by explicit index, never positional
    # recycling: a single vectorised st_difference() call can silently drop
    # empty-result rows, desyncing a positionally-assigned result
    frac_outside <- .compute_frac_outside(x1, x2, idx_out, seg, len_total, poly_u, crs)

    g_ij <- as.numeric(rowsum(frac_outside, pair_rep)) / qsq   # aggregate q^2 sub-pairs -> 1 value per piece-pair
    list(edges = edges, g_ij = g_ij)
}
