## =========================================================================
## 6. Running all indices over many polygons: independent rows, or one
## weighted collection
## =========================================================================
## byrow = TRUE (default): each row indexed independently. Triangulates
## once per row (shared between the CDT-based indices) and can
## parallelise at the row level (`parallel_rows`) - the only parallelism
## this package offers.
##
## byrow = FALSE: every row is a weighted sub-polygon of one overall
## shape, st_union() of the rows that aren't holes (weight 0/NA). Returns
## one row, so the caller supplies `id`. No separate rows to parallelise
## across in this mode - `parallel_rows = TRUE` errors rather than
## silently doing nothing.
##
## `which` selects a subset of the 13 available indices (default "all").
## Seven (convexity, moment_of_inertia, moment_isotropy, directional_balance,
## span, radial_concentration, depth) need a CDT mesh; six (hull_ratio,
## polsby_popper, width_length_ratio, reock, detour, exchange) are classic
## boundary/hull/bbox scores that don't. Requesting only the latter skips
## triangulation (and, for byrow = FALSE, the area-weighted overlay)
## entirely - a real cost saving, not just a display filter.

#' Filters `dots` down to the names a function actually accepts, so one
#' shared `...` can feed several functions with only partly-overlapping
#' arguments (e.g. `deterministic`/`n_quad`/`n_lines`, shared by
#' convexity_index() and span_index(), vs `weight`, accepted by both of
#' those plus radial_concentration_index()) without do.call() erroring on
#' an unused argument.
#' @param fn the target function
#' @param dots a named list
#' @return the subset of `dots` whose names are in `names(formals(fn))`
#' @noRd
.dots_for <- function(fn, dots) dots[names(dots) %in% names(formals(fn))]

#' Draws ONE pool of points to share across convexity_index()/span_index()/
#' radial_concentration_index()'s Monte Carlo modes when more than one is
#' requested together with deterministic = FALSE - the Monte Carlo
#' analogue of sharing one `prep` triangulation across indices, avoiding
#' 2-3 independent draws of the same size from what's ultimately the same
#' weighted density. Uses exactly the same sampling method each engine
#' would use drawing on its own (st_sample() unweighted, top-up loop
#' included; .sample_weighted_points() weighted), so a shared draw is
#' indistinguishable in distribution from any one of them sampling
#' independently - the only difference is paying for one draw instead of
#' up to three.
#' @param tri sf data frame of CDT triangles, with an `area` column
#' @param poly_u the (multi)polygon `tri` triangulates
#' @param weight optional numeric vector, one entry per triangle; NULL
#'   samples uniformly
#' @param n number of points to draw
#' @return n x 2 coordinate matrix, or NULL if sampling failed (caller
#'   falls back to each engine drawing independently, same as if this
#'   function had never run)
#' @noRd
.draw_shared_mc_points <- function(tri, poly_u, weight, n) {
    if (!is.null(weight)) {
        return(.sample_weighted_points(tri, weight, n))
    }
    pts <- st_sample(poly_u, size = n, type = "random", exact = TRUE)
    tries <- 0
    while (length(pts) < n && tries < 10) {
        pts <- c(pts, st_sample(poly_u, size = n, type = "random", exact = TRUE))
        tries <- tries + 1
    }
    if (length(pts) < n) return(NULL)
    pts <- pts[seq_len(n)]
    st_coordinates(pts)[, 1:2, drop = FALSE]
}

## canonical order - also what "all" expands to, and what output columns/
## vector names appear in regardless of the order `which` was given in
.ALL_SHAPE_INDICES <- c("convexity", "moment_of_inertia", "moment_isotropy", "directional_balance",
                        "span", "radial_concentration", "depth", "hull_ratio", "polsby_popper",
                        "width_length_ratio", "reock", "detour", "exchange")

.MESH_INDICES <- c("convexity", "moment_of_inertia", "moment_isotropy", "directional_balance",
                   "span", "radial_concentration", "depth")

## the six classic metrics have no weighted form at all (see
## classical-metrics.R's own file header for why a rearrangement-style
## weighted reference doesn't apply to any of them)
.CLASSIC_INDICES <- c("hull_ratio", "polsby_popper", "width_length_ratio", "reock", "detour", "exchange")

#' Resolves shape_indices()/shape_indices_sf()'s `which` argument.
#' @param which "all", or a character vector naming a subset of
#'   .ALL_SHAPE_INDICES
#' @return character vector, the requested subset in canonical order
#' @noRd
.resolve_which <- function(which) {
    if (identical(which, "all")) return(.ALL_SHAPE_INDICES)
    if (!is.character(which)) {
        stop("`which` must be \"all\" or a character vector of index names; got ",
             class(which)[1], ".")
    }
    unknown <- setdiff(which, .ALL_SHAPE_INDICES)
    if (length(unknown) > 0) {
        stop("Unknown index name(s) in `which`: ", paste(unknown, collapse = ", "),
             ". Valid choices: ", paste(.ALL_SHAPE_INDICES, collapse = ", "), ", or \"all\".")
    }
    .ALL_SHAPE_INDICES[.ALL_SHAPE_INDICES %in% which]
}

#' All indices (or a chosen subset) for a single (multi)polygon
#'
#' Shares one triangulation when any CDT-based index is requested.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @param which "all" (default), or a character vector naming a subset of
#'   these thirteen values - each listed here with the function it actually
#'   calls, since the `which` string and the function name aren't always
#'   identical:
#'
#'   * `"convexity"` - [convexity_index()]
#'   * `"moment_of_inertia"` - [moment_of_inertia_index()]
#'   * `"moment_isotropy"` - [moment_isotropy_index()]
#'   * `"directional_balance"` - [directional_balance_index()]
#'   * `"span"` - [span_index()]
#'   * `"radial_concentration"` - [radial_concentration_index()]
#'   * `"depth"` - [depth_index()]
#'   * `"hull_ratio"` - [hull_ratio_index()]
#'   * `"polsby_popper"` - [polsby_popper_index()]
#'   * `"width_length_ratio"` - [width_length_ratio_index()]
#'   * `"reock"` - [reock_index()]
#'   * `"detour"` - [detour_index()]
#'   * `"exchange"` - [exchange_index()]
#'
#'   The first seven need a CDT mesh; the last six are classic
#'   boundary/hull/bounding-box metrics that don't - requesting only those
#'   six skips triangulation entirely. An unrecognised name in `which`
#'   errors immediately, listing all thirteen valid values.
#' @param deterministic_max_tri if set, meshes with more triangles than
#'   this use `deterministic = FALSE` automatically - real-world polygons
#'   can triangulate into thousands of CDT triangles, where the
#'   exhaustive/full-subdivision methods are minutes-to-hours per polygon.
#'   Can only push `deterministic` from TRUE/default to FALSE, never
#'   override an explicit `deterministic = FALSE` back to TRUE. Shared by
#'   convexity_index()/span_index()/radial_concentration_index()/
#'   directional_balance_index()/depth_index() alike; has no effect on the
#'   six classic metrics, which have no mesh to switch at all.
#' @param simplify_tolerance passed to `prepare_polygon()` - see its own
#'   doc. Has no effect when `which` requests only the six classic
#'   metrics, since no mesh is built for those to begin with.
#' @param ... passed to convexity_index(), span_index(),
#'   radial_concentration_index(), directional_balance_index(), and
#'   depth_index() (e.g. `weight` is shared by all five, and also
#'   forwarded to moment_of_inertia_index()/moment_isotropy_index(), so
#'   every weighted index stays consistently weighted;
#'   `deterministic`/`n_lines` are shared by all five mesh-based Monte
#'   Carlo indices - `n_lines` sets every one of their sample counts with
#'   one argument, even though
#'   span_index()/radial_concentration_index()/directional_balance_index()/
#'   depth_index() never actually build line geometry; `n_quad` is
#'   convexity_index()/span_index() only, since the other three mesh
#'   indices' deterministic modes have no quadrature refinement of that
#'   kind to control).
#' @return named numeric vector, one entry per requested index, in
#'   canonical order
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' wake <- nc[nc$NAME == "Wake", ]
#' shape_indices(wake)
#'
#' # a subset - skips triangulation entirely, since none of these need a mesh
#' shape_indices(wake, which = c("hull_ratio", "polsby_popper", "reock"))
#'
#' # deterministic_max_tri forces the cheaper Monte Carlo estimator once
#' # the mesh exceeds it, regardless of how simple the polygon actually is -
#' # here set deliberately low (5) just to demonstrate the switch on a
#' # small county
#' shape_indices(wake, deterministic_max_tri = 5, n_lines = 2000, seed = 1)
#'
#' # simplify_tolerance trades a small amount of boundary detail (here 50m)
#' # for a smaller mesh - passed straight through to prepare_polygon()
#' shape_indices(wake, simplify_tolerance = 50)
#' @export
shape_indices <- function(poly, which = "all", deterministic_max_tri = NULL,
                           simplify_tolerance = NULL, ...) {
    which <- .resolve_which(which)
    needs_mesh <- any(.MESH_INDICES %in% which)
    prep <- if (needs_mesh) prepare_polygon(poly, simplify_tolerance = simplify_tolerance) else NULL
    poly_clean <- if (!is.null(prep)) prep$poly else NULL

    dots <- list(...)
    # isFALSE(), not !is.null(): deterministic_max_tri only pushes TRUE ->
    # FALSE, never overrides an explicit deterministic = FALSE back to TRUE
    if (!is.null(deterministic_max_tri) && needs_mesh && !isFALSE(dots$deterministic)) {
        n_tri <- if (is.null(prep$tri)) 0 else nrow(prep$tri)
        dots$deterministic <- n_tri <= deterministic_max_tri
    }

    # share one drawn point sample across every requested Monte Carlo mesh
    # index instead of each drawing independently (see
    # .draw_shared_mc_points()'s own comments) - only worth it once there's
    # more than one such index actually requested, and only possible once
    # there's a real mesh to sample from
    n_tri_mesh <- if (is.null(prep$tri)) 0 else nrow(prep$tri)
    mc_indices <- intersect(which, c("convexity", "span", "radial_concentration", "directional_balance", "depth"))
    if (n_tri_mesh > 0 && isFALSE(dots$deterministic) && length(mc_indices) > 1) {
        n_lines_shared <- dots$n_lines %||% 3000
        if (is.function(n_lines_shared)) n_lines_shared <- max(1L, round(n_lines_shared(n_tri_mesh)))
        if (!is.null(dots$seed)) set.seed(dots$seed)
        shared_points <- .draw_shared_mc_points(prep$tri, prep$poly, dots$weight, 2 * n_lines_shared)
        if (!is.null(shared_points)) dots$points <- shared_points
    }

    out <- c()
    if ("convexity" %in% which) {
        out["convexity"] <- do.call(convexity_index, c(list(poly = poly, prep = prep), .dots_for(convexity_index, dots)))$index
    }
    if ("moment_of_inertia" %in% which) {
        out["moment_of_inertia"] <- moment_of_inertia_index(poly, prep = prep, weight = dots$weight)$index
    }
    if ("moment_isotropy" %in% which) {
        out["moment_isotropy"] <- moment_isotropy_index(poly, prep = prep, weight = dots$weight)$index
    }
    if ("directional_balance" %in% which) {
        out["directional_balance"] <- do.call(directional_balance_index, c(list(poly = poly, prep = prep), .dots_for(directional_balance_index, dots)))$index
    }
    if ("span" %in% which) {
        out["span"] <- do.call(span_index, c(list(poly = poly, prep = prep), .dots_for(span_index, dots)))$index
    }
    if ("radial_concentration" %in% which) {
        out["radial_concentration"] <- do.call(radial_concentration_index, c(list(poly = poly, prep = prep), .dots_for(radial_concentration_index, dots)))$index
    }
    if ("depth" %in% which) {
        out["depth"] <- do.call(depth_index, c(list(poly = poly, prep = prep), .dots_for(depth_index, dots)))$index
    }
    # poly_clean, not poly: the classic metrics have no `prep` arg, so they
    # need the already-cleaned geometry explicitly (when it exists) or
    # every requested index could disagree about what it's measuring
    if ("hull_ratio" %in% which) out["hull_ratio"] <- hull_ratio_index(poly_clean %||% poly)$index
    if ("polsby_popper" %in% which) out["polsby_popper"] <- polsby_popper_index(poly_clean %||% poly)$index
    if ("width_length_ratio" %in% which) out["width_length_ratio"] <- width_length_ratio_index(poly_clean %||% poly)$index
    if ("reock" %in% which) out["reock"] <- reock_index(poly_clean %||% poly)$index
    if ("detour" %in% which) out["detour"] <- detour_index(poly_clean %||% poly)$index
    if ("exchange" %in% which) out["exchange"] <- exchange_index(poly_clean %||% poly)$index

    out[which]
}

## byrow = FALSE support: .weighted_mesh() below triangulates the union of
## non-hole rows ONCE (ties triangle count to the union's boundary, not
## row count), so a triangle can straddle several rows; weight is then
## assigned by area-weighted overlay against the original row boundaries -
## the TOTAL is conserved exactly, but a triangle straddling rows of
## different density gets their average baked in uniformly, which is an
## approximation at the per-triangle level. When weights is genuinely
## supplied, .shape_indices_sf_grouped() tries the exact constrained-
## triangulation mesh first instead (constrained-triangulation.R), where
## no triangle can ever straddle a row boundary at all; .weighted_mesh()
## remains the automatic fallback for weights = NULL (density is then
## uniform, so this approximation is a no-op) and for input that mesh
## can't safely handle (see its own file header).

#' Resolves `weights` and excludes hole rows (weight 0/NA), unioning
#' what's left - the common first step for shape_indices_sf(byrow =
#' FALSE), split out from the actual triangulation so a `which` request
#' naming only the classic (non-mesh) metrics can skip triangulating
#' entirely.
#'
#' `simplify_tolerance`, if given, simplifies `poly_u` ONLY - never the
#' returned `geoms` (the original per-row boundaries), which the weight
#' overlay (.weighted_mesh()) still needs at full precision. Simplifying
#' each row independently BEFORE union is a different, worse operation:
#' adjacent rows sharing an edge get that edge simplified differently on
#' each side, breaking what used to be a clean tessellation into thousands
#' of sliver gaps - which inflates poly_u's total boundary edge count
#' rather than shrinking it (verified: a real 50,424-row dataset that
#' unions to ONE clean ~50K-vertex polygon produces 1,851 disjoint sliver
#' parts instead, after independently simplifying every row first). Doing
#' it here, once, after the union already dissolved every shared internal
#' edge, has no such failure mode.
#' @param x an sf data frame, already projected
#' @param weights NULL, a column-name string, or a numeric vector - see
#'   .resolve_row_weights()
#' @param simplify_tolerance optional st_simplify() tolerance, applied to
#'   poly_u only - see above. Always in `x`'s own CRS unit at this point,
#'   NOT necessarily metres - see prepare_polygon()'s own doc for why
#'   (same caveat, same warning if that unit isn't metres).
#' @return list(geoms, w_row, poly_u, raw_total) - geoms/w_row are the
#'   kept (non-hole) rows' geometry (unsimplified) and resolved weight,
#'   poly_u is st_union() of them (simplified if requested), raw_total is
#'   sum(weights) before normalisation (NA/hole rows excluded)
#' @noRd
.resolve_union <- function(x, weights, simplify_tolerance = NULL) {
    n_row <- nrow(x)
    geoms <- .make_valid_warn(st_geometry(x))
    w_row <- .resolve_row_weights(x, weights, n_row, geoms)

    is_hole <- is.na(w_row) | w_row == 0
    if (any(is_hole)) {
        warning(sprintf("%d of %d rows have zero or NA weight; ", sum(is_hole), n_row),
                "treating them as holes - excluded from the triangulated union entirely ",
                "(not just zero-weighted area). Pass weights without zeros/NA for those ",
                "rows if you want that area included with zero weight instead.")
    }
    keep <- !is_hole
    if (!any(keep)) {
        stop("Every row has zero or NA weight; nothing left to triangulate.")
    }

    poly_u <- st_union(geoms[keep])
    if (!is.null(simplify_tolerance)) {
        .warn_if_nonmetric_tolerance(poly_u, simplify_tolerance)
        poly_u <- st_simplify(poly_u, dTolerance = simplify_tolerance, preserveTopology = TRUE)
        poly_u <- .make_valid_warn(poly_u)
    }

    list(geoms = geoms[keep], w_row = w_row[keep], poly_u = poly_u,
         raw_total = sum(w_row, na.rm = TRUE))
}

#' Builds the combined triangle mesh + weights for shape_indices_sf(byrow
#' = FALSE) via area-weighted overlay: triangulate poly_u once, then
#' st_intersection() each triangle against the row boundaries and
#' area-weight-average their density onto it. Used unconditionally when
#' weights = NULL (density is then uniform, so there's nothing to smear
#' regardless of mesh coarseness) and as the automatic fallback, with a
#' warning, when weights is genuinely supplied but the exact constrained-
#' triangulation mesh (.constrained_weighted_mesh(), constrained-
#' triangulation.R) can't be built safely - see that file's own header
#' for why this overlay can smear density across triangle boundaries, and
#' .shape_indices_sf_grouped() for exactly when each path is chosen. See
#' .resolve_union() for the hole-exclusion/simplification step this
#' builds on.
#' @param x an sf data frame, already projected
#' @param weights NULL, a column-name string, or a numeric vector - see
#'   .resolve_row_weights()
#' @param simplify_tolerance see .resolve_union() - applied to poly_u
#'   before triangulating; the weight overlay below still uses the
#'   original, unsimplified row geometries, so a triangle right at a
#'   simplified boundary can capture slightly more or less of the true
#'   row-level area than one from the unsimplified mesh would, bounded by
#'   the tolerance and typically negligible relative to the whole shape
#' @return list(pieces, poly_u, raw_total) - pieces is an sf with one row
#'   per triangle of poly_u plus a `weight` column summing to 1, poly_u is
#'   st_union() of the non-hole rows, raw_total is sum(weights) before
#'   normalisation (NA/hole rows excluded)
#' @noRd
.weighted_mesh <- function(x, weights, simplify_tolerance = NULL) {
    ru <- .resolve_union(x, weights, simplify_tolerance)
    row_area    <- as.numeric(st_area(ru$geoms))
    row_density <- ru$w_row / row_area

    tri <- cdt_triangles(ru$poly_u)
    if (is.null(tri) || nrow(tri) == 0) {
        stop("Union of kept rows triangulated to no triangles; nothing to compute.")
    }

    # area-weighted overlay: intersect the coarse mesh against the fine row
    # boundaries, weight each piece by its row's density, sum per triangle
    rows_kept <- st_sf(density = row_density, geometry = ru$geoms)
    tri_sfc   <- st_sf(tri_id = seq_len(nrow(tri)), geometry = st_geometry(tri))
    ov <- suppressWarnings(st_intersection(tri_sfc, rows_kept))
    ov_weight     <- as.numeric(st_area(ov)) * ov$density
    tri_weight_by_id <- tapply(ov_weight, ov$tri_id, sum)
    tri_weight <- numeric(nrow(tri))
    tri_weight[as.integer(names(tri_weight_by_id))] <- as.numeric(tri_weight_by_id)

    tri$weight <- .normalize_weight(tri_weight)

    list(pieces = tri, poly_u = ru$poly_u, raw_total = ru$raw_total)
}

#' Coerces a value that's supposed to already be numeric weights - errors
#' instead of silently corrupting the data the way as.numeric() would for
#' a factor or a non-numeric character vector.
#' @param v the candidate weight vector
#' @param what short description of `v`, for the error message
#' @param code a valid R expression (as a string) evaluating to `v`, used
#'   in the fix-it snippet - defaults to `what`, which is only valid code
#'   when `what` is itself a bare variable name
#' @return `v` coerced to numeric
#' @noRd
.coerce_weight_vector <- function(v, what, code = what) {
    if (is.factor(v)) {
        stop(what, " is a factor; as.numeric() on a factor silently returns its integer ",
             "level codes, not the labelled values. Convert explicitly first, e.g. ",
             "as.numeric(as.character(", code, ")).")
    }
    if (!is.numeric(v)) {
        stop(what, " must be numeric, got ", class(v)[1], ".")
    }
    as.numeric(v)
}

#' Resolves shape_indices_sf(byrow = FALSE)'s `weights` argument to one
#' validated numeric value per row of `x`.
#' @param x an sf data frame (used to look up a named weight column)
#' @param weights NULL (each row's own area), a column-name string, or a
#'   numeric vector (one entry per row)
#' @param n_row expected length (nrow(x))
#' @param geoms `x`'s already-validated/projected geometry
#' @return numeric vector of length `n_row`, non-negative, no Inf/-Inf - NA
#'   is allowed (treated as a hole)
#' @noRd
.resolve_row_weights <- function(x, weights, n_row, geoms) {
    w_row <- if (is.null(weights)) {
        as.numeric(st_area(geoms))
    } else if (is.character(weights) && length(weights) == 1) {
        x_attr <- sf::st_drop_geometry(x)
        if (!weights %in% names(x_attr)) {
            # checked against attribute columns, not names(x) - the
            # geometry column name is technically "in" names(x) too
            stop("weights = \"", weights, "\" is not a column of x", if (weights %in% names(x)) {
                paste0(" - it's x's geometry column (", weights, "), not an attribute column, so ",
                       "there's no numeric data there to use as a weight")
            }, ".")
        }
        .coerce_weight_vector(x_attr[[weights]], sprintf("column \"%s\"", weights),
                               code = sprintf('x[["%s"]]', weights))
    } else {
        .coerce_weight_vector(weights, "weights")
    }
    if (length(w_row) != n_row) {
        stop("weights must have one entry per row of x (", n_row, "), got ", length(w_row), ".")
    }
    if (any(!is.na(w_row) & !is.finite(w_row))) {
        stop("weights must be finite (no Inf/-Inf).")
    }
    if (any(w_row < 0, na.rm = TRUE)) {
        stop("weights must be non-negative.")
    }
    w_row
}

#' All indices (or a chosen subset) for every feature in an sf data frame.
#' @param x an sf data frame
#' @param byrow if TRUE (default), index each row separately (`weights`
#'   and `id` must be NULL). If FALSE, treat every row as a weighted
#'   sub-polygon of one overall shape (`st_union(x)`) and return a
#'   single-row sf.
#' @param which "all" (default), or a character vector naming a subset of
#'   these thirteen values - each listed here with the function it actually
#'   calls, since the `which` string and the function name aren't always
#'   identical:
#'
#'   * `"convexity"` - [convexity_index()]
#'   * `"moment_of_inertia"` - [moment_of_inertia_index()]
#'   * `"moment_isotropy"` - [moment_isotropy_index()]
#'   * `"directional_balance"` - [directional_balance_index()]
#'   * `"span"` - [span_index()]
#'   * `"radial_concentration"` - [radial_concentration_index()]
#'   * `"depth"` - [depth_index()]
#'   * `"hull_ratio"` - [hull_ratio_index()]
#'   * `"polsby_popper"` - [polsby_popper_index()]
#'   * `"width_length_ratio"` - [width_length_ratio_index()]
#'   * `"reock"` - [reock_index()]
#'   * `"detour"` - [detour_index()]
#'   * `"exchange"` - [exchange_index()]
#'
#'   The first seven need a CDT mesh; the last six are classic
#'   boundary/hull/bounding-box metrics that don't - requesting only those
#'   six skips triangulation entirely, in both `byrow` modes. An
#'   unrecognised name in `which` errors immediately, listing all thirteen
#'   valid values.
#' @param weights only used when `byrow = FALSE`. NULL (default) weights
#'   each row by its own area; otherwise a numeric vector or column name.
#'   Weight 0/NA excludes that row as a HOLE (with a warning) rather than
#'   zero-weighting it - since this changes `poly_u` itself, even
#'   `hull_ratio_index` can differ between calls that only differ in `weights`.
#'   When rows genuinely differ in weight, the combined mesh is built by
#'   an exact constrained triangulation (every kept row's own boundary
#'   supplied as a triangulation constraint, so no triangle can ever
#'   straddle a row boundary and density allocation is exact) rather than
#'   triangulating the union once and area-weight-averaging row density
#'   onto whatever it overlaps. Falls back to that coarser overlay
#'   automatically, with a warning, if kept rows genuinely overlap or
#'   triangulation otherwise fails; `weights = NULL` never needs either
#'   approach, since every row then weights by its own area and density
#'   is uniform everywhere.
#' @param id required when `byrow = FALSE`: a single scalar value
#'   identifying this group in the collapsed single-row result.
#' @param parallel_rows only meaningful when `byrow = TRUE` (errors if
#'   TRUE with `byrow = FALSE`). If TRUE (default) and furrr is
#'   installed, rows run across the active `future::plan()` - pass
#'   `deterministic_max_tri` too on real-world polygons first, or
#'   unbounded per-row cost plus parallelism turns "slow" into "OOM crash".
#' @param ... passed to shape_indices() when `byrow = TRUE`. When
#'   `byrow = FALSE`: `n_lines`, `seed` are read the same way and set
#'   every mesh-based Monte Carlo index's sample count together; `n_quad`
#'   likewise, but only affects convexity/span (radial_concentration's,
#'   directional_balance's, and depth's deterministic modes have no
#'   quadrature refinement of that kind to control); `deterministic` picks
#'   the method directly (default TRUE), shared by convexity, span,
#'   radial_concentration, directional_balance, and depth alike;
#'   `deterministic_max_tri` can push `deterministic` from TRUE to FALSE
#'   based on the combined triangle count, but never overrides an
#'   explicit `deterministic = FALSE`; on the exact constrained mesh
#'   above, convexity's/span's deterministic mode is additionally capped
#'   by a safety ceiling derived from memory actually available on this
#'   machine right now (that mesh can be far larger than the union-only
#'   one, since it's built at the rows' own boundary resolution) - left
#'   at its default, a mesh above that ceiling falls back to the Monte
#'   Carlo estimator silently; an explicit `deterministic_max_tri` that
#'   would still exceed the safe ceiling errors instead, rather than
#'   either attempting the computation or silently overriding what was
#'   asked for; `simplify_tolerance` simplifies
#'   `st_union(x)` itself before triangulating - do not simplify the rows
#'   yourself beforehand instead: adjacent rows sharing an edge get that
#'   edge simplified differently on each side, fragmenting the union into
#'   slivers (see `vignette("a-basic-usage")`'s simplify_tolerance
#'   section).
#' @return if `byrow = TRUE`: x with one new `<name>_index` column per
#'   requested index. If `byrow = FALSE`: a single-row sf with id, those
#'   same columns, and total_weight.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#'
#' # byrow = TRUE (default): index each county independently
#' res <- shape_indices_sf(nc[1:5, ])
#' res$convexity_index
#'
#' # a subset, requested the same way as shape_indices()
#' shape_indices_sf(nc[1:5, ], which = c("hull_ratio", "reock"))
#'
#' # byrow = TRUE, with deterministic_max_tri forcing the Monte Carlo
#' # estimator once a row's own mesh exceeds it - `...` passes n_lines/seed
#' # through to shape_indices() for every row
#' res_rli <- shape_indices_sf(nc[1:5, ], byrow = TRUE, deterministic_max_tri = 5,
#'                              n_lines = 2000, seed = 1)
#' res_rli$convexity_index
#'
#' # byrow = FALSE: treat several adjacent counties as one weighted shape -
#' # Wake, Durham, Orange and Chatham (the Research Triangle) are contiguous
#' triangle <- nc[nc$NAME %in% c("Wake", "Durham", "Orange", "Chatham"), ]
#' shape_indices_sf(triangle, byrow = FALSE, id = "triangle_by_area")
#' shape_indices_sf(triangle, byrow = FALSE, weights = "BIR74",
#'                   id = "triangle_by_births")
#'
#' # weight 0/NA treats a row as a hole - excluded from the shape itself,
#' # not just zero-weighted - e.g. dropping Chatham (the least urban/most
#' # rural of the four) out of the footprint entirely
#' w <- ifelse(triangle$NAME == "Chatham", 0, triangle$BIR74)
#' shape_indices_sf(triangle, byrow = FALSE, weights = w, id = "triangle_minus_chatham")
#'
#' # excluding an INTERIOR row (rather than one on the edge, as above) can
#' # punch an actual hole through the union, not just trim the boundary. A
#' # 3x3 grid of unit squares is a solid convex block when every cell is
#' # kept; excluding only the centre cell turns it into a square ring (CI
#' # 1 -> 0.86, hull_ratio_index 1 -> 0.89 = 8/9) - two calls differing only in
#' # `weights` can disagree on hull_ratio_index despite hull_ratio_index having no
#' # weighted form of its own.
#' cell <- function(cx, cy) sf::st_polygon(list(rbind(
#'   c(cx - 0.5, cy - 0.5), c(cx + 0.5, cy - 0.5),
#'   c(cx + 0.5, cy + 0.5), c(cx - 0.5, cy + 0.5), c(cx - 0.5, cy - 0.5))))
#' ctr <- expand.grid(x = -1:1, y = -1:1)
#' grid9 <- sf::st_sf(pop = ifelse(ctr$x == 0 & ctr$y == 0, 0, 10),
#'                     geometry = sf::st_sfc(mapply(cell, ctr$x, ctr$y, SIMPLIFY = FALSE),
#'                                            crs = 3857))
#' shape_indices_sf(grid9, byrow = FALSE, id = "solid_block")            # weights = NULL: keeps all 9
#' shape_indices_sf(grid9, byrow = FALSE, weights = "pop", id = "ring")  # centre cell excluded
#' @export
shape_indices_sf <- function(x, byrow = TRUE, which = "all", weights = NULL, id = NULL,
                              parallel_rows = TRUE, ...) {
    which <- .resolve_which(which)
    if (byrow) {
        if (!is.null(weights)) {
            stop("weights is only used when byrow = FALSE; with byrow = TRUE each row ",
                 "is indexed on its own, unweighted.")
        }
        if (!is.null(id)) {
            stop("id is only used when byrow = FALSE; with byrow = TRUE each row of x ",
                 "already identifies itself (it's returned as its own row), there's ",
                 "nothing extra to attach.")
        }
        return(.shape_indices_sf_byrow(x, which, parallel_rows, ...))
    }
    if (!missing(parallel_rows) && parallel_rows) {
        stop("parallel_rows = TRUE has no meaning when byrow = FALSE - there are no ",
             "separate rows to distribute (byrow = FALSE is a single combined ",
             "computation over the whole collection). Leave parallel_rows unset, or ",
             "pass parallel_rows = FALSE explicitly.")
    }
    .shape_indices_sf_grouped(x, which, weights, id, ...)
}

#' byrow = TRUE path: one row of x -> one row out, each indexed
#' independently and unweighted.
#' @param x an sf data frame
#' @param which see shape_indices_sf()'s own doc
#' @param parallel_rows see shape_indices_sf()'s own doc
#' @param ... passed to shape_indices() for every row
#' @return x with one new `<name>_index` column per requested index
#' @noRd
.shape_indices_sf_byrow <- function(x, which, parallel_rows, ...) {
    geoms <- st_geometry(x)
    dots  <- c(list(...), list(which = which))

    # geoms[i], not lapply(geoms, .) directly - the latter strips to a bare
    # sfg with no CRS of its own
    compute_one <- function(i) do.call(shape_indices, c(list(poly = geoms[i]), dots))
    idx <- seq_along(geoms)

    use_parallel <- parallel_rows && requireNamespace("furrr", quietly = TRUE) &&
        requireNamespace("future", quietly = TRUE)

    if (use_parallel) {
        if (inherits(future::plan(), "sequential")) {
            warning("parallel_rows = TRUE but no parallel future::plan() is active ",
                    "(still on the default sequential plan) - running in order. ",
                    "Call future::plan(future::multisession, workers = ...) first ",
                    "for actual multi-core speedup.")
        }
        res <- furrr::future_map(idx, compute_one,
                                 .options = furrr::furrr_options(scheduling = 1))
    } else {
        res <- lapply(idx, compute_one)
    }

    res_m <- do.call(rbind, res)
    # dynamic, not hardcoded per index - `which` can be any subset
    for (nm in colnames(res_m)) x[[paste0(nm, "_index")]] <- res_m[, nm]
    x
}

#' byrow = FALSE path: every row of x -> one weighted sub-polygon of a
#' single overall shape, st_union(x). See .weighted_mesh()/.resolve_union()
#' for the mesh-building step this dispatches on.
#' @param x an sf data frame
#' @param which see shape_indices_sf()'s own doc
#' @param weights see shape_indices_sf()'s own doc
#' @param id see shape_indices_sf()'s own doc
#' @param ... n_quad, n_lines, deterministic, deterministic_max_tri, seed,
#'   simplify_tolerance - see shape_indices_sf()'s own doc
#' @return single-row sf: id, one `<name>_index` column per requested
#'   index, total_weight, geometry = st_union(x)
#' @noRd
.shape_indices_sf_grouped <- function(x, which, weights, id, ...) {
    if (is.null(id)) {
        stop("byrow = FALSE collapses every row of x into a single result, so there is ",
             "no way afterwards to tell which group of rows it came from unless you say ",
             "so up front. Construct an id for this group yourself (e.g. the ",
             "metro/cluster/group value these rows share) and pass it as `id`.")
    }
    if (length(id) != 1) {
        stop("id must be a single value identifying this group, got length ", length(id), ".")
    }

    dots <- list(...)
    n_quad                <- dots$n_quad %||% 3
    n_lines               <- dots$n_lines %||% 3000
    deterministic_max_tri <- dots$deterministic_max_tri
    deterministic         <- dots$deterministic %||% TRUE
    seed                  <- dots$seed
    simplify_tolerance    <- dots$simplify_tolerance

    # project the WHOLE collection to one shared CRS before triangulating
    # or unioning rows against each other
    x <- .ensure_projected(x)

    needs_mesh <- any(.MESH_INDICES %in% which)
    array_mesh <- NULL
    if (needs_mesh) {
        # the exact constrained-triangulation mesh (see constrained-
        # triangulation.R's own file header) only fixes a defect that
        # exists in the first place because DIFFERENT rows can have
        # DIFFERENT density - with weights = NULL every row's weight is
        # its own area (.resolve_row_weights()), so density is uniform
        # everywhere and .weighted_mesh()'s coarse overlay can't smear
        # anything; only attempt it when weights is genuinely supplied
        if (!is.null(weights)) {
            cwm <- .constrained_weighted_mesh(x, weights, simplify_tolerance)
            if (isTRUE(cwm$ok)) {
                array_mesh <- cwm
            } else {
                warning(
                    "Falling back to the coarse (row-boundary-unaware) weighted mesh: ",
                    cwm$reason, ". Density may be smeared across triangle boundaries near ",
                    "rows with very different densities - see .weighted_mesh()'s own ",
                    "documentation for what this means for the affected indices."
                )
            }
        }

        if (is.null(array_mesh)) {
            wm        <- .weighted_mesh(x, weights, simplify_tolerance)
            pieces    <- wm$pieces
            poly_u    <- wm$poly_u
            raw_total <- wm$raw_total
            n_tri     <- nrow(pieces)
        } else {
            poly_u    <- array_mesh$poly_u
            raw_total <- array_mesh$raw_total
            n_tri     <- nrow(array_mesh$T)
        }

        # point-cloud indices (moment_of_inertia/moment_isotropy always;
        # radial_concentration/directional_balance's own deterministic
        # mode) are O(n) regardless of mesh size - same one-directional
        # deterministic_max_tri precedence as shape_indices(), unaffected
        # by how large the exact constrained mesh gets
        use_deterministic <- if (!is.null(deterministic_max_tri) && deterministic) {
            n_tri <= deterministic_max_tri
        } else {
            deterministic
        }

        # convexity/span's deterministic mode is O(n^2) in triangle pairs.
        # On the fallback/unweighted (sf pieces) mesh this is exactly
        # today's existing behaviour - deterministic_max_tri is the only
        # guard, same as use_deterministic above. The exact constrained
        # mesh can be orders of magnitude larger, so it's additionally
        # capped by a memory-derived safety ceiling - see
        # .safe_deterministic_tri_ceiling()'s own comments for why a fixed
        # constant is the wrong shape for that. Below the ceiling, nothing
        # changes; above it, default settings silently and safely fall
        # back to Monte Carlo (no user burden), while an explicit
        # deterministic_max_tri that would have permitted deterministic
        # mode here fails gracefully instead of either attempting the
        # computation or silently overriding what was actually asked for.
        if (is.null(array_mesh)) {
            use_deterministic_pairwise <- use_deterministic
        } else if (!deterministic) {
            use_deterministic_pairwise <- FALSE
        } else {
            safe_ceiling <- min(.safe_deterministic_tri_ceiling(n_quad, "convexity"),
                                 .safe_deterministic_tri_ceiling(n_quad, "span"))
            if (n_tri <= safe_ceiling) {
                use_deterministic_pairwise <- if (!is.null(deterministic_max_tri)) {
                    n_tri <= deterministic_max_tri
                } else {
                    TRUE
                }
            } else if (is.null(deterministic_max_tri)) {
                use_deterministic_pairwise <- FALSE
            } else if (n_tri <= deterministic_max_tri) {
                stop(sprintf(
                    paste("deterministic_max_tri = %d permits convexity/span's deterministic",
                          "(exhaustive, O(n^2)-in-triangle-pairs) mode at this mesh's %d",
                          "triangles, but that exceeds the memory-safe ceiling (%d triangles at",
                          "n_quad = %d) estimated from memory actually available on this",
                          "machine right now. Lower deterministic_max_tri, lower n_quad, set",
                          "deterministic = FALSE to use the random-line/random-pair estimator",
                          "instead, or leave deterministic_max_tri unset to fall back to it",
                          "automatically."),
                    deterministic_max_tri, n_tri, safe_ceiling, n_quad))
            } else {
                use_deterministic_pairwise <- FALSE
            }
        }

        # share one drawn point sample across every requested Monte Carlo
        # mesh index instead of each drawing independently (see
        # .draw_shared_mc_points()'s own comments) - only applies on the
        # sf-pieces mesh, where every mesh index still shares one
        # determinism flag exactly as before
        mc_indices <- intersect(which, c("convexity", "span", "radial_concentration", "directional_balance", "depth"))
        shared_points <- NULL
        if (is.null(array_mesh) && !use_deterministic && length(mc_indices) > 1) {
            if (!is.null(seed)) set.seed(seed)
            shared_points <- .draw_shared_mc_points(pieces, poly_u, pieces$weight, 2 * n_lines)
        }

        # array-mesh Monte Carlo point draws, split the same way the
        # determinism flags above are: convexity/span (pairwise, 2*n_lines
        # points each pair needs) and radial_concentration/
        # directional_balance (point-cloud, n_lines points) can disagree
        # on determinism now that only the former is ceiling-gated, so
        # they can no longer always share one draw the way the sf-pieces
        # path does
        array_pairwise_points   <- NULL
        array_pointcloud_points <- NULL
        array_pieces_lazy       <- NULL
        if (!is.null(array_mesh)) {
            if (any(c("convexity", "span") %in% which) && !use_deterministic_pairwise) {
                if (!is.null(seed)) set.seed(seed)
                array_pairwise_points <- .sample_weighted_points_array(
                    array_mesh$P, array_mesh$T, array_mesh$tri_weight, 2 * n_lines)
            }
            if (any(c("radial_concentration", "directional_balance", "depth") %in% which) && !use_deterministic) {
                if (!is.null(seed)) set.seed(seed)
                array_pointcloud_points <- .sample_weighted_points_array(
                    array_mesh$P, array_mesh$T, array_mesh$tri_weight, n_lines)
            }
            # built once, lazily, and reused by both convexity and span
            # when both are requested and both landed in deterministic
            # mode - see .array_mesh_to_sf_pieces()'s own comments for why
            # this is only ever attempted below the safety ceiling, where
            # it's cheap
            if (use_deterministic_pairwise && any(c("convexity", "span") %in% which)) {
                array_pieces_lazy <- .array_mesh_to_sf_pieces(
                    array_mesh$P, array_mesh$T, array_mesh$tri_area, array_mesh$crs)
            }
        }
    } else {
        ru        <- .resolve_union(x, weights, simplify_tolerance)
        poly_u    <- ru$poly_u
        raw_total <- ru$raw_total
    }

    # after weight resolution succeeds (not before - an invalid `weights`
    # should just error, not also print this first)
    requested_classic <- .CLASSIC_INDICES[.CLASSIC_INDICES %in% which]
    if (!is.null(weights) && length(requested_classic) > 0) {
        warning(
            paste(requested_classic, collapse = "/"), " don't use `weights` in their own ",
            "formula at all - only whether a row's weight is exactly 0/NA (excluded as a ",
            "hole, which changes poly_u itself) has any effect. The magnitude of a nonzero ",
            "weight (5 vs 5000) is completely invisible to any of the six classic metrics; ",
            "if you meant them to reflect the weighting itself, they won't."
        )
    }

    out <- list(id = id)
    if ("convexity" %in% which) {
        out$convexity_index <- if (!is.null(array_mesh)) {
            if (use_deterministic_pairwise) {
                .mesh_convexity_index(array_pieces_lazy, poly_u, n_quad = n_quad, plot = FALSE,
                                       weight = array_mesh$tri_weight)$index
            } else {
                .random_line_index(poly_u, n_lines = n_lines,
                                    prep = list(poly = poly_u, tri = data.frame(area = array_mesh$tri_area)),
                                    seed = seed, plot = FALSE, weight = array_mesh$tri_weight,
                                    points = array_pairwise_points)$index
            }
        } else if (use_deterministic) {
            .mesh_convexity_index(pieces, poly_u, n_quad = n_quad, plot = FALSE,
                                   weight = pieces$weight)$index
        } else {
            .random_line_index(poly_u, n_lines = n_lines, prep = list(poly = poly_u, tri = pieces),
                                seed = seed, plot = FALSE, weight = pieces$weight, points = shared_points)$index
        }
    }
    # shared between moment_of_inertia and moment_isotropy, computed at most
    # once (they read different fields off the same core result) since
    # both are in the default "all" request together
    moi_core <- NULL
    if (any(c("moment_of_inertia", "moment_isotropy") %in% which)) {
        moi_core <- if (!is.null(array_mesh)) {
            .moment_of_inertia_core_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                           array_mesh$tri_weight / array_mesh$tri_area, array_mesh$crs)
        } else {
            .moment_of_inertia_core(st_geometry(pieces), pieces$area, pieces$weight / pieces$area, poly_u)
        }
    }
    if ("moment_of_inertia" %in% which) {
        out$moment_of_inertia_index <- moi_core$index
    }
    if ("moment_isotropy" %in% which) {
        eig <- .principal_moments(moi_core$Ixx, moi_core$Iyy, moi_core$Ixy)
        out$moment_isotropy_index <- eig$lambda_min / eig$lambda_max
    }
    if ("directional_balance" %in% which) {
        out$directional_balance_index <- if (!is.null(array_mesh)) {
            if (use_deterministic) {
                .mesh_directional_balance_index_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                                        array_mesh$tri_weight, array_mesh$crs)$index
            } else {
                .random_point_directional_balance_index_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                                                array_mesh$tri_weight, array_pointcloud_points,
                                                                array_mesh$crs)$index
            }
        } else if (use_deterministic) {
            .mesh_directional_balance_index(pieces, weight = pieces$weight)$index
        } else {
            .random_point_directional_balance_index(poly_u, n_lines = n_lines, prep = list(poly = poly_u, tri = pieces),
                                                      seed = seed, weight = pieces$weight, points = shared_points)$index
        }
    }
    if ("span" %in% which) {
        out$span_index <- if (!is.null(array_mesh)) {
            if (use_deterministic_pairwise) {
                .mesh_span_index(array_pieces_lazy, n_quad = n_quad, weight = array_mesh$tri_weight)$index
            } else {
                .random_pair_span_index(poly_u, n_lines = n_lines,
                                         prep = list(poly = poly_u, tri = data.frame(area = array_mesh$tri_area)),
                                         seed = seed, weight = array_mesh$tri_weight,
                                         points = array_pairwise_points)$index
            }
        } else if (use_deterministic) {
            .mesh_span_index(pieces, n_quad = n_quad, weight = pieces$weight)$index
        } else {
            .random_pair_span_index(poly_u, n_lines = n_lines, prep = list(poly = poly_u, tri = pieces),
                                     seed = seed, weight = pieces$weight, points = shared_points)$index
        }
    }
    if ("radial_concentration" %in% which) {
        out$radial_concentration_index <- if (!is.null(array_mesh)) {
            if (use_deterministic) {
                .mesh_radial_concentration_index_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                                         array_mesh$tri_weight, array_mesh$crs)$index
            } else {
                .random_point_radial_index_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                                  array_mesh$tri_weight, array_pointcloud_points,
                                                  array_mesh$crs)$index
            }
        } else if (use_deterministic) {
            .mesh_radial_concentration_index(pieces, weight = pieces$weight)$index
        } else {
            .random_point_radial_index(poly_u, n_lines = n_lines, prep = list(poly = poly_u, tri = pieces),
                                        seed = seed, weight = pieces$weight, points = shared_points)$index
        }
    }
    if ("depth" %in% which) {
        out$depth_index <- if (!is.null(array_mesh)) {
            if (use_deterministic) {
                .mesh_depth_index_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                         array_mesh$tri_weight, poly_u)$index
            } else {
                .random_point_depth_index_array(array_mesh$P, array_mesh$T, array_mesh$tri_area,
                                                 array_mesh$tri_weight, array_pointcloud_points,
                                                 poly_u)$index
            }
        } else if (use_deterministic) {
            .mesh_depth_index(pieces, poly_u, weight = pieces$weight)$index
        } else {
            .random_point_depth_index(poly_u, n_lines = n_lines, prep = list(poly = poly_u, tri = pieces),
                                       seed = seed, weight = pieces$weight, points = shared_points)$index
        }
    }
    if ("hull_ratio" %in% which) out$hull_ratio_index <- hull_ratio_index(poly_u)$index
    if ("polsby_popper" %in% which) out$polsby_popper_index <- polsby_popper_index(poly_u)$index
    if ("width_length_ratio" %in% which) out$width_length_ratio_index <- width_length_ratio_index(poly_u)$index
    if ("reock" %in% which) out$reock_index <- reock_index(poly_u)$index
    if ("detour" %in% which) out$detour_index <- detour_index(poly_u)$index
    if ("exchange" %in% which) out$exchange_index <- exchange_index(poly_u)$index

    out$total_weight <- raw_total
    out$geometry <- poly_u
    do.call(st_sf, out)
}
