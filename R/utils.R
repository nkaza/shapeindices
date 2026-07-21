## -- small shared utilities ------------------------------------------------

# b if a is NULL, else a - avoids repeated is.null(x) defaulting
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Normalise a weight vector to sum to 1. Every weighted index in this
#' package is scale-invariant in weight, so this only fixes the reported
#' scale, not the index itself.
#' @param weight numeric vector, not yet normalised
#' @return `weight / sum(weight)`
#' @noRd
.normalize_weight <- function(weight) {
    if (!is.numeric(weight)) {
        extra <- if (is.factor(weight)) {
            paste(" Factors silently coerce to their integer level codes, not their labels -",
                  "convert explicitly first, e.g. as.numeric(as.character(weight)).")
        } else ""
        stop("weight must be numeric, got ", class(weight)[1], ".", extra)
    }
    if (anyNA(weight)) {
        stop("weight contains NA; every triangle/row needs a value.")
    }
    if (!all(is.finite(weight))) {
        stop("weight must be finite (no Inf/-Inf).")
    }
    # checked before the sum: e.g. c(-5, 10) sums positive but has no valid
    # interpretation as a weight
    if (any(weight < 0)) {
        stop("weight must be non-negative (a negative weight has no meaningful ",
             "interpretation as an area/importance weight, even if the total ",
             "still sums to something positive).")
    }
    s <- sum(weight)
    if (s <= 0) {
        stop("weight must sum to a positive total (all-zero weights aren't meaningful).")
    }
    weight / s
}

#' Auto-project geographic (lon/lat) input to a local equal-area CRS, since
#' this package's math assumes planar coordinates. Warns (not errors) on
#' missing CRS; errors on empty/degenerate geometry.
#' @param x sfc/sfg (multi)polygon, or a full sf data frame
#' @return `x`, projected if it was geographic
#' @noRd
.ensure_projected <- function(x) {
    geom <- st_geometry(x)
    if (length(geom) == 0 || anyNA(st_bbox(geom)) || all(st_is_empty(geom))) {
        stop("Geometry is empty or has NA/degenerate coordinates (bounding box: ",
             paste(signif(st_bbox(geom), 6), collapse = ", "),
             "); nothing meaningful to compute.")
    }
    crs <- st_crs(geom)
    if (is.na(crs)) {
        warning("Geometry has no CRS set; proceeding as if it's already planar. ",
                "If these are actually longitude/latitude coordinates, set the ",
                "CRS first (st_set_crs()) so this can auto-project correctly - ",
                "otherwise every length/area/convexity computation below will ",
                "be wrong.")
        return(x)
    }
    if (!isTRUE(st_is_longlat(crs))) {
        return(x)   # already planar
    }
    ctr <- suppressWarnings(st_coordinates(st_centroid(st_union(geom))))[1, 1:2]
    local_crs <- sprintf("+proj=laea +lat_0=%.10f +lon_0=%.10f +datum=WGS84 +units=m +no_defs",
                          ctr[2], ctr[1])
    message("Input is in geographic (lon/lat) coordinates; auto-projecting to a local ",
            "azimuthal-equal-area CRS centred on the data (lat_0 = ", round(ctr[2], 4),
            ", lon_0 = ", round(ctr[1], 4), ") before computing - pass already-projected ",
            "data instead if you need a specific CRS.")
    st_transform(x, crs = local_crs)
}

#' Warns if `simplify_tolerance` is being applied against a CRS whose
#' linear unit isn't metres - `st_simplify()`'s `dTolerance` is always in
#' the geometry's OWN CRS units, and `.ensure_projected()` only guarantees
#' metres when it performs its own auto-projection (geographic input);
#' already-projected input keeps whatever unit its CRS uses (many US State
#' Plane CRSs are US survey feet, not metres), silently changing what a
#' tolerance value like 50 actually means. No-op if `simplify_tolerance`
#' is NULL (nothing to warn about) or the unit genuinely is metres.
#' @param geom sfc/sfg, already the CRS simplification will actually run in
#' @param simplify_tolerance the value about to be passed as `dTolerance`,
#'   or NULL
#' @return invisible NULL
#' @noRd
.warn_if_nonmetric_tolerance <- function(geom, simplify_tolerance) {
    if (is.null(simplify_tolerance)) return(invisible(NULL))
    crs <- st_crs(geom)
    if (is.na(crs)) return(invisible(NULL))   # .ensure_projected() already warned about this
    unit <- crs$units_gdal
    if (is.null(unit) || tolower(unit) %in% c("metre", "meter")) return(invisible(NULL))
    warning(
        "simplify_tolerance is interpreted in this geometry's own CRS unit, which is \"",
        unit, "\", not metres - a value of ", simplify_tolerance, " means ", simplify_tolerance,
        " ", unit, "(s), not ", simplify_tolerance, " metres. This only happens for input ",
        "that was already in a projected CRS (this package only forces metres when it ",
        "auto-projects geographic lon/lat input itself); reproject to a metric CRS first ",
        "if you want simplify_tolerance to mean metres here."
    )
    invisible(NULL)
}

#' Exact mass centroid of a triangle mesh: mass-weighted average of each
#' triangle's own centroid. Exact, not an approximation - a uniform-density
#' triangle's centroid IS its own centre of mass, so a mass-weighted average
#' across triangles gives the true centroid of the whole piecewise-uniform-
#' density shape. Shared by .moment_of_inertia_core() (which then computes
#' Ixx/Iyy/Ixy about this point) and directional_balance_index() (which
#' needs only this, not the full second-moment tensor).
#' @param tri_geom sfc of triangle POLYGON geometries
#' @param mass numeric vector, each triangle's own mass (rho * area;
#'   physical area itself when unweighted)
#' @return numeric(2), c(x, y)
#' @noRd
.mass_centroid <- function(tri_geom, mass) {
    tri_centroids <- vapply(tri_geom, function(t) {
        colMeans(st_coordinates(t)[1:3, 1:2, drop = FALSE])
    }, numeric(2))
    c(sum(mass * tri_centroids[1, ]), sum(mass * tri_centroids[2, ])) / sum(mass)
}

#' Area-adaptive subdivision depth for a set of triangles: how many levels
#' of medial subdivision (each level quarters a triangle's area) each one
#' needs so that every triangle ends up with roughly the same FINAL
#' sub-triangle size, rather than every triangle getting the same fixed
#' depth regardless of its own size. `target_area` is pinned to the set's
#' own largest triangle at `max_depth`, so that triangle's result is
#' identical to what a fixed depth of `max_depth` would have given (no
#' regression there) - every smaller triangle needs proportionally less.
#' Shared by radial_concentration_index() (`.radial_point_cloud()`) and
#' the public `subdivide_mesh()`.
#' @param tri_area numeric vector, each triangle's own area
#' @param max_depth subdivision depth ceiling
#' @return integer vector, same length as tri_area, each in `[0, max_depth]`
#' @noRd
.adaptive_tri_depth <- function(tri_area, max_depth) {
    target_area <- max(tri_area) / 4^max_depth
    pmax(0L, pmin(max_depth, as.integer(ceiling(log(tri_area / target_area) / log(4)))))
}

#' Best-effort estimate of currently-available system memory, in MB - used
#' only to size a processing chunk, never to guarantee anything, so a wrong
#' answer just makes chunks more/less conservative rather than incorrect.
#' Reads /proc/meminfo's MemAvailable on Linux, vm_stat's free pages on
#' macOS; falls back to a fixed conservative constant everywhere else
#' (including if either parse fails) - safe because UNDER-estimating only
#' makes chunks smaller (slower, never unsafe), while OVER-estimating could
#' reintroduce whatever memory problem the caller is trying to avoid.
#' @return numeric, estimated available memory in MB
#' @noRd
.available_memory_mb <- function() {
    fallback_mb <- 2048
    tryCatch({
        sys <- Sys.info()[["sysname"]]
        if (identical(sys, "Linux") && file.exists("/proc/meminfo")) {
            meminfo <- readLines("/proc/meminfo", n = 5L, warn = FALSE)
            avail <- grep("^MemAvailable:", meminfo, value = TRUE)
            kb <- if (length(avail) == 1) {
                as.numeric(regmatches(avail, regexpr("[0-9]+", avail)))
            } else NA_real_
            if (is.finite(kb) && kb > 0) return(kb / 1024)
            fallback_mb
        } else if (identical(sys, "Darwin")) {
            vmstat <- system("vm_stat", intern = TRUE, ignore.stderr = TRUE)
            page_line <- grep("page size of", vmstat, value = TRUE)
            page_bytes <- if (length(page_line) == 1) {
                as.numeric(regmatches(page_line, regexpr("[0-9]+", page_line)))
            } else NA_real_
            if (!is.finite(page_bytes) || page_bytes <= 0) page_bytes <- 4096
            free_line <- grep("^Pages free:", vmstat, value = TRUE)
            pages <- if (length(free_line) == 1) {
                as.numeric(regmatches(free_line, regexpr("[0-9]+", free_line)))
            } else NA_real_
            if (is.finite(pages) && pages > 0) return(pages * page_bytes / 1024^2)
            fallback_mb
        } else {
            fallback_mb
        }
    }, error = function(e) fallback_mb, warning = function(w) fallback_mb)
}

#' How many candidate lines to test against all `E` boundary edges in one
#' vectorised pass (see .frac_outside_vectorised()'s own header for why this
#' needs bounding at all): sized so the batch's own intermediate arrays fit
#' within a fraction of currently-available memory, regardless of how large
#' `E` or the total candidate-line count gets.
#' @param E number of boundary edges being tested against
#' @param bytes_per_pair conservative bytes needed per (line, edge) pair
#'   across every intermediate vector built for one batch (not user-facing,
#'   a deliberately generous estimate covering the ~13 double + 2 int/
#'   logical vectors .frac_outside_vectorised() builds, plus R's own
#'   copy-on-modify overhead during arithmetic)
#' @param mem_fraction fraction of estimated available memory to budget for
#'   one batch - conservative, since other live objects (the mesh, the
#'   polygon itself) already occupy some of what's "available"
#' @return integer >= 1
#' @noRd
.choose_line_chunk_size <- function(E, bytes_per_pair = 200, mem_fraction = 0.2) {
    if (E <= 0) return(1L)
    budget_bytes <- .available_memory_mb() * 1024^2 * mem_fraction
    max(1L, as.integer(floor(budget_bytes / (E * bytes_per_pair))))
}

#' Warn before repairing an invalid/non-simple polygon (e.g. a
#' self-intersecting ring) - st_make_valid() can change the shape's
#' topology, not just clean up numerical noise.
#' @param poly sfc (MULTI)POLYGON, one or more features, already non-empty
#' @return `st_make_valid(poly)`
#' @noRd
.make_valid_warn <- function(poly) {
    invalid <- !st_is_valid(poly)
    invalid[is.na(invalid)] <- TRUE   # NA from st_is_valid() treated as invalid
    if (any(invalid)) {
        n <- sum(invalid)
        warning(
            if (length(poly) == 1) {
                "Geometry is not a valid, simple polygon (e.g. a self-intersecting ring). "
            } else {
                sprintf("%d of %d geometries are not valid, simple polygons (e.g. a self-intersecting ring). ",
                        n, length(poly))
            },
            "Repairing with sf::st_make_valid() before computing - this can change the ",
            "shape itself (a self-intersecting ring can split into separate polygons ",
            "touching at a point), not just clean up numerical noise, so the index below ",
            "may not describe the polygon you thought you passed in. Call ",
            "sf::st_make_valid() yourself first if you want to inspect the repaired ",
            "geometry before it's used."
        )
    }
    st_make_valid(poly)
}
