## =========================================================================
## 8. Classic compactness scores. No
## triangulation needed for any of these, just area/perimeter/hull/bbox
## ==========================================================================
##
## Six widely-used, purely boundary/hull-based compactness scores,
## unrelated to (and much cheaper than) the CDT-based indices elsewhere in
## this package. All are in (0, 1], 1 = a circle (or, for
## width_length_ratio_index, a square bounding rectangle), larger = more
## compact - the same "larger is more compact" convention as every other
## index here.
##
##   hull_ratio_index          = area(polygon) / area(convex hull)
##   polsby_popper_index       = 4*pi*area / perimeter^2
##   width_length_ratio_index  = min(mbr side1, mbr side2) / max(...), mbr =
##                               sf::st_minimum_rotated_rectangle(polygon) -
##                               the minimum-AREA bounding rectangle at any
##                               rotation (GEOS-backed), not the axis-aligned
##                               bounding box - see that function's own doc
##                               for why the axis-aligned version this used
##                               to be was a real problem, not a style choice
##   reock_index               = area(polygon) / area(minimum bounding circle)
##   detour_index              = perimeter(equal-area circle) / perimeter(convex hull)
##   exchange_index            = area(polygon INTERSECT equal-area circle at centroid) / area(polygon)
##
## WEIGHT: none of the six depend on how mass is distributed WITHIN the
## polygon - area, perimeter, convex hull and minimum bounding circle/
## rectangle are all properties of the boundary/extent alone, unchanged by
## any rearrangement of interior density. That's exactly why MOI/span/
## radial_concentration's "concentric rings" weighted reference doesn't
## carry over here: rearranging density never changes area, perimeter or
## hull, so there's nothing for a rearrangement argument to act on.
##
## Substituting a weighted quantity for "area" while leaving perimeter^2
## or MBC-area untouched was considered and rejected for polsby_popper_index()/
## reock_index() specifically: their (0, 1] range comes from both sides of the
## ratio being tied together by the SAME theorem (isoperimetric inequality
## for Polsby-Popper; area(polygon) <= area(MBC) trivially, since the MBC
## contains the whole polygon, for Reock). A weight has no natural
## relationship to perimeter or MBC at all, so nothing keeps the ratio
## bounded once one side is swapped out - it becomes arbitrarily sensitive
## to whatever scale `weight` happens to be on (population counts vs.
## square metres, say), not a meaningful compactness score any more. Same
## reasoning rules out a weighted detour_index(): the hull's perimeter has no
## natural relationship to a per-triangle weight either. exchange_index() has no
## natural weighted form either, for a related but distinct reason: it's
## an area-overlap fraction, and there's nothing for a per-triangle weight
## to attach to beyond which triangle falls inside the reference circle.
##
## All six are therefore left genuinely unweighted, no `weight` argument
## at all - the only way `weights` can move any of them, via
## shape_indices_sf(byrow = FALSE), is the same indirect route as any
## other index: a weight of exactly 0/NA excludes that row as a hole
## BEFORE any index runs, which changes poly_u itself. The magnitude of a
## nonzero weight has zero effect on any of the six - see
## .shape_indices_sf_grouped()'s own warning about this.
##
## detour_index() = ratio of the equal-area circle's perimeter to the polygon's
## own convex hull's perimeter (Angel, Parent & Civco 2010's Detour Index -
## a measure of how hard a shape is to circumnavigate as an obstacle).
## Provably in (0, 1]: hull_area >= poly_area (hull contains the polygon),
## so by the isoperimetric inequality applied to the hull itself,
## perimeter(hull) >= circle_perimeter(hull_area) >= circle_perimeter(poly_area)
## - the same chaining-inequality pattern used elsewhere in this package
## (e.g. radial_concentration_index() bounding its own centroid-anchored
## cousin). Robust to holes and multi-part dispersal by construction: the
## convex hull is taken over the union of all parts, exactly like
## hull_ratio_index()/reock_index() already do, and a hole never reaches the hull at
## all (same known blind spot hull_ratio_index() already has, not a new one).
##
## exchange_index() = share of the polygon's own area that falls inside the
## equal-area circle centred at its centroid (Angel, Parent & Civco 2010's
## Exchange Index). Bounded in [0, 1] trivially (an intersection's area
## can't exceed either operand's). KNOWN PATHOLOGY, kept anyway: for a
## multi-part polygon, the circle is centred at the OVERALL centroid,
## which for well-separated parts sits in the empty gap between them - once
## the circle's radius is smaller than the distance to the nearest part,
## the intersection is exactly zero, not just low. Verified empirically:
## two unit squares a couple of units apart already score under 0.05, and
## the index hits exactly 0 once they're far enough apart - Angel et al.
## themselves note the same failure mode for real districts split by water
## ("does not account for separation by bodies of water"). Included
## despite this, matching width_length_ratio_index()'s own precedent of shipping
## a real, useful, but documented-blind-spot metric rather than omitting
## it - see width_length_ratio_index()'s own "KNOWN LIMITATION" section.
##
## reock_index needs a minimum enclosing circle (MEC), computed by a
## deterministic (no shuffling, so no RNG entanglement with the caller's
## own seed) iterative Welzl-style algorithm on the polygon's own convex
## hull vertices - the MEC of a point set depends only on its hull, so
## this is exact, not an approximation, despite only using O(n^3) worst
## case (fine at the vertex counts a convex hull actually has).

## -- minimum enclosing circle, for reock_index() ---------------------------------

#' @noRd
.circle_from_2 <- function(p1, p2) {
    list(center = (p1 + p2) / 2, r = sqrt(sum((p1 - p2)^2)) / 2)
}

#' Circumcircle of 3 points, falling back to the widest pairwise circle if
#' they're (near-)collinear.
#' @noRd
.circle_from_3 <- function(p1, p2, p3) {
    ax <- p1[1]; ay <- p1[2]; bx <- p2[1]; by <- p2[2]; cx <- p3[1]; cy <- p3[2]
    d <- 2 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by))
    if (abs(d) < 1e-9) {
        cands <- list(.circle_from_2(p1, p2), .circle_from_2(p2, p3), .circle_from_2(p1, p3))
        return(cands[[which.max(vapply(cands, `[[`, numeric(1), "r"))]])
    }
    ux <- ((ax^2 + ay^2) * (by - cy) + (bx^2 + by^2) * (cy - ay) + (cx^2 + cy^2) * (ay - by)) / d
    uy <- ((ax^2 + ay^2) * (cx - bx) + (bx^2 + by^2) * (ax - cx) + (cx^2 + cy^2) * (bx - ax)) / d
    center <- c(ux, uy)
    list(center = center, r = sqrt(sum((p1 - center)^2)))
}

#' @noRd
.point_in_circle <- function(p, circle, eps = 1e-7) {
    sqrt(sum((p - circle$center)^2)) <= circle$r + eps
}

#' Minimum enclosing circle of a point set - the classic iterative
#' Welzl-style incremental algorithm, without the usual random shuffle (so
#' it never touches the caller's RNG state); worst case O(n^3), fine for
#' the vertex counts a convex hull actually has.
#' @param pts an Nx2 coordinate matrix
#' @return list(center, r)
#' @noRd
.min_enclosing_circle <- function(pts) {
    n <- nrow(pts)
    if (n == 1) return(list(center = pts[1, ], r = 0))

    circle <- .circle_from_2(pts[1, ], pts[2, ])
    for (i in seq_len(n)) {
        if (.point_in_circle(pts[i, ], circle)) next
        circle <- list(center = pts[i, ], r = 0)
        for (j in seq_len(i - 1)) {
            if (.point_in_circle(pts[j, ], circle)) next
            circle <- .circle_from_2(pts[i, ], pts[j, ])
            for (k in seq_len(j - 1)) {
                if (.point_in_circle(pts[k, ], circle)) next
                circle <- .circle_from_3(pts[i, ], pts[j, ], pts[k, ])
            }
        }
    }
    circle
}

## -- the six scores -----------------------------------------------------

#' Convex hull area ratio of a (multi)polygon
#'
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(index, area, hull_area, hull). index is in (0, 1], 1 = the
#'   polygon is convex (equals its own convex hull).
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' hull_ratio_index(nc[nc$NAME == "Wake", ])$index
#' hull_ratio_index(nc[nc$NAME == "Dare", ])$index
#' @export
hull_ratio_index <- function(poly) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    area <- sum(as.numeric(st_area(poly)))
    hull <- st_convex_hull(st_union(poly))
    hull_area <- as.numeric(st_area(hull))
    index <- if (hull_area > 0) area / hull_area else NA_real_
    list(index = index, area = area, hull_area = hull_area, hull = hull)
}

#' Polsby-Popper compactness score of a (multi)polygon
#'
#' 4*pi*area/perimeter^2, in (0, 1], 1 = a circle. Perimeter is the full
#' boundary length (outer ring plus any holes).
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(index, area, perimeter)
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' polsby_popper_index(nc[nc$NAME == "Wake", ])$index
#' polsby_popper_index(nc[nc$NAME == "Dare", ])$index
#' @export
polsby_popper_index <- function(poly) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    area <- sum(as.numeric(st_area(poly)))
    perimeter <- sum(as.numeric(st_length(st_boundary(poly))))
    index <- if (perimeter > 0) (4 * pi * area) / perimeter^2 else NA_real_
    list(index = index, area = area, perimeter = perimeter)
}

#' Width-length ratio of a (multi)polygon's minimum-area bounding rectangle
#'
#' The shorter of the bounding rectangle's two side lengths over the
#' longer (`length` is the longer one, matching everyday usage, not the
#' shorter one), in (0, 1], 1 = a square bounding rectangle. Uses the
#' MINIMUM-AREA rectangle at any rotation (`sf::st_minimum_rotated_rectangle()`,
#' GEOS-backed - no new package dependency), not the axis-aligned bounding
#' box: an earlier version of this function used the axis-aligned box,
#' which meant the same shape could score anywhere from its true ratio up
#' to a spurious 1 depending purely on which way it happened to be drawn
#' relative to the coordinate axes - orientation is not a property of the
#' shape itself, so an index built on it shouldn't depend on it either.
#' Verified directly: a fixed 2:1 rectangle rotated from 0 to 90 degrees
#' now returns the same 0.5 throughout, rather than swinging up to 1.0
#' at 45 degrees the way the axis-aligned version did.
#'
#' KNOWN LIMITATION: this looks only at the bounding rectangle, so it's
#' blind to both holes and multi-part dispersal in ways
#' `hull_ratio_index`/`reock_index` are not. A solid square and the same
#' square with a large hole punched through it score identically (a hole
#' never extends past the outer boundary, so it never moves the bounding
#' rectangle). Likewise, several small disconnected pieces sitting near
#' the corners of a roughly square overall extent can score close to 1
#' despite being scattered, not compact - the bounding rectangle is taken
#' over every part combined, independent of how much of it is actually
#' filled.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(index, length, width)
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' width_length_ratio_index(nc[nc$NAME == "Wake", ])$index
#' width_length_ratio_index(nc[nc$NAME == "Dare", ])$index
#' @export
width_length_ratio_index <- function(poly) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    mbr <- st_minimum_rotated_rectangle(poly)
    xy <- st_coordinates(mbr)[, 1:2, drop = FALSE]
    side1 <- sqrt(sum((xy[2, ] - xy[1, ])^2))
    side2 <- sqrt(sum((xy[3, ] - xy[2, ])^2))
    width  <- min(side1, side2)
    length <- max(side1, side2)
    index  <- if (length > 0) width / length else NA_real_
    list(index = index, length = length, width = width)
}

#' Reock compactness score of a (multi)polygon
#'
#' area/area(minimum bounding circle), in (0, 1], 1 = a circle.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(index, area, mbc_area, mbc). `mbc` is the minimum bounding
#'   circle itself, as an sfc POLYGON, for plotting.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' reock_index(nc[nc$NAME == "Wake", ])$index
#' reock_index(nc[nc$NAME == "Dare", ])$index
#' @export
reock_index <- function(poly) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    area <- sum(as.numeric(st_area(poly)))
    hull <- st_convex_hull(st_union(poly))
    pts  <- unique(st_coordinates(hull)[, 1:2, drop = FALSE])
    mec  <- .min_enclosing_circle(pts)
    mbc_area <- pi * mec$r^2
    index <- if (mbc_area > 0) area / mbc_area else NA_real_
    mbc <- st_buffer(st_sfc(st_point(mec$center), crs = st_crs(poly)), dist = mec$r, nQuadSegs = 90)
    list(index = index, area = area, mbc_area = mbc_area, mbc = mbc)
}

#' Detour compactness score of a (multi)polygon
#'
#' The ratio of the equal-area circle's perimeter to the perimeter of the
#' polygon's own convex hull, in (0, 1], 1 = a circle (whose hull is
#' itself). Angel, Parent & Civco (2010) introduce this to measure how
#' hard a shape is to circumnavigate as an obstacle - a smooth, rounded
#' hull scores high regardless of how convoluted the actual boundary
#' inside that hull is, since only the hull's own perimeter enters the
#' ratio.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(index, area, hull_perimeter, hull). `hull` is the convex
#'   hull itself, as an sfc POLYGON, for plotting.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' detour_index(nc[nc$NAME == "Wake", ])$index
#' detour_index(nc[nc$NAME == "Dare", ])$index
#' @export
detour_index <- function(poly) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    area <- sum(as.numeric(st_area(poly)))
    hull <- st_convex_hull(st_union(poly))
    hull_perimeter <- as.numeric(st_length(st_boundary(hull)))
    circle_perimeter <- 2 * sqrt(pi * area)
    index <- if (hull_perimeter > 0) circle_perimeter / hull_perimeter else NA_real_
    list(index = index, area = area, hull_perimeter = hull_perimeter, hull = hull)
}

#' Exchange compactness score of a (multi)polygon
#'
#' The share of the polygon's own area that falls inside the equal-area
#' circle centred at its centroid, in `[0, 1]`, 1 = a circle. Angel,
#' Parent & Civco (2010) introduce this as a natural metric for
#' gerrymandering: a district that "reaches out" to grab distant voters,
#' or excludes nearby ones, has most of its area outside that circle.
#'
#' KNOWN LIMITATION: for a multi-part polygon, the reference circle is
#' centred at the *overall* centroid, which for well-separated parts sits
#' in the empty space between them. Once the circle's radius is smaller
#' than the distance to the nearest part, the intersection - and so the
#' index - is exactly 0, not just low; two individually-compact pieces a
#' few units apart can already score under 0.05. Angel et al. themselves
#' note the same failure mode on real districts split by water. Included
#' despite this, the same way `width_length_ratio_index()` ships a real,
#' useful, but documented-blind-spot metric rather than omitting it.
#' @param poly a single sfg/sfc (MULTI)POLYGON
#' @return list(index, area, circle_area, circle). `circle` is the
#'   equal-area reference circle itself, as an sfc POLYGON, for plotting.
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' exchange_index(nc[nc$NAME == "Wake", ])$index
#' exchange_index(nc[nc$NAME == "Dare", ])$index
#' @export
exchange_index <- function(poly) {
    poly <- st_geometry(poly)
    poly <- .ensure_projected(poly)
    poly <- .make_valid_warn(poly)
    poly_u <- st_union(poly)
    area <- sum(as.numeric(st_area(poly)))
    cen  <- st_centroid(poly_u)
    r    <- sqrt(area / pi)
    circle <- st_buffer(cen, dist = r, nQuadSegs = 90)
    inter  <- suppressWarnings(st_intersection(poly_u, circle))
    inter_area <- if (length(inter) == 0) 0 else sum(as.numeric(st_area(inter)))
    index <- if (area > 0) inter_area / area else NA_real_
    list(index = index, area = area, circle_area = as.numeric(st_area(circle)), circle = circle)
}

## -- vectorised _sf() wrappers ----------------------------------------

#' hull_ratio_index() for every row of an sf data frame
#'
#' Each row indexed independently (see shape_indices_sf() for a version
#' that also computes the other indices).
#' @param x an sf data frame
#' @return `x` with one new column, hull_ratio_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- hull_ratio_index_sf(nc[1:5, ])
#' res$hull_ratio_index
#' @export
hull_ratio_index_sf <- function(x) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) hull_ratio_index(geoms[i]))
    x$hull_ratio_index <- vapply(res, function(r) r$index, numeric(1))
    x
}

#' polsby_popper_index() for every row of an sf data frame
#'
#' @param x an sf data frame
#' @return `x` with one new column, polsby_popper_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- polsby_popper_index_sf(nc[1:5, ])
#' res$polsby_popper_index
#' @export
polsby_popper_index_sf <- function(x) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) polsby_popper_index(geoms[i]))
    x$polsby_popper_index <- vapply(res, function(r) r$index, numeric(1))
    x
}

#' width_length_ratio_index() for every row of an sf data frame
#'
#' @param x an sf data frame
#' @return `x` with one new column, width_length_ratio_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- width_length_ratio_index_sf(nc[1:5, ])
#' res$width_length_ratio_index
#' @export
width_length_ratio_index_sf <- function(x) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) width_length_ratio_index(geoms[i]))
    x$width_length_ratio_index <- vapply(res, function(r) r$index, numeric(1))
    x
}

#' reock_index() for every row of an sf data frame
#'
#' @param x an sf data frame
#' @return `x` with one new column, reock_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- reock_index_sf(nc[1:5, ])
#' res$reock_index
#' @export
reock_index_sf <- function(x) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) reock_index(geoms[i]))
    x$reock_index <- vapply(res, function(r) r$index, numeric(1))
    x
}

#' detour_index() for every row of an sf data frame
#'
#' @param x an sf data frame
#' @return `x` with one new column, detour_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- detour_index_sf(nc[1:5, ])
#' res$detour_index
#' @export
detour_index_sf <- function(x) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) detour_index(geoms[i]))
    x$detour_index <- vapply(res, function(r) r$index, numeric(1))
    x
}

#' exchange_index() for every row of an sf data frame
#'
#' @param x an sf data frame
#' @return `x` with one new column, exchange_index
#' @examples
#' nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
#' res <- exchange_index_sf(nc[1:5, ])
#' res$exchange_index
#' @export
exchange_index_sf <- function(x) {
    geoms <- st_geometry(x)
    res <- lapply(seq_along(geoms), function(i) exchange_index(geoms[i]))
    x$exchange_index <- vapply(res, function(r) r$index, numeric(1))
    x
}
