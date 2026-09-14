# Drawing a hypergraph: every hyperedge is a smooth pebble around its members
# -- the rounded convex hull with its corners low-pass filtered away -- in a
# translucent fill parted from its neighbours by a white seam.
# The default layout places the star expansion (the bipartite graph of nodes
# and hyperedges), which gives every hyperedge its own position among its
# members: blobs fan out instead of piling onto shared nodes, and each
# hyperedge label has a natural home. The force-directed layout is computed
# here (Fruchterman-Reingold), the ring layout by cograph; ggplot2 draws.

# Layout of the star expansion, or of the clique projection, as two tables in
# one shared [0, 1] frame: `nodes` (node, x, y) and `edges` (hyperedge, x, y).
# Under a projection layout a hyperedge has no vertex of its own and sits at
# the centroid of its members. Both axes are scaled by the same factor so the
# layout keeps its shape.
.thg_positions <- function(hg, layout, seed, center = NULL) {
  inc <- as.matrix(hg$incidence != 0)
  n <- nrow(inc)
  m <- ncol(inc)
  star <- identical(layout, "bipartite")
  centred <- hg$nodes %in% center
  if (is.data.frame(layout)) {
    if (!all(c("node", "x", "y") %in% names(layout))) {
      .thg_bad_input("a `layout` data.frame needs `node`, `x` and `y` columns")
    }
    missing <- setdiff(hg$nodes, as.character(layout$node))
    if (length(missing)) {
      .thg_bad_input(sprintf("`layout` lacks coordinates for %d node(s)",
                             length(missing)))
    }
    idx <- match(hg$nodes, as.character(layout$node))
    xy <- cbind(as.numeric(layout$x[idx]), as.numeric(layout$y[idx]))
  } else if (star) {
    # vertex ids, not names: a node and a hyperedge may share a label
    adjacency <- matrix(0, n + m, n + m,
                        dimnames = rep(list(paste0("v", seq_len(n + m))), 2L))
    adjacency[seq_len(n), n + seq_len(m)] <- inc
    adjacency <- adjacency + t(adjacency)
    xy <- .thg_layout_fr(adjacency, seed = seed,
                         pinned = c(centred, rep(FALSE, m)))
  } else {
    layout <- match.arg(layout, c("spring", "circle"))
    projection <- as.matrix(hg_project(hg, method = "clique",
                                       weighted = FALSE, what = "matrix"))
    xy <- if (identical(layout, "spring")) {
      .thg_layout_fr((projection != 0) * 1, seed = seed, pinned = centred)
    } else if (any(centred)) {
      # the centred nodes on a small inner ring, the rest on the outer ring
      .thg_ring(centred, radius = 0.25) + .thg_ring(!centred)
    } else {
      pos <- cograph::layout_circle(cograph::as_cograph(projection,
                                                        directed = FALSE))
      cbind(as.numeric(pos$x), as.numeric(pos$y))
    }
  }
  if (!star) {
    members <- sweep(inc, 2L, pmax(colSums(inc), 1), "/")
    xy <- rbind(xy, t(members) %*% xy)
  }
  lower <- apply(xy, 2L, min)
  span <- max(apply(xy, 2L, function(z) diff(range(z))))
  if (span < sqrt(.Machine$double.eps)) span <- 1
  xy <- sweep(xy, 2L, lower) / span
  list(
    nodes = data.frame(node = hg$nodes, x = xy[seq_len(n), 1L],
                       y = xy[seq_len(n), 2L], stringsAsFactors = FALSE),
    edges = data.frame(hyperedge = colnames(hg$incidence) %||%
                         paste0("h", seq_len(m)),
                       x = xy[n + seq_len(m), 1L], y = xy[n + seq_len(m), 2L],
                       stringsAsFactors = FALSE)
  )
}

# Fruchterman-Reingold placement of an undirected 0/1 adjacency matrix. Every
# pair repels with force k^2 / d and every edge attracts with d^2 / k (k = 1),
# each vertex moving at most the current temperature, which cools linearly.
# One sweep is vectorised over all pairs; the sweeps run through Reduce().
# Chosen over cograph::layout_spring() because it keeps non-members out of
# the pebbles: on three test hypergraphs, 10 seeds each, it averaged 0.3, 0.3
# and 0 non-member nodes inside a pebble against 5, 1.8 and 10. The seed is
# scoped to this call; the caller's random stream is restored. `pinned`
# vertices are held on a small ring at the origin and never move; the rest of
# the layout is then wrapped around them by .thg_surround().
.thg_layout_fr <- function(adjacency, seed = 1L, n_iter = 500L,
                           pinned = rep(FALSE, nrow(adjacency))) {
  n <- nrow(adjacency)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) {
    assign(".Random.seed", old_seed, envir = globalenv())
  } else {
    rm(".Random.seed", envir = globalenv())
  }, add = TRUE)
  set.seed(seed)
  side <- sqrt(n)
  start <- list(x = stats::runif(n, -side / 2, side / 2),
                y = stats::runif(n, -side / 2, side / 2))
  hold <- .thg_ring(pinned, radius = side / 4)
  start$x[pinned] <- hold[pinned, 1L]
  start$y[pinned] <- hold[pinned, 2L]
  temperatures <- seq(side / 10, side / 1000, length.out = n_iter)
  sweep_once <- function(state, temperature) {
    dx <- outer(state$x, state$x, "-")
    dy <- outer(state$y, state$y, "-")
    d <- pmax(sqrt(dx^2 + dy^2), 0.01)
    # both forces divided by d, so multiplying by (dx, dy) gives the vector
    force <- 1 / d^2 - adjacency * d
    diag(force) <- 0
    move_x <- rowSums(dx * force)
    move_y <- rowSums(dy * force)
    magnitude <- pmax(sqrt(move_x^2 + move_y^2), 1e-12)
    capped <- pmin(magnitude, temperature) * !pinned
    list(x = state$x + move_x / magnitude * capped,
         y = state$y + move_y / magnitude * capped)
  }
  placed <- Reduce(sweep_once, temperatures, start)
  xy <- cbind(placed$x, placed$y)
  if (any(pinned)) xy <- .thg_surround(xy, pinned)
  xy
}

# Wrap the free vertices around the pinned ring. Pinning alone holds the ring
# still but lets the rest settle on one side, which leaves the "centre" at the
# rim; gravity toward the ring did not fix that (2 of 10 layouts at strength
# 1). Instead each free vertex keeps its distance from the ring and its
# angular order, the angles are spread evenly around the full circle, and no
# free vertex sits closer than 1.5 times the ring's radius. The ring radius
# (a quarter of the layout side) and that floor gave the fewest non-members
# inside pebbles on two step-resolution hypergraphs; centring still costs
# some, because convex pebbles around the steps must cover the middle.
.thg_surround <- function(xy, pinned) {
  free <- !pinned
  angle <- atan2(xy[free, 2L], xy[free, 1L])
  ring_radius <- max(sqrt(rowSums(xy[pinned, , drop = FALSE]^2)),
                     sqrt(.Machine$double.eps))
  radius <- pmax(sqrt(xy[free, 1L]^2 + xy[free, 2L]^2), 1.5 * ring_radius)
  even <- min(angle) +
    2 * pi * (rank(angle, ties.method = "first") - 1) / sum(free)
  xy[free, 1L] <- radius * cos(even)
  xy[free, 2L] <- radius * sin(even)
  xy
}

# Points evenly spaced on a circle for the flagged entries of `flag`, zero
# elsewhere: an n x 2 matrix, so rings for disjoint flags add up.
.thg_ring <- function(flag, radius = 1) {
  k <- sum(flag)
  angle <- 2 * pi * (seq_len(k) - 1) / max(k, 1L) + pi / 2
  out <- matrix(0, length(flag), 2L)
  if (k == 1L) return(out)
  out[flag, 1L] <- radius * cos(angle)
  out[flag, 2L] <- radius * sin(angle)
  out
}

# Node coordinates in [0, 1]^2 from a layout name or a user table.
.thg_layout <- function(hg, layout, seed) {
  .thg_positions(hg, layout, seed)$nodes
}

# A rounded hull: the convex hull of the member points each widened by a
# circle of `radius` (the Minkowski sum of the hull and a disc), traced as a
# closed polygon. Two members give a stadium, one member a disc.
.thg_hull <- function(x, y, radius, n_arc = 36L) {
  theta <- seq(0, 2 * pi, length.out = n_arc + 1L)[-1L]
  px <- as.vector(outer(x, radius * cos(theta), "+"))
  py <- as.vector(outer(y, radius * sin(theta), "+"))
  ring <- grDevices::chull(px, py)
  data.frame(x = px[ring], y = py[ring])
}

# A pebble: the rounded hull with its corners smoothed away. The outline is
# resampled evenly by arc length and read as one complex signal x + iy; each
# harmonic f is damped by exp(-(f / detail)^2 / 2). A Gaussian taper, unlike a
# hard cut-off, adds no ripple. Smoothing a convex curve keeps it convex but
# pulls its corners in, so the curve is then pushed out along its normal just
# far enough that every member keeps `clearance * radius` of room; for a convex
# curve that offset stays smooth and raises every inside distance by exactly
# the offset. `detail = Inf` returns the rounded hull unsmoothed.
.thg_pebble <- function(x, y, radius, detail, clearance = 0.6,
                        n_points = 256L) {
  ring <- .thg_hull(x, y, radius, n_arc = 48L)
  if (!is.finite(detail)) return(ring)
  closed_x <- c(ring$x, ring$x[1L])
  closed_y <- c(ring$y, ring$y[1L])
  arc <- c(0, cumsum(sqrt(diff(closed_x)^2 + diff(closed_y)^2)))
  at <- seq(0, arc[length(arc)], length.out = n_points + 1L)[-(n_points + 1L)]
  signal <- complex(real = stats::approx(arc, closed_x, at)$y,
                    imaginary = stats::approx(arc, closed_y, at)$y)
  half <- n_points %/% 2L
  freq <- c(0:half, -((half - 1L):1L))
  spectrum <- stats::fft(signal) * exp(-(freq / detail)^2 / 2)
  curve <- stats::fft(spectrum, inverse = TRUE) / n_points
  tangent <- stats::fft(spectrum * 1i * freq, inverse = TRUE) / n_points
  normal <- -1i * tangent / Mod(tangent)
  # Re(Conj(a) * b) is the dot product: orient the normals outward
  if (sum(Re(Conj(normal) * (curve - mean(curve)))) < 0) normal <- -normal
  members <- complex(real = x, imaginary = y)
  # sampled offsets are only nearly exact, so re-measure; two passes suffice
  # in practice and the bound keeps the loop finite
  for (pass in seq_len(3L)) {
    deficit <- clearance * radius - .thg_inside_distance(curve, normal, members)
    if (max(deficit) <= 0) break
    curve <- curve + max(deficit) * normal
  }
  data.frame(x = Re(curve), y = Im(curve))
}

# Signed distance from each point to a closed convex outline given as complex
# vertices with outward normals: positive inside, negative outside. The
# distance is to the outline's edges, not its vertices, which would overstate
# it.
.thg_inside_distance <- function(curve, normal, points) {
  start <- curve
  edge <- c(curve[-1L], curve[1L]) - start
  vapply(points, function(p) {
    along <- pmin(1, pmax(0, Re(Conj(edge) * (p - start)) / Mod(edge)^2))
    gap <- min(Mod(p - (start + along * edge)))
    if (all(Re(Conj(normal) * (p - curve)) <= 0)) gap else -gap
  }, numeric(1L))
}

# Is the point (px, py) inside the closed convex outline (x, y)? True when it
# lies on the same side of every edge, whichever way the outline winds.
.thg_in_convex <- function(x, y, px, py) {
  ex <- c(x[-1L], x[1L]) - x
  ey <- c(y[-1L], y[1L]) - y
  side <- ex * (py - y) - ey * (px - x)
  all(side >= 0) || all(side <= 0)
}

# One value per hyperedge from a selector: NULL, "size", an edge_data
# column, a named vector keyed by hyperedge, or a vector of length m.
.thg_edge_aesthetic <- function(hg, by, arg) {
  if (is.null(by)) return(NULL)
  edge_names <- colnames(hg$incidence)
  m <- hg$n_hyperedges
  if (is.character(by) && length(by) == 1L) {
    if (identical(by, "size")) return(lengths(hg$hyperedges))
    if (!is.null(hg$edge_data) && by %in% names(hg$edge_data)) {
      values <- hg$edge_data[[by]]
      return(values[match(edge_names, as.character(hg$edge_data$edge))])
    }
    if (m != 1L) {
      .thg_bad_input(sprintf(
        "`%s` must be \"size\", a column of the edge metadata, or one value per hyperedge",
        arg
      ))
    }
  }
  if (!is.null(names(by))) {
    missing <- setdiff(edge_names, names(by))
    if (length(missing)) {
      .thg_bad_input(sprintf("`%s` has no value for %d hyperedge(s)", arg,
                             length(missing)))
    }
    return(unname(by[edge_names]))
  }
  if (length(by) != m) {
    .thg_bad_input(sprintf("`%s` must have one value per hyperedge (%d)", arg, m))
  }
  by
}

# Colours for a fill selector: a named palette over the levels of a
# discrete variable, or a ramp over the range of a numeric one.
.thg_fill_colours <- function(values) {
  if (is.numeric(values)) {
    rng <- range(values)
    unit <- if (diff(rng) > 0) (values - rng[1L]) / diff(rng) else rep(0.5, length(values))
    ramp <- grDevices::colorRamp(.thg_okabe_ito_ramp())
    rgb <- ramp(unit)
    return(list(colours = grDevices::rgb(rgb[, 1L], rgb[, 2L], rgb[, 3L],
                                         maxColorValue = 255),
                palette = NULL))
  }
  levels_fill <- if (is.factor(values)) levels(values) else
    sort(unique(as.character(values)))
  palette <- if (length(levels_fill) <= length(.thg_okabe_ito)) {
    .thg_okabe_ito[seq_along(levels_fill)]
  } else {
    grDevices::colorRampPalette(.thg_okabe_ito_ramp())(length(levels_fill))
  }
  names(palette) <- levels_fill
  list(colours = unname(palette[as.character(values)]), palette = palette)
}

.thg_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

#' Plot a hypergraph with hyperedges as pebbles
#'
#' Draws every hyperedge as a smooth pebble around its members -- the convex
#' hull of the members, widened by `padding` and with its corners smoothed
#' away -- in a translucent fill parted from its neighbours by a white seam.
#' Every member is guaranteed to lie inside its pebble. The default layout places the bipartite star expansion
#' (nodes and hyperedges together), which gives each hyperedge a position of
#' its own among its members; hulls fan out instead of piling onto shared
#' nodes, and hyperedge labels sit at those positions. Hulls can be coloured
#' by hyperedge size, by a column of the edge metadata, or by any value
#' supplied per hyperedge, and outlined with a line type by the same kinds of
#' selector. A hyperedge with a single member is drawn as its node alone.
#'
#' @param x A `net_hypergraph` with at least one hyperedge of two or more
#'   members. Draw a part of a large hypergraph by passing [hg_subset()]
#'   first.
#' @param layout `"bipartite"` (default; a Fruchterman-Reingold layout of the
#'   star expansion, nodes and hyperedges placed together), `"spring"` (the
#'   same on the clique projection), `"circle"`, or a data.frame with `node`,
#'   `x` and `y` columns to reuse node coordinates across panels (hyperedges
#'   then sit at the centroid of their members).
#' @param center Node names to place in the middle of the picture, such as
#'   the outcome states every hyperedge ends in. Under `"bipartite"` and
#'   `"spring"` they are held on a small ring at the centre while the rest of
#'   the layout arranges around them; under `"circle"` they sit inside the
#'   ring of the other nodes. Names that are not nodes of this hypergraph are
#'   skipped, so one vector can serve several subsets; an error is raised
#'   only when none of them is. `NULL` (default) places every node freely.
#'   Not used with a `layout` table.
#' @param seed Seed for the force-directed layouts (default `1`); the
#'   caller's random number stream is left untouched.
#' @param color_by Colour of the hulls: `NULL` (one colour), `"size"`
#'   (hyperedge cardinality, sequential scale), the name of a column in the
#'   edge metadata (`x$edge_data`), a vector named by hyperedge, or a vector
#'   with one value per hyperedge. Character or factor values get the
#'   Okabe-Ito palette; numeric values a sequential scale built from it.
#' @param linetype_by Outline line type of the hulls, resolved like
#'   `color_by`; discrete only.
#' @param dismantled Draw one panel per hyperedge instead of one figure with
#'   every hull overlaid (default `FALSE`). All panels share the layout, so
#'   positions are comparable; each panel greys out the non-members and is
#'   titled with its hyperedge name. `edge_labels` does not apply.
#' @param ncol Columns in the panel grid when `dismantled = TRUE`
#'   (default: roughly square).
#' @param edge_labels Name the hyperedges on the figure. `FALSE` (default)
#'   writes nothing, `TRUE` writes the hyperedge names, or pass a vector named
#'   by hyperedge to write something else. Under the bipartite layout each
#'   label sits at the hyperedge's own position; under the projection layouts,
#'   just outside the member furthest from the centre, clear of the crowded
#'   overlap.
#' @param edge_label_size Text size for `edge_labels` (default `3`).
#' @param labels `TRUE` (default) writes the node names, `FALSE` writes
#'   none, and a character vector named by node replaces the names shown.
#' @param label_size,node_size Text and point sizes.
#' @param detail How much of the hull's shape the smoothing keeps: the width,
#'   in harmonics, of the Gaussian low-pass applied to the outline. The
#'   default `5` gives rounded, tapering petals; `3` is rounder still, `8`
#'   stays close to the members, and `Inf` draws the unsmoothed rounded hull.
#' @param outline Outline colour: `"white"` (default) draws seams that part
#'   overlapping hyperedges; `"fill"` outlines each hyperedge in its own
#'   colour (pair it with a lower `alpha`, e.g. `0.15`); any other colour is
#'   used as given.
#' @param alpha Fill transparency (default `0.45`); the outline is always
#'   drawn at full strength.
#' @param linewidth Width of the outlines (default `1.1`).
#' @param padding Room around the member positions as a fraction of the layout
#'   extent (default `0.045`).
#' @param legend_title Legend title for `color_by`; defaults to the
#'   selector's name.
#' @param ... Unused; for S3 consistency.
#' @return A ggplot object (with `dismantled = TRUE`, one facet per
#'   hyperedge).
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#'
#'   Fruchterman, T. M. J., & Reingold, E. M. (1991). Graph drawing by
#'   force-directed placement. *Software: Practice and Experience*, 21(11),
#'   1129-1164. \doi{10.1002/spe.4380211102}
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "b", "c", "d", "d", "e", "f"),
#'   event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e4")
#' )
#' hg <- group_hypergraph(dat, "member", "event")
#' plot(hg, color_by = "size", edge_labels = TRUE)
#' plot(hg, detail = Inf, outline = "fill", alpha = 0.15)
#' plot(hg, center = c("b", "c"))
#' plot(hg, dismantled = TRUE)
#' @export
plot.net_hypergraph <- function(x, layout = c("bipartite", "spring", "circle"),
                                center = NULL, seed = 1L, color_by = NULL,
                                linetype_by = NULL,
                                labels = TRUE, label_size = 3,
                                edge_labels = FALSE, edge_label_size = 3,
                                dismantled = FALSE, ncol = NULL,
                                node_size = 2.5, detail = 5,
                                outline = "white", alpha = 0.45,
                                linewidth = 1.1, padding = 0.045,
                                legend_title = NULL, ...) {
  .thg_check_hg(x)
  stopifnot(
    "`alpha` must be a number in (0, 1]" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha <= 1,
    "`padding` must be a positive number" =
      is.numeric(padding) && length(padding) == 1L && padding > 0,
    "`linewidth` must be a non-negative number" =
      is.numeric(linewidth) && length(linewidth) == 1L && linewidth >= 0,
    "`detail` must be a positive number (Inf for no smoothing)" =
      is.numeric(detail) && length(detail) == 1L && !is.na(detail) && detail > 0,
    "`outline` must be a single colour or \"fill\"" =
      is.character(outline) && length(outline) == 1L && !is.na(outline)
  )
  if (x$n_nodes == 0L) .thg_bad_input("the hypergraph has no nodes to draw")
  drawable <- lengths(x$hyperedges) >= 2L
  if (!any(drawable)) {
    .thg_bad_input("nothing to draw: no hyperedge has two or more members")
  }
  if (isTRUE(dismantled) && !isFALSE(edge_labels)) {
    .thg_bad_input(
      "`edge_labels` does not apply when `dismantled = TRUE`: each panel is already titled with its hyperedge name"
    )
  }
  if (!is.data.frame(layout)) layout <- match.arg(layout)
  if (!is.null(center)) {
    if (!is.character(center) || anyNA(center)) {
      .thg_bad_input("`center` must be a character vector of node names")
    }
    # names absent from this hypergraph are skipped, so one `center` vector
    # serves every subset drawn from the same data; none matching is a typo
    if (!any(center %in% x$nodes)) {
      .thg_bad_input(sprintf("none of the `center` names is a node of the hypergraph: %s",
                             paste(center, collapse = ", ")))
    }
    if (is.data.frame(layout)) {
      .thg_bad_input("`center` cannot be combined with a `layout` table, which fixes every position already")
    }
  }
  fill <- .thg_edge_aesthetic(x, color_by, "color_by")
  ltype <- .thg_edge_aesthetic(x, linetype_by, "linetype_by")
  if (!is.null(ltype) && is.numeric(ltype)) ltype <- as.character(ltype)

  shown <- if (isTRUE(labels)) {
    x$nodes
  } else if (isFALSE(labels)) {
    rep("", x$n_nodes)
  } else {
    if (is.null(names(labels))) .thg_bad_input("`labels` must be TRUE, FALSE or a vector named by node")
    replacement <- unname(labels[x$nodes])
    ifelse(is.na(replacement), x$nodes, as.character(replacement))
  }

  pos <- .thg_positions(x, layout, seed, center)
  edge_names <- pos$edges$hyperedge
  drawn <- which(drawable)
  # largest hyperedges first, so a small pebble is never buried under a big one
  drawn <- drawn[order(-lengths(x$hyperedges)[drawn], drawn)]
  drawn_names <- edge_names[drawn]

  fill_map <- if (is.null(fill)) NULL else .thg_fill_colours(fill)
  edge_colours <- if (is.null(fill_map)) {
    rep(.thg_okabe_ito[[4L]], x$n_hyperedges)
  } else {
    fill_map$colours
  }
  ltype_levels <- if (is.null(ltype)) NULL else
    (if (is.factor(ltype)) levels(ltype) else sort(unique(as.character(ltype))))
  ltype_map <- if (is.null(ltype)) NULL else
    stats::setNames(rep_len(.thg_linetypes, length(ltype_levels)), ltype_levels)
  if (!is.null(fill_map)) {
    attr(fill_map, "title") <- legend_title %||%
      (if (is.character(color_by) && length(color_by) == 1L) color_by else "value")
  }
  if (!is.null(ltype_map)) {
    # one variable driving both gets one title, so ggplot2 merges the legends
    attr(ltype_map, "title") <- if (identical(linetype_by, color_by)) {
      attr(fill_map, "title")
    } else if (is.character(linetype_by) && length(linetype_by) == 1L) {
      linetype_by
    } else {
      "line type"
    }
  }

  hulls <- do.call(rbind, lapply(drawn, function(k) {
    members <- x$hyperedges[[k]]
    ring <- .thg_pebble(pos$nodes$x[members], pos$nodes$y[members], padding,
                        detail)
    ring$hyperedge <- edge_names[k]
    ring$colour <- edge_colours[k]
    ring$value <- if (is.null(fill_map)) NA else fill[k]
    ring$linetype <- if (is.null(ltype)) "solid" else as.character(ltype[k])
    ring
  }))
  hulls$hyperedge <- factor(hulls$hyperedge, levels = drawn_names)
  if (!is.null(fill_map$palette)) {
    hulls$value <- factor(as.character(hulls$value),
                          levels = names(fill_map$palette))
  }

  node_data <- data.frame(x = pos$nodes$x, y = pos$nodes$y, label = shown,
                          stringsAsFactors = FALSE)
  node_ink <- "#2B2B2B"

  if (isTRUE(dismantled)) {
    # every panel repeats the full node set; members are inked, the rest grey
    panel_nodes <- do.call(rbind, lapply(drawn, function(k) {
      member <- seq_len(x$n_nodes) %in% x$hyperedges[[k]]
      data.frame(node_data, hyperedge = edge_names[k], member = member,
                 stringsAsFactors = FALSE)
    }))
    panel_nodes$hyperedge <- factor(panel_nodes$hyperedge, levels = drawn_names)
    # separate layers, not a per-row colour vector: faceting reorders rows,
    # and a constant vector does not travel with them
    members_only <- panel_nodes[panel_nodes$member, , drop = FALSE]
    p <- .thg_hull_layers(ggplot2::ggplot(), hulls, fill_map, ltype_map, alpha,
                          linewidth, outline) +
      ggplot2::geom_point(
        data = panel_nodes[!panel_nodes$member, , drop = FALSE],
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        colour = "#C8C8C8", size = node_size * 0.7
      ) +
      ggplot2::geom_point(
        data = members_only,
        mapping = ggplot2::aes(x = .data$x, y = .data$y),
        shape = 21, fill = node_ink, colour = "white", stroke = 0.5,
        size = node_size
      ) +
      ggplot2::geom_text(
        data = members_only,
        mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
        size = label_size, vjust = -1.1, colour = node_ink
      ) +
      ggplot2::facet_wrap(
        ggplot2::vars(.data$hyperedge),
        ncol = as.integer(ncol %||% ceiling(sqrt(length(drawn))))
      )
    return(.thg_hull_theme(p) +
             ggplot2::theme(strip.text = ggplot2::element_text(
               size = edge_label_size * 3, face = "bold"
             )))
  }

  p <- .thg_hull_layers(ggplot2::ggplot(), hulls, fill_map, ltype_map, alpha,
                        linewidth, outline) +
    ggplot2::geom_point(
      data = node_data, mapping = ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, fill = node_ink, colour = "white", stroke = 0.5,
      size = node_size
    )
  if (!isFALSE(labels)) {
    p <- p + ggplot2::geom_text(
      data = node_data,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = label_size, vjust = -1.1, colour = node_ink
    )
  }

  if (!isFALSE(edge_labels)) {
    edge_shown <- if (isTRUE(edge_labels)) {
      drawn_names
    } else {
      if (is.null(names(edge_labels))) {
        .thg_bad_input("`edge_labels` must be TRUE, FALSE or a vector named by hyperedge")
      }
      replaced <- unname(edge_labels[drawn_names])
      ifelse(is.na(replaced), drawn_names, as.character(replaced))
    }
    anchor_xy <- if (identical(layout, "bipartite")) {
      # the hyperedge's own position, unless the layout left it outside its
      # pebble (a two-member hyperedge's vertex sits off the line between its
      # members); then the centroid of the members, which is always inside
      do.call(rbind, lapply(drawn, function(k) {
        own <- c(pos$edges$x[k], pos$edges$y[k])
        ring <- hulls[hulls$hyperedge == edge_names[k], , drop = FALSE]
        if (.thg_in_convex(ring$x, ring$y, own[1L], own[2L])) return(own)
        members <- x$hyperedges[[k]]
        c(mean(pos$nodes$x[members]), mean(pos$nodes$y[members]))
      }))
    } else {
      # no hyperedge vertex: anchor beyond the member furthest from the
      # middle, where the hull is on a free edge rather than the crowded core
      centre <- c(mean(pos$nodes$x), mean(pos$nodes$y))
      do.call(rbind, lapply(drawn, function(k) {
        m <- cbind(pos$nodes$x, pos$nodes$y)[x$hyperedges[[k]], , drop = FALSE]
        away <- sqrt((m[, 1L] - centre[1L])^2 + (m[, 2L] - centre[2L])^2)
        tip <- m[which.max(away), ]
        v <- tip - centre
        len <- sqrt(sum(v^2))
        if (!is.finite(len) || len == 0) tip else tip + 2.5 * padding * v / len
      }))
    }
    edge_label_data <- data.frame(x = anchor_xy[, 1L], y = anchor_xy[, 2L],
                                  label = edge_shown, stringsAsFactors = FALSE)
    p <- p + ggplot2::geom_label(
      data = edge_label_data,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label),
      size = edge_label_size, fontface = "bold", fill = "white",
      border.colour = edge_colours[drawn], text.colour = node_ink,
      linewidth = 0.4, label.padding = grid::unit(0.15, "lines"),
      inherit.aes = FALSE
    )
  }
  .thg_hull_theme(p)
}

# The shape layer with its fill, outline and line-type scales. With
# `outline = "fill"` the outline shares the fill mapping, so the legend key
# looks like the shapes; otherwise the outline is one fixed colour.
.thg_hull_layers <- function(p, hulls, fill_map, ltype_map, alpha, linewidth,
                             outline) {
  mapped <- !is.null(fill_map)
  own_colour <- identical(outline, "fill")
  aesthetics <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$hyperedge,
                             fill = .data$value)
  if (!mapped) aesthetics$fill <- quote(.data$colour)
  if (own_colour) aesthetics$colour <- aesthetics$fill
  if (!is.null(ltype_map)) aesthetics$linetype <- quote(.data$linetype)
  layer_args <- list(data = hulls, mapping = aesthetics, alpha = alpha,
                     linewidth = linewidth)
  if (!own_colour) layer_args$colour <- outline
  p <- p + do.call(ggplot2::geom_polygon, layer_args)
  scaled <- if (own_colour) c("fill", "colour") else "fill"
  if (!mapped) {
    p <- p + ggplot2::scale_fill_identity(aesthetics = scaled)
  }
  title <- attr(fill_map, "title")
  if (mapped) {
    p <- p + if (is.null(fill_map$palette)) {
      ggplot2::scale_fill_gradientn(colours = .thg_okabe_ito_ramp(),
                                    name = title, aesthetics = scaled,
                                    breaks = .thg_count_breaks(hulls$value))
    } else {
      ggplot2::scale_fill_manual(values = fill_map$palette, name = title,
                                 drop = FALSE, aesthetics = scaled)
    }
  }
  if (!is.null(ltype_map)) {
    p <- p + ggplot2::scale_linetype_manual(values = ltype_map, drop = FALSE,
                                            name = attr(ltype_map, "title")) +
      if (!identical(attr(ltype_map, "title"), title)) {
        ggplot2::guides(linetype = ggplot2::guide_legend(
          override.aes = list(fill = NA, colour = "black", alpha = 1)
        ))
      }
  }
  p
}

# Legend breaks for a numeric colour scale: whole numbers when every value is
# a whole number (hyperedge sizes, counts), so a legend never offers "2.25
# members"; ggplot2's own breaks otherwise.
.thg_count_breaks <- function(values) {
  whole <- all(abs(values - round(values)) < sqrt(.Machine$double.eps))
  if (!whole) return(ggplot2::waiver())
  function(limits) {
    candidates <- unique(round(pretty(limits)))
    candidates[candidates >= limits[1L] & candidates <= limits[2L]]
  }
}

.thg_hull_theme <- function(p) {
  p + ggplot2::coord_equal(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = 0.08)) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(legend.position = "right",
                   plot.margin = ggplot2::margin(8, 8, 8, 8))
}

# Sequential ramp built from Okabe-Ito colours: light yellow through green to
# blue, matching the yellow-to-blue ordering of the paper's cardinality scale.
.thg_okabe_ito_ramp <- function() {
  c("#F0E442", "#009E73", "#0072B2")
}
