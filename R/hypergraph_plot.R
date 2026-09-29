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
# layout keeps its shape. `padding` is the pebble reach the packing has to
# leave between disconnected components.
.thg_positions <- function(hg, layout, seed, center = NULL, padding = 0.045) {
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
    xy <- .thg_layout_packed(adjacency, seed = seed,
                             pinned = c(centred, rep(FALSE, m)),
                             padding = padding)
  } else {
    layout <- match.arg(layout, c("spring", "circle"))
    projection <- as.matrix(hg_project(hg, method = "clique",
                                       weighted = FALSE, what = "matrix"))
    xy <- if (identical(layout, "spring")) {
      .thg_layout_packed((projection != 0) * 1, seed = seed, pinned = centred,
                         padding = padding)
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

# Force-directed placement, component by component. Fruchterman-Reingold
# repels every pair but attracts only along edges, so two disconnected
# components feel repulsion and nothing else: they drift apart for as long as
# the temperature allows, and because the frame is then scaled by its widest
# span, the islands set the scale and the connected structure is squeezed into
# a corner of it. Laying each component out on its own and packing the boxes
# afterwards keeps every component at its own natural size and spends the
# frame on structure instead of on empty space between islands. One component
# is the ordinary case and goes straight to the placement.
.thg_layout_packed <- function(adjacency, seed = 1L,
                               pinned = rep(FALSE, nrow(adjacency)),
                               padding = 0.045) {
  component <- .thg_components(adjacency)
  if (max(component) == 1L) {
    return(.thg_layout_fr(adjacency, seed = seed, pinned = pinned))
  }
  parts <- lapply(seq_len(max(component)), function(k) {
    keep <- component == k
    xy <- if (sum(keep) == 1L) {
      matrix(0, 1L, 2L)
    } else {
      .thg_layout_fr(adjacency[keep, keep, drop = FALSE], seed = seed,
                     pinned = pinned[keep])
    }
    sweep(xy, 2L, apply(xy, 2L, min))
  })
  size <- t(vapply(parts, function(p) apply(p, 2L, max), numeric(2L)))
  # The gap between components has to outrun the pebbles, whose reach is
  # `padding` of the finished extent -- but that extent is itself what the gap
  # sets, so the two are settled by iteration from one ideal edge length. It
  # converges from below: a wider gap only widens the extent sub-linearly.
  gap <- Reduce(function(g, pass) {
    extent <- max(apply(.thg_pack(parts, size, g, component), 2L,
                        function(z) diff(range(z))))
    max(1, 2.5 * padding * extent)
  }, seq_len(3L), init = 1)
  .thg_pack(parts, size, gap, component)
}

# Place the component layouts at their packed offsets and restore the original
# vertex order. `parts` is in component order, and each part holds its
# component's vertices in increasing original index, which is exactly how
# split() groups them.
.thg_pack <- function(parts, size, gap, component) {
  offset <- .thg_shelf_offsets(size, gap)
  stacked <- do.call(rbind, Map(function(p, k) sweep(p, 2L, offset[k, ], "+"),
                                parts, seq_along(parts)))
  rows <- unlist(split(seq_along(component), component), use.names = FALSE)
  stacked[order(rows), , drop = FALSE]
}

# Connected components of an undirected adjacency matrix as labels 1..K.
# Label propagation: each vertex repeatedly takes the smallest label in its
# closed neighbourhood, settling after as many passes as the graph is wide.
# Graph traversal is the justified exception to the no-loop rule.
.thg_components <- function(adjacency) {
  neighbourhood <- adjacency != 0
  diag(neighbourhood) <- TRUE
  label <- seq_len(nrow(adjacency))
  repeat {
    settled <- apply(neighbourhood, 1L, function(near) min(label[near]))
    if (identical(settled, label)) break
    label <- settled
  }
  match(label, sort(unique(label)))
}

# Shelf packing of the component boxes into a roughly square frame: boxes are
# laid tallest first along a shelf until the next one would overrun the side
# of a square of the same total area, then a new shelf starts above. Shelves
# are centred horizontally and their boxes centred vertically, so the packing
# reads as one figure rather than a ragged staircase. A square frame is what
# makes coord_equal() spend the whole panel on the picture. Returns one (x, y)
# offset per box, in the order the boxes were given.
.thg_shelf_offsets <- function(size, gap) {
  wide <- size[, 1L] + gap
  high <- size[, 2L] + gap
  tall <- order(-high, -wide)
  target <- max(sqrt(sum(wide * high)), max(wide))
  # a running total is sequential by nature, so the shelf a box lands on is
  # carried forward by Reduce() rather than recomputed
  running <- Reduce(function(state, w) {
    if (state[["used"]] > 0 && state[["used"]] + w > target) {
      c(shelf = state[["shelf"]] + 1, used = w)
    } else {
      c(shelf = state[["shelf"]], used = state[["used"]] + w)
    }
  }, wide[tall], init = c(shelf = 1, used = 0), accumulate = TRUE)[-1L]
  shelf <- vapply(running, function(s) s[["shelf"]], numeric(1L))
  index <- match(shelf, sort(unique(shelf)))
  shelf_height <- vapply(split(high[tall], index), max, numeric(1L))
  left <- stats::ave(wide[tall], index, FUN = function(w) cumsum(w) - w)
  shelf_width <- stats::ave(wide[tall], index, FUN = sum)
  offset <- matrix(0, nrow(size), 2L)
  offset[tall, 1L] <- left + (target - shelf_width) / 2 + gap / 2
  offset[tall, 2L] <- cumsum(c(0, shelf_height))[index] +
    (shelf_height[index] - high[tall]) / 2 + gap / 2
  offset
}

# Distance beyond which two vertices stop repelling, in ideal edge lengths
# (k = 1 here, the layout being spread over an area of n). Unbounded
# repulsion leaves a weakly attached cluster no equilibrium: it is pushed by
# every vertex in the graph and pulled back by one edge, so it settles roughly
# n^(1/3) ideal lengths out and the hyperedge bridging it stretches across the
# page. Cutting repulsion off beyond a few ideal lengths is Fruchterman and
# Reingold's own grid variant, whose whole point is that distant pairs
# contribute nothing. On a core-plus-satellite fixture over 30 seeds, 2.5
# brought the satellite from 0.84 of the frame away to 0.70, widened the core
# from 0.31 of the frame to 0.43 and narrowed the bridging pebble from 0.43 to
# 0.37; it also left the fewest non-members inside a pebble of the values
# tried (0.20 per layout against 0.32 uncut). See tmp/measure_layout.R.
.thg_repulsion_cutoff <- 2.5

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
                           pinned = rep(FALSE, nrow(adjacency)),
                           cutoff = .thg_repulsion_cutoff) {
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
    force <- (d < cutoff) / d^2 - adjacency * d
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
  # sampled offsets are only nearly exact, so re-measure and push again, at
  # most three times (two suffice in practice); a pass with no deficit left
  # changes nothing
  curve <- Reduce(function(outline, pass) {
    deficit <- max(clearance * radius -
                     .thg_inside_distance(outline, normal, members))
    if (deficit <= 0) outline else outline + deficit * normal
  }, seq_len(3L), curve)
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
#' A hypergraph whose hyperedges carry one numeric attribute -- such as a
#' `trials` count that [group_hypergraph()] kept because it is constant
#' within each group -- is drawn by that attribute with no further argument:
#' it colours the pebbles, writes a title box with the count beside each
#' one, and names the unit of the legends. Of several numeric attributes,
#' the one that varies between hyperedges is read when only one does (a
#' membership share of 1 on every hyperedge is passed over); otherwise none
#' is. `color_by`, `titles` and `unit` override it.
#'
#' Node labels are bold with a white halo, the legends sit below the plot
#' (a colour bar 2.5 cm per key, discrete legends at most four keys to a
#' row) and the plot has wide margins (40, 130, 30 and 130 pt) so title
#' boxes outside the pebbles are not cut. The look of
#' earlier versions is partly available through arguments: `alpha = 0.45`,
#' `label_size = 3` and `pieces = "packed"`; the legend position and margins
#' through `ggplot2::theme()` added to the result.
#'
#' A hypergraph that falls into several disconnected pieces is laid out one
#' piece at a time and the pieces are then packed into a roughly square frame.
#' A force-directed layout has no force at all between two disconnected
#' components -- they repel and nothing pulls back -- so laying the whole
#' hypergraph out at once lets a stray hyperedge drift to the edge of the
#' picture and set the scale for everything else, leaving the connected
#' structure a speck in the middle. Packing keeps every component at its own
#' size and the frame spent on structure.
#'
#' @param x A `net_hypergraph` with at least one hyperedge of two or more
#'   members. Draw a part of a large hypergraph by passing [hg_subset()]
#'   first.
#' @param layout `"bipartite"` (default; a Fruchterman-Reingold layout of the
#'   star expansion, nodes and hyperedges placed together), `"spring"` (the
#'   same on the clique projection), `"circle"`, or a data.frame with `node`,
#'   `x` and `y` columns to reuse node coordinates across panels (hyperedges
#'   then sit at the centroid of their members). The two force-directed
#'   layouts place each connected component separately and pack the results;
#'   a `layout` table is used exactly as given.
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
#' @param color_by Colour of the hulls: `NULL` (the hyperedges' count
#'   attribute described above, else one colour), `"size"`
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
#' @param edge_label_size Text size for `edge_labels` (default `3`), and for
#'   the panel titles when `dismantled = TRUE`. A title is never cut to fit
#'   its panel -- a clipped hyperedge name reads as a different hyperedge --
#'   so lower this, or `ncol`, if long names collide.
#' @param labels `TRUE` (default) writes the node names, `FALSE` writes
#'   none, and a character vector named by node replaces the names shown.
#' @param label_size,node_size Text and point sizes (defaults `4.2` and
#'   `2.5`).
#' @param detail How much of the hull's shape the smoothing keeps: the width,
#'   in harmonics, of the Gaussian low-pass applied to the outline. The
#'   default `5` gives rounded, tapering petals; `3` is rounder still, `8`
#'   stays close to the members, and `Inf` draws the unsmoothed rounded hull.
#' @param outline Outline colour: `"white"` (default) draws seams that part
#'   overlapping hyperedges; `"fill"` outlines each hyperedge in its own
#'   colour (pair it with a lower `alpha`, e.g. `0.15`); any other colour is
#'   used as given.
#' @param alpha Fill transparency (default `0.5`); the outline is always
#'   drawn at full strength.
#' @param linewidth Width of the outlines (default `1.1`).
#' @param padding Room around the member positions as a fraction of the layout
#'   extent (default `0.045`). It also sets how far apart the layout packs
#'   disconnected components, so that two components' pebbles never touch and
#'   imply a member they do not share.
#' @param legend_title Legend title for `color_by`; defaults to the
#'   selector's name.
#' @param node_sizes Size every node by a value, such as how often the node
#'   occurs: a data.frame with a node column and a value column, or a numeric
#'   vector named by node. The columns are `node` and `value` when those names
#'   are present; otherwise the one character (or factor) column holds the
#'   node names and the one numeric column the values, so a table such as
#'   `data.frame(state, trials)` is read as it is. Every node of the
#'   hypergraph needs a non-negative value; names that are not nodes are
#'   skipped, so one table serves every [hg_subset()]. Without `direction`
#'   the nodes are points whose area follows the value
#'   ([ggplot2::scale_size_area()], largest 10 mm), with a size legend, and
#'   each label sits just above its dot. With `direction` they are circles in
#'   data units (the figure has a fixed 1:1 aspect): the largest has a radius
#'   of 4% of the layout extent and none is smaller than 30% of that radius,
#'   so the rarest nodes stay visible, and a caption says what area and
#'   triangle show. `NULL` (default) draws plain points.
#' @param direction Mark on each circle where the walk goes next: a
#'   data.frame with columns `from`, `to` and a weight -- `weight`, or else
#'   the one other numeric column, such as `trials` (how often `to` follows
#'   `from`; repeated pairs are summed). Each node gets a triangle pointing at
#'   the other node that follows it with the largest weight, ties broken by
#'   node name; a node that nothing else follows (a sink, or only itself)
#'   gets no triangle. Rows naming a node that is not in the hypergraph are
#'   skipped. Requires `node_sizes`.
#' @param arrow_style `"inside"` (default; a triangle inscribed in the
#'   circle, its apex on the rim facing the next node) or `"outside"` (the
#'   circle drawn out into a tip beyond its rim). An outside tip that points
#'   upward moves that node's label below the circle.
#' @param node_fill Colour of the circles drawn by `node_sizes` with
#'   `direction` (default black).
#' @param arrow_fill Colour of the `direction` triangles. `NULL` (default) is
#'   a dark slate (`"#1F3A44"`) for `"inside"` and `node_fill` for
#'   `"outside"`.
#' @param transitions Draw the moves between nodes as curved arrows: a
#'   data.frame with `from`, `to` and a weight, read like `direction`. Line
#'   width follows the weight; a node followed by itself gets a small loop
#'   below it. Rows naming a node that is not in the hypergraph are skipped.
#'   The width legend is titled "<Unit> making the transition" when `unit` is
#'   given and "Transition weight" otherwise.
#' @param size_title What the node area shows (e.g. `"trials with the
#'   event"`): the caption reads "Circle area: <size_title>." and the point
#'   size legend is titled with it. `NULL` (default) uses "<unit> with the
#'   event" when there is a unit, and otherwise says the area is
#'   proportional to the event value.
#' @param titles Write a title box beside every hyperedge, just outside its
#'   pebble on the side facing away from the centre of the figure. `NULL`
#'   (default) writes the hyperedges' count attribute described above (not
#'   when `dismantled = TRUE`), else none; `FALSE` writes
#'   none; `TRUE` writes the hyperedge names; a selector
#'   resolved like `color_by` (an edge-metadata column, a vector named by
#'   hyperedge, or one value per hyperedge) must be numeric and adds a second
#'   line with that number and `unit` (e.g. "3,667 trials"). A hyperedge with
#'   a single member gets its box beside that node. Boxes are placed in
#'   order of decreasing value (name breaks ties) and a box that would sit on
#'   an earlier one moves below it.
#' @param title_prefix Text written before each hyperedge name in its title
#'   box (default `""`), such as `"Mostly "`.
#' @param notes A further line for some title boxes: a character vector named
#'   by hyperedge, or a data.frame with a `note` column and the hyperedge
#'   names in a `hyperedge`, `edge` or `group` column. Hyperedges without a
#'   note get none; names that are not hyperedges are skipped. Requires
#'   `titles`.
#' @param unit The word for what the values count (e.g. `"trials"`,
#'   `"steps"`). It titles the `color_by` legend (capitalised, unless
#'   `legend_title` is given), follows the number in each title box, and
#'   names the node area ("<unit> with the event") and the transition
#'   widths. `NULL` (default) is the unit the hypergraph records
#'   (`"sequences"` for [group_hypergraph()] of a clustering), else the
#'   name of the hyperedges' count attribute described above, else those
#'   texts stay generic.
#' @param pieces How the bipartite layout places a hypergraph that falls into
#'   disconnected pieces: `"row"` (default; every piece is laid out in a
#'   frame of its own, as if drawn alone, and the frames are set side by
#'   side 1.4 frame widths apart, in order of the piece's most prominent
#'   hyperedge -- see `titles` -- with a piece of one node at the middle of
#'   its frame) or `"packed"` (pieces are laid out one by one and packed into
#'   a square frame). Only the bipartite layout has pieces; `"row"` given
#'   with another layout is an error.
#' @param title_gap Distance of each title box beyond its pebble, in layout
#'   units (the layout spans 0 to 1; default `0.06`).
#' @param group Draw one group only: the name of a value of the hyperedges'
#'   `group` attribute, such as `"Cluster 1"` of a [group_hypergraph()] built
#'   from a clustering of sequences. Only that group's hyperedges and their
#'   members are drawn; for clustered sequences each node is sized by the
#'   sequences of the group containing it (unless `node_sizes` is given),
#'   each set is named by the states it adds to those every drawn set
#'   shares (listed in the subtitle), the figure is titled with the group
#'   and its number of sequences, `titles` defaults to `FALSE` and the
#'   legends are stacked. A
#'   name that is not a group, or a hypergraph without a `group` attribute,
#'   raises `hypernets_bad_input`. `NULL` (default) draws every hyperedge.
#' @param ... Unused; for S3 consistency.
#' @return A ggplot object (with `dismantled = TRUE`, one facet per
#'   hyperedge). With `node_sizes` and `direction` the nodes are polygon
#'   layers in data units and the caption names what circle area and
#'   triangle show; with `node_sizes` alone they are points with a size
#'   legend. Title boxes are the last layer.
#' @section Conditions:
#' `hypernets_bad_input` for an invalid selector or layout, a `node_sizes`
#' that misses a node or holds a negative value, a `direction` or
#' `transitions` table without `from`, `to` and one non-negative numeric
#' weight, a `direction` without `node_sizes`, `titles` that are not
#' numeric, `notes` without titles, an invalid `unit`, `title_prefix` or
#' `title_gap`, `pieces = "row"` with a layout other than
#' `"bipartite"`, or node overlays or titles combined with
#' `dismantled = TRUE`.
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
#'
#' # disconnected pieces are laid out separately and packed
#' apart <- group_hypergraph(
#'   data.frame(member = c("a", "b", "c", "b", "c", "d", "x", "y", "z"),
#'              event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e3")),
#'   "member", "event"
#' )
#' plot(apart, color_by = "size")
#'
#' # circles sized by how often each node occurs, pointing where it goes next
#' counts <- c(a = 40, b = 25, c = 12, d = 30, e = 5, f = 2)
#' moves <- data.frame(from = c("a", "b", "c", "d", "d", "e"),
#'                     to = c("b", "c", "a", "e", "b", "f"),
#'                     weight = c(20, 9, 4, 8, 3, 1))
#' plot(hg, node_sizes = counts, direction = moves,
#'      size_title = "occurrences of the node")
#' plot(hg, node_sizes = counts, direction = moves, transitions = moves,
#'      arrow_style = "outside")
#'
#' # groups of events with their trial counts: the count colours the pebbles
#' # and writes a title box per group; points sized by how many trials
#' # contain each event
#' members <- data.frame(
#'   group = c("Hint", "Hint", "Hint", "Question", "Question", "Question"),
#'   state = c("Wrong", "Hint", "Retry", "Wrong", "Question", "Retry"),
#'   trials = c(120, 120, 120, 80, 80, 80)
#' )
#' trial_groups <- group_hypergraph(members, actor = "state", group = "group")
#' event_trials <- data.frame(state = c("Wrong", "Hint", "Question", "Retry"),
#'                            trials = c(200, 120, 80, 190))
#' plot(trial_groups, node_sizes = event_trials)
#' @export
plot.net_hypergraph <- function(x, layout = c("bipartite", "spring", "circle"),
                                center = NULL, seed = 1L, color_by = NULL,
                                linetype_by = NULL,
                                labels = TRUE, label_size = 4.2,
                                edge_labels = FALSE, edge_label_size = 3,
                                dismantled = FALSE, ncol = NULL,
                                node_size = 2.5, detail = 5,
                                outline = "white", alpha = 0.5,
                                linewidth = 1.1, padding = 0.045,
                                legend_title = NULL, node_sizes = NULL,
                                direction = NULL,
                                arrow_style = c("inside", "outside"),
                                node_fill = "#000000", arrow_fill = NULL,
                                transitions = NULL, size_title = NULL,
                                titles = NULL, title_prefix = "",
                                notes = NULL, unit = NULL, pieces = NULL,
                                title_gap = 0.06, group = NULL, ...) {
  .thg_check_hg(x)
  group_title <- NULL
  if (!is.null(group)) {
    selected <- .thg_select_group(x, group)
    x <- selected$hypergraph
    node_sizes <- node_sizes %||% selected$node_sizes
    # one group's sets are told apart by their members and the subtitle, as
    # in a table of the group's sets; title boxes would crowd the pebbles
    titles <- titles %||% FALSE
    group_title <- selected[c("title", "subtitle")]
  }
  arrow_style <- match.arg(arrow_style)
  explicit_pieces <- !is.null(pieces)
  pieces <- match.arg(pieces %||% "row", c("packed", "row"))
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
  overlaid <- !is.null(node_sizes) || !is.null(direction) ||
    !is.null(transitions)
  if (isTRUE(dismantled) && overlaid) {
    .thg_bad_input(
      "`node_sizes`, `direction` and `transitions` do not apply when `dismantled = TRUE`"
    )
  }
  if (!is.null(direction) && is.null(node_sizes)) {
    .thg_bad_input("`direction` needs `node_sizes`: the triangles sit inside the circles")
  }
  .thg_check_string(node_fill, "node_fill", colour = TRUE)
  .thg_check_string(arrow_fill, "arrow_fill", null_ok = TRUE, colour = TRUE)
  .thg_check_string(size_title, "size_title", null_ok = TRUE)
  .thg_check_string(unit, "unit", null_ok = TRUE)
  .thg_check_string(title_prefix, "title_prefix")
  if (!is.numeric(title_gap) || length(title_gap) != 1L || !is.finite(title_gap) ||
      title_gap < 0) {
    .thg_bad_input("`title_gap` must be one non-negative number")
  }
  # A numeric hyperedge attribute (such as a count of trials that
  # group_hypergraph() kept because it is constant within each group) is
  # what the pebbles show unless the call says otherwise: it colours them,
  # writes the title boxes and names the unit. Of several, an attribute with
  # one value on every hyperedge (a membership share of 1 throughout) tells
  # them apart no better than none and gives way to the one that varies;
  # otherwise the call has to choose.
  numeric_attributes <- Filter(function(column) {
    is.numeric(x$edge_data[[column]]) && all(is.finite(x$edge_data[[column]]))
  }, setdiff(names(x$edge_data), "edge"))
  varying <- Filter(function(column) length(unique(x$edge_data[[column]])) > 1L,
                    numeric_attributes)
  count <- if (length(numeric_attributes) == 1L) numeric_attributes else
    if (length(varying) == 1L) varying
  color_by <- color_by %||% count
  unit <- unit %||% x$params$unit %||% count
  titles <- titles %||% (if (is.null(count) || isTRUE(dismantled)) FALSE else count)
  titled <- !isFALSE(titles)
  if (!is.null(notes) && !titled) {
    .thg_bad_input("`notes` are written into the title boxes: give `titles` too")
  }
  if (isTRUE(dismantled) && titled) {
    .thg_bad_input("`titles` do not apply when `dismantled = TRUE`: each panel is already titled with its hyperedge name")
  }
  title_values <- if (is.logical(titles) && length(titles) == 1L && !is.na(titles)) {
    NULL
  } else {
    values <- .thg_edge_aesthetic(x, titles, "titles")
    if (!is.numeric(values) || any(!is.finite(values))) {
      .thg_bad_input("`titles` must be TRUE, FALSE or a selector of finite numbers (a count per hyperedge)")
    }
    values
  }
  title_notes <- if (is.null(notes)) NULL else
    .thg_keyed_values(notes, colnames(x$incidence), "notes",
                      key_columns = c("hyperedge", "edge", "group"),
                      value_column = "note", type = "character")
  if (explicit_pieces && identical(pieces, "row") &&
      !identical(layout, "bipartite")) {
    .thg_bad_input("`pieces = \"row\"` applies to the bipartite layout only")
  }
  node_value <- if (is.null(node_sizes)) NULL else
    .thg_keyed_values(node_sizes, x$nodes, "node_sizes", key_columns = "node",
                      value_column = "value")
  direction <- if (is.null(direction)) NULL else
    .thg_move_table(direction, "direction", x$nodes)
  transitions <- if (is.null(transitions)) NULL else
    .thg_move_table(transitions, "transitions", x$nodes)
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

  # the order in which hyperedges claim room: title boxes, and the pieces of
  # a row layout, go by decreasing title value (else numeric colour value),
  # then by name
  edge_names <- colnames(x$incidence) %||% paste0("h", seq_len(x$n_hyperedges))
  prominence <- title_values %||% (if (is.numeric(fill)) fill)
  precedence <- if (is.null(prominence)) order(edge_names) else
    order(-prominence, edge_names)
  # decided before a row of pieces turns the layout into a table: every
  # hyperedge of a bipartite layout has a position inside its pebble (its own
  # vertex, or its members' centroid), where its label goes
  bipartite <- identical(layout, "bipartite")
  if (identical(pieces, "row") && bipartite) {
    row <- .thg_row_layout(x, seed, center, padding, precedence)
    if (!is.null(row)) {
      layout <- row
      center <- NULL
    }
  }
  pos <- .thg_positions(x, layout, seed, center, padding)
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
      (if (!is.null(unit)) tools::toTitleCase(unit)) %||%
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
    # strip.clip = "off": a panel title is a hyperedge name, and ggplot2
    # otherwise cuts it to the panel's width, which silently turns
    # "32006L0123" into "2006L012" -- an identifier that reads as a different
    # hyperedge. Letting a long name overflow its strip is visible and the
    # reader can act on it; a truncated one is not.
    return(.thg_hull_theme(p) +
             ggplot2::theme(strip.clip = "off",
                            strip.text = ggplot2::element_text(
                              size = edge_label_size * 3, face = "bold"
                            )))
  }

  p <- .thg_hull_layers(ggplot2::ggplot(), hulls, fill_map, ltype_map, alpha,
                        linewidth, outline)
  circles <- !is.null(node_value) && !is.null(direction)
  geometry <- if (circles) {
    .thg_node_geometry(node_data, x$nodes, node_value, direction, arrow_style)
  }
  label_data <- .thg_label_positions(node_data, node_value, geometry)
  p <- p + .thg_node_layers(
    node_data, x$nodes, label_data, value = node_value, geometry = geometry,
    transitions = transitions, style = arrow_style, node_fill = node_fill,
    arrow_fill = arrow_fill, labels = labels, label_size = label_size,
    node_size = node_size,
    area_what = size_title %||% (if (!is.null(unit)) paste(unit, "with the event")),
    unit = unit
  )

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
    anchor_xy <- if (bipartite) {
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
  if (titled) {
    boxes <- .thg_title_boxes(hulls, pos, x$hyperedges, edge_names, precedence,
                              title_values, title_prefix, unit, title_notes,
                              gap = title_gap)
    p <- p + ggplot2::geom_label(
      data = boxes,
      mapping = ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                             hjust = .data$hjust),
      inherit.aes = FALSE, size = 3.4, fontface = "bold", lineheight = 0.9,
      linewidth = 0.3, fill = "white"
    )
  }
  # legends below the plot, wide margins so title boxes outside the pebbles
  # are not cut; three legends do not fit on one row
  p <- .thg_hull_theme(p) +
    ggplot2::theme(plot.margin = ggplot2::margin(40, 130, 30, 130),
                   legend.position = "bottom",
                   legend.key.width = grid::unit(2.5, "cm"))
  if (!is.null(transitions)) p <- p + ggplot2::theme(legend.box = "vertical")
  if (!is.null(group_title)) {
    # two legends side by side overrun a one-group figure; stack them
    p <- p + ggplot2::labs(title = group_title$title,
                           subtitle = group_title$subtitle) +
      ggplot2::theme(legend.box = "vertical")
  }
  p
}

# One group of a hypergraph whose hyperedges carry a `group` attribute: its
# hyperedges, its nodes, and -- for clustered sequences -- the node sizes of
# that group (sequences containing the state) and a title naming it.
.thg_select_group <- function(x, group) {
  if (!is.character(group) || length(group) != 1L || is.na(group)) {
    .thg_bad_input("`group` must be one group name")
  }
  groups <- x$edge_data$group
  if (is.null(groups)) {
    .thg_bad_input("`group` needs hyperedges with a `group` attribute, such as group_hypergraph() of a clustering")
  }
  if (!group %in% groups) {
    .thg_bad_input(sprintf("`group` \"%s\" is not a group of this hypergraph; groups: %s",
                           group, paste(unique(groups), collapse = ", ")))
  }
  hypergraph <- hg_subset(x, where = list(group = group))
  counts <- x$state_counts
  node_sizes <- NULL
  title <- group
  subtitle <- NULL
  if (!is.null(counts)) {
    # a set is named by what it adds to the states every drawn set shares,
    # so the title boxes stay short; the shared states go in the subtitle
    members <- lapply(hypergraph$hyperedges, \(idx) hypergraph$nodes[idx])
    core <- Reduce(intersect, members)
    added <- vapply(members, \(m) {
      extra <- setdiff(m, core)
      if (length(extra)) paste(extra, collapse = " + ") else "(shared states only)"
    }, character(1L))
    added <- make.unique(added, sep = " ")
    old <- colnames(hypergraph$incidence)
    colnames(hypergraph$incidence) <- added
    hypergraph$edge_data$edge <- added[match(hypergraph$edge_data$edge, old)]
    # hyperedges in name order, as group_hypergraph() orders a table of the
    # group's sets, so the seeded layout is the one that table would get
    by_name <- match(sort(added), added)
    hypergraph <- .thg_rebuild(hypergraph,
                               hypergraph$incidence[, by_name, drop = FALSE],
                               by_name)
    if (length(core)) {
      subtitle <- sprintf("Shared by every set: %s", paste(core, collapse = ", "))
    }
    in_group <- counts$group == group
    node_sizes <- stats::setNames(counts$count[in_group], counts$node[in_group])
    n <- x$group_sizes$sequences[match(group, x$group_sizes$group)]
    title <- sprintf("%s (%s %s)", group,
                     formatC(n, big.mark = ",", format = "d"),
                     x$params$unit %||% "sequences")
  }
  list(hypergraph = hypergraph, node_sizes = node_sizes, title = title,
       subtitle = subtitle)
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
      # the wide keys of the theme are for the colour bar: a discrete legend
      # keeps the default key width and at most four keys to a row, so a
      # legend of many levels stays narrow
      ggplot2::scale_fill_manual(
        values = fill_map$palette, name = title, drop = FALSE,
        aesthetics = scaled,
        guide = ggplot2::guide_legend(
          nrow = ceiling(length(fill_map$palette) / 4), byrow = TRUE,
          theme = ggplot2::theme(legend.key.width = grid::unit(1.2, "lines"))
        )
      )
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

# ---- node overlays: sized circles, direction triangles, transitions ---------
# Ported from the pipeline helpers the package author drew these figures with
# (plot_blobs(), .direction_node_layers(), .transition_arrows()), rebuilt as
# ordinary layers so no caller edits a finished plot's layers by index.

# Stop with hypernets_bad_input unless `x` is one non-missing string -- one
# grDevices can read as a colour when `colour = TRUE` -- or, with
# `null_ok = TRUE`, NULL.
.thg_check_string <- function(x, arg, null_ok = FALSE, colour = FALSE) {
  if (null_ok && is.null(x)) return(invisible(NULL))
  ok <- is.character(x) && length(x) == 1L && !is.na(x) &&
    (!colour ||
       tryCatch(is.matrix(grDevices::col2rgb(x)), error = function(e) FALSE))
  if (!ok) {
    .thg_bad_input(sprintf("`%s` must be %s%s", arg,
                           if (null_ok) "NULL or " else "",
                           if (colour) "a single colour" else "a single string"))
  }
  invisible(x)
}

# Numbers as text with a thousands separator, never in scientific notation.
# Legend labels are each as short as they can be written ("12,000", "2.5");
# `count = TRUE`, for title boxes, writes whole numbers as integers ("3,667")
# and others to three significant digits.
.thg_comma <- function(x, count = FALSE) {
  vapply(x, function(v) {
    if (is.na(v)) return(NA_character_)
    if (!count) return(format(v, big.mark = ",", scientific = FALSE, trim = TRUE))
    if (abs(v - round(v)) < sqrt(.Machine$double.eps)) {
      formatC(round(v), big.mark = ",", format = "d")
    } else {
      formatC(v, digits = 3L, format = "fg", big.mark = ",")
    }
  }, character(1L))
}

# The one column of `data` whose values satisfy `test`, excluding `skip`:
# how a two-column table such as (state, trials) is read without naming its
# columns. NULL when there is not exactly one.
.thg_only_column <- function(data, test, skip = character()) {
  hits <- Filter(function(column) test(data[[column]]),
                 setdiff(names(data), skip))
  if (length(hits) == 1L) hits else NULL
}

# Values keyed by name, returned in the order of `keys` (NA where a key has
# none; the first value given for a name counts), from a vector named by key
# or a data.frame. The data.frame's key column is the first of `key_columns`
# it has -- else, for a numeric `type`, its one text column -- and its value
# column is `value_column`, else, for a numeric `type`, its one numeric
# column. A numeric table (`node_sizes`) must also name every key once, with
# finite non-negative values, at least one positive. Names that are not keys
# are skipped, so one table serves every subset.
.thg_keyed_values <- function(tab, keys, arg, key_columns, value_column,
                              type = c("numeric", "character")) {
  type <- match.arg(type)
  numeric <- identical(type, "numeric")
  is_type <- if (numeric) is.numeric else is.character
  if (is.data.frame(tab)) {
    key <- intersect(key_columns, names(tab))[1L]
    if (is.na(key) && numeric) {
      key <- .thg_only_column(tab, function(v) is.character(v) || is.factor(v)) %||%
        NA_character_
    }
    value <- if (value_column %in% names(tab)) value_column else if (numeric) {
      .thg_only_column(tab, is.numeric, skip = key)
    }
    if (is.na(key) || is.null(value)) {
      .thg_bad_input(sprintf(
        "a `%s` data.frame needs a `%s` column and the names in a %s column%s",
        arg, value_column,
        paste(sprintf("`%s`", key_columns), collapse = " or "),
        if (numeric) ", or exactly one text column (the names) and one numeric column (the values)" else ""))
    }
    names_given <- as.character(tab[[key]])
    values <- tab[[value]]
  } else {
    if (!is_type(tab) || is.null(names(tab))) {
      .thg_bad_input(sprintf("`%s` must be a data.frame or a %s vector named by %s",
                             arg, type, key_columns[[1L]]))
    }
    names_given <- names(tab)
    values <- unname(tab)
  }
  if (!is_type(values)) .thg_bad_input(sprintf("`%s` values must be %s", arg, type))
  if (!numeric) return(values[match(keys, names_given)])
  if (anyDuplicated(names_given)) .thg_bad_input(sprintf("`%s` names a node more than once", arg))
  missing <- setdiff(keys, names_given)
  if (length(missing)) {
    .thg_bad_input(sprintf("`%s` has no value for %d node(s): %s", arg,
                           length(missing),
                           paste(utils::head(missing, 5L), collapse = ", ")))
  }
  out <- values[match(keys, names_given)]
  if (any(!is.finite(out)) || any(out < 0)) {
    .thg_bad_input(sprintf("`%s` values must be finite and non-negative", arg))
  }
  if (max(out) <= 0) .thg_bad_input(sprintf("`%s` needs at least one positive value", arg))
  out
}

# A from/to/weight table (the weight named `weight`, or the one other
# numeric column, such as `trials`) restricted to pairs of drawn nodes,
# repeated pairs summed, zero weights dropped (a pair that never happens is
# no move), in a deterministic order.
.thg_move_table <- function(tab, arg, nodes) {
  weight_column <- if (!is.data.frame(tab)) NULL else if ("weight" %in% names(tab)) {
    "weight"
  } else {
    .thg_only_column(tab, is.numeric, skip = c("from", "to"))
  }
  if (!is.data.frame(tab) || !all(c("from", "to") %in% names(tab)) ||
      is.null(weight_column)) {
    .thg_bad_input(sprintf(
      "`%s` must be a data.frame with `from`, `to` and a weight: a `weight` column, or exactly one other numeric column",
      arg))
  }
  weight <- tab[[weight_column]]
  if (!is.numeric(weight) || any(!is.finite(weight)) || any(weight < 0)) {
    .thg_bad_input(sprintf("`%s$weight` must be finite non-negative numbers", arg))
  }
  from <- as.character(tab$from)
  to <- as.character(tab$to)
  if (anyNA(from) || anyNA(to)) {
    .thg_bad_input(sprintf("`%s` has missing `from` or `to` names", arg))
  }
  inside <- from %in% nodes & to %in% nodes
  if (length(from) && !any(inside)) {
    .thg_bad_input(sprintf("no row of `%s` joins two nodes of the hypergraph", arg))
  }
  keep <- inside & weight > 0
  if (!any(keep)) {
    return(data.frame(from = character(), to = character(), weight = numeric()))
  }
  out <- stats::aggregate(weight ~ from + to,
                          data = data.frame(from = from[keep], to = to[keep],
                                            weight = weight[keep]),
                          FUN = sum)
  out <- out[order(out$from, out$to), , drop = FALSE]
  rownames(out) <- NULL
  out
}

# Circle radius and heading for every node. Radius: area proportional to
# the value, the largest 4% of the layout extent, none below 30% of that.
# Heading: towards the other node that follows with the largest weight, ties
# broken by name; NA for a node nothing else follows. `below`: an outside
# tip points up, so that node's label goes under its circle.
.thg_node_geometry <- function(node_data, nodes, value, direction,
                               style = "inside") {
  span <- max(diff(range(node_data$y)), diff(range(node_data$x)), 1e-9)
  r_max <- 0.04 * span
  radius <- pmax(r_max * sqrt(value / max(value)), 0.3 * r_max)
  target <- rep(NA_integer_, length(nodes))
  if (!is.null(direction) && nrow(direction)) {
    moves <- direction[direction$from != direction$to, , drop = FALSE]
    moves <- moves[order(moves$from, -moves$weight, moves$to), , drop = FALSE]
    best <- moves[!duplicated(moves$from), , drop = FALSE]
    target <- match(best$to[match(nodes, best$from)], nodes)
  }
  heading <- atan2(node_data$y[target] - node_data$y,
                   node_data$x[target] - node_data$x)
  below <- identical(style, "outside") & !is.na(heading) & sin(heading) > 0.5
  data.frame(node = nodes, x = node_data$x, y = node_data$y, radius = radius,
             target = nodes[target], heading = heading, below = below,
             stringsAsFactors = FALSE)
}

# Where each node label is written: above its node (-1.1 label heights), just
# above its dot when points are sized by a value (the largest dot, 10 mm
# across, has a radius of about 1.2 label heights at label size 4.2, plus a
# small gap), and just above its circle -- or just below when an outside
# tip points up -- when circles are drawn.
.thg_label_positions <- function(node_data, value = NULL, geometry = NULL) {
  y <- node_data$y
  vjust <- rep(-1.1, nrow(node_data))
  if (!is.null(geometry)) {
    y <- node_data$y + ifelse(geometry$below, -1, 1) * geometry$radius
    vjust <- ifelse(geometry$below, 1.35, -0.35)
  } else if (!is.null(value)) {
    vjust <- -0.4 - 1.8 * sqrt(value / max(value))
  }
  data.frame(x = node_data$x, y = y, vjust = vjust, label = node_data$label,
             stringsAsFactors = FALSE)
}

# Everything drawn on the nodes, as a list of layers, scales and labels to
# add to the pebbles, in drawing order:
#   * `transitions`: curved arrows whose width follows the weight. Each end
#     is pulled back from its node (past its circle, when circles are drawn)
#     so the head stays visible; a node followed by itself gets a small loop
#     below it, or above it when its label sits below.
#   * the nodes: with `geometry` (node_sizes plus direction) a disc per node
#     and a triangle turned onto its heading -- the triangle given in its own
#     frame, +x towards the next node, in radii of its circle -- and a caption
#     saying what area and triangle show; with `value` alone, points whose
#     area follows the value (each label just above its own dot: the largest
#     dot, 10 mm across, has a radius of about 1.2 label heights at label
#     size 4.2, plus a small gap); otherwise plain points.
#   * the labels, bold with a halo as cograph's .add_text_with_halo():
#     stamped eight times in white at points on a small circle around it
#     (alpha 0.6), then once on top.
.thg_node_layers <- function(node_data, nodes, label_data, value = NULL,
                             geometry = NULL,
                             transitions = NULL, style = "inside",
                             node_fill = "#000000", arrow_fill = NULL,
                             labels = TRUE, label_size = 4.2, node_size = 2.5,
                             area_what = NULL, unit = NULL) {
  ink <- "#2B2B2B"
  arrows <- list()
  if (!is.null(transitions)) {
    radius <- geometry$radius %||% rep(0, length(nodes))
    flip <- geometry$below %||% rep(FALSE, length(nodes))
    span <- max(diff(range(node_data$y)), diff(range(node_data$x)), 1e-9)
    gap <- 0.035 * span
    i <- match(transitions$from, nodes)
    j <- match(transitions$to, nodes)
    d <- data.frame(transitions, x = node_data$x[i], y = node_data$y[i],
                    xend = node_data$x[j], yend = node_data$y[j],
                    r_from = radius[i], r_to = radius[j], flip = flip[i])
    head <- ggplot2::arrow(length = grid::unit(3, "mm"), type = "closed")
    links <- d[d$from != d$to, , drop = FALSE]
    len <- pmax(sqrt((links$xend - links$x)^2 + (links$yend - links$y)^2), 1e-9)
    pull_from <- pmin(pmax(gap, links$r_from + 0.01 * span) / len, 0.45)
    pull_to <- pmin(pmax(gap, links$r_to + 0.01 * span) / len, 0.45)
    links$x0 <- links$x + (links$xend - links$x) * pull_from
    links$y0 <- links$y + (links$yend - links$y) * pull_from
    links$x1 <- links$xend - (links$xend - links$x) * pull_to
    links$y1 <- links$yend - (links$yend - links$y) * pull_to
    loops <- d[d$from == d$to, , drop = FALSE]
    # swapping the ends turns the bulge of the curve to the other side
    side <- ifelse(loops$flip, -1, 1)
    loops$half <- side * pmax(gap, loops$r_from)
    loops$ly <- loops$y - side * (loops$r_from + 0.4 * gap)
    width_title <- if (is.null(unit)) "Transition weight" else
      paste(tools::toTitleCase(unit), "making the transition")
    arrows <- list(
      if (nrow(links)) ggplot2::geom_curve(
        data = links,
        mapping = ggplot2::aes(x = .data$x0, y = .data$y0, xend = .data$x1,
                               yend = .data$y1, linewidth = .data$weight),
        curvature = 0.15, arrow = head, colour = "#333333", alpha = 0.75,
        inherit.aes = FALSE
      ),
      if (nrow(loops)) ggplot2::geom_curve(
        data = loops,
        mapping = ggplot2::aes(x = .data$x - .data$half, y = .data$ly,
                               xend = .data$x + .data$half, yend = .data$ly,
                               linewidth = .data$weight),
        curvature = 1.6, arrow = head, colour = "#333333", alpha = 0.75,
        inherit.aes = FALSE
      ),
      ggplot2::scale_linewidth(range = c(0.3, 2.6), name = width_title,
                               labels = .thg_comma)
    )
  }
  shape_mapping <- ggplot2::aes(x = .data$x, y = .data$y, group = .data$id)
  caption <- list()
  if (!is.null(geometry)) {
    a <- seq(0, 2 * pi, length.out = 64L)
    discs <- do.call(rbind, lapply(seq_len(nrow(geometry)), function(i) {
      data.frame(id = i, x = geometry$x[i] + geometry$radius[i] * cos(a),
                 y = geometry$y[i] + geometry$radius[i] * sin(a))
    }))
    tip <- if (identical(style, "inside")) {
      list(x = c(0.95, -0.2, -0.2), y = c(0, 0.86, -0.86))
    } else {
      list(x = c(1.55, 0.55, 0.55), y = c(0, 0.62, -0.62))
    }
    triangles <- do.call(rbind, lapply(which(!is.na(geometry$heading)), function(i) {
      tx <- geometry$radius[i] * tip$x
      ty <- geometry$radius[i] * tip$y
      h <- geometry$heading[i]
      data.frame(id = i, x = geometry$x[i] + tx * cos(h) - ty * sin(h),
                 y = geometry$y[i] + tx * sin(h) + ty * cos(h))
    }))
    marks <- list(
      ggplot2::geom_polygon(data = discs, mapping = shape_mapping,
                            fill = node_fill, colour = NA, inherit.aes = FALSE),
      if (!is.null(triangles)) ggplot2::geom_polygon(
        data = triangles, mapping = shape_mapping,
        fill = arrow_fill %||% (if (identical(style, "inside")) "#1F3A44" else node_fill),
        colour = NA, inherit.aes = FALSE
      )
    )
    caption <- list(ggplot2::labs(caption = sprintf(
      "Circle area: %s. Triangle: points to the event that most often follows it.",
      area_what %||% "proportional to the event value")))
  } else if (!is.null(value)) {
    marks <- list(
      ggplot2::geom_point(
        data = data.frame(x = node_data$x, y = node_data$y, value = value),
        mapping = ggplot2::aes(x = .data$x, y = .data$y, size = .data$value),
        shape = 21, fill = ink, colour = "white", stroke = 0.5
      ),
      ggplot2::scale_size_area(
        max_size = 10, labels = .thg_comma,
        name = if (is.null(area_what)) {
          "Event value"
        } else {
          paste0(toupper(substr(area_what, 1L, 1L)), substring(area_what, 2L))
        }
      )
    )
  } else {
    marks <- list(ggplot2::geom_point(
      data = node_data, mapping = ggplot2::aes(x = .data$x, y = .data$y),
      shape = 21, fill = ink, colour = "white", stroke = 0.5, size = node_size
    ))
  }
  text <- if (isFALSE(labels)) {
    list()
  } else {
    halo_width <- 0.003 * max(diff(range(label_data$y)), 1e-9)
    label_mapping <- ggplot2::aes(x = .data$x, y = .data$y, label = .data$label,
                                  vjust = .data$vjust)
    stamps <- lapply(seq(0, 2 * pi, length.out = 9L)[-9L], function(a) {
      shifted <- label_data
      shifted$x <- shifted$x + halo_width * cos(a)
      shifted$y <- shifted$y + halo_width * sin(a)
      ggplot2::geom_text(data = shifted, mapping = label_mapping,
                         colour = "white", fontface = "bold",
                         size = label_size, alpha = 0.6, inherit.aes = FALSE)
    })
    c(stamps, list(ggplot2::geom_text(data = label_data, mapping = label_mapping,
                                      colour = ink, fontface = "bold",
                                      size = label_size, inherit.aes = FALSE)))
  }
  c(arrows, marks, text, caption)
}

# ---- title boxes and side-by-side pieces -------------------------------------
# Ported from the pipeline helper plot_blobs() (see the node overlays above):
# the same placement rule and constants, so the package draws its figures.

# One title box per hyperedge, in `precedence` order: the name (after
# `prefix`), then the value and `unit`, then the note, each on its own line.
# A box sits `gap` layout units beyond the pebble's furthest reach along the
# direction from the figure's centre (the middle of the pebbles' extent) to
# the pebble's centre, left-, centre- or right-justified by the side it is
# on; a one-member hyperedge, which has no pebble, is placed from its node.
# A box within a quarter of the figure's height across and 7% of it
# vertically of an earlier box moves to 8% below the lowest such box.
.thg_title_boxes <- function(hulls, pos, hyperedges, edge_names, precedence,
                             values, prefix, unit, notes, gap = 0.06) {
  centre <- c(mean(range(hulls$x)), mean(range(hulls$y)))
  boxes <- do.call(rbind, lapply(precedence, function(k) {
    outline <- hulls[hulls$hyperedge == edge_names[k], c("x", "y")]
    if (!nrow(outline)) outline <- pos$nodes[hyperedges[[k]], c("x", "y")]
    heading <- c(mean(outline$x), mean(outline$y)) - centre
    length_heading <- sqrt(sum(heading^2))
    # a pebble centred on the figure has no side of its own: write above it
    heading <- if (length_heading > 0) heading / length_heading else c(0, 1)
    reach <- max((outline$x - centre[1L]) * heading[1L] +
                   (outline$y - centre[2L]) * heading[2L])
    at <- centre + (reach + gap) * heading
    count <- if (is.null(values)) "" else
      sprintf("\n%s%s", .thg_comma(values[k], count = TRUE),
              if (is.null(unit)) "" else paste0(" ", unit))
    note <- if (is.null(notes) || is.na(notes[k])) "" else paste0("\n", notes[k])
    data.frame(x = at[1L], y = at[2L],
               hjust = 0.5 - 0.5 * sign(round(heading[1L], 1L)),
               label = paste0(prefix, edge_names[k], count, note),
               stringsAsFactors = FALSE)
  }))
  spread <- diff(range(hulls$y))
  # each box depends on where the earlier ones ended up, so the positions are
  # carried forward by Reduce()
  boxes$y <- Reduce(function(ys, i) {
    earlier <- seq_len(i - 1L)
    clash <- abs(boxes$x[earlier] - boxes$x[i]) < 0.25 * spread &
      abs(ys[earlier] - ys[i]) < 0.07 * spread
    if (any(clash)) ys[i] <- min(ys[earlier][clash]) - 0.08 * spread
    ys
  }, seq_len(nrow(boxes))[-1L], boxes$y)
  rownames(boxes) <- NULL
  boxes
}

# Node positions for `pieces = "row"`: every connected piece (hyperedges
# joined by a shared node) is laid out on its own with the bipartite layout,
# in its own unit frame, exactly as it would be drawn alone, and piece i is
# shifted 1.4 (i - 1) to the right. Pieces are ordered by their earliest
# hyperedge in `precedence`; a piece of one node sits at (0.5, 0.5) of its
# frame. NULL for a connected hypergraph, which needs no row.
.thg_row_layout <- function(hg, seed, center, padding, precedence) {
  inc <- as.matrix(hg$incidence != 0)
  piece <- .thg_components(crossprod(inc) * 1)
  if (max(piece) == 1L) return(NULL)
  earliest <- vapply(seq_len(max(piece)), function(k) {
    min(match(which(piece == k), precedence))
  }, integer(1L))
  do.call(rbind, lapply(seq_along(earliest), function(i) {
    edges <- which(piece == order(earliest)[i])
    rows <- which(rowSums(inc[, edges, drop = FALSE]) > 0)
    shift <- (i - 1L) * 1.4
    if (length(rows) == 1L) {
      return(data.frame(node = hg$nodes[rows], x = 0.5 + shift, y = 0.5,
                        stringsAsFactors = FALSE))
    }
    part <- list(incidence = hg$incidence[rows, edges, drop = FALSE],
                 nodes = hg$nodes[rows])
    at <- .thg_positions(part, "bipartite", seed, center, padding)$nodes
    data.frame(node = at$node, x = at$x + shift, y = at$y,
               stringsAsFactors = FALSE)
  }))
}
