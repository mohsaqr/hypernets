# Growth of a temporal hypergraph over its event times, and the two small
# result classes shared by the descriptive layer: a time series table
# (`honets_series`) and an empirical distribution table
# (`honets_distribution`). Both are plain data.frames with a plot method.

#' Growth of a temporal hypergraph over time
#'
#' Counts, on a measurement grid, the nodes and hyperedges present in the
#' temporal hypergraph: the descriptive time series of Coupette et al.
#' (2024, Figures 4a and 5b). Hyperedges are counted twice, as they are
#' (`n_edges`, the multi-hypergraph) and as distinct member sets
#' (`n_edges_distinct`, the binary hypergraph), so the two representations of
#' the paper's Table 2 can be read against each other over time.
#'
#' `mode = "active"` counts what each window measures -- an interval
#' hyperedge active in it, a contact hyperedge occurring in it -- and adds
#' the `_cumulative` columns, the hyperedges begun by the end of the window,
#' the static aggregate that a point-aggregation model would report.
#' `mode = "cumulative"` reports only that growing view, and then counts a
#' node from its own entry time when the constructor received a node table
#' with a `start` column, otherwise from the first hyperedge that contains
#' it.
#'
#' @inheritParams hypergraph_snapshots
#' @param x A [temporal_hypergraph()].
#' @param components Also report the connectivity of the measured hypergraph
#'   at each time: the number of connected components, the share of its nodes
#'   in the largest, and its diameter (longest shortest path between two
#'   nodes sharing a chain of hyperedges)? Default `FALSE`; the computation
#'   is a breadth-first search per node and time.
#' @return A base data.frame of class `honets_series`, one row per window,
#'   with `time` (the window start on the hypergraph's clock), `n_nodes`,
#'   `n_edges`, `n_edges_distinct`, `n_memberships`, and for
#'   `mode = "active"` also `n_nodes_cumulative`, `n_edges_cumulative` and
#'   `n_edges_distinct_cumulative`. With `components = TRUE`, further
#'   `n_components`, `largest_component` (share of active nodes) and
#'   `diameter`. `plot()` draws each count against time, on the calendar
#'   when the hypergraph has one.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' seats <- data.frame(
#'   case = rep(c("A", "B", "C"), each = 3),
#'   arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
#'   constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3)
#' )
#' thg <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
#'                            start = "constituted", end = "concluded")
#' growth <- hg_growth(thg)
#' growth
#' plot(growth)
#' # a regular grid: every two steps, each window two steps wide
#' hg_growth(thg, step = 2)
#' @export
hg_growth <- function(x, start = NULL, end = NULL, step = NULL, window = NULL,
                      at = NULL, mode = c("active", "cumulative"),
                      components = FALSE) {
  .thg_check_temporal(x)
  mode <- .thg_check_mode(mode, "hg_growth")
  stopifnot(
    "`components` must be TRUE or FALSE" =
      is.logical(components) && length(components) == 1L && !is.na(components)
  )
  grid <- .thg_grid(x, start, end, step, window, at)
  at <- grid$times
  ed <- x$edge_data
  mem <- x$memberships
  edge_index <- match(mem$edge, ed$edge)
  # Signature of each hyperedge's member set, for the distinct-edge count.
  members_by_edge <- split(mem$member, factor(edge_index, levels = seq_len(nrow(ed))))
  signature <- vapply(members_by_edge, function(v) paste(sort(unique(v)), collapse = "\r"),
                      character(1L))
  size <- lengths(members_by_edge)
  # Node entry: the node table's own start when given, else first membership.
  first_membership <- tapply(ed$start[edge_index], mem$member, min)
  node_start <- first_membership[x$nodes]
  if (!is.null(x$node_data) && any(!is.na(x$node_data$start))) {
    given <- !is.na(x$node_data$start)
    node_start[given] <- x$node_data$start[given]
  }
  count <- function(keep) {
    c(n_nodes = length(unique(mem$member[keep[edge_index]])),
      n_edges = sum(keep),
      n_edges_distinct = length(unique(signature[keep])),
      n_memberships = sum(size[keep]))
  }
  window_end <- function(t) t + grid$window

  if (identical(mode, "cumulative")) {
    start_first_signature <- tapply(ed$start, signature, min)
    begun <- t(vapply(at, function(t) {
      count(.thg_edges_in_window(x, t, grid$window, "cumulative", grid$closed))
    }, numeric(4L)))
    counts <- data.frame(
      time = at,
      n_nodes = vapply(at, function(t) sum(node_start <= window_end(t), na.rm = TRUE),
                       numeric(1L)),
      n_edges = begun[, "n_edges"],
      n_edges_distinct = begun[, "n_edges_distinct"],
      n_memberships = begun[, "n_memberships"]
    )
  } else {
    active <- t(vapply(at, function(t) {
      count(.thg_edges_in_window(x, t, grid$window, "active", grid$closed))
    }, numeric(4L)))
    begun <- t(vapply(at, function(t) {
      count(.thg_edges_in_window(x, t, grid$window, "cumulative", grid$closed))
    }, numeric(4L)))
    counts <- data.frame(
      time = at,
      n_nodes = active[, "n_nodes"], n_edges = active[, "n_edges"],
      n_edges_distinct = active[, "n_edges_distinct"],
      n_memberships = active[, "n_memberships"],
      n_nodes_cumulative = begun[, "n_nodes"],
      n_edges_cumulative = begun[, "n_edges"],
      n_edges_distinct_cumulative = begun[, "n_edges_distinct"]
    )
  }
  integer_columns <- setdiff(names(counts), "time")
  counts[integer_columns] <- lapply(counts[integer_columns], as.integer)

  if (components) {
    connectivity <- t(vapply(at, function(t) {
      snap <- .thg_snapshot_at(x, t, grid$window, mode, TRUE, grid$closed)
      .thg_connectivity(snap)
    }, numeric(3L)))
    counts$n_components <- as.integer(connectivity[, 1L])
    counts$largest_component <- connectivity[, 2L]
    counts$diameter <- as.integer(connectivity[, 3L])
  }
  rownames(counts) <- NULL
  .thg_series(counts, x)
}

# A series or distribution table remembers the clock it was measured on, so
# its plot can label the axis with dates.
.thg_series <- function(out, x, class = "honets_series") {
  class(out) <- c(class, "data.frame")
  attr(out, "time_unit") <- x$time_unit
  attr(out, "origin") <- x$origin
  out
}

#' @rdname hg_growth
#' @export
hypergraph_growth <- hg_growth

# Connected components of a static hypergraph through shared hyperedges:
# their number, the share of nodes in the largest, and its diameter.
.thg_connectivity <- function(hg) {
  n <- hg$n_nodes
  if (n == 0L) return(c(n_components = 0, largest = NA_real_, diameter = NA_real_))
  b <- .thg_binary(hg$incidence)
  adjacency <- Matrix::tcrossprod(b) > 0
  # Frontier expansion is many small matrix-vector products; below a few
  # thousand nodes a dense matrix is far cheaper than a sparse one.
  if (n <= 2000L) adjacency <- as.matrix(adjacency)
  membership <- .thg_components(adjacency)
  sizes <- tabulate(membership)
  largest <- which.max(sizes)
  inside <- which(membership == largest)
  diameter <- if (length(inside) <= 1L) 0 else
    .thg_diameter(adjacency[inside, inside, drop = FALSE])
  c(n_components = length(sizes), largest = max(sizes) / n, diameter = diameter)
}

# One row per connected component of a static hypergraph.
.thg_component_table <- function(hg) {
  empty <- data.frame(component = integer(), n_nodes = integer(),
                      n_edges = integer(), share = numeric(),
                      diameter = integer())
  if (hg$n_nodes == 0L) return(empty)
  b <- .thg_binary(hg$incidence)
  adjacency <- Matrix::tcrossprod(b) > 0
  if (hg$n_nodes <= 2000L) adjacency <- as.matrix(adjacency)
  membership <- .thg_components(adjacency)
  # Every member of a hyperedge lies in one component, so its first member
  # names the hyperedge's component.
  first_member <- vapply(hg$hyperedges, function(v) if (length(v)) v[[1L]] else NA_integer_,
                         integer(1L))
  edge_component <- membership[first_member]
  sizes <- tabulate(membership)
  ord <- order(-sizes, seq_along(sizes))
  rows <- lapply(ord, function(k) {
    inside <- which(membership == k)
    diameter <- if (length(inside) <= 1L) 0L else
      as.integer(.thg_diameter(adjacency[inside, inside, drop = FALSE]))
    data.frame(component = k, n_nodes = sizes[k],
               n_edges = sum(edge_component == k, na.rm = TRUE),
               share = sizes[k] / hg$n_nodes, diameter = diameter)
  })
  out <- do.call(rbind, rows)
  out$component <- match(out$component, ord)
  rownames(out) <- NULL
  out
}

# Component labels by repeated frontier expansion over a logical adjacency.
.thg_components <- function(adjacency) {
  n <- nrow(adjacency)
  label <- integer(n)
  current <- 0L
  # Each pass seeds one unlabelled node and floods its component; the number
  # of passes is the number of components, sequential by nature.
  while (any(label == 0L)) {
    current <- current + 1L
    reached <- logical(n)
    reached[which(label == 0L)[1L]] <- TRUE
    repeat {
      nxt <- reached | (as.numeric(adjacency %*% reached) > 0)
      if (identical(nxt, reached)) break
      reached <- nxt
    }
    label[reached] <- current
  }
  label
}

# Diameter of a connected adjacency: the smallest k for which every node
# reaches every other within k steps, found by multiplying the reachability
# matrix by (A + I) until it is full. Equivalent to the largest breadth-first
# eccentricity, without a per-node loop.
.thg_diameter <- function(adjacency) {
  n <- nrow(adjacency)
  if (n <= 1L) return(0L)
  step <- adjacency * 1
  diag(step) <- 1
  reach <- step
  k <- 1L
  # bounded by n - 1 multiplications; a disconnected input stops when the
  # reachability matrix no longer changes
  while (any(reach == 0)) {
    nxt <- (reach %*% step > 0) * 1
    if (isTRUE(all.equal(nxt, reach, check.attributes = FALSE))) return(NA_integer_)
    reach <- nxt
    k <- k + 1L
    if (k > n) break
  }
  k
}

#' @rdname hg_growth
#' @param x A `honets_series` table (for `plot`).
#' @param columns Which numeric columns to draw; default all but `time`.
#' @param facets Draw one panel per column (default `TRUE`, because node
#'   and hyperedge counts live on different scales) or all columns in one
#'   panel with a legend.
#' @param ... Unused; for S3 consistency.
#' @return For `plot`, a ggplot object.
#' @export
plot.honets_series <- function(x, columns = NULL, facets = TRUE, ...) {
  d <- as.data.frame(x)
  numeric_columns <- setdiff(names(d)[vapply(d, is.numeric, logical(1L))], "time")
  columns <- columns %||% numeric_columns
  unknown <- setdiff(columns, numeric_columns)
  if (length(unknown)) {
    .thg_bad_input(sprintf("unknown series column(s): %s",
                           paste(unknown, collapse = ", ")))
  }
  origin <- attr(x, "origin")
  time_unit <- attr(x, "time_unit")
  axis_time <- if (!is.null(origin) && inherits(origin, "POSIXt")) {
    .thg_calendar(d$time, origin, time_unit)
  } else {
    d$time
  }
  long <- do.call(rbind, lapply(columns, function(column) {
    data.frame(time = axis_time, measure = column, value = as.numeric(d[[column]]),
               stringsAsFactors = FALSE)
  }))
  long$measure <- factor(long$measure, levels = columns)
  palette <- rep_len(.thg_okabe_ito, length(columns))
  p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$time, y = .data$value,
                                          colour = .data$measure,
                                          linetype = .data$measure)) +
    ggplot2::geom_step(linewidth = 0.7) +
    ggplot2::scale_colour_manual(values = palette, name = NULL) +
    ggplot2::scale_linetype_manual(values = rep_len(c("solid", "dashed", "dotted",
                                                      "dotdash", "longdash",
                                                      "twodash"),
                                                    length(columns)),
                                   name = NULL) +
    ggplot2::labs(x = if (!is.null(origin) && inherits(origin, "POSIXt")) "date" else
                    if (is.null(time_unit) || identical(time_unit, "step")) "time" else
                      time_unit,
                  y = NULL) +
    ggplot2::theme_minimal(base_size = 12)
  if (isTRUE(facets)) {
    p + ggplot2::facet_wrap(~ measure, scales = "free_y") +
      ggplot2::theme(legend.position = "none")
  } else {
    p + ggplot2::theme(legend.position = "top")
  }
}

#' Empirical distribution of a hypergraph measure
#'
#' The table returned by [hg_edges()] and [hg_measures()] with
#' `what = "distribution"`: one row per distinct observed value with its
#' count, proportion and complementary cumulative share, the CCDF that the
#' descriptive figures of Coupette et al. (2024) draw on a logarithmic axis.
#' A grouping column (`time`, `s`) present in the table gives one curve per
#' group.
#'
#' @param x A `honets_distribution` table.
#' @param log Draw the CCDF on a logarithmic y axis (default `TRUE`).
#' @param ... Unused; for S3 consistency.
#' @return A ggplot object.
#' @keywords internal
#' @export
plot.honets_distribution <- function(x, log = TRUE, ...) {
  d <- as.data.frame(x)
  group <- intersect(c("time", "s", "series"), names(d))
  d$series <- if (length(group)) {
    factor(as.character(d[[group[[1L]]]]), levels = unique(as.character(d[[group[[1L]]]])))
  } else {
    factor("all")
  }
  n_series <- nlevels(d$series)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$value, y = .data$ccdf,
                                       colour = .data$series,
                                       linetype = .data$series)) +
    ggplot2::geom_step(direction = "vh", linewidth = 0.7) +
    ggplot2::geom_point(size = 1.2) +
    ggplot2::scale_colour_manual(values = rep_len(.thg_okabe_ito, n_series),
                                 name = if (length(group)) group[[1L]] else NULL) +
    ggplot2::scale_linetype_manual(values = rep_len(c("solid", "dashed", "dotted",
                                                      "dotdash", "longdash",
                                                      "twodash"), n_series),
                                   name = if (length(group)) group[[1L]] else NULL) +
    ggplot2::labs(x = attr(x, "measure") %||% "value", y = "P(X >= value)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = if (n_series > 1L) "right" else "none")
  if (isTRUE(log)) p <- p + ggplot2::scale_y_log10()
  p
}
