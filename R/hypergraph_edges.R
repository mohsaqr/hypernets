# The hyperedge-level tier. Everything else in the package is keyed on
# vertices; these are the measures that describe the hyperedges themselves.

#' Hyperedge-level measures, as a tidy table
#'
#' One row per hyperedge. `hg_measures(what = "edges")` reports cardinality
#' alone; this adds the incident weight, how many other hyperedges a
#' hyperedge meets, and how large its vertex neighbourhood is — the
#' hyperedge-side counterparts of degree and neighbourhood size, and the
#' descriptive layer used to characterise higher-order structure (Coupette
#' et al. 2024). A temporal hypergraph is evaluated snapshot by snapshot,
#' with a leading `time` column, so the distribution of hyperedge sizes at
#' several dates (the paper's Figure 4b) or the neighbourhood size of a
#' tribunal over time (Figure 5c) come from one call.
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], any hypernets
#'   `net_hypergraph`, or a [temporal_hypergraph()].
#' @param what `"edges"` (default) for one row per hyperedge,
#'   `"distribution"` for the empirical distribution of `measure` across
#'   hyperedges, or `"summary"` for its mean, standard deviation and
#'   quartiles.
#' @param measure Which column `what = "distribution"` and
#'   `what = "summary"` describe: `"size"` (default), `"weight"`,
#'   `"n_incident_edges"` or `"n_neighbors"`.
#' @param s Minimum number of shared vertices for another hyperedge to count
#'   as incident. The default `1` is ordinary incidence; larger values match
#'   the thresholds used by [hg_edge_centrality()]. Several values give one
#'   block of rows each, with an `s` column.
#' @param start,end,step,window,at Measurement grid passed to
#'   [hypergraph_snapshots()] when `hg` is temporal: the bounds of the
#'   period, how often to look, how much time each look covers, or the
#'   instants themselves.
#' @param snapshot_mode,multiedges Snapshot `mode` (`"active"` or
#'   `"cumulative"`) and multi-edge handling passed to
#'   [hypergraph_snapshots()] when `hg` is temporal.
#' @return With `what = "edges"`, a base data.frame with one row per
#'   hyperedge and columns `edge` (name), `size` (integer, vertices it
#'   contains), `weight` (numeric, its incidence weights summed),
#'   `n_incident_edges` (integer, other hyperedges sharing at least `s`
#'   vertices) and `n_neighbors` (integer, vertices adjacent to a member
#'   without being one). With `what = "distribution"`, a `hypernets_distribution`
#'   table with one row per distinct observed value of `measure`, ascending,
#'   with columns `value`, `n`, `proportion` and `ccdf` — the complementary
#'   cumulative distribution \eqn{P(X \ge value)}, so the first row's `ccdf`
#'   is always 1; its `plot()` draws the CCDF. With `what = "summary"`, one
#'   row with `n_edges`, `mean`, `sd`, `min`, `q25`, `median`, `q75` and
#'   `max`. Several `s` values add an `s` column; temporal input adds a
#'   leading `time` column, and the summary is then a `hypernets_series` whose
#'   `plot()` draws each statistic against time.
#' @references
#' Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
#' *Philosophical Transactions of the Royal Society A*, 382(2270), 20230141.
#' \doi{10.1098/rsta.2023.0141}
#' @seealso [hg_measures()] for vertex-level measures, [hg_line_graph()] for
#'   the graph whose degrees `n_incident_edges` reports.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
#'                         c = "stars and salt"))
#' hg_edges(hg)
#' hg_edges(hg, what = "distribution")
#' hg_edges(hg, what = "summary", measure = "n_neighbors")
#' @export
hg_edges <- function(hg, what = c("edges", "distribution", "summary"),
                     measure = c("size", "weight", "n_incident_edges",
                                 "n_neighbors"), s = 1L, start = NULL,
                     end = NULL, step = NULL, window = NULL, at = NULL,
                     snapshot_mode = c("active", "cumulative"),
                     multiedges = TRUE) {
  what <- match.arg(what)
  measure <- match.arg(measure)
  snapshot_mode <- .thg_check_mode(snapshot_mode, "hg_edges", "snapshot_mode")
  if (!is.numeric(s) || length(s) < 1L || any(!is.finite(s)) || any(s < 1) ||
      any(abs(s - round(s)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`s` must contain positive whole numbers")
  }
  s <- as.integer(round(s))

  if (inherits(hg, "net_temporal_hypergraph")) {
    snaps <- hypergraph_snapshots(hg, start = start, end = end, step = step,
                                  window = window, mode = snapshot_mode, at = at,
                                  multiedges = multiedges)
    rows <- lapply(seq_along(snaps), function(i) {
      ans <- hg_edges(snaps[[i]], what = what, measure = measure, s = s)
      if (nrow(ans) == 0L) return(NULL)
      time <- snaps[[i]]$params$at %||% names(snaps)[i]
      data.frame(time = rep(time, nrow(ans)), as.data.frame(ans),
                 row.names = NULL, stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    out <- .thg_edges_class(out, what, measure)
    attr(out, "time_unit") <- hg$time_unit
    attr(out, "origin") <- hg$origin
    return(out)
  }

  .thg_check_hg(hg)
  # `n_neighbors` costs two products over the node-by-node adjacency, which
  # densifies on a hub-heavy network; compute it only when it is asked for.
  blocks <- .thg_edge_tables(hg, s, neighbors = identical(what, "edges") ||
                                                identical(measure, "n_neighbors"))
  if (length(s) > 1L) {
    blocks <- lapply(seq_along(s), function(i) {
      data.frame(s = rep(s[i], nrow(blocks[[i]])), blocks[[i]],
                 stringsAsFactors = FALSE)
    })
  }
  edges <- do.call(rbind, blocks)
  rownames(edges) <- NULL
  if (identical(what, "edges")) return(edges)
  groups <- if (length(s) > 1L) split(edges, edges$s) else list(edges)
  out <- do.call(rbind, lapply(seq_along(groups), function(i) {
    values <- groups[[i]][[measure]]
    tab <- if (identical(what, "distribution")) .thg_distribution(values) else
      .thg_summary_row(values)
    if (length(s) > 1L) tab <- data.frame(s = rep(s[i], nrow(tab)), tab)
    tab
  }))
  rownames(out) <- NULL
  .thg_edges_class(out, what, measure)
}

#' @rdname hg_edges
#' @export
hypergraph_edges <- hg_edges

# One row per hyperedge for one intersection threshold `s`.
# One edge table per value of `s`. Everything except `n_incident_edges` is
# invariant in `s`, so the products are taken once and only the threshold is
# swept -- an `s` sweep used to recompute them all per value.
.thg_edge_tables <- function(hg, s, neighbors = TRUE) {
  incidence <- hg$incidence
  b <- .thg_binary(incidence)
  size <- as.integer(Matrix::colSums(b))
  overlap <- crossprod(b)
  diag(overlap) <- 0
  weight <- as.numeric(Matrix::colSums(incidence))

  n_neighbors <- if (isTRUE(neighbors)) {
    # A vertex is a neighbour of hyperedge e when it shares some other
    # hyperedge with a member of e without being a member itself. Members of a
    # hyperedge of size >= 2 are adjacent to each other, so they always appear
    # in the reach and are subtracted back out.
    adjacency <- tcrossprod(b)
    diag(adjacency) <- 0
    reach <- as.integer(Matrix::colSums((adjacency %*% b) > 0))
    as.integer(reach - ifelse(size >= 2L, size, 0L))
  } else {
    rep(NA_integer_, length(size))
  }

  lapply(s, function(ss) {
    data.frame(
      edge = colnames(incidence),
      size = size,
      weight = weight,
      n_incident_edges = as.integer(Matrix::colSums(overlap >= ss)),
      n_neighbors = n_neighbors,
      row.names = NULL
    )
  })
}

.thg_edge_table <- function(hg, s) {
  .thg_edge_tables(hg, s, neighbors = TRUE)[[1L]]
}

.thg_edges_class <- function(out, what, measure) {
  if (identical(what, "distribution")) {
    class(out) <- c("hypernets_distribution", "data.frame")
    attr(out, "measure") <- measure
  } else if (identical(what, "summary") && "time" %in% names(out)) {
    class(out) <- c("hypernets_series", "data.frame")
  }
  out
}

# Empirical distribution of a numeric vector: one row per distinct value,
# ascending, with the complementary cumulative P(X >= value).
.thg_distribution <- function(x) {
  counts <- table(x)
  value <- as.numeric(names(counts))
  n <- as.integer(counts)
  data.frame(
    value = value,
    n = n,
    proportion = n / length(x),
    ccdf = rev(cumsum(rev(n))) / length(x),
    row.names = NULL
  )
}

# Mean, spread and quartiles of a numeric vector as one row.
.thg_summary_row <- function(x) {
  q <- stats::quantile(x, c(0.25, 0.5, 0.75), names = FALSE)
  data.frame(
    n_edges = length(x), mean = mean(x), sd = stats::sd(x),
    min = min(x), q25 = q[1L], median = q[2L], q75 = q[3L], max = max(x),
    row.names = NULL
  )
}
