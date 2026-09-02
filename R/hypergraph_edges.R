# The hyperedge-level tier. Everything else in the package is keyed on
# vertices; these are the measures that describe the hyperedges themselves.

#' Hyperedge-level measures, as a tidy table
#'
#' One row per hyperedge. `hg_measures(what = "edges")` reports cardinality
#' alone; this adds the incident weight, how many other hyperedges a
#' hyperedge meets, and how large its vertex neighbourhood is — the
#' hyperedge-side counterparts of degree and neighbourhood size, and the
#' descriptive layer used to characterise higher-order structure (Coupette
#' et al. 2024).
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any honets
#'   `net_hypergraph`.
#' @param what `"edges"` (default) for one row per hyperedge, or
#'   `"distribution"` for the empirical distribution of `measure` across
#'   hyperedges.
#' @param measure Which column `what = "distribution"` summarises: `"size"`
#'   (default), `"weight"`, `"n_incident_edges"` or `"n_neighbors"`.
#' @param s Minimum number of shared vertices for another hyperedge to count
#'   as incident. The default `1` is ordinary incidence; larger values match
#'   the thresholds used by [hg_edge_centrality()].
#' @return With `what = "edges"`, a base data.frame with one row per
#'   hyperedge and columns `edge` (name), `size` (integer, vertices it
#'   contains), `weight` (numeric, its incidence weights summed),
#'   `n_incident_edges` (integer, other hyperedges sharing at least one
#'   vertex) and `n_neighbors` (integer, vertices adjacent to a member
#'   without being one). With `what = "distribution"`, one row per distinct
#'   observed value of `measure`, ascending, with columns `value`, `n`,
#'   `proportion` and `ccdf` — the complementary cumulative distribution
#'   \eqn{P(X \ge value)}, so the first row's `ccdf` is always 1.
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
#' @export
hg_edges <- function(hg, what = c("edges", "distribution"),
                     measure = c("size", "weight", "n_incident_edges",
                                 "n_neighbors"), s = 1L) {
  .thg_check_hg(hg)
  what <- match.arg(what)
  measure <- match.arg(measure)
  if (length(s) != 1L || !is.numeric(s) || !is.finite(s) || s < 1 ||
      abs(s - round(s)) > sqrt(.Machine$double.eps)) {
    .thg_bad_input("`s` must be one positive whole number")
  }
  s <- as.integer(s)
  incidence <- hg$incidence
  b <- .thg_binary(incidence)
  size <- as.integer(Matrix::colSums(b))

  overlap <- crossprod(b)
  diag(overlap) <- 0
  n_incident <- as.integer(Matrix::colSums(overlap >= s))

  # A vertex is a neighbour of hyperedge e when it shares some other
  # hyperedge with a member of e without being a member itself. Members of a
  # hyperedge of size >= 2 are adjacent to each other, so they always appear
  # in the reach and are subtracted back out.
  adjacency <- tcrossprod(b)
  diag(adjacency) <- 0
  reach <- as.integer(Matrix::colSums((adjacency %*% b) > 0))
  n_neighbors <- reach - ifelse(size >= 2L, size, 0L)

  edges <- data.frame(
    edge = colnames(incidence),
    size = size,
    weight = as.numeric(Matrix::colSums(incidence)),
    n_incident_edges = n_incident,
    n_neighbors = as.integer(n_neighbors),
    row.names = NULL
  )
  if (identical(what, "edges")) return(edges)
  .thg_distribution(edges[[measure]])
}

#' @rdname hg_edges
#' @export
hypergraph_edges <- hg_edges

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
