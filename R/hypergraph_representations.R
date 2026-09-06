# The representation table of Coupette et al. (2024, Table 2): the same data
# read as a binary graph, a multi-graph, a binary hypergraph and a
# multi-hypergraph, with the edge count and degree statistics of each.

#' Graph and hypergraph representations of the same data
#'
#' Reports, side by side, the four ways Coupette et al. (2024, Table 2) model
#' one dataset: the hypergraph as it is (`mh`, repeated member sets counted
#' separately), the binary hypergraph of distinct member sets (`bh`), and the
#' graph the hypergraph reduces to, with (`mg`) or without (`bg`) edge
#' multiplicities. The graph is the clique expansion for collaboration data
#' (every pair of members of a tribunal is an edge) or the citation graph
#' for data with hyperedge sources (the citing decision is joined to every
#' decision in the block). The degree statistics show how much the choice
#' of representation changes what a node looks like.
#'
#' @param hg A static `net_hypergraph`, typically a snapshot of a
#'   [temporal_hypergraph()].
#' @param graph `"clique"` (default) or `"citation"`, the graph projection
#'   compared; see [hg_project()]. `"citation"` needs hyperedge sources.
#' @param edge_source Hyperedge sources for `graph = "citation"`, as in
#'   [hg_project()]; omitted when the hypergraph carries them.
#' @return A base data.frame with one row per representation, in the order
#'   `bg`, `mg`, `bh`, `mh`, and columns `representation`, `type` (`"graph"`
#'   or `"hypergraph"`), `n_nodes`, `n_edges`, `mean_degree`, `median_degree`.
#'   A directed citation graph adds `median_out_degree` and
#'   `median_in_degree`, and its `mean_degree` is the mean out-degree, which
#'   equals the mean in-degree. Hypergraph degree is the number of hyperedges
#'   containing a node; every node of the hypergraph counts, including nodes
#'   in no hyperedge.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' tribunals <- data.frame(
#'   member = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a1", "a2"),
#'   case = rep(c("A", "B", "C"), each = 3)
#' )
#' hg <- group_hypergraph(tribunals, "member", "case")
#' hg_representations(hg)
#' @export
hg_representations <- function(hg, graph = c("clique", "citation"),
                               edge_source = NULL) {
  .thg_check_hg(hg)
  graph <- match.arg(graph)
  b <- .thg_binary(hg$incidence)
  n <- hg$n_nodes

  multi_degree <- as.numeric(Matrix::rowSums(b))
  collapsed <- .thg_collapse_duplicate_edges(hg)
  binary_degree <- as.numeric(Matrix::rowSums(.thg_binary(collapsed$incidence)))
  hyper_rows <- data.frame(
    representation = c("bh", "mh"), type = "hypergraph", n_nodes = n,
    n_edges = c(collapsed$n_hyperedges, hg$n_hyperedges),
    mean_degree = c(mean(binary_degree), mean(multi_degree)),
    median_degree = c(stats::median(binary_degree), stats::median(multi_degree)),
    stringsAsFactors = FALSE
  )

  if (identical(graph, "clique")) {
    weighted <- hg_project(hg, method = "clique", weighted = FALSE,
                           what = "matrix")
    binary <- (weighted != 0) * 1
    graph_rows <- data.frame(
      representation = c("bg", "mg"), type = "graph", n_nodes = n,
      n_edges = c(sum(binary[upper.tri(binary)]), sum(weighted[upper.tri(weighted)])),
      mean_degree = c(mean(Matrix::rowSums(binary)), mean(Matrix::rowSums(weighted))),
      median_degree = c(stats::median(as.numeric(Matrix::rowSums(binary))),
                        stats::median(as.numeric(Matrix::rowSums(weighted)))),
      stringsAsFactors = FALSE
    )
  } else {
    directed <- lapply(c("collapse", "count"), function(dup) {
      hg_project(hg, method = "citation", what = "matrix",
                 duplicate_edges = dup, edge_source = edge_source,
                 directed = TRUE)
    })
    out_degree <- lapply(directed, function(w) as.numeric(Matrix::rowSums(w)))
    in_degree <- lapply(directed, function(w) as.numeric(Matrix::colSums(w)))
    graph_rows <- data.frame(
      representation = c("bg", "mg"), type = "graph",
      n_nodes = nrow(directed[[1L]]),
      n_edges = vapply(directed, function(w) sum(w), numeric(1L)),
      mean_degree = vapply(out_degree, mean, numeric(1L)),
      median_degree = vapply(seq_along(directed), function(i) {
        stats::median(out_degree[[i]] + in_degree[[i]])
      }, numeric(1L)),
      median_out_degree = vapply(out_degree, stats::median, numeric(1L)),
      median_in_degree = vapply(in_degree, stats::median, numeric(1L)),
      stringsAsFactors = FALSE
    )
    hyper_rows$median_out_degree <- NA_real_
    hyper_rows$median_in_degree <- NA_real_
  }
  out <- rbind(graph_rows, hyper_rows)
  out$n_edges <- as.integer(round(out$n_edges))
  rownames(out) <- NULL
  out
}

#' @rdname hg_representations
#' @export
hypergraph_representations <- hg_representations
