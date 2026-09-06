# ---- Hypergraph -> pairwise network projection (HON-7) -------------------
# Standard clique expansion: for each hyperedge, every pair of members
# contributes a unit (or weighted) tally to the resulting pairwise network.
# Closes the I/O cycle with build_hypergraph() / group_hypergraph().

#' Clique expansion of a hypergraph
#'
#' Projects a [net_hypergraph][build_hypergraph] to a standard pairwise
#' netobject (the *clique expansion* - also called the
#' "downgrade" of a hypergraph to a dyadic graph). Each hyperedge of size
#' k contributes 1 (or its weight) to every pair of its members. The
#' resulting edge weight `W[i, j]` equals the number of hyperedges
#' containing both `i` and `j` (binary incidence) or the sum of incidence
#' products (weighted incidence).
#'
#' @param hg A `net_hypergraph` object as returned by [build_hypergraph()]
#'   or [group_hypergraph()].
#' @param weighted Logical. If `TRUE` (default), use the hypergraph's
#'   incidence values directly (so weighted hypergraphs from
#'   [group_hypergraph()] produce weighted projections). If `FALSE`,
#'   binarise the incidence first so `W[i, j]` is just the count of
#'   shared hyperedges.
#'
#' @return A `netobject` (also `cograph_network`) with `method =
#'   "clique_expansion"`, undirected, with weighted symmetric adjacency
#'   `W = incidence %*% t(incidence)` and zero diagonal.
#'
#' @details
#' The clique expansion is the standard "loss-y but lossless-on-pairwise"
#' projection: it preserves *which pairs co-occurred* and *how often* but
#' discards the higher-order grouping. Comparing `clique_expansion(hg)` to
#' a directly-estimated pairwise network (e.g. via \code{Nestimate::cooccurrence()} on
#' the same data) quantifies how much information was carried by the
#' hyperedge structure.
#'
#' Computed in one BLAS call via `tcrossprod(incidence)`; runs in
#' `O(n_nodes^2 * n_hyperedges)` time, fast for typical sizes.
#'
#' Closes the I/O cycle: event data -> [group_hypergraph()] ->
#' `clique_expansion()` -> any function that accepts a `netobject`
#' (centrality, bootstrap, clustering, plotting via cograph).
#'
#' @seealso [build_hypergraph()], [group_hypergraph()], \code{Nestimate::build_network()}.
#'
#' @examples
#' df <- data.frame(
#'   person  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
#' )
#' hg  <- group_hypergraph(df, actor = "person", group = "session")
#' net <- clique_expansion(hg)
#' net
#'
#' @references
#' Tian, Y., & Zafarani, R. (2024). Higher-order network analysis methods.
#' \emph{SIGKDD Explorations} 26(1), Section 5.1.5.
#'
#' @note (experimental) Validated against `tcrossprod(incidence)` with zero
#'   diagonal. No external R package exposes clique expansion as a primitive;
#'   the implementation is a direct one-line restatement of the definition.
#'
#' @export
clique_expansion <- function(hg, weighted = TRUE) {
  stopifnot(
    inherits(hg, "net_hypergraph"),
    is.logical(weighted), length(weighted) == 1L
  )

  if (hg$n_hyperedges == 0L) {
    weights <- matrix(0, hg$n_nodes, hg$n_nodes,
                      dimnames = list(hg$nodes, hg$nodes))
  } else {
    inc <- hg$incidence
    if (!weighted) {
      inc <- (inc > 0) * 1L
    }
    storage.mode(inc) <- "double"
    weights <- tcrossprod(inc)
    diag(weights) <- 0
    rownames(weights) <- colnames(weights) <- hg$nodes
  }

  net <- .wrap_netobject(weights,
                         method   = "clique_expansion",
                         directed = FALSE)
  net$params <- list(
    source                       = "clique_expansion",
    weighted                     = weighted,
    n_hyperedges                 = hg$n_hyperedges,
    hypergraph_size_distribution = hg$size_distribution
  )
  net
}
