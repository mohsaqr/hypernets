# Dual hypergraph: swap the vertex and hyperedge roles. For a bag text
# hypergraph this is exactly the nodes = "doc"/"word" flip, available on a
# fitted object without re-tokenizing.

#' Dual of a hypergraph
#'
#' Returns the dual hypergraph: every hyperedge becomes a vertex and every
#' vertex becomes a hyperedge, with the transposed weighted incidence. For a
#' bag-construction [text_hypergraph()], the dual is identical to rebuilding
#' with the opposite `nodes` orientation (tested), so document-level and
#' word-level analyses can share one constructed object. Duals of windowed
#' and kNN hypergraphs are returned as plain `net_hypergraph` objects (their
#' corpus bookkeeping does not transpose meaningfully).
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any honets
#'   `net_hypergraph`.
#' @return A hypergraph whose incidence is the transpose of `hg`'s: a
#'   `text_hypergraph` with flipped `nodes` for bag constructions, otherwise
#'   a `net_hypergraph`. Accepted by all `hg_*` verbs.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
#' dual <- dual_hypergraph(hg)
#' dual
#' hg_measures(dual, what = "edges")
#' @export
dual_hypergraph <- function(hg) {
  .thg_check_hg(hg)
  if (.thg_is_sparse(hg)) {
    triplet <- methods::as(hg$incidence, "TsparseMatrix")
    long <- data.frame(
      vertex = colnames(hg$incidence)[triplet@j + 1L],
      edge = rownames(hg$incidence)[triplet@i + 1L],
      w = triplet@x
    )
    return(.thg_sparse_bipartite(long, actor = "vertex", group = "edge",
                                 weight = "w"))
  }
  nz <- which(hg$incidence != 0, arr.ind = TRUE)
  long <- data.frame(
    vertex = colnames(hg$incidence)[nz[, "col"]],
    edge = rownames(hg$incidence)[nz[, "row"]],
    w = as.numeric(hg$incidence[nz])
  )
  dual <- group_hypergraph(long, actor = "vertex",
                                      group = "edge", weight = "w")
  if (inherits(hg, "text_hypergraph") &&
      identical(hg$text$construction, "bag")) {
    dual$text <- hg$text
    dual$text$nodes <- if (identical(hg$text$nodes, "doc")) "word" else "doc"
    class(dual) <- c("text_hypergraph", class(dual))
  }
  dual
}
