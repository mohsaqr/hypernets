#' Low-Dimensional Hypergraph Embedding
#'
#' Returns the node coordinates computed by honets' existing hypergraph
#' spectral or symmetric-NMF engine without exposing the incidental k-means
#' assignments produced by [hg_cluster()]. This is the direct embedding verb;
#' `hg_embed()` and `hypergraph_embed()` are identical names for it.
#'
#' @param hg Any honets `net_hypergraph`.
#' @param dimensions Number of embedding coordinates, between 2 and
#'   `n_nodes - 1`.
#' @param type Laplacian type: `"zhou"` or `"random_walk"`.
#' @param method `"spectral"` (default) or Hayashi et al.'s
#'   `"symnmf"`. SymNMF requires a dense incidence matrix.
#' @param seed Optional initialization seed.
#' @param nstart Number of k-means starts for the shared engine; for spectral
#'   embeddings this affects only the discarded cluster assignment, not the
#'   returned coordinates. For SymNMF it is the number of factor restarts.
#' @param max_iter,tol SymNMF iteration controls.
#' @return A data frame with `node`, stationary mass `pi`, and
#'   `dim1` through `dim<dimensions>`.
#' @examples
#' h <- group_hypergraph(
#'   data.frame(node = c("a", "b", "c", "b", "c", "d"),
#'              edge = rep(c("e1", "e2"), each = 3)),
#'   "node", "edge"
#' )
#' hg_embed(h, dimensions = 2)
#' @export
hg_embed <- function(hg, dimensions = 2L,
                     type = c("zhou", "random_walk"),
                     method = c("spectral", "symnmf"), seed = NULL,
                     nstart = 25L, max_iter = 500L, tol = 1e-6) {
  type <- match.arg(type)
  method <- match.arg(method)
  if (!is.numeric(dimensions) || length(dimensions) != 1L ||
      !is.finite(dimensions) || dimensions != as.integer(dimensions)) {
    stop("`dimensions` must be one whole number.", call. = FALSE)
  }
  fit <- hg_cluster(
    hg, k = as.integer(dimensions), type = type, seed = seed,
    nstart = nstart, what = "embedding", algorithm = method,
    max_iter = max_iter, tol = tol
  )
  coordinates <- grep("^dim[0-9]+$", names(fit), value = TRUE)
  out <- fit[, c("node", "pi", coordinates), drop = FALSE]
  attr(out, "method") <- method
  attr(out, "laplacian") <- type
  out
}

#' @rdname hg_embed
#' @export
hypergraph_embed <- hg_embed
