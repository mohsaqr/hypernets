# kNN embedding hypergraph: each item and its k nearest neighbors in an
# embedding space form one hyperedge, weighted by cosine similarity. The
# tier-3 -> tier-2 bridge: contrastive embeddings feed the explicit spectral
# machinery instead of a black-box UMAP+HDBSCAN pipeline.

#' Build a k-nearest-neighbor hypergraph from embeddings
#'
#' Every row of `embeddings` becomes one hyperedge containing the row itself
#' and its `k` nearest neighbors by cosine similarity (ties broken
#' deterministically by row name). With `weight = "cosine"` the incidence
#' entries are the similarities to the hyperedge's center (1 for the center
#' itself); `weight = "binary"` gives the unweighted membership hypergraph
#' (the `HyperG::knn_hypergraph()` construction, which this generalizes).
#'
#' Cosine similarity and Euclidean distance order neighbors identically on
#' L2-normalized embeddings, so normalized encoder output (the default of
#' `sbert::encode()`) gives the conventional kNN structure.
#'
#' @param embeddings Numeric matrix, one row per item, with unique non-empty
#'   rownames (the item IDs). No missing values; no all-zero rows.
#' @param k Number of neighbors per hyperedge (between 1 and
#'   `nrow(embeddings) - 1`).
#' @param weight `"cosine"` (default) or `"binary"`.
#'
#' @return A `net_hypergraph` (from [group_hypergraph()]) with one
#'   hyperedge per item, each of size `k + 1`, plus a `knn` field recording
#'   `k` and the weighting. Accepted by [hg_measures()], [hg_centrality()],
#'   [hg_cluster()], and [hg_classify()].
#'
#' @section Conditions: Raises `honets_bad_input` for broken contracts, and
#'   `honets_nonpositive_similarity` when `weight = "cosine"` selects a neighbor
#'   with similarity <= 0 -- a non-positive incidence weight would invalidate
#'   the random-walk machinery downstream; use `weight = "binary"` or a
#'   smaller `k` instead.
#'
#' @examples
#' set.seed(1)
#' emb <- matrix(rnorm(20, mean = 3), nrow = 5,
#'               dimnames = list(paste0("doc_", 1:5), NULL))
#' hg <- knn_hypergraph(emb, k = 2)
#' hg_measures(hg, what = "edges")
#' @export
knn_hypergraph <- function(embeddings, k, weight = c("cosine", "binary")) {
  weight <- match.arg(weight)
  stopifnot(
    "`embeddings` must be a numeric matrix" =
      is.matrix(embeddings) && is.numeric(embeddings),
    "`embeddings` must have at least two rows" = nrow(embeddings) >= 2L,
    "`embeddings` must not contain missing values" = !anyNA(embeddings),
    "`k` must be a single integer between 1 and nrow(embeddings) - 1" =
      length(k) == 1L && is.finite(k) && k >= 1 && k <= nrow(embeddings) - 1L
  )
  ids <- rownames(embeddings)
  if (is.null(ids) || anyNA(ids) || anyDuplicated(ids) > 0L ||
      !all(nzchar(ids))) {
    stop(errorCondition(
      "`embeddings` must carry unique, non-empty rownames (the item IDs)",
      class = "honets_bad_input", call = NULL
    ))
  }
  norms <- sqrt(rowSums(embeddings^2))
  if (any(norms < sqrt(.Machine$double.eps))) {
    stop(errorCondition(
      "`embeddings` contains all-zero rows; cosine similarity is undefined",
      class = "honets_bad_input", call = NULL
    ))
  }

  sims <- tcrossprod(embeddings / norms)
  diag(sims) <- 1

  # k+1 members per edge: the center plus its k most similar rows,
  # similarity descending, ties broken by ID for determinism.
  members <- lapply(seq_len(nrow(sims)), \(i) {
    ranked <- order(-sims[i, ], ids)
    c(i, ranked[ranked != i][seq_len(k)])
  })
  member_idx <- unlist(members, use.names = FALSE)
  center_idx <- rep(seq_along(ids), each = k + 1L)
  w <- if (identical(weight, "cosine")) {
    sims[cbind(center_idx, member_idx)]
  } else {
    rep(1, length(member_idx))
  }
  if (any(w <= 0)) {
    stop(errorCondition(
      paste(
        "a selected neighbor has cosine similarity <= 0, which would put a",
        "non-positive weight in the incidence; use weight = \"binary\" or a",
        "smaller `k`"
      ),
      class = "honets_nonpositive_similarity", call = NULL
    ))
  }

  long <- data.frame(
    item = ids[member_idx],
    edge = ids[center_idx],
    w = w
  )
  hg <- group_hypergraph(long, actor = "item", group = "edge",
                                    weight = "w")
  hg$knn <- list(k = as.integer(k), weight = weight)
  hg
}
