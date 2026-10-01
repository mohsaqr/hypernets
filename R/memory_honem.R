# ---- hg_get() for net_honem ----
#
# honem() returns the embedding object of the HONEM builder unchanged; its
# tables are read here.

#' Tables of a higher-order network embedding
#'
#' @param x A `net_honem` object from [honem()].
#' @param what `"embeddings"` (default) for the node coordinates, or
#'   `"variance"` for the per-dimension singular values.
#' @param ... Unused.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return A data.frame. For `what = "embeddings"`, one row per higher-order
#'   node: `node` followed by one column per embedding dimension (`dim1`,
#'   `dim2`, ...). For `what = "variance"`, one row per dimension: `dimension`,
#'   `singular_value`, `proportion` (that dimension's share of the total
#'   squared singular value).
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' emb <- honem(hon(seqs, max_order = 2), dim = 2)
#' hg_get(emb)
#' hg_get(emb, what = "variance")
#' @export
hg_get.net_honem <- function(x, what = c("embeddings", "variance"), ...,
                             top = NULL) {
  what <- match.arg(what)
  if (what == "variance") {
    return(.ho_top(data.frame(
      dimension      = seq_along(x$singular_values),
      singular_value = as.numeric(x$singular_values),
      proportion     = as.numeric(x$singular_values^2 /
                                    sum(x$singular_values^2)),
      stringsAsFactors = FALSE
    ), top))
  }
  emb <- x$embeddings
  out <- data.frame(node = x$nodes, stringsAsFactors = FALSE)
  emb_df <- as.data.frame(emb, stringsAsFactors = FALSE)
  names(emb_df) <- paste0("dim", seq_len(ncol(emb)))
  out <- cbind(out, emb_df)
  rownames(out) <- NULL
  .ho_top(out, top)
}
