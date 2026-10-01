# ---- hg_get() for net_mogen ----
#
# mogen() returns the multi-order generative model unchanged; its tables are
# read here. "paths" replaces a raw-data path counter: the k-step path counts
# are the counts the fitted order-(k - 1) layer already holds, so they are
# read off the model rather than recomputed from the sequences.

#' Tables of a multi-order generative model
#'
#' @param x A `net_mogen` object from [mogen()].
#' @param what `"orders"` (default) for the per-order model-selection table,
#'   `"transitions"` for the fitted transition probabilities of one layer,
#'   `"paths"` for the frequency of every observed path of `k` states, or
#'   `"pathways"` for the layer's transitions written as pathways
#'   (`"a b -> c"`), the form [simplicial()] turns into simplices.
#' @param ... For `what = "pathways"`, `min_prob` (drop transitions below
#'   this probability). Otherwise unused.
#' @param order Integer or `NULL`. For `what = "transitions"` and
#'   `"pathways"`: which layer to read. Default `NULL` uses the selected
#'   order.
#' @param k Integer vector, for `what = "paths"`: the path length or
#'   lengths in states, each between `2` and the model's highest order plus
#'   one, e.g. `k = 3` or `k = 2:4`. Default `2` (single transitions).
#' @param min_count Integer, for `what = "transitions"` and `"pathways"`:
#'   drop transitions observed fewer times. Default `1`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return A data.frame. For `what = "orders"`, one row per candidate order:
#'   `order`, `log_likelihood`, `aic`, `bic`, `df`, `layer_df`, and
#'   `optimal` (logical, `TRUE` for the selected order). For
#'   `what = "transitions"`, one row per transition: `path`, `count`,
#'   `probability`, `from`, `to`. For `what = "paths"`, one row per distinct
#'   path of each length in `k` (`top` applies within each length), grouped
#'   by length and most frequent first (ties alphabetical): `k`, `path`,
#'   `count`, `proportion` (share of all `k`-state paths, rounded to 4
#'   decimals). For `what = "pathways"`, one row per transition with the
#'   column `pathway`.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' mg <- mogen(seqs, max_order = 2)
#' hg_get(mg)
#' hg_get(mg, what = "transitions", top = 5)
#' hg_get(mg, what = "paths", k = 3)
#' @export
hg_get.net_mogen <- function(x, what = c("orders", "transitions", "paths",
                                         "pathways"), ...,
                             order = NULL, k = 2L, min_count = 1L,
                             top = NULL) {
  what <- match.arg(what)
  if (what == "transitions") {
    return(.ho_top(Nestimate::mogen_transitions(
      x, order = order %||% x$optimal_order, min_count = min_count), top))
  }
  if (what == "pathways") {
    return(.ho_top(.ho_pathway_table(Nestimate::pathways(
      x, order = order, min_count = min_count, ...)), top))
  }
  if (what == "paths") {
    stopifnot("`k` must be a non-empty vector of path lengths" =
                is.numeric(k) && length(k) >= 1L && !anyNA(k))
    tables <- lapply(unique(k), \(kk) {
      counts <- .ho_top(.hmg_path_counts(x, kk), top)
      data.frame(k = rep(as.integer(kk), nrow(counts)), counts,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, tables)
    rownames(out) <- NULL
    return(out)
  }
  out <- data.frame(
    order          = as.integer(x$orders),
    log_likelihood = as.numeric(x$log_likelihood),
    aic            = as.numeric(x$aic),
    bic            = as.numeric(x$bic),
    df             = as.integer(x$dof),
    layer_df       = as.integer(x$layer_dof),
    stringsAsFactors = FALSE
  )
  out$optimal <- out$order == x$optimal_order
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Path counts of length k read off a fitted multi-order model
#'
#' The order-(k - 1) layer of a MOGen model counts every transition between
#' overlapping (k - 1)-grams, i.e. every observed path of k states. The table
#' is laid out as a raw-data k-gram count: most frequent first, ties in
#' sorted path order, proportion of all k-state paths rounded to 4 decimals.
#'
#' @param x A `net_mogen`.
#' @param k Integer path length in states.
#' @return data.frame: path, count, proportion.
#' @noRd
.hmg_path_counts <- function(x, k) {
  max_k <- max(x$orders) + 1L
  if (!is.numeric(k) || length(k) != 1L || !is.finite(k) || k != round(k) ||
      k < 2L || k > max_k) {
    stop(errorCondition(
      sprintf(paste0("`k` must be a single whole number from 2 to %d (the ",
                     "model's highest order plus one); refit mogen() with a ",
                     "larger `max_order` for longer paths."), max_k),
      class = "hypernets_bad_input", call = NULL))
  }
  tr <- Nestimate::mogen_transitions(x, order = as.integer(k) - 1L,
                                     min_count = 1L)
  if (!nrow(tr)) {
    return(data.frame(path = character(0), count = integer(0),
                      proportion = numeric(0), stringsAsFactors = FALSE))
  }
  counts <- tapply(as.integer(tr$count), tr$path, sum)
  counts <- counts[sort(names(counts))]
  counts <- counts[order(as.integer(counts), decreasing = TRUE,
                         method = "radix")]
  out <- data.frame(path = names(counts), count = as.integer(counts),
                    proportion = round(as.numeric(counts) / sum(counts), 4),
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  out
}
