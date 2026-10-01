# ---- hg_get() for net_path_dependence ----
#
# memory() wraps the Nestimate estimator and returns its object
# unchanged; its per-context table is read here.

#' Per-context table of a path-dependence analysis
#'
#' @param x A `net_path_dependence` object from [memory()].
#' @param what `"contexts"`, the only table.
#' @param ... Additional arguments (ignored).
#' @param min_count Integer or `NULL`. Keep only contexts observed at least
#'   this many times.
#' @param flips `NULL` (default, every context), `TRUE` for only the contexts
#'   whose most probable next state differs from the first-order prediction
#'   (`flips` is `TRUE`), or `FALSE` for only those where it agrees.
#' @param sort_by `"kl"` (default) sorts by the order-k vs order-1
#'   divergence, largest first; `NULL` keeps the construction order.
#' @return A data.frame, one row per context: `context`, `count` (times
#'   observed), `entropy_first_order` and `entropy_order_k` (entropy of the
#'   next state given the last state and given the whole context),
#'   `entropy_drop` (their difference), `kl` (Kullback-Leibler divergence of
#'   the order-k from the first-order distribution), `top_first_order` and
#'   `top_order_k` (the most probable next state under each), and `flips`
#'   (`TRUE` when the two differ).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' # one sequence per row
#' wide <- data.frame(t1 = c("a", "x", "a", "x"), t2 = c("b", "b", "b", "b"),
#'                    t3 = c("c", "d", "c", "d"), t4 = c("a", "x", "a", "x"),
#'                    t5 = c("b", "b", "b", "b"), t6 = c("c", "d", "c", "d"))
#' pd <- memory(wide, order = 2, min_count = 1)
#' hg_get(pd)
#' hg_get(pd, flips = TRUE)
#' @export
hg_get.net_path_dependence <- function(x, what = "contexts", ...,
                                       min_count = NULL, flips = NULL,
                                       sort_by = "kl", top = NULL) {
  match.arg(what, "contexts")
  out <- .ho_rename(x$contexts,
                    c(n = "count", H_order1 = "entropy_first_order",
                      H_orderk = "entropy_order_k", H_drop = "entropy_drop",
                      KL = "kl", top_o1 = "top_first_order",
                      top_ok = "top_order_k"))
  if (!is.null(min_count)) {
    stopifnot("`min_count` must be a single integer >= 1" =
                is.numeric(min_count) && length(min_count) == 1L &&
                min_count >= 1)
    out <- out[out$count >= min_count, , drop = FALSE]
  }
  if (!is.null(flips)) {
    stopifnot("`flips` must be NULL, TRUE or FALSE" =
                is.logical(flips) && length(flips) == 1L && !is.na(flips))
    out <- out[out$flips %in% flips, , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, "kl")
    out <- out[order(-out[[sort_by]], out$context), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
