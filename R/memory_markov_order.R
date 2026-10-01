# ---- hg_get() for the Markov order test ----
#
# markov_order() wraps the Nestimate estimator and returns its object
# unchanged (a net_markov_order, or a net_markov_order_group for a
# netobject_group); their tables are read here.

#' Tables of a grouped Markov order test
#'
#' @description
#' The per-group test tables stacked into one data.frame, with a `group`
#' column naming the group each row came from. This is the supported way to
#' read a `net_markov_order_group`; the print method shows each group in turn
#' but a comparison across groups needs one table.
#'
#' @param x A `net_markov_order_group`, as returned by
#'   [markov_order()] on a `netobject_group`.
#' @param what `"orders"` (default) for one row per group and order, or
#'   `"null"` for the permutation null draws, one row per group, order and
#'   replicate.
#' @param ... Unused.
#' @param top Keep only the first `top` rows, applied last.
#' @return A base `data.frame`. For `what = "orders"`, one row per group and
#'   order, with `group` first and then the columns of
#'   [hg_get.net_markov_order()]. For `what = "null"`, one row per
#'   group, order and replicate with columns `group`, `order`, `replicate`,
#'   `g2`. Returns a zero-row frame with those columns when the group is
#'   empty.
#' @examples
#' seqs <- list(c("a", "b", "a", "c"), c("b", "a", "c", "a"))
#' fit <- markov_order(seqs, max_order = 2L, n_perm = 20L, seed = 1L)
#' hg_get(fit)
#' @export
hg_get.net_markov_order_group <- function(x, what = c("orders", "null"),
                                          ..., top = NULL) {
  what <- match.arg(what)
  groups <- names(x) %||% as.character(seq_along(x))
  parts <- Map(function(fit, label) {
    tab <- hg_get(fit, what = what)
    if (nrow(tab) == 0L) return(NULL)
    cbind(group = rep.int(label, nrow(tab)), tab, stringsAsFactors = FALSE)
  }, x, groups)
  out <- do.call(rbind, parts)
  if (is.null(out)) {
    out <- if (identical(what, "orders")) {
      data.frame(group = character(0), stringsAsFactors = FALSE)
    } else {
      data.frame(group = character(0), order = integer(0),
                 replicate = integer(0), g2 = numeric(0),
                 stringsAsFactors = FALSE)
    }
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Tables of a Markov order test
#'
#' @param x A `net_markov_order` object from [markov_order()].
#' @param ... Additional arguments (ignored).
#' @param what `"orders"` (default) for the per-order test table, or
#'   `"null"` for the permutation null distribution behind it.
#' @return A data.frame. For `what = "orders"`, one row per candidate order:
#'   `order`, `log_likelihood`, `aic`, `bic`, `df`, `g2` (the
#'   likelihood-ratio statistic against the order below), `p_value` (its
#'   permutation p-value), `p_asymptotic` (its chi-squared p-value) and
#'   `significant`. For `what = "null"`, one row per
#'   permutation replicate: `order`, `replicate`, `g2` (the replicate's
#'   likelihood-ratio statistic under the null).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' seqs <- list(c("a", "b", "a", "c", "a", "b"), c("b", "a", "c", "a", "b", "a"))
#' fit <- markov_order(seqs, max_order = 2L, n_perm = 20L, seed = 1L)
#' hg_get(fit)
#' hg_get(fit, what = "null", top = 5)
#' @export
hg_get.net_markov_order <- function(x, what = c("orders", "null"), ...,
                                    top = NULL) {
  what <- match.arg(what)
  if (what == "orders") {
    out <- .ho_rename(x$test_table,
                      c(loglik = "log_likelihood", AIC = "aic", BIC = "bic",
                        p_permutation = "p_value"))
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  nulls <- x$permutation_null
  if (!length(nulls)) {
    return(data.frame(order = integer(0), replicate = integer(0),
                      g2 = numeric(0), stringsAsFactors = FALSE))
  }
  # Null draws are named by order ("order_1", "order_2", ...); otherwise
  # fall back to position.
  numbered <- !is.null(names(nulls)) &&
    all(grepl("^(order_)?[0-9]+$", names(nulls)))
  ord <- if (numbered) {
    as.integer(sub("^order_", "", names(nulls)))
  } else {
    seq_along(nulls)
  }
  out <- do.call(rbind, Map(function(v, k) {
    data.frame(order = rep.int(k, length(v)),
               replicate = seq_along(v),
               g2 = as.numeric(v),
               stringsAsFactors = FALSE)
  }, nulls, ord))
  rownames(out) <- NULL
  .ho_top(out, top)
}
