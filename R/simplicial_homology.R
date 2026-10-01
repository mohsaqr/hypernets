# ---- hg_get() for persistent_homology ----
#
# hg_homology() wraps the Nestimate estimator and returns its object
# unchanged; its tables are read here.

#' Tables of a persistent homology result
#'
#' @param x A `persistent_homology` object from
#'   [hg_homology()].
#' @param ... Additional arguments (ignored).
#' @param what `"persistence"` (default) for the persistence diagram, or
#'   `"betti"` for the Betti curves across the filtration.
#' @param dimension Integer or `NULL`. Keep only this homological dimension.
#' @param sort_by `NULL` (construction order, default) or `"persistence"` -
#'   longest-lived features first. Only for `what = "persistence"`.
#' @return A data.frame. For `what = "persistence"`, one row per feature:
#'   `dimension`, `birth`, `death`, `persistence`. For `what = "betti"`, one
#'   row per (threshold, dimension): `threshold`, `dimension`, `betti`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' set.seed(1)
#' w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
#' w <- (w + t(w)) / 2
#' diag(w) <- 0
#' ph <- hg_homology(w, n_steps = 10)
#' hg_get(ph, sort_by = "persistence")
#' hg_get(ph, what = "betti", dimension = 0)
#' @export
hg_get.persistent_homology <- function(
    x, what = c("persistence", "betti"), ..., dimension = NULL,
    sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  out <- if (what == "betti") x$betti_curve else x$persistence
  if (!is.null(dimension)) {
    stopifnot("`dimension` must be a single integer >= 0" =
                is.numeric(dimension) && length(dimension) == 1L &&
                dimension >= 0)
    out <- out[out$dimension == as.integer(dimension), , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, "persistence")
    if (what == "betti") {
      stop(errorCondition(
        "`sort_by` applies only to what = \"persistence\"",
        class = "hypernets_bad_input", call = NULL))
    }
    out <- out[order(-out$persistence, out$dimension, out$birth), ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
