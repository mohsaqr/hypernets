# ---- hg_get() for net_markov_stability ----
#
# hg_markov_stability() wraps the Nestimate estimator and returns its object
# unchanged: list(stability = <rounded per-state table>, mpt = <net_mpt:
# matrix, stationary, return_times, states>). The "passage_time" and
# "stationary" tables are read from the unrounded `mpt` component; the
# "states" table is the estimator's, which rounds its columns (4 / 2
# decimals).

#' Tables of a Markov stability analysis
#'
#' @param x A `net_markov_stability` object from
#'   [hg_markov_stability()].
#' @param ... Ignored.
#' @param what Which table: `"states"` (default, one row per state),
#'   `"passage_time"` (one row per ordered pair of states) or
#'   `"stationary"` (one row per state, distribution only).
#' @param sort_by `NULL` (construction order, the default) or a column name
#'   of the selected table.
#' @param decreasing Logical. Sort `sort_by` from largest to smallest
#'   (`TRUE`, the default) or smallest to largest (`FALSE`). Ties are broken
#'   by the row key, ascending in both directions: `state`, or `from` then
#'   `to` for the passage table.
#' @param from,to `NULL` (default, every state) or a character vector of
#'   states. Restrict the `what = "passage_time"` table to passages leaving
#'   the `from` states and arriving at the `to` states. Only that table has
#'   `from`/`to` columns; using either with another `what`, or naming a state
#'   the chain does not have, raises `hypernets_bad_input`.
#' @param top Integer or `NULL`. Keep only the first `top` rows, applied
#'   last: after `what`, `from`/`to` and `sort_by`.
#' @return A base `data.frame`. For `what = "states"`, one row per state with
#'   `state`, `persistence`, `stationary`, `return_time`,
#'   `sojourn_time`, `avg_time_to_others` and `avg_time_from_others`, as
#'   the estimator stores them (probabilities rounded to 4 decimals, times to 2). For
#'   `what = "passage_time"`, one row per ordered pair with `from`, `to` and
#'   `steps`. For `what = "stationary"`, one row per state with `state`,
#'   `stationary` and `return_time`; the passage and stationary tables
#'   are unrounded.
#' @examples
#' P <- matrix(c(0.5, 0.3, 0.2,
#'               0.2, 0.5, 0.3,
#'               0.1, 0.4, 0.5), 3, 3, byrow = TRUE,
#'             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' ms <- hg_markov_stability(P)
#' hg_get(ms, sort_by = "stationary")
#' hg_get(ms, what = "passage_time", from = "A")
#' hg_get(ms, what = "stationary")
#' @export
hg_get.net_markov_stability <- function(x, what = c("states",
                                                    "passage_time",
                                                    "stationary"), ...,
                                        sort_by = NULL,
                                               decreasing = TRUE,
                                               from = NULL, to = NULL,
                                               top = NULL) {
  what <- match.arg(what)
  stopifnot(
    "`decreasing` must be a single TRUE or FALSE" =
      is.logical(decreasing) && length(decreasing) == 1L && !is.na(decreasing)
  )
  if (!identical(what, "passage_time") && (!is.null(from) || !is.null(to))) {
    stop(errorCondition(
      sprintf(paste0("`from` and `to` filter the passage table; ",
                     "what = \"%s\" has no from/to columns. ",
                     "Use what = \"passage_time\"."), what),
      class = "hypernets_bad_input", call = NULL))
  }
  mpt <- x$mpt
  .hms_check_states(from, "from", mpt$states)
  .hms_check_states(to, "to", mpt$states)

  out <- switch(
    what,
    states = .ho_rename(x$stability, c(stationary_prob = "stationary")),
    stationary = data.frame(
      state           = mpt$states,
      stationary      = unname(mpt$stationary),
      return_time     = unname(mpt$return_times),
      stringsAsFactors = FALSE
    ),
    passage_time = {
      M <- mpt$matrix
      data.frame(
        from  = rep(rownames(M), times = ncol(M)),
        to    = rep(colnames(M), each  = nrow(M)),
        steps = as.vector(M),
        stringsAsFactors = FALSE
      )
    }
  )
  if (!is.null(from)) out <- out[out$from %in% from, , drop = FALSE]
  if (!is.null(to))   out <- out[out$to %in% to, , drop = FALSE]

  if (!is.null(sort_by)) {
    stopifnot(
      "`sort_by` must be a single column name of the returned table" =
        is.character(sort_by) && length(sort_by) == 1L &&
        sort_by %in% names(out)
    )
    # Primary key signed so one ascending order() serves both directions;
    # xtfrm() makes a character column sortable the same way. The row key
    # breaks ties deterministically and stays ascending either way.
    primary <- xtfrm(out[[sort_by]])
    if (decreasing) primary <- -primary
    keys <- if (identical(what, "passage_time")) {
      list(out$from, out$to)
    } else {
      list(out$state)
    }
    out <- out[do.call(order, c(list(primary), keys)), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Validate a from/to state filter against the chain's states
#'
#' @param states `NULL` or the character vector the caller passed.
#' @param arg Argument name, for the message.
#' @param known Character vector of the chain's states.
#' @return `NULL`, invisibly; raises `hypernets_bad_input` on bad input.
#' @noRd
.hms_check_states <- function(states, arg, known) {
  if (is.null(states)) return(invisible(NULL))
  if (!is.character(states) || !length(states) || anyNA(states)) {
    stop(errorCondition(
      sprintf("`%s` must be a non-empty character vector of states.", arg),
      class = "hypernets_bad_input", call = NULL))
  }
  unknown <- setdiff(states, known)
  if (length(unknown)) {
    stop(errorCondition(
      sprintf("`%s` names state(s) the chain does not have: %s.", arg,
              paste(unknown, collapse = ", ")),
      class = "hypernets_bad_input", call = NULL))
  }
  invisible(NULL)
}
