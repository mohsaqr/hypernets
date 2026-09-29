# ---- Markov stability: dwell, recurrence and first passage ---------------
#
# Per-state stability of the random walk carried by a transition matrix.
# Ported verbatim (mathematics and floating-point pathway) from Nestimate
# 0.9.x `R/markov.R` (`markov_stability()` + the `.mpt_*` helpers it shares
# with `passage_time()`); the JS sibling tnaj carries the same method as
# `markovStability`.
#
# What hypernets adds is the INPUT, not the arithmetic. `build_hon()` returns
# a row-stochastic matrix over HIGHER-ORDER states ("c -> b" is the state b
# carrying the memory of c), so the same Kemeny-Snell machinery applied to
# that matrix describes the higher-order random walk: how long the process
# dwells in a memory state, how often it returns to one, how many steps it
# takes to reach one from another. Nestimate's first-order netobject cannot
# express that question.
#
# Deliberate departures from Nestimate, both name/presentation only:
#   * the table is stored at FULL PRECISION and rounded by print(). Nestimate
#     stores round(pi, 4); on a higher-order network with thousands of states
#     a stationary probability is ~1e-4 and four decimals would destroy it.
#     `test-equiv-markov-stability.R` proves identity by comparing the
#     unrounded quantities bit-for-bit AND the table after applying
#     Nestimate's own rounding.
#   * the mean-first-passage matrix is a field of this object rather than a
#     separate exported `net_mpt` class, reached with
#     `as.data.frame(x, what = "passage_time")` per the taxonomy.
# No computed value differs.

# ---------------------------------------------------------------------------
# Internal: input coercion
# ---------------------------------------------------------------------------

#' Pull a transition matrix out of whatever the caller supplied
#'
#' @param x A numeric matrix, a model object carrying `$weights`
#'   (`net_hon`, `net_hypa`, `net_mogen`, `cograph_network`, `netobject`,
#'   `tna`), or sequence data accepted by [build_hon()].
#' @return A square numeric matrix with row and column names.
#' @noRd
.hms_extract_P <- function(x) {
  if (is.matrix(x) && is.numeric(x)) return(x)

  if (inherits(x, "net_hon") || inherits(x, "net_hypa") ||
      inherits(x, "net_mogen") || inherits(x, "netobject") ||
      inherits(x, "cograph_network") || inherits(x, "tna")) {
    P <- x$weights
    if (is.null(P) || !is.matrix(P) || !is.numeric(P)) {
      stop(errorCondition(
        paste0("`x` is a ", class(x)[1L], " with no numeric weight matrix; ",
               "markov_stability() needs one."),
        class = "hypernets_bad_input", call = NULL))
    }
    return(P)
  }

  # Sequence input: the first-order chain is build_hon(max_order = 1), whose
  # weight matrix is the relative transition matrix.
  if (is.data.frame(x) || is.list(x)) {
    return(build_hon(x, max_order = 1L, min_freq = 1L)$weights)
  }

  stop(errorCondition(
    paste0("`x` must be a numeric transition matrix, a net_hon / net_hypa / ",
           "net_mogen / cograph_network / netobject / tna object, or ",
           "sequence data accepted by build_hon(); got ", class(x)[1L], "."),
    class = "hypernets_bad_input", call = NULL))
}

#' Check and, when permitted, rescale the rows of a transition matrix
#'
#' @param P Square numeric matrix.
#' @param state_names Character vector of state names.
#' @param normalize Logical. Rescale rows that do not sum to 1?
#' @return `P`, with rows summing to 1.
#' @noRd
.hms_normalize_rows <- function(P, state_names, normalize = TRUE) {
  row_sums <- rowSums(P)

  zero_rows <- which(row_sums <= 0)
  if (length(zero_rows)) {
    stop(errorCondition(
      paste0("Transition matrix has zero-sum row(s) for state(s): ",
             paste(state_names[zero_rows], collapse = ", "),
             ". These states have no outgoing transitions, so the chain is ",
             "not ergodic and mean first passage times are undefined. ",
             "Remove the state(s) or supply a different transition matrix; ",
             "for a pruned higher-order network, lower `min_freq` in ",
             "build_hon()."),
      class = c("hypernets_not_ergodic", "hypernets_bad_input"), call = NULL))
  }

  if (any(abs(row_sums - 1) > 1e-6)) {
    if (!normalize) {
      stop(errorCondition(
        "Transition matrix rows must sum to 1 when `normalize = FALSE`.",
        class = "hypernets_bad_input", call = NULL))
    }
    warning(warningCondition(
      "Rows do not sum to 1; normalizing.",
      class = "hypernets_renormalized", call = NULL))
    P <- P / row_sums
  }

  P
}

#' Stationary distribution of a row-stochastic matrix
#'
#' Left eigenvector of `P` for eigenvalue 1, taken as the eigenvector of
#' `t(P)` whose eigenvalue is closest to 1, made positive and normalized.
#'
#' @param P Square row-stochastic numeric matrix.
#' @return Numeric vector summing to 1, unnamed.
#' @noRd
.hms_stationary <- function(P) {
  ev   <- eigen(t(P))
  idx  <- which.min(abs(Re(ev$values) - 1))
  stat <- abs(Re(ev$vectors[, idx]))
  stat / sum(stat)
}

#' Is the support digraph of a transition matrix strongly connected?
#'
#' Reachability by repeated boolean squaring of `I + A`: after
#' `ceiling(log2(n))` squarings the matrix holds every path of length `< n`,
#' so an all-TRUE result is strong connectivity. The squarings run through
#' Reduce(): `log2(n)` rounds, not `n`.
#'
#' @param P Square numeric matrix.
#' @return Logical scalar.
#' @noRd
.hms_irreducible <- function(P) {
  n <- nrow(P)
  if (n == 1L) return(TRUE)
  reach <- (diag(n) + (P != 0)) != 0
  steps <- max(1L, ceiling(log2(n)))
  # repeated squaring: log2(n) rounds, each doubling the path length
  reach <- Reduce(function(r, i) (r %*% r) != 0, seq_len(steps), reach)
  all(reach)
}

#' Mean first passage times by the Kemeny-Snell fundamental matrix
#'
#' `Z = (I - P + Pi)^-1`, `M[i, j] = (Z[j, j] - Z[i, j]) / pi[j]`, with the
#' diagonal replaced by the mean recurrence time `1 / pi[i]`.
#'
#' @param P Square row-stochastic numeric matrix.
#' @param stat Stationary distribution of `P`.
#' @return Numeric matrix of expected steps from row state to column state.
#' @noRd
.hms_passage <- function(P, stat) {
  n     <- nrow(P)
  PI    <- matrix(stat, nrow = n, ncol = n, byrow = TRUE)
  Z     <- solve(diag(n) - P + PI)
  Zd    <- matrix(diag(Z), nrow = n, ncol = n, byrow = TRUE)
  pimat <- matrix(stat,   nrow = n, ncol = n, byrow = TRUE)
  M     <- (Zd - Z) / pimat
  diag(M) <- 1 / stat
  M
}


# ---------------------------------------------------------------------------
# markov_stability()
# ---------------------------------------------------------------------------

#' Markov Stability of a First- or Higher-Order Chain
#'
#' @description
#' Per-state stability of the random walk carried by a transition matrix:
#' how persistent a state is, how much of the walk's time it holds, how long
#' the process dwells in it before leaving, how long until it comes back, and
#' how many steps separate it from every other state.
#'
#' Given a `net_hon` the states are **higher-order** states — `"c -> b"` is
#' the state `b` reached from `c` — so the same quantities describe the
#' higher-order random walk rather than the memoryless one. `build_hon()`
#' returns a row-stochastic matrix over those states, which is exactly the
#' object this verb consumes.
#'
#' @param x A `net_hon` (from [build_hon()]), `net_hypa`, `net_mogen`,
#'   `cograph_network`, `netobject` or `tna` object carrying a weight matrix;
#'   a row-stochastic numeric transition matrix; or sequence data accepted by
#'   [build_hon()] (a wide `data.frame`, one trajectory per row, or a list of
#'   character vectors), in which case the first-order chain
#'   `build_hon(x, max_order = 1)` is used.
#' @param normalize Logical. Rescale rows that do not sum to 1 (with a
#'   `hypernets_renormalized` warning)? Default `TRUE`. `FALSE` makes a
#'   non-stochastic matrix an error instead.
#'
#' @return An object of class `"net_markov_stability"`. Reach it with
#'   [as.data.frame()]:
#'   \describe{
#'     \item{`as.data.frame(x)`}{one row per state, columns `state`,
#'       `persistence` (\eqn{P_{ii}}), `stationary_prob` (\eqn{\pi_i}),
#'       `return_time` (\eqn{1/\pi_i}), `sojourn_time`
#'       (\eqn{1/(1 - P_{ii})}), `avg_time_to_others` (mean first passage
#'       leaving the state) and `avg_time_from_others` (mean first passage
#'       arriving at it).}
#'     \item{`as.data.frame(x, what = "passage_time")`}{one row per ordered
#'       pair of states, columns `from`, `to`, `steps`; the diagonal carries
#'       the mean recurrence time.}
#'     \item{`as.data.frame(x, what = "stationary")`}{one row per state,
#'       columns `state`, `stationary_prob`, `return_time`.}
#'   }
#'   `summary()` names the attractor and the stickiest state, `print()`
#'   shows the table, and `plot()` draws a stability landscape, a per-state
#'   profile, a first-passage heatmap or the transition network (see
#'   [plot.net_markov_stability()]). The row-normalised transition matrix
#'   the quantities were computed from is kept for the network view.
#'
#' @details
#' \strong{Persistence} is the self-transition probability \eqn{P_{ii}}.
#'
#' \strong{Sojourn time} is the expected number of consecutive steps spent in
#' a state before leaving it, \eqn{1/(1 - P_{ii})} — the mean of the
#' geometric dwell distribution. A state with \eqn{P_{ii} = 1} is absorbing
#' and its sojourn time is `Inf`.
#'
#' \strong{Return time} is the mean recurrence time \eqn{1/\pi_i}, the
#' expected number of steps between successive visits.
#'
#' \strong{avg_time_to_others} / \strong{avg_time_from_others} are the
#' off-diagonal row and column means of the mean-first-passage matrix: how
#' far the state is from the rest of the chain, and how accessible it is from
#' it (a small value marks an attractor).
#'
#' Mean first passage times use the Kemeny-Snell fundamental matrix
#' \deqn{Z = (I - P + \Pi)^{-1}, \qquad
#'       M_{ij} = \frac{Z_{jj} - Z_{ij}}{\pi_j},}
#' with \eqn{\Pi_{ij} = \pi_j}. This requires an **irreducible** chain: a
#' reducible one has no unique stationary distribution and
#' \eqn{(I - P + \Pi)} is singular. That is not a rare corner for a
#' higher-order network — pruning a `net_hon` with `min_freq` can leave an
#' absorbing or unreachable higher-order state — so irreducibility is checked
#' up front and reported as `hypernets_not_ergodic` rather than surfacing as
#' a LAPACK singularity.
#'
#' \strong{Heavy pruning.} `min_freq` in [build_hon()] drops rare
#' transitions, and on a small corpus that can leave a higher-order state
#' with no outgoing transition at all, or split the chain into parts that
#' cannot reach each other. On the AI turns of `ai_long`,
#' `build_hon(ai_long, action = "code", actor = "session_id",
#' time = "order_in_session", max_order = 2, min_freq = 50)` leaves the
#' state `Ask` with no outgoing transition, and `min_freq = 20` leaves a
#' reducible chain; `markov_stability()` refuses both with
#' `hypernets_not_ergodic`. Lower `min_freq` (`min_freq = 10` gives an
#' irreducible chain there) rather than analysing a chain whose passage
#' times are undefined.
#'
#' @section Conditions:
#' Raises `hypernets_bad_input` for an unsupported `x`, a non-square or
#' single-state matrix, missing values, and for rows that do not sum to 1
#' under `normalize = FALSE`. Raises `hypernets_not_ergodic` (which also
#' inherits `hypernets_bad_input`) when a state has no outgoing transitions
#' or the chain is reducible. Warns `hypernets_renormalized` when rows are
#' rescaled.
#'
#' @seealso [build_hon()] for the higher-order chain, [markov_order_test()]
#'   for how far the memory reaches, [hon_centrality()] for path-based
#'   centralities of the same network.
#'
#' @references
#' Kemeny, J. G. and Snell, J. L. (1976). \emph{Finite Markov Chains}.
#' Springer-Verlag.
#'
#' Rosvall, M., Esquivel, A. V., Lancichinetti, A., West, J. D. and
#' Lambiotte, R. (2014). Memory in network flows and its effects on spreading
#' dynamics and community detection. \emph{Nature Communications}, 5, 4630.
#' \doi{10.1038/ncomms5630}
#'
#' Scholtes, I., Wider, N., Pfitzner, R., Garas, A., Tessone, C. J. and
#' Schweitzer, F. (2014). Causality-driven slow-down and speed-up of
#' diffusion in non-Markovian temporal networks. \emph{Nature
#' Communications}, 5, 5024. \doi{10.1038/ncomms6024}
#'
#' @examples
#' # A first-order chain, straight from a transition matrix
#' P <- matrix(c(0.7, 0.2, 0.1,
#'               0.3, 0.5, 0.2,
#'               0.2, 0.3, 0.5),
#'             nrow = 3, byrow = TRUE,
#'             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' markov_stability(P)
#'
#' # The higher-order walk: states that carry their memory
#' hon <- build_hon(split(human_long$code, human_long$session_id),
#'                  max_order = 2L, min_freq = 50L)
#' ms <- markov_stability(hon)
#' as.data.frame(ms, sort_by = "stationary_prob", top = 5L)
#' as.data.frame(ms, what = "passage_time", sort_by = "steps", top = 5L)
#' # the least persistent states, and the quickest passages out of "Specify"
#' as.data.frame(ms, sort_by = "persistence", decreasing = FALSE, top = 5L)
#' as.data.frame(ms, what = "passage_time", from = "Specify",
#'               sort_by = "steps", decreasing = FALSE, top = 5L)
#'
#' @export
markov_stability <- function(x, normalize = TRUE) {
  stopifnot(
    "`normalize` must be a single TRUE or FALSE" =
      is.logical(normalize) && length(normalize) == 1L && !is.na(normalize)
  )

  P <- .hms_extract_P(x)
  if (nrow(P) != ncol(P)) {
    stop(errorCondition(
      sprintf("`x` must give a square transition matrix; got %d x %d.",
              nrow(P), ncol(P)),
      class = "hypernets_bad_input", call = NULL))
  }
  if (anyNA(P)) {
    stop(errorCondition(
      "Transition matrix has missing values.",
      class = "hypernets_bad_input", call = NULL))
  }
  if (nrow(P) < 2L) {
    stop(errorCondition(
      paste0("A Markov stability analysis needs at least 2 states; `x` has ",
             nrow(P), "."),
      class = "hypernets_bad_input", call = NULL))
  }

  state_names <- colnames(P)
  if (is.null(state_names)) {
    state_names <- paste0("S", seq_len(nrow(P)))
    colnames(P) <- rownames(P) <- state_names
  }

  P <- .hms_normalize_rows(P, state_names, normalize = normalize)

  # Kemeny-Snell needs an irreducible chain: a reducible one has no unique
  # stationary distribution and (I - P + Pi) is singular, so solve() would
  # fail with a LAPACK message that says nothing about the data. Higher-order
  # networks hit this routinely -- a pruned HON can leave an absorbing
  # higher-order state -- so the diagnosis is worth making explicit.
  if (!.hms_irreducible(P)) {
    stop(errorCondition(
      paste0("The chain is reducible: its ", nrow(P), " states are not all ",
             "mutually reachable, so the stationary distribution is not ",
             "unique and mean first passage times are undefined. For a ",
             "higher-order network this usually means pruning left an ",
             "absorbing or unreachable higher-order state; change `min_freq` ",
             "or `max_order` in build_hon(), or analyze one recurrent ",
             "component."),
      class = c("hypernets_not_ergodic", "hypernets_bad_input"), call = NULL))
  }

  stat <- .hms_stationary(P)
  names(stat) <- state_names
  # Defensive: Perron-Frobenius guarantees a strictly positive stationary
  # vector for the irreducible chain checked above, so this fires only on a
  # numerical underflow in eigen() on a very large or badly scaled chain.
  if (any(stat <= 0)) {
    warning(warningCondition(
      paste0("Stationary probabilities underflowed to zero; the passage ",
             "times below are not trustworthy."),
      class = "hypernets_not_ergodic", call = NULL))
    stat <- pmax(stat, .Machine$double.eps)
    stat <- stat / sum(stat)
  }

  M <- .hms_passage(P, stat)
  dimnames(M) <- list(state_names, state_names)

  n           <- nrow(P)
  persistence <- diag(P)
  sojourn     <- ifelse(persistence >= 1 - .Machine$double.eps,
                        Inf, 1 / (1 - persistence))
  avg_to   <- vapply(seq_len(n), function(i) mean(M[i, -i]), numeric(1L))
  avg_from <- vapply(seq_len(n), function(i) mean(M[-i, i]), numeric(1L))

  stability <- data.frame(
    state                = state_names,
    persistence          = unname(persistence),
    stationary_prob      = unname(stat),
    return_time          = unname(1 / stat),
    sojourn_time         = unname(sojourn),
    avg_time_to_others   = avg_to,
    avg_time_from_others = avg_from,
    stringsAsFactors     = FALSE
  )
  rownames(stability) <- NULL

  structure(
    list(
      stability    = stability,
      passage_time = M,
      transition   = P,
      stationary   = stat,
      return_times = 1 / stat,
      states       = state_names,
      n_states     = n,
      normalize    = normalize
    ),
    class = "net_markov_stability"
  )
}


# ---------------------------------------------------------------------------
# Accessor
# ---------------------------------------------------------------------------

#' Tidy accessor for net_markov_stability
#'
#' @param x A `net_markov_stability` object from [markov_stability()].
#' @param row.names Ignored, present for S3 consistency.
#' @param optional Ignored, present for S3 consistency.
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
#'   `state`, `persistence`, `stationary_prob`, `return_time`,
#'   `sojourn_time`, `avg_time_to_others` and `avg_time_from_others`. For
#'   `what = "passage_time"`, one row per ordered pair with `from`, `to` and
#'   `steps`. For `what = "stationary"`, one row per state with `state`,
#'   `stationary_prob` and `return_time`.
#' @inherit markov_stability examples
#' @export
as.data.frame.net_markov_stability <- function(x, row.names = NULL,
                                               optional = FALSE, ...,
                                               what = c("states",
                                                        "passage_time",
                                                        "stationary"),
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
  .hms_check_states(from, "from", x$states)
  .hms_check_states(to, "to", x$states)

  out <- switch(
    what,
    states = x$stability,
    stationary = data.frame(
      state           = x$states,
      stationary_prob = unname(x$stationary),
      return_time     = unname(x$return_times),
      stringsAsFactors = FALSE
    ),
    passage_time = {
      M <- x$passage_time
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


# ---------------------------------------------------------------------------
# print / summary / plot
# ---------------------------------------------------------------------------

#' Print Method for net_markov_stability
#'
#' @param x A `net_markov_stability` object.
#' @param digits Integer. Decimal places for the probability columns;
#'   time columns are shown at two fewer. Default `4`.
#' @param top Integer or `NULL`. Show only the `top` states with the largest
#'   stationary probability. Default `NULL` (all states).
#' @param ... Ignored.
#' @return `x` invisibly.
#' @inherit markov_stability examples
#' @export
print.net_markov_stability <- function(x, digits = 4L, top = NULL, ...) {
  cat(sprintf("Markov Stability -- %d states\n\n", x$n_states))
  d <- as.data.frame(x, sort_by = if (is.null(top)) NULL else "stationary_prob",
                     top = top)
  d$persistence     <- round(d$persistence, digits)
  d$stationary_prob <- round(d$stationary_prob, digits)
  time_cols <- c("return_time", "sojourn_time",
                 "avg_time_to_others", "avg_time_from_others")
  d[time_cols] <- lapply(d[time_cols], round, max(digits - 2L, 0L))
  print(d, row.names = FALSE)
  if (!is.null(top) && top < x$n_states) {
    cat(sprintf("\n... %d further states.\n", x$n_states - as.integer(top)))
  }
  cat("\nas.data.frame(x) for the full table,",
      "as.data.frame(x, what = \"passage_time\") for first passage times.\n")
  invisible(x)
}

#' Summary Method for net_markov_stability
#'
#' @param object A `net_markov_stability` object.
#' @param ... Ignored.
#' @return The per-state `data.frame` (as [as.data.frame()] returns it),
#'   with the attractor, the stickiest state and the chain's Kemeny constant
#'   attached as attributes. Returned **invisibly**: `summary(x)` prints the
#'   summary; assign the result to keep the table.
#' @inherit markov_stability examples
#' @export
summary.net_markov_stability <- function(object, ...) {
  d <- as.data.frame(object)
  attractor <- d$state[which.min(d$avg_time_from_others)]
  stickiest <- d$state[which.max(d$sojourn_time)]
  # Kemeny constant: expected steps to a state drawn from the stationary
  # distribution, the same from every starting state in an ergodic chain.
  kemeny <- sum(object$passage_time[1L, ] * object$stationary) -
    object$passage_time[1L, 1L] * object$stationary[1L]

  cat(sprintf("Markov Stability -- %d states\n\n", object$n_states))
  cat(sprintf("Most accessible state (attractor): %s (%.2f steps from others)\n",
              attractor, min(d$avg_time_from_others)))
  cat(sprintf("Most persistent state (stickiest): %s (%.2f steps of dwell)\n",
              stickiest, max(d$sojourn_time)))
  cat(sprintf("Kemeny constant (mixing):          %.2f steps\n",
              unname(kemeny)))
  cat(sprintf("Stationary probability range:      %.4g to %.4g\n",
              min(d$stationary_prob), max(d$stationary_prob)))

  attr(d, "attractor") <- attractor
  attr(d, "stickiest") <- stickiest
  attr(d, "kemeny")    <- unname(kemeny)
  invisible(d)
}

#' Plot Method for net_markov_stability
#'
#' @description
#' Four views of the same chain. Three delegate to cograph's own plotting
#' functions, which own the layout, labels, legend and theme; only the
#' landscape, which cograph has no counterpart for, is drawn here.
#'
#' \describe{
#'   \item{`what = "landscape"` (default)}{\emph{Which states matter, and how?}
#'     One point per state: stationary share (x, log scale) against
#'     persistence (y). Dashed lines at the even share \eqn{1/n} and at the
#'     mean persistence split the plane into four named quadrants:
#'     \strong{hubs} (common and sticky), \strong{relays} (common, quickly
#'     left), \strong{traps} (rare but sticky) and \strong{transients} (rare,
#'     quickly left). A state is common when its share is strictly above
#'     \eqn{1/n} and sticky when its persistence is strictly above the mean.
#'     Point size is the mean stay (sojourn time); shape separates memory
#'     states (`"a -> b"`) from first-order states. Every state is labelled
#'     up to 15 states, otherwise the most common and most persistent are.
#'     Drawn with ggplot2: cograph has no scatter or quadrant view.}
#'   \item{`what = "states"`}{\emph{How does each state compare with the
#'     rest?} [cograph::plot_centrality()] on the per-state table, one panel
#'     per metric, states ordered by stationary share. The `k` largest values
#'     in each panel are highlighted, where `k` is the number of states whose
#'     share exceeds the even share \eqn{1/n}.}
#'   \item{`what = "passage_time"`}{\emph{How far apart are the states?}
#'     [cograph::plot_heatmap()] on the mean-first-passage matrix, rows and
#'     columns in stationary-share order, values printed for up to 12
#'     states. The diagonal is the mean return time \eqn{1/\pi_i}. The
#'     subtitle names the fastest and the slowest off-diagonal passage.}
#'   \item{`what = "network"`}{\emph{What does the walk look like?}
#'     [cograph::plot_tna()], cograph's TNA network view, on the transition
#'     matrix: node size from the stationary share, the node's pie (donut)
#'     filled to its persistence, edges the transition probabilities (below
#'     0.05 hidden). By default only the 10 most common states are drawn; the
#'     title says how many were left out.}
#' }
#'
#' `top` restricts every view to the `top` states with the largest
#' stationary share. The landscape's reference lines are computed over the
#' whole chain, so a state's quadrant does not change with `top`.
#'
#' @param x A `net_markov_stability` object.
#' @param what Which view: `"landscape"` (default), `"states"`,
#'   `"passage_time"` or `"network"`. Supplying `metrics` without `what`
#'   selects `"states"`, the view `metrics` belongs to.
#' @param metrics Character vector, any subset of `"persistence"`,
#'   `"stationary_prob"`, `"return_time"`, `"sojourn_time"`,
#'   `"avg_time_to_others"`, `"avg_time_from_others"`: the panels of the
#'   `"states"` view. Default: all six. Other views refuse it.
#' @param top Integer or `NULL`. Show only the `top` states with the largest
#'   stationary share. Default `NULL`: every state, except the network view,
#'   which draws at most 10.
#' @param from,to `NULL` or character vectors of states: restrict the
#'   `"passage_time"` heatmap to passages leaving `from` / arriving at `to`
#'   (applied together with `top`). Other views refuse them.
#' @param ... Passed to the cograph function behind the view
#'   ([cograph::plot_centrality()], [cograph::plot_heatmap()] or
#'   [cograph::plot_tna()]), overriding the defaults set here. Ignored by
#'   the landscape.
#' @return
#' \describe{
#'   \item{`"landscape"`}{a `ggplot`; its data has one row per state with
#'     `state`, `stationary_prob`, `persistence`, `sojourn_time`, `memory`
#'     (`"memory state"` / `"first-order state"`) and `quadrant`.}
#'   \item{`"states"`}{the `ggplot` from [cograph::plot_centrality()]; its
#'     data has one row per state x metric.}
#'   \item{`"passage_time"`}{the `ggplot` from [cograph::plot_heatmap()]; its
#'     data has one row per ordered pair shown (`row` = from, `col` = to,
#'     `value` = mean steps).}
#'   \item{`"network"`}{`x`, invisibly. [cograph::plot_tna()] draws with
#'     base graphics and has no ggplot to return.}
#' }
#' @section Conditions:
#' Raises `hypernets_bad_input` for an unknown `what`, for `metrics` with a
#' view other than `"states"`, for `from`/`to` with a view other than
#' `"passage_time"`, and for `from`/`to` naming unknown states.
#' @examples
#' P <- matrix(c(0.7, 0.2, 0.1,
#'               0.3, 0.5, 0.2,
#'               0.2, 0.3, 0.5),
#'             nrow = 3, byrow = TRUE,
#'             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' ms <- markov_stability(P)
#' plot(ms)                           # landscape
#' plot(ms, what = "states")          # cograph::plot_centrality()
#' plot(ms, what = "passage_time")    # cograph::plot_heatmap()
#' plot(ms, what = "network")         # cograph::plot_tna()
#'
#' # A second-order chain: memory states are drawn as triangles
#' hon <- build_hon(split(human_long$code, human_long$session_id),
#'                  max_order = 2L, min_freq = 50L)
#' plot(markov_stability(hon), top = 10L)
#' @export
plot.net_markov_stability <- function(x,
                                      what = c("landscape", "states",
                                               "passage_time", "network"),
                                      metrics = c("persistence",
                                                  "stationary_prob",
                                                  "return_time",
                                                  "sojourn_time",
                                                  "avg_time_to_others",
                                                  "avg_time_from_others"),
                                      top = NULL,
                                      from = NULL, to = NULL,
                                      ...) {
  views <- c("landscape", "states", "passage_time", "network")
  metrics_given <- !missing(metrics)
  if (missing(what)) {
    what <- if (metrics_given) "states" else "landscape"
  }
  if (!is.character(what) || length(what) != 1L || !(what %in% views)) {
    stop(errorCondition(
      sprintf("`what` must be one of %s.",
              paste0("\"", views, "\"", collapse = ", ")),
      class = "hypernets_bad_input", call = NULL))
  }
  if (metrics_given && !identical(what, "states")) {
    stop(errorCondition(
      sprintf("`metrics` chooses the panels of what = \"states\"; what = \"%s\" has none.",
              what),
      class = "hypernets_bad_input", call = NULL))
  }
  if ((!is.null(from) || !is.null(to)) && !identical(what, "passage_time")) {
    stop(errorCondition(
      sprintf("`from` and `to` restrict what = \"passage_time\"; what = \"%s\" does not take them.",
              what),
      class = "hypernets_bad_input", call = NULL))
  }
  switch(what,
         landscape    = .hms_plot_landscape(x, top),
         states       = .hms_plot_states(x, metrics, top, ...),
         passage_time = .hms_plot_passage(x, top, from, to, ...),
         network      = .hms_plot_network(x, top, ...))
}


# ---- plot helpers ----------------------------------------------------------

#' Which state labels carry memory (higher-order states)
#' @param states Character vector of state labels.
#' @return Logical vector.
#' @noRd
.hms_is_memory <- function(states) grepl(" -> ", states, fixed = TRUE)

#' Chain order and state counts, for plot titles
#' @param x A `net_markov_stability`.
#' @return A list with `order_word` and `counts` strings.
#' @noRd
.hms_describe_chain <- function(x) {
  memory <- .hms_is_memory(x$states)
  order <- max(lengths(strsplit(x$states, " -> ", fixed = TRUE)))
  order_word <- switch(as.character(order), `1` = "first-order",
                       `2` = "second-order", `3` = "third-order",
                       sprintf("order-%d", order))
  list(
    order_word = order_word,
    counts = if (any(memory)) {
      sprintf("%d states (%d carry memory)", x$n_states, sum(memory))
    } else {
      sprintf("%d states", x$n_states)
    }
  )
}

#' States in stationary-share order, cut to `top`
#' @param x A `net_markov_stability`.
#' @param top Integer or `NULL`.
#' @return The states data.frame.
#' @noRd
.hms_plot_table <- function(x, top) {
  as.data.frame(x, sort_by = "stationary_prob", top = top)
}

#' Quadrant of each state in the stability landscape
#'
#' Common = stationary share strictly above the even share `1/n`; sticky =
#' persistence strictly above the chain's mean persistence. Strict
#' inequalities keep a chain of identical states out of the upper quadrants.
#' @param share,persistence Numeric vectors.
#' @param n_states Number of states in the whole chain.
#' @param mean_persistence Mean persistence over the whole chain.
#' @return Factor with levels hub, relay, trap, transient.
#' @noRd
.hms_quadrant <- function(share, persistence, n_states, mean_persistence) {
  common <- share > 1 / n_states
  sticky <- persistence > mean_persistence
  q <- ifelse(common, ifelse(sticky, "hub", "relay"),
              ifelse(sticky, "trap", "transient"))
  factor(q, levels = c("hub", "relay", "trap", "transient"))
}

.hms_quadrant_palette <- c(hub = "#D55E00", relay = "#E69F00",
                           trap = "#0072B2", transient = "#999999")

#' Stability landscape: share x persistence, four quadrants
#'
#' The one view drawn here: cograph has no scatter or quadrant plot, and its
#' `theme_cograph_*()` objects are network themes, not ggplot themes.
#' @noRd
.hms_plot_landscape <- function(x, top) {
  d <- .hms_plot_table(x, top)
  n_states <- x$n_states
  even_share <- 1 / n_states
  mean_persistence <- mean(x$stability$persistence)
  d$memory <- factor(ifelse(.hms_is_memory(d$state), "memory state",
                            "first-order state"),
                     levels = c("first-order state", "memory state"))
  d$quadrant <- .hms_quadrant(d$stationary_prob, d$persistence, n_states,
                              mean_persistence)
  # an absorbing state has Inf sojourn; size it as the largest finite stay
  finite_stay <- d$sojourn_time[is.finite(d$sojourn_time)]
  d$stay <- pmin(d$sojourn_time, if (length(finite_stay)) max(finite_stay) else 1)

  # every label when few; otherwise the most common and the most persistent.
  # Rows are in share order, so check_overlap keeps the more common label.
  labelled <- nrow(d) <= 15L |
    seq_len(nrow(d)) <= 6L |
    (rank(-d$persistence, ties.method = "first") <= 4L & d$persistence > 0)
  label_df <- d[labelled, , drop = FALSE]

  # explicit limits: a log axis cannot place a corner label at -Inf
  x_range <- range(c(d$stationary_prob, even_share))
  x_lim <- exp(log(x_range) +
                 c(-1, 1) * max(diff(log(x_range)), log(4)) * 0.12)
  y_range <- range(c(d$persistence, mean_persistence))
  y_pad <- max(diff(y_range), 0.05) * 0.12
  y_lim <- c(y_range[1L] - 2.5 * y_pad, y_range[2L] + 3.5 * y_pad)
  corners <- data.frame(
    x = x_lim[c(2L, 2L, 1L, 1L)], y = y_lim[c(2L, 1L, 2L, 1L)],
    hjust = c(1, 1, 0, 0), vjust = c(1, 0, 1, 0),
    label = c("HUBS\ncommon and sticky", "RELAYS\ncommon, quickly left",
              "TRAPS\nrare but sticky", "TRANSIENTS\nrare, quickly left"),
    quadrant = factor(names(.hms_quadrant_palette),
                      levels = names(.hms_quadrant_palette)))
  chain <- .hms_describe_chain(x)

  ggplot2::ggplot(d, ggplot2::aes(x = .data$stationary_prob,
                                  y = .data$persistence)) +
    ggplot2::geom_vline(xintercept = even_share, linetype = "dashed",
                        colour = "#999999") +
    ggplot2::geom_hline(yintercept = mean_persistence, linetype = "dashed",
                        colour = "#999999") +
    ggplot2::geom_text(data = corners,
                       ggplot2::aes(x = .data$x, y = .data$y,
                                    label = .data$label,
                                    colour = .data$quadrant,
                                    hjust = .data$hjust, vjust = .data$vjust),
                       size = 3.2, fontface = "bold", lineheight = 0.9,
                       inherit.aes = FALSE) +
    ggplot2::geom_point(ggplot2::aes(size = .data$stay, shape = .data$memory,
                                     fill = .data$quadrant),
                        colour = "#000000", stroke = 0.5) +
    ggplot2::geom_text(data = label_df, ggplot2::aes(label = .data$state),
                       vjust = 2, size = 3.1, check_overlap = TRUE) +
    ggplot2::scale_x_log10(
      limits = x_lim, expand = ggplot2::expansion(mult = 0.01),
      labels = function(v) sprintf("%g%%", signif(100 * v, 2))) +
    ggplot2::scale_y_continuous(
      limits = y_lim, expand = ggplot2::expansion(mult = 0.01),
      breaks = function(l) pretty(pmax(l, 0))) +
    ggplot2::scale_shape_manual(values = c("first-order state" = 21,
                                           "memory state" = 24),
                                name = NULL) +
    ggplot2::scale_fill_manual(values = .hms_quadrant_palette,
                               guide = "none") +
    ggplot2::scale_colour_manual(values = .hms_quadrant_palette,
                                 guide = "none") +
    ggplot2::scale_size_continuous(range = c(2.5, 8),
                                   name = "Mean stay\n(steps)") +
    ggplot2::labs(
      x = sprintf("Stationary share (log scale; dashed = even share 1/%d)",
                  n_states),
      y = "Persistence (probability of repeating)",
      title = sprintf("Stability landscape of a %s chain", chain$order_word),
      subtitle = paste0(
        "Right of the dashed line: visited more than an even share. ",
        "Above it: more likely than average to repeat itself.\n",
        "Bigger point = longer stay. ", chain$counts,
        if (nrow(d) < n_states) sprintf("; the %d most common shown", nrow(d)),
        ".")) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   plot.title = ggplot2::element_text(face = "bold"),
                   plot.subtitle = ggplot2::element_text(size = 10))
}

#' Per-state metrics through cograph::plot_centrality()
#' @noRd
.hms_plot_states <- function(x, metrics, top, ...) {
  metrics <- match.arg(metrics,
                       c("persistence", "stationary_prob", "return_time",
                         "sojourn_time", "avg_time_to_others",
                         "avg_time_from_others"),
                       several.ok = TRUE)
  d <- .hms_plot_table(x, top)
  table <- data.frame(node = d$state, d[metrics], check.names = FALSE)
  chain <- .hms_describe_chain(x)
  defaults <- list(
    x = table, style = "dot",
    order_by = if ("stationary_prob" %in% metrics) "stationary_prob" else NULL,
    highlight = sum(x$stability$stationary_prob > 1 / x$n_states),
    title = sprintf("State profile of a %s chain", chain$order_word),
    subtitle = sprintf("%s%s; highlighted: the largest values in each panel",
                       chain$counts,
                       if (nrow(d) < x$n_states) sprintf(", %d shown", nrow(d)) else ""))
  do.call(cograph::plot_centrality, utils::modifyList(defaults, list(...)))
}

#' Mean first-passage heatmap through cograph::plot_heatmap()
#' @noRd
.hms_plot_passage <- function(x, top, from, to, ...) {
  .hms_check_states(from, "from", x$states)
  .hms_check_states(to, "to", x$states)
  keep <- .hms_plot_table(x, top)$state
  rows <- if (is.null(from)) keep else keep[keep %in% from]
  cols <- if (is.null(to)) keep else keep[keep %in% to]
  if (!length(rows) || !length(cols)) {
    stop(errorCondition(
      "No passages left to draw: `from`/`to` name no state among the `top` shown.",
      class = "hypernets_bad_input", call = NULL))
  }
  M <- x$passage_time[rows, cols, drop = FALSE]
  pairs <- as.data.frame(x, what = "passage_time", from = rows, to = cols)
  off <- pairs[pairs$from != pairs$to, , drop = FALSE]
  extremes <- if (nrow(off)) {
    fast <- off[which.min(off$steps), ]
    slow <- off[which.max(off$steps), ]
    sprintf("\nFastest: %s to %s, %.1f steps. Slowest: %s to %s, %.1f steps.",
            fast$from, fast$to, fast$steps, slow$from, slow$to, slow$steps)
  } else {
    ""
  }
  chain <- .hms_describe_chain(x)
  defaults <- list(
    x = M, colors = c("#FFFFFF", "#4A6FE3"), legend_title = "Mean steps",
    show_values = max(dim(M)) <= 12L, value_digits = 1L,
    xlab = "to", ylab = "from",
    title = sprintf("Mean first-passage time in a %s chain",
                    chain$order_word),
    subtitle = paste0("Steps from the row state to the column state; ",
                      "the diagonal is the mean return time.", extremes))
  do.call(cograph::plot_heatmap, utils::modifyList(defaults, list(...)))
}

#' Transition network through cograph::plot_tna()
#' @noRd
.hms_plot_network <- function(x, top, ...) {
  P <- x$transition
  if (is.null(P)) {
    stop(errorCondition(
      paste0("This net_markov_stability object predates hypernets 0.5.1 and ",
             "carries no transition matrix; recompute it with ",
             "markov_stability()."),
      class = "hypernets_bad_input", call = NULL))
  }
  d <- .hms_plot_table(x, top %||% min(x$n_states, 10L))
  hidden <- x$n_states - nrow(d)
  chain <- .hms_describe_chain(x)
  defaults <- list(
    x = P[d$state, d$state, drop = FALSE],
    vsize = 3 + 5 * sqrt(d$stationary_prob / max(d$stationary_prob)),
    pie = d$persistence, minimum = 0.05, label_size = 0.7,
    title = if (hidden > 0L) {
      sprintf("Top %d of %d states by share (%d hidden)",
              nrow(d), x$n_states, hidden)
    } else {
      sprintf("%s chain, %d states", chain$order_word, x$n_states)
    })
  do.call(cograph::plot_tna, utils::modifyList(defaults, list(...)))
  invisible(x)
}
