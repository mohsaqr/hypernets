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
             "Remove the state(s) or supply a different transition matrix."),
      class = "hypernets_bad_input", call = NULL))
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
#' so an all-TRUE result is strong connectivity. Squaring is the justified
#' loop here — it runs `log2(n)` times, not `n` times.
#'
#' @param P Square numeric matrix.
#' @return Logical scalar.
#' @noRd
.hms_irreducible <- function(P) {
  n <- nrow(P)
  if (n == 1L) return(TRUE)
  reach <- (diag(n) + (P != 0)) != 0
  steps <- max(1L, ceiling(log2(n)))
  # repeated squaring: log2(n) iterations, each doubling the path length
  for (i in seq_len(steps)) {
    reach <- (reach %*% reach) != 0
    if (all(reach)) return(TRUE)
  }
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
#'   `summary()` names the attractor and the stickiest state, `plot()` draws
#'   the metrics as a faceted bar chart, and `print()` shows the table.
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
#' @section Conditions:
#' Raises `hypernets_bad_input` for an unsupported `x`, a non-square or
#' single-state matrix, missing values, a state with no outgoing transitions,
#' and for rows that do not sum to 1 under `normalize = FALSE`. Raises
#' `hypernets_not_ergodic` (which also inherits `hypernets_bad_input`) when
#' the chain is reducible. Warns `hypernets_renormalized` when rows are
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
#'   of the selected table; sorts decreasing.
#' @param top Integer or `NULL`. Keep only the first `top` rows, applied
#'   after `what` and `sort_by`.
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
                                               top = NULL) {
  what <- match.arg(what)
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

  if (!is.null(sort_by)) {
    stopifnot(
      "`sort_by` must be a single column name of the returned table" =
        is.character(sort_by) && length(sort_by) == 1L &&
        sort_by %in% names(out)
    )
    out <- out[order(out[[sort_by]], decreasing = TRUE), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
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
#' One horizontal bar panel per requested metric, states ordered by
#' stationary probability. Colour is the metric (Okabe-Ito, colour-blind
#' safe) and is redundant with the facet label, so no distinction rests on
#' colour alone.
#'
#' @param x A `net_markov_stability` object.
#' @param metrics Character vector, any subset of `"persistence"`,
#'   `"stationary_prob"`, `"return_time"`, `"sojourn_time"`,
#'   `"avg_time_to_others"`, `"avg_time_from_others"`. Default: all six.
#' @param top Integer or `NULL`. Plot only the `top` states with the largest
#'   stationary probability. Default `NULL` (all states).
#' @param ... Ignored.
#' @return A `ggplot` object, faceted by metric.
#' @inherit markov_stability examples
#' @export
plot.net_markov_stability <- function(x,
                                      metrics = c("persistence",
                                                  "stationary_prob",
                                                  "return_time",
                                                  "sojourn_time",
                                                  "avg_time_to_others",
                                                  "avg_time_from_others"),
                                      top = NULL,
                                      ...) {
  metrics <- match.arg(metrics, several.ok = TRUE)
  d <- as.data.frame(x, sort_by = "stationary_prob", top = top)

  labels <- c(persistence          = "Persistence",
              stationary_prob      = "Stationary probability",
              return_time          = "Return time (steps)",
              sojourn_time         = "Sojourn time (steps)",
              avg_time_to_others   = "Mean steps to others",
              avg_time_from_others = "Mean steps from others")
  pal <- c("#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#0072B2", "#D55E00")
  names(pal) <- labels[c("persistence", "stationary_prob", "return_time",
                         "sojourn_time", "avg_time_to_others",
                         "avg_time_from_others")]

  plot_df <- do.call(rbind, lapply(metrics, function(m) {
    data.frame(state  = d$state,
               metric = unname(labels[m]),
               value  = d[[m]],
               stringsAsFactors = FALSE)
  }))
  # an absorbing state has Inf sojourn time; drop it from the bars rather
  # than let ggplot fail on a non-finite value, and say so in the subtitle
  n_inf <- sum(!is.finite(plot_df$value))
  plot_df <- plot_df[is.finite(plot_df$value), , drop = FALSE]
  plot_df$state  <- factor(plot_df$state, levels = rev(d$state))
  plot_df$metric <- factor(plot_df$metric, levels = unname(labels[metrics]))

  ggplot2::ggplot(plot_df,
                  ggplot2::aes(x = .data$state, y = .data$value,
                               fill = .data$metric)) +
    ggplot2::geom_col(width = 0.7, show.legend = FALSE) +
    ggplot2::scale_fill_manual(values = pal) +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_x", ncol = 2L) +
    ggplot2::labs(
      x = NULL, y = NULL,
      title = "Markov stability",
      subtitle = if (n_inf > 0L) {
        sprintf("%d states; %d non-finite value(s) omitted (absorbing state)",
                nrow(d), n_inf)
      } else {
        sprintf("%d states, ordered by stationary probability", nrow(d))
      }) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.minor = ggplot2::element_blank(),
      strip.text       = ggplot2::element_text(face = "bold"),
      plot.title       = ggplot2::element_text(face = "bold")
    )
}
