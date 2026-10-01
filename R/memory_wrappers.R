# ---- Memory-network constructors and tests ----
#
# hon(), honem(), mogen(), hypa() on sequences, markov_order(),
# memory() and hg_markov_stability() compute nothing: each reads
# its input through .ho_sequence_input() (long / wide / list / model
# objects), calls the Nestimate estimator, and returns Nestimate's object
# unchanged. Every number is therefore Nestimate's; identity with the direct
# call is asserted in tests/testthat/test-memory_wrappers.R.

#' Sequence input shared by the sequence verbs
#'
#' @name sequence-input
#' @section Sequence input:
#' Every verb that reads sequences reads them the same way:
#' \describe{
#'   \item{Long event table}{a data.frame with one event per row. `action`
#'     names the state column, `actor` the column (or columns) that
#'     identify whose events they are, `session` the column (or columns)
#'     that split an actor's events into sessions, and `time` the column
#'     that orders them: timestamps in any common date-time format, Unix
#'     times, or positions; row order when `NULL`. A column named
#'     `action`, `time`, `session` or `session_id` (in any case) is used
#'     when its argument is `NULL`; `session = FALSE` switches the session
#'     detection off. When `time` is given, a gap of more than
#'     `time_threshold` seconds between consecutive events starts a new
#'     sequence (default `900`, fifteen minutes); `time_threshold = FALSE`
#'     keeps each actor or session in one sequence however long the gaps.
#'     `timezone` is the time zone of timestamps that carry none (default
#'     `"UTC"`). Events with a missing actor or session raise
#'     `hypernets_bad_input`; a missing action stays in its sequence as a
#'     gap. Without `actor` all events form one sequence, announced by the
#'     message `hypernets_single_sequence`.}
#'   \item{Wide data.frame or character matrix}{one sequence per row;
#'     trailing `NA`s end a sequence.}
#'   \item{List}{one character vector per sequence.}
#'   \item{Model object}{a `netobject`, `netobject_group`, `tna` or
#'     `cograph_network` that carries its sequence data.}
#' }
#' These are the conventions of the tna family of packages, and the same
#' call builds the same sequences there.
#'
#' A data.frame with columns named like an event table (`code`, `state`,
#' `user`, `timestamp`, ...) that is passed without `action =` and has no
#' `action` column raises `hypernets_long_format` (a
#' `hypernets_bad_input`): read as wide, its actor ids and times would
#' silently become states. `actor`, `time` or `session` without an action
#' column raises `hypernets_bad_input`.
#' @keywords internal
NULL

#' Higher-order network (HON) from sequences
#'
#' Builds a higher-order network: a state is split into memory nodes
#' (`"a -> b"`: state `b` reached from `a`) wherever the next-state
#' distribution depends on where the sequence came from, as decided by the
#' BuildHON / BuildHON+ rule extraction.
#'
#' @param data Sequences in any form described in [sequence-input]: a long
#'   event table (with `action`), a wide data.frame, a list of vectors, or a
#'   model object carrying its sequences.
#' @param max_order Integer. Highest order considered. Default `5`.
#' @param min_freq Integer. Transitions observed fewer times are treated as
#'   zero. Default `1`.
#' @param collapse_repeats Logical. Collapse adjacent repeated states before
#'   extraction. Default `FALSE`.
#' @param method `"hon+"` (default, parameter-free BuildHON+) or `"hon"`
#'   (the original eager BuildHON).
#' @param action,actor,time,session,time_threshold,timezone Long-format
#'   arguments (see [sequence-input]); leave the column names `NULL` for
#'   wide, list or model input.
#' @return A `net_hon` object (also a `cograph_network`, so it plots
#'   directly). Read its tables with [hg_get()]: the rules
#'   (`what = "rules"`, the default), the memory nodes (`"nodes"`) and the
#'   higher-order pathways (`"pathways"`).
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances}, 2(5),
#' e1600028. \doi{10.1126/sciadv.1600028}
#'
#' Saebi, M., Xu, J., Kaplan, L. M., Ribeiro, B., & Chawla, N. V. (2020).
#' Efficient modeling of higher-order dependencies in networks: from
#' algorithm to application for anomaly detection. \emph{EPJ Data Science},
#' 9, 15. \doi{10.1140/epjds/s13688-020-00233-y}
#' @param group `NULL` (default) for one network, or the grouping of the
#'   sequences for a group model: the name of a column of `data` (long or
#'   wide format), or a vector with one label per sequence (per row of a
#'   wide data frame, per element of a list). A group model holds one
#'   network per group and the data of each group; it is the input of
#'   [hg_compare()] and can be passed to [hg_bootstrap()] and [hg_get()].
#' @seealso [honem()], [hg_centrality()], [hg_communities()],
#'   [hg_bootstrap()], [hg_compare()]
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' net <- hon(seqs, max_order = 2)
#' hg_get(net)
#'
#' # the same network from a long event table
#' events <- data.frame(
#'   who  = rep(c("s1", "s2", "s3", "s4"), each = 6),
#'   step = rep(1:6, 4),
#'   what = unlist(seqs)
#' )
#' hg_get(hon(events, action = "what", actor = "who", time = "step",
#'            max_order = 2), what = "nodes")
#' @export
hon <- function(data, max_order = 5L, min_freq = 1L,
                collapse_repeats = FALSE, method = "hon+",
                action = NULL, actor = NULL, time = NULL, session = NULL,
                time_threshold = 900, timezone = "UTC", group = NULL) {
  if (!is.null(group)) {
    return(.hon_group(data, group = group, max_order = max_order,
                      min_freq = min_freq,
                      collapse_repeats = collapse_repeats, method = method,
                      action = action, actor = actor, time = time,
                      session = session, time_threshold = time_threshold,
                      timezone = timezone))
  }
  data <- .ho_sequence_input(data, action = action, actor = actor,
                             time = time, session = session,
                             time_threshold = time_threshold,
                             timezone = timezone)
  .ho_result(Nestimate::build_hon(data, max_order = max_order,
                                  min_freq = min_freq,
                                  collapse_repeats = collapse_repeats,
                                  method = method))
}

#' Embedding of a higher-order network (HONEM)
#'
#' Low-dimensional node embedding of a higher-order network that preserves
#' its higher-order dependencies: exponentially decaying powers of the
#' network's transition matrix, followed by a truncated singular value
#' decomposition.
#'
#' @param hon A `net_hon` from [hon()], or a square weighted adjacency
#'   matrix.
#' @param dim Integer. Embedding dimension, capped at the number of nodes
#'   minus one. Default `32`.
#' @param max_power Integer. Longest walk length in the neighbourhood
#'   matrix. Default `10`.
#' @return A `net_honem` object. Read it with [hg_get()]: the node
#'   coordinates (`what = "embeddings"`, the default) or the variance per
#'   dimension (`"variance"`).
#' @references
#' Saebi, M., Ciampaglia, G. L., Kaplan, L. M., & Chawla, N. V. (2020).
#' HONEM: Learning embedding for higher order networks. \emph{Big Data},
#' 8(4), 255--269. \doi{10.1089/big.2019.0169}
#' @seealso [hon()]
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' emb <- honem(hon(seqs, max_order = 2), dim = 2)
#' hg_get(emb)
#' @export
honem <- function(hon, dim = 32L, max_power = 10L) {
  .ho_result(Nestimate::build_honem(.ho_unresult(hon), dim = dim,
                                    max_power = max_power))
}

#' Multi-order generative model (MOGen) of sequences
#'
#' Fits Markov chains of every order from 0 to `max_order` as layers of one
#' multi-order model (higher-order De Bruijn graphs) and selects the order
#' that best balances fit and parsimony.
#'
#' @inheritParams hon
#' @param max_order Integer. Highest Markov order tested. Default `5`;
#'   capped at the longest sequence minus one, with a message.
#' @param criterion `"aic"` (default), `"bic"` or `"lrt"` (likelihood-ratio
#'   test).
#' @param lrt_alpha Significance level of the likelihood-ratio test.
#'   Default `0.01`.
#' @return A `net_mogen` object (also a `cograph_network`). Read it with
#'   [hg_get()]: the order-selection table (`what = "orders"`, the default),
#'   the fitted transitions of a layer (`"transitions"`), the path counts
#'   (`"paths"`, with `k =`) and the pathways (`"pathways"`).
#' @references
#' Scholtes, I. (2017). When is a network a network? Multi-order graphical
#' model selection in pathways and temporal networks. \emph{Proceedings of
#' KDD 2017}, 1037--1046. \doi{10.1145/3097983.3098145}
#'
#' Gote, C., Casiraghi, G., Schweitzer, F., & Scholtes, I. (2023). Predicting
#' variable-length paths in networked systems using multi-order generative
#' models. \emph{Applied Network Science}, 8, 68.
#' \doi{10.1007/s41109-023-00596-0}
#' @seealso [markov_order()], [memory()]
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' mg <- mogen(seqs, max_order = 2)
#' hg_get(mg)
#' hg_get(mg, what = "paths", k = 3)
#' @export
mogen <- function(data, max_order = 5L, criterion = c("aic", "bic", "lrt"),
                  lrt_alpha = 0.01,
                  action = NULL, actor = NULL, time = NULL, session = NULL,
                  time_threshold = 900, timezone = "UTC") {
  data <- .ho_sequence_input(data, action = action, actor = actor,
                             time = time, session = session,
                             time_threshold = time_threshold,
                             timezone = timezone)
  .ho_result(Nestimate::build_mogen(data, max_order = max_order,
                                    criterion = criterion,
                                    lrt_alpha = lrt_alpha))
}

#' Test the Markov order of sequences
#'
#' At each order k the likelihood-ratio statistic G2 tests whether the k-th
#' state back adds information about the next state given the k - 1 most
#' recent ones; its null distribution comes from an exact within-context
#' permutation. The order is raised while the test rejects.
#'
#' @inheritParams hon
#' @param data Sequences in any form described in [sequence-input]; a
#'   `netobject_group` gives one test per group.
#' @param max_order Integer. Highest order tested. Default `3`.
#' @param n_perm Integer. Permutations per order. Default `500`.
#' @param alpha Significance level of the sequential selection. Default
#'   `0.05`.
#' @param parallel Logical. Run the permutations with
#'   `parallel::mclapply()` (not on Windows). Default `FALSE`.
#' @param n_cores Integer. Cores when `parallel = TRUE`. Default `2`.
#' @param seed Integer or `NULL`. Seed for the permutations.
#' @return A `net_markov_order` object (a `net_markov_order_group` for a
#'   `netobject_group`). Read it with [hg_get()]: the per-order test table
#'   (`what = "orders"`, the default) or the permutation null draws
#'   (`"null"`).
#' @references
#' Anderson, T. W., & Goodman, L. A. (1957). Statistical inference about
#' Markov chains. \emph{The Annals of Mathematical Statistics}, 28(1),
#' 89--110. \doi{10.1214/aoms/1177707039}
#'
#' Good, P. (2005). \emph{Permutation, Parametric and Bootstrap Tests of
#' Hypotheses} (3rd ed.). Springer.
#' @seealso [mogen()], [memory()]
#' @examples
#' seqs <- list(c("a", "b", "a", "c", "a", "b"), c("b", "a", "c", "a", "b", "a"))
#' fit <- markov_order(seqs, max_order = 2L, n_perm = 20L, seed = 1L)
#' hg_get(fit)
#' @export
markov_order <- function(data, max_order = 3L, n_perm = 500L,
                            alpha = 0.05, parallel = FALSE, n_cores = 2L,
                            seed = NULL, action = NULL, actor = NULL,
                            time = NULL, session = NULL,
                            time_threshold = 900, timezone = "UTC") {
  data <- .ho_sequence_input(data, action = action, actor = actor,
                             time = time, session = session,
                             time_threshold = time_threshold,
                             timezone = timezone,
                             models = "decode")
  .ho_result(Nestimate::markov_order_test(data, max_order = max_order,
                                          n_perm = n_perm, alpha = alpha,
                                          parallel = parallel,
                                          n_cores = n_cores, seed = seed))
}

#' Per-context path dependence
#'
#' For every context of `order - 1` states, compares the observed
#' next-state distribution with the first-order prediction from the most
#' recent state alone: the Kullback-Leibler divergence and entropy drop
#' say how much the longer history adds, and `flips` marks the contexts
#' where it changes the most likely next state.
#'
#' @inheritParams hon
#' @param order Integer >= 2. Length of the conditioning context plus one:
#'   `2` compares two-step memory with one-step. Default `2`.
#' @param min_count Integer. Contexts observed fewer times are dropped.
#'   Default `5`.
#' @param base Logarithm base of the entropies and divergences. Default `2`
#'   (bits).
#' @return A `net_path_dependence` object. Read it with [hg_get()]: one row
#'   per context (`what = "contexts"`), sorted by divergence, largest first.
#' @references
#' Kullback, S., & Leibler, R. A. (1951). On information and sufficiency.
#' \emph{The Annals of Mathematical Statistics}, 22(1), 79--86.
#' \doi{10.1214/aoms/1177729694}
#'
#' Cover, T. M., & Thomas, J. A. (2006). \emph{Elements of Information
#' Theory} (2nd ed.). Wiley.
#' @seealso [markov_order()], [mogen()]
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' pd <- memory(seqs, order = 2, min_count = 1)
#' hg_get(pd)
#' @export
memory <- function(data, order = 2L, min_count = 5L, base = 2,
                               action = NULL, actor = NULL, time = NULL,
                               session = NULL,
                               time_threshold = 900, timezone = "UTC") {
  data <- .ho_sequence_input(data, action = action, actor = actor,
                             time = time, session = session,
                             time_threshold = time_threshold,
                             timezone = timezone,
                             lists = "wide", models = "decode")
  .ho_result(Nestimate::path_dependence(data, order = order,
                                        min_count = min_count, base = base))
}

#' Markov stability of the states of a transition network
#'
#' Per-state persistence (self-transition probability), stationary
#' probability, return time, sojourn time and mean first-passage times to
#' and from the other states of a first-order Markov chain.
#'
#' @inheritParams hon
#' @param data A row-stochastic transition matrix, a network object
#'   (`netobject`, `tna`, `cograph_network`, `netobject_group`) whose
#'   weights are the transition matrix, or sequences in any form described
#'   in [sequence-input], from which the first-order transition
#'   probabilities are estimated.
#' @param normalize Logical. Normalise the rows to sum to one. Default
#'   `TRUE`.
#' @return A `net_markov_stability` object. Read it with [hg_get()]: one row
#'   per state (`what = "states"`, the default), the mean first-passage
#'   times (`"passage_time"`) or the stationary distribution
#'   (`"stationary"`).
#' @section Conditions:
#' The quantities are defined only for an irreducible chain, one in which
#' every state can reach every other. A chain that is not irreducible (a
#' transient or absorbing state, or several closed classes, as a pruned
#' higher-order network can have) raises `hypernets_not_ergodic`, a
#' `hypernets_bad_input`, naming those states.
#' @references
#' Kemeny, J. G., & Snell, J. L. (1976). \emph{Finite Markov Chains}.
#' Springer.
#' @seealso [markov_order()]
#' @examples
#' P <- matrix(c(0.5, 0.3, 0.2,
#'               0.2, 0.5, 0.3,
#'               0.1, 0.4, 0.5), 3, 3, byrow = TRUE,
#'             dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
#' ms <- hg_markov_stability(P)
#' hg_get(ms, sort_by = "stationary")
#' hg_get(ms, what = "passage_time", from = "A")
#' @export
hg_markov_stability <- function(data, normalize = TRUE,
                                action = NULL, actor = NULL, time = NULL,
                                session = NULL,
                                time_threshold = 900, timezone = "UTC") {
  data <- .ho_sequence_input(data, action = action, actor = actor,
                             time = time, session = session,
                             time_threshold = time_threshold,
                             timezone = timezone,
                             lists = "wide")
  data <- .ho_unresult(data)
  .hms_check_irreducible(.hms_chain_structure(data))
  .ho_result(Nestimate::markov_stability(data, normalize = normalize))
}

# The chain's structure; a state with no outgoing transition (a dead end a
# pruned network can leave) is reported by the estimator with a plain error,
# raised here as hypernets_not_ergodic. Every other error passes unchanged.
.hms_chain_structure <- function(data) {
  tryCatch(Nestimate::chain_structure(data), error = function(e) {
    message <- conditionMessage(e)
    if (!grepl("no outgoing transitions", message, fixed = TRUE)) stop(e)
    dead_ends <- sub(".*for state\\(s\\): ([^.]*)\\..*", "\\1", message)
    stop(errorCondition(sprintf(paste0(
      "the chain is not irreducible (no outgoing transition from: %s), so ",
      "its stationary distribution, return times and passage times are not ",
      "defined. Lower `min_freq` or `max_order` of the network."),
      dead_ends), class = c("hypernets_not_ergodic", "hypernets_bad_input"),
      call = NULL))
  })
}

# The stationary distribution, return and passage times exist only for an
# irreducible chain (Kemeny & Snell 1976); a chain that is not is refused
# with the states that break it.
.hms_check_irreducible <- function(structure) {
  if (isTRUE(structure$is_irreducible)) return(invisible(NULL))
  name_some <- function(states) {
    shown <- utils::head(states, 5L)
    paste0(paste(shown, collapse = ", "),
           if (length(states) > 5L) sprintf(" and %d more", length(states) - 5L))
  }
  kinds <- structure$classification
  transient <- names(kinds)[kinds == "transient"]
  absorbing <- structure$absorbing_states
  parts <- c(
    sprintf("%d recurrent classes", length(structure$recurrent_classes)),
    if (length(transient)) sprintf("transient: %s", name_some(transient)),
    if (length(absorbing)) sprintf("absorbing: %s", name_some(absorbing)))
  stop(errorCondition(sprintf(paste0(
    "the chain is not irreducible (%s), so its stationary distribution, ",
    "return times and passage times are not defined. Lower `min_freq` or ",
    "`max_order` of the network, or analyse one recurrent class."),
    paste(parts, collapse = "; ")),
    class = c("hypernets_not_ergodic", "hypernets_bad_input"), call = NULL))
}
