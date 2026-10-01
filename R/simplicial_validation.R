# ---- Statistically validated window sets ----
#
# simplicial(type = "window", validate = TRUE) keeps a set of actions only
# when it occurs in more windows than a null model predicts. Two nulls:
#
# "swap" (default): swap randomization of the windows-by-actions incidence
#   (Gionis et al. 2007), the checkerboard swap chain hg_null_test() uses
#   (Gotelli 2000). Every window keeps its number of actions and every action
#   its number of windows; only which actions share a window is randomized.
# "hypergeometric": the statistically validated hypergraph of Musciotto,
#   Battiston & Mantegna (2021). Sets of n actions are tested against the
#   windows with exactly n actions, and each action falls into them
#   independently of the others. That independence ignores that a window
#   holds exactly n actions, so on a small set of actions the null expects
#   more co-occurrence than the windows can hold.

# Distribution of the size X of the intersection of length(degrees) random
# subsets, of sizes `degrees`, of `total` items, on 0..min(degrees): the
# hierarchical convolution of hypergeometric intersections of Musciotto et
# al. (2021), Eq. 3. The intersection of the first j subsets meets subset
# j + 1 in a hypergeometric number of items; the support of the j-th
# intersection is 0..min(degrees[1..j]).
.svh_distribution <- function(degrees, total) {
  caps <- cummin(degrees)
  first <- stats::dhyper(0:caps[2L], degrees[1L], total - degrees[1L],
                         degrees[2L])
  Reduce(\(d, j) {
    step <- outer(0:caps[j - 1L], 0:caps[j],
                  \(x, y) stats::dhyper(y, x, total - x, degrees[j]))
    as.vector(d %*% step)
  }, seq_along(degrees)[-(1:2)], first)
}

# Upper tail P(X >= n_obs) of that distribution.
.svh_upper_tail <- function(n_obs, degrees, total) {
  dist <- .svh_distribution(degrees, total)
  support <- seq_along(dist) - 1L
  min(1, sum(dist[support >= n_obs]))
}

# Windows-by-actions incidence: one row per window (a set seen in k windows
# gives k rows), 1 where the action is in the window.
.svh_incidence <- function(sets, counts, n_nodes) {
  rows <- rep(seq_along(sets), counts)
  m <- matrix(0L, length(rows), n_nodes)
  m[cbind(rep(seq_along(rows), lengths(sets[rows])),
          as.integer(unlist(sets[rows])))] <- 1L
  m
}

# One key per row of a binary matrix.
.svh_row_keys <- function(m) {
  do.call(paste0, as.data.frame(m))
}

# Counts of the `targets` rows (keys) in `n_null` swap-randomized copies of
# the incidence `m`: a sequential chain with burn-in 10 * nnz swap attempts
# and nnz attempts between samples, as in hg_null_test(). Returns a
# length(targets) x n_null matrix.
.svh_swap_counts <- function(m, targets, n_null) {
  nnz <- sum(m)
  state <- .thg_swap_chain(m, attempts = 10L * nnz)
  vapply(seq_len(n_null), \(i) {
    state <<- .thg_swap_chain(state, attempts = nnz)
    tabulate(match(.svh_row_keys(state), targets), nbins = length(targets))
  }, numeric(length(targets)))
}

# Validation of the window sets. `sets` are integer node indices, `counts`
# the number of windows with exactly that set of actions. Sets of one action
# are not tested. p-values are corrected with the Benjamini-Hochberg false
# discovery rate within each set size, over all possible sets of that size:
# choose(n_active, n), n_active being the actions in windows of size n for
# the hypergeometric null (Musciotto et al. 2021, Methods) and the actions
# in any window for the swap null. Returns the table and the validated flag
# of every set, in the order of `sets`.
.sc_validate_windows <- function(sets, counts, nodes, alpha,
                                 null = c("swap", "hypergeometric"),
                                 n_null = NULL) {
  null <- match.arg(null)
  sizes <- lengths(sets)
  tested <- which(sizes >= 2L)
  max_size <- max(sizes)
  # windows of each size, and each action's number of windows of that size
  total <- vapply(sizes, \(n) sum(counts[sizes == n]), numeric(1L))
  degree_of <- lapply(seq_len(max_size), \(n) {
    idx <- which(sizes == n)
    tabulate(rep(as.integer(unlist(sets[idx])), rep(counts[idx], each = n)),
             nbins = length(nodes))
  })
  if (identical(null, "hypergeometric")) {
    n_tests <- vapply(seq_len(max_size),
                      \(n) choose(sum(degree_of[[n]] > 0), n), numeric(1L))
    moments <- vapply(tested, \(i) {
      dist <- .svh_distribution(degree_of[[sizes[i]]][sets[[i]]], total[i])
      support <- seq_along(dist) - 1L
      mean_x <- sum(support * dist)
      c(mean_x, sqrt(max(0, sum(support^2 * dist) - mean_x^2)),
        min(1, sum(dist[support >= counts[i]])))
    }, numeric(3L))
    # the null needs sparse co-occurrence: most possible sets rare or absent,
    # so each tested set is expected in less than one window. When the
    # median tested set is expected in a window or more, the null's
    # independence overstates co-occurrence by a factor the observed counts
    # rarely exceed.
    by_size <- split(seq_along(tested), sizes[tested])
    dense <- vapply(by_size, \(j) stats::median(moments[1L, j]) >= 1,
                    logical(1L))
    if (any(dense)) {
      inflation <- vapply(by_size[dense], \(j) {
        sum(moments[1L, j]) / total[tested[j[1L]]]
      }, numeric(1L))
      median_expected <- vapply(by_size[dense], \(j) {
        stats::median(moments[1L, j])
      }, numeric(1L))
      warning(warningCondition(sprintf(paste0(
        "the hypergeometric null (Musciotto et al. 2021) needs sparse ",
        "co-occurrence, each set expected in less than one window. Here the ",
        "median set is expected in %s windows and the expected counts sum ",
        "to %s times the windows (sizes %s), so few sets can pass. ",
        "null = \"swap\" keeps the number of actions in each window fixed."),
        paste(signif(median_expected, 3), collapse = ", "),
        paste(signif(inflation, 3), collapse = ", "),
        paste(names(by_size)[dense], collapse = ", ")),
        class = "hypernets_dense_cooccurrence", call = NULL))
    }
  } else {
    n_active <- sum(Reduce(`+`, degree_of) > 0)
    n_tests <- choose(n_active, seq_len(max_size))
    m <- .svh_incidence(sets, counts, length(nodes))
    target_rows <- .svh_incidence(sets[tested], rep(1L, length(tested)),
                                  length(nodes))
    # a set passes alone only if the smallest possible p-value clears the
    # first Benjamini-Hochberg step of its size: 1 / (n_null + 1) <=
    # alpha / n_tests
    sizes_tested <- sort(unique(sizes[tested]))
    needed <- if (length(tested) > 0L) {
      as.integer(ceiling(max(n_tests[sizes_tested]) / alpha))
    } else {
      0L
    }
    n_null <- as.integer(n_null %||% max(999L, needed))
    short <- sizes_tested[1 / (n_null + 1) > alpha / n_tests[sizes_tested]]
    if (length(short) > 0L) {
      warning(warningCondition(sprintf(paste0(
        "with n_null = %d the smallest p-value is %.3g, above the ",
        "Benjamini-Hochberg threshold alpha / n_tests for sets of size %s; ",
        "a set can then pass only together with many others. ",
        "n_null = %d resolves every size."),
        n_null, 1 / (n_null + 1), paste(short, collapse = ", "), needed),
        class = "hypernets_low_resolution", call = NULL))
    }
    draws <- .svh_swap_counts(m, .svh_row_keys(target_rows), n_null)
    draws <- matrix(draws, nrow = length(tested))
    exceed <- rowSums(draws >= counts[tested])
    # permutation p-value, never zero (Phipson & Smyth 2010)
    moments <- rbind(rowMeans(draws), apply(draws, 1L, stats::sd),
                     (1 + exceed) / (n_null + 1))
  }
  moments <- matrix(moments, nrow = 3L)
  expected <- moments[1L, ]
  null_sd <- moments[2L, ]
  p_value <- moments[3L, ]
  p_adj <- p_value
  if (length(tested) > 0L) {
    p_adj <- unsplit(lapply(split(seq_along(tested), sizes[tested]), \(j) {
      stats::p.adjust(p_value[j], method = "BH",
                      n = n_tests[sizes[tested[j[1L]]]])
    }), sizes[tested])
  }
  observed <- counts[tested]
  z <- rep(NA_real_, length(tested))
  spread <- null_sd > 0
  z[spread] <- (observed[spread] - expected[spread]) / null_sd[spread]
  table <- data.frame(
    members = vapply(sets[tested],
                     \(s) paste(nodes[sort(s)], collapse = ", "),
                 character(1L)),
    size = as.integer(sizes[tested]),
    count = as.integer(observed),
    n_windows = as.integer(total[tested]),
    expected = expected,
    z = z,
    p_value = p_value,
    p_adj = p_adj,
    n_tests = n_tests[sizes[tested]],
    significant = p_adj <= alpha,
    stringsAsFactors = FALSE
  )
  table <- table[order(table$size, table$p_value, -table$z, table$members), ,
                 drop = FALSE]
  rownames(table) <- NULL
  validated <- logical(length(sets))
  validated[tested] <- p_adj <= alpha
  list(table = table, validated = validated,
       n_null = if (identical(null, "swap")) n_null)
}
