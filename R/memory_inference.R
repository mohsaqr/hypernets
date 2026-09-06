# ---- Inference for HON rules: bootstrap CIs + two-sample comparison ------
#
# Both verbs follow the family bootstrap discipline: per-sequence
# observation counts are precomputed ONCE (.hi_seq_counts), and every
# bootstrap / permutation replicate is a weighted aggregation of those
# counts (.hi_count_env) followed by the counts-based rule extraction
# (.hon_extract_rules_count) - trajectories are never re-scanned. All
# randomness is drawn serially up front, so parallel and serial runs give
# identical results under the same seed.
#
# Rule extraction in replicates uses the eager counts-based path
# (method = "hon"); the shipped test suite pins that "hon" and "hon+"
# produce identical networks, so the results apply to build_hon() output
# of either method.

# ---------------------------------------------------------------------------
# Shared input parsing (wide / list / model objects / long format)
# ---------------------------------------------------------------------------

#' Parse inference input into a list of character trajectories
#'
#' Long format (action given): one trajectory per actor, ordered by time
#' (row order within actor when time is NULL). Otherwise delegates to the
#' package input contract (.coerce_sequence_input + .hon_parse_input).
#'
#' @param data Sequence data (wide data.frame, list, tna, netobject, or
#'   long data.frame with `action`).
#' @param action,actor,time Long-format column names or NULL.
#' @param collapse_repeats Logical, as in build_hon().
#' @return List of character trajectories (length >= 1).
#' @noRd
.hi_parse <- function(data, action = NULL, actor = NULL, time = NULL,
                      collapse_repeats = FALSE) {
  if (!is.null(action)) {
    stopifnot(
      "`data` must be a data.frame when `action` is given" =
        is.data.frame(data),
      "`action` must name a column of `data`" =
        is.character(action) && length(action) == 1L &&
        action %in% names(data),
      "`actor` must be NULL or name a column of `data`" =
        is.null(actor) ||
        (is.character(actor) && length(actor) == 1L && actor %in% names(data)),
      "`time` must be NULL or name a column of `data`" =
        is.null(time) ||
        (is.character(time) && length(time) == 1L && time %in% names(data))
    )
    a <- as.character(data[[action]])
    g <- if (is.null(actor)) rep("sequence_1", length(a)) else
      as.character(data[[actor]])
    keep <- !is.na(g)
    a <- a[keep]
    g <- g[keep]
    o <- if (is.null(time)) order(g) else order(g, data[[time]][keep])
    trajectories <- split(a[o], g[o])
  } else {
    stopifnot(
      "`actor` requires `action`" = is.null(actor),
      "`time` requires `action`"  = is.null(time)
    )
    data <- .coerce_sequence_input(data)
    trajectories <- .hon_parse_input(data,
                                     collapse_repeats = FALSE)
  }
  if (isTRUE(collapse_repeats)) {
    trajectories <- lapply(trajectories, function(traj) {
      if (length(traj) <= 1L) return(traj)
      traj[c(TRUE, traj[-1L] != traj[-length(traj)])]
    })
  }
  trajectories <- trajectories[lengths(trajectories) >= 2L]
  if (length(trajectories) == 0L) {
    stop("No valid trajectories (each must have at least 2 states).",
         call. = FALSE)
  }
  trajectories
}

# ---------------------------------------------------------------------------
# Per-sequence counts and weighted aggregation
# ---------------------------------------------------------------------------

#' Precompute per-sequence observation counts as a tidy triplet table
#'
#' @param trajectories List of character trajectories.
#' @param max_order Integer.
#' @return data.frame: seq (integer), source_key (encoded context),
#'   target (state), n (integer count).
#' @noRd
.hi_seq_counts <- function(trajectories, max_order) {
  per <- lapply(seq_along(trajectories), function(i) {
    env <- .hon_build_observations(trajectories[i], max_order)
    keys <- ls(env)
    if (length(keys) == 0L) return(NULL)
    counts <- lapply(keys, function(k) env[[k]])
    lens <- vapply(counts, length, integer(1L))
    data.frame(
      seq        = rep(i, sum(lens)),
      source_key = rep(keys, lens),
      target     = unlist(lapply(counts, names), use.names = FALSE),
      n          = unlist(counts, use.names = FALSE),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, per)
}

#' Aggregate per-sequence counts under sequence weights into a count env
#'
#' Equivalent to .hon_build_observations() on the weighted multiset of
#' trajectories (counts are additive), without re-scanning them. Returns a
#' fresh environment on every call because .hon_build_distributions()
#' mutates it in place.
#'
#' @param sc Triplet data.frame from .hi_seq_counts().
#' @param w Integer weights, one per sequence (0 = not drawn).
#' @return Environment: source_key -> named integer vector of counts.
#' @noRd
.hi_count_env <- function(sc, w) {
  wk <- w[sc$seq]
  keep <- wk > 0L
  if (!any(keep)) {
    return(new.env(hash = TRUE, parent = emptyenv()))
  }
  n2 <- sc$n[keep] * wk[keep]
  pair <- paste(sc$source_key[keep], sc$target[keep], sep = "\x03")
  agg <- rowsum(n2, pair, reorder = FALSE)
  parts <- strsplit(rownames(agg), "\x03", fixed = TRUE)
  src <- vapply(parts, `[[`, character(1L), 1L)
  tgt <- vapply(parts, `[[`, character(1L), 2L)
  vals <- as.integer(agg[, 1L])
  lst <- lapply(split(seq_along(src), src), function(idx) {
    stats::setNames(vals[idx], tgt[idx])
  })
  list2env(lst, hash = TRUE, parent = emptyenv())
}

#' Tidy rule table from an extraction result
#'
#' @param extraction list(rules, count) from .hon_extract_rules_count().
#' @return data.frame, one row per rule edge: from, to, order, count,
#'   probability - ordered by (order, from, to).
#' @noRd
.hi_rule_table <- function(extraction) {
  rules <- extraction$rules
  count <- extraction$count
  keys <- ls(rules)
  rows <- lapply(keys, function(k) {
    p <- rules[[k]]
    p <- p[p > 0]
    if (length(p) == 0L) return(NULL)
    data.frame(
      source_key = k,
      from   = .hon_sequence_to_node(.hon_decode(k)),
      to     = names(p),
      order  = .hon_key_len(k),
      count  = as.integer(count[[k]][names(p)]),
      probability = as.numeric(p),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  out <- out[order(out$order, out$from, out$to), , drop = FALSE]
  rownames(out) <- NULL
  out
}

#' Conditional probabilities of fixed (context, target) pairs under a
#' distribution environment
#'
#' @param distr Environment from .hon_build_distributions().
#' @param source_key,target Parallel character vectors.
#' @return Numeric vector: probability; 0 when the context is observed but
#'   the target is not retained; NA when the context is unobserved.
#' @noRd
.hi_lookup_prob <- function(distr, source_key, target) {
  vapply(seq_along(source_key), function(i) {
    p <- distr[[source_key[i]]]
    if (is.null(p)) return(NA_real_)
    out <- p[target[i]]
    if (is.na(out)) 0 else as.numeric(out)
  }, numeric(1L))
}

#' Run a replicate function serially or via parallel::mclapply
#' @noRd
.hi_apply <- function(n, fn, parallel, n_cores) {
  if (isTRUE(parallel) && .Platform$OS.type != "windows") {
    parallel::mclapply(seq_len(n), fn, mc.cores = n_cores)
  } else {
    lapply(seq_len(n), fn)
  }
}

# ---------------------------------------------------------------------------
# bootstrap_hon
# ---------------------------------------------------------------------------

#' Bootstrap inference for higher-order network rules
#'
#' Nonparametric bootstrap over sequences for the rules of a higher-order
#' network (see [build_hon()]): sequences are resampled with replacement,
#' and for every rule edge of the observed network the replicate
#' distribution yields a percentile confidence interval for its
#' conditional probability and a *support* - the fraction of replicates in
#' which the rule's context is extracted as a rule at all. Support close
#' to 1 for a higher-order context means the order elevation is stable
#' under resampling, not an artifact of the particular sample; first-order
#' contexts have support 1 by construction.
#'
#' Per-sequence observation counts are computed once; every replicate is a
#' weighted aggregation of those counts followed by counts-based rule
#' extraction (never a re-scan of the data). All randomness is drawn
#' before any parallel work, so `parallel = TRUE` reproduces the serial
#' result under the same `seed`.
#'
#' @param data Sequence data: wide data.frame (one sequence per row), list
#'   of vectors, `tna`/`netobject` model objects, or a long data.frame
#'   together with `action` (and optionally `actor`, `time`).
#' @param n_boot Integer >= 2. Bootstrap replicates. Default `500`.
#' @param level Confidence level in (0, 1). Default `0.95` (percentile
#'   interval).
#' @param max_order,min_freq,collapse_repeats As in [build_hon()].
#' @param action,actor,time Long-format column names: `action` holds the
#'   categorical state, `actor` groups events into sequences, `time`
#'   orders them within an actor (row order when `NULL`). Leave `NULL`
#'   for wide/list input.
#' @param parallel Logical. Use `parallel::mclapply` for the replicates
#'   (not on Windows). Results are identical to the serial run.
#' @param n_cores Integer. Cores when `parallel = TRUE`.
#' @param seed Optional integer seed.
#'
#' @return An object of class `net_hon_boot`: a list with `edges` (the
#'   tidy inference table, one row per rule edge of the observed network:
#'   `from`, `to`, `order`, `count`, `probability`, `ci_lower`,
#'   `ci_upper`, `support`, `n_boot_used`), `n_boot`, `level`,
#'   `max_order`, `min_freq`, `n_trajectories`, and `seed`. Has `print`,
#'   `summary`, `plot`, and `as.data.frame` methods; `as.data.frame()`
#'   returns the inference table (optionally filtered with
#'   `min_support =` or restricted with `order_min =`).
#'
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances} 2(5),
#' e1600028. \doi{10.1126/sciadv.1600028}
#'
#' Efron, B., & Tibshirani, R. J. (1993). \emph{An Introduction to the
#' Bootstrap}. Chapman & Hall.
#'
#' @examples
#' hg_seqs <- list(
#'   c("a", "b", "c", "a", "b", "c"),
#'   c("x", "b", "d", "x", "b", "d"),
#'   c("a", "b", "c", "a", "b", "c"),
#'   c("x", "b", "d", "x", "b", "d")
#' )
#' bs <- bootstrap_hon(hg_seqs, n_boot = 50, max_order = 2, seed = 1)
#' bs
#' rules <- as.data.frame(bs)
#' head(rules)
#'
#' @seealso [build_hon()], [compare_hon()], [markov_order_test()]
#'
#' @export
bootstrap_hon <- function(data, n_boot = 500L, level = 0.95,
                          max_order = 5L, min_freq = 1L,
                          collapse_repeats = FALSE,
                          action = NULL, actor = NULL, time = NULL,
                          parallel = FALSE, n_cores = 2L, seed = NULL) {
  stopifnot(
    "`n_boot` must be a single integer >= 2" =
      is.numeric(n_boot) && length(n_boot) == 1L && is.finite(n_boot) &&
      n_boot >= 2,
    "`level` must be a single number in (0, 1)" =
      is.numeric(level) && length(level) == 1L && is.finite(level) &&
      level > 0 && level < 1,
    "`max_order` must be >= 1" = is.numeric(max_order) && max_order >= 1,
    "`min_freq` must be >= 1" = is.numeric(min_freq) && min_freq >= 1
  )
  n_boot <- as.integer(n_boot)
  max_order <- as.integer(max_order)
  min_freq <- as.integer(min_freq)

  trajectories <- .hi_parse(data, action = action, actor = actor,
                            time = time,
                            collapse_repeats = collapse_repeats)
  n_seq <- length(trajectories)
  if (n_seq < 2L) {
    stop("Bootstrap over sequences needs at least 2 sequences.",
         call. = FALSE)
  }

  sc <- .hi_seq_counts(trajectories, max_order)

  # Observed fit from the same counts functional the replicates use
  observed <- .hi_rule_table(
    .hon_extract_rules_count(.hi_count_env(sc, rep(1L, n_seq)),
                             max_order, min_freq))

  # All randomness up front (serial), so parallel == serial under a seed
  if (!is.null(seed)) set.seed(as.integer(seed))
  draws <- matrix(sample.int(n_seq, n_seq * n_boot, replace = TRUE),
                  nrow = n_seq, ncol = n_boot)

  one_rep <- function(b) {
    w <- tabulate(draws[, b], nbins = n_seq)
    env <- .hi_count_env(sc, w)
    extraction <- .hon_extract_rules_count(env, max_order, min_freq)
    # distributions were built inside extraction; rebuild lookup from the
    # rules env for extracted contexts and the distr env for probabilities
    distr <- .hon_build_distributions(env, min_freq)
    list(
      prob = .hi_lookup_prob(distr, observed$source_key, observed$to),
      is_rule = observed$source_key %in% ls(extraction$rules)
    )
  }
  reps <- .hi_apply(n_boot, one_rep, parallel, n_cores)

  prob_mat <- do.call(cbind, lapply(reps, `[[`, "prob"))
  rule_mat <- do.call(cbind, lapply(reps, `[[`, "is_rule"))

  alpha <- 1 - level
  ci <- apply(prob_mat, 1L, stats::quantile,
              probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
  n_used <- rowSums(!is.na(prob_mat))

  edges <- data.frame(
    from        = observed$from,
    to          = observed$to,
    order       = observed$order,
    count       = observed$count,
    probability = observed$probability,
    ci_lower    = as.numeric(ci[1L, ]),
    ci_upper    = as.numeric(ci[2L, ]),
    support     = rowMeans(rule_mat),
    n_boot_used = as.integer(n_used),
    stringsAsFactors = FALSE
  )

  structure(
    list(
      edges = edges,
      n_boot = n_boot,
      level = level,
      max_order = max_order,
      min_freq = min_freq,
      n_trajectories = n_seq,
      seed = seed
    ),
    class = "net_hon_boot"
  )
}

# ---------------------------------------------------------------------------
# compare_hon
# ---------------------------------------------------------------------------

#' Two-sample permutation comparison of higher-order network rules
#'
#' Tests whether two cohorts of sequences differ in their higher-order
#' rule probabilities. The rule set is extracted from the pooled data (see
#' [build_hon()]); for every pooled rule edge the statistic is the
#' absolute difference of the two cohorts' conditional probabilities, and
#' its null distribution comes from permuting cohort labels over
#' sequences. Per-edge p-values are Benjamini-Hochberg adjusted; a global
#' test aggregates the edge differences weighted by pooled counts.
#'
#' As in [bootstrap_hon()], per-sequence counts are precomputed once and
#' every permutation is a weighted aggregation; permutations are drawn
#' before any parallel work, so `parallel = TRUE` reproduces the serial
#' result under the same `seed`. Edges whose context is unobserved in a
#' cohort under some permutation contribute only their valid permutations
#' (`n_perm_used`).
#'
#' @param x,y The two cohorts of sequence data, each in any input format
#'   accepted by [bootstrap_hon()] (wide data.frame, list,
#'   `tna`/`netobject`, or long data.frame with `action`/`actor`/`time`).
#' @param n_perm Integer >= 2. Label permutations. Default `1000`.
#' @param alpha Significance level for the `significant` flag on the
#'   adjusted p-values. Default `0.05`.
#' @param max_order,min_freq,collapse_repeats As in [build_hon()].
#' @param action,actor,time Long-format column names applied to both `x`
#'   and `y`; `NULL` for wide/list input.
#' @param names Character vector of length 2 naming the cohorts in the
#'   output (default `c("x", "y")`).
#' @param parallel,n_cores,seed As in [bootstrap_hon()].
#'
#' @return An object of class `net_hon_compare`: a list with `edges` (one
#'   row per pooled rule edge: `from`, `to`, `order`, `count`,
#'   `count_x`/`count_y` and `prob_x`/`prob_y` (columns named after
#'   `names`), `diff` (prob difference, first minus second), `p_value`,
#'   `p_adj` (BH), `significant`, `n_perm_used`), `global` (list:
#'   `statistic` - the pooled-count-weighted mean absolute difference -
#'   and `p_value`), `names`, `n_perm`, `alpha`, `max_order`, `min_freq`,
#'   `n_trajectories` (per cohort), and `seed`. Has `print`, `summary`,
#'   `plot`, and `as.data.frame` methods; `as.data.frame()` returns the
#'   edge table (`significant = TRUE` restricts it).
#'
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances} 2(5),
#' e1600028. \doi{10.1126/sciadv.1600028}
#'
#' Good, P. (2005). \emph{Permutation, Parametric and Bootstrap Tests of
#' Hypotheses} (3rd ed.). Springer.
#'
#' @examples
#' first_order  <- replicate(6, sample(c("a", "b", "c"), 12, replace = TRUE),
#'                           simplify = FALSE)
#' second_order <- replicate(6, rep(c("a", "b", "c", "b"), 3),
#'                           simplify = FALSE)
#' cmp <- compare_hon(first_order, second_order, n_perm = 99,
#'                    max_order = 2, seed = 1)
#' cmp
#' head(as.data.frame(cmp))
#'
#' @seealso [bootstrap_hon()], [build_hon()], [markov_order_test()]
#'
#' @export
compare_hon <- function(x, y, n_perm = 1000L, alpha = 0.05,
                        max_order = 5L, min_freq = 1L,
                        collapse_repeats = FALSE,
                        action = NULL, actor = NULL, time = NULL,
                        names = c("x", "y"),
                        parallel = FALSE, n_cores = 2L, seed = NULL) {
  stopifnot(
    "`n_perm` must be a single integer >= 2" =
      is.numeric(n_perm) && length(n_perm) == 1L && is.finite(n_perm) &&
      n_perm >= 2,
    "`alpha` must be a single number in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
      alpha > 0 && alpha < 1,
    "`names` must be two distinct cohort labels" =
      is.character(names) && length(names) == 2L && !anyNA(names) &&
      names[1L] != names[2L],
    "`max_order` must be >= 1" = is.numeric(max_order) && max_order >= 1,
    "`min_freq` must be >= 1" = is.numeric(min_freq) && min_freq >= 1
  )
  n_perm <- as.integer(n_perm)
  max_order <- as.integer(max_order)
  min_freq <- as.integer(min_freq)

  tr_x <- .hi_parse(x, action = action, actor = actor, time = time,
                    collapse_repeats = collapse_repeats)
  tr_y <- .hi_parse(y, action = action, actor = actor, time = time,
                    collapse_repeats = collapse_repeats)
  n_x <- length(tr_x)
  n_y <- length(tr_y)
  if (n_x < 2L || n_y < 2L) {
    stop("Each cohort needs at least 2 sequences.", call. = FALSE)
  }

  trajectories <- c(tr_x, tr_y)
  n_seq <- n_x + n_y
  is_x_obs <- c(rep(TRUE, n_x), rep(FALSE, n_y))

  sc <- .hi_seq_counts(trajectories, max_order)

  # Pooled rule set defines the edges under test
  pooled <- .hi_rule_table(
    .hon_extract_rules_count(.hi_count_env(sc, rep(1L, n_seq)),
                             max_order, min_freq))

  probs_for <- function(is_x) {
    distr_x <- .hon_build_distributions(
      .hi_count_env(sc, as.integer(is_x)), min_freq)
    distr_y <- .hon_build_distributions(
      .hi_count_env(sc, as.integer(!is_x)), min_freq)
    cbind(.hi_lookup_prob(distr_x, pooled$source_key, pooled$to),
          .hi_lookup_prob(distr_y, pooled$source_key, pooled$to))
  }

  obs_probs <- probs_for(is_x_obs)
  obs_diff <- obs_probs[, 1L] - obs_probs[, 2L]
  w_edge <- pooled$count / sum(pooled$count)
  obs_global <- sum(w_edge * abs(obs_diff), na.rm = TRUE)

  count_x_env <- .hi_count_env(sc, as.integer(is_x_obs))
  count_y_env <- .hi_count_env(sc, as.integer(!is_x_obs))
  count_of <- function(env) {
    vapply(seq_len(nrow(pooled)), function(i) {
      counts <- env[[pooled$source_key[i]]]
      if (is.null(counts)) return(0L)  # context absent from this cohort
      v <- counts[pooled$to[i]]
      if (is.na(v)) 0L else as.integer(v)
    }, integer(1L))
  }

  # All permutations drawn up front (serial), so parallel == serial
  if (!is.null(seed)) set.seed(as.integer(seed))
  perms <- replicate(n_perm, sample(is_x_obs))

  one_perm <- function(p) {
    pp <- probs_for(perms[, p])
    list(diff = abs(pp[, 1L] - pp[, 2L]),
         global = sum(w_edge * abs(pp[, 1L] - pp[, 2L]), na.rm = TRUE))
  }
  reps <- .hi_apply(n_perm, one_perm, parallel, n_cores)

  diff_mat <- do.call(cbind, lapply(reps, `[[`, "diff"))
  global_null <- vapply(reps, `[[`, numeric(1L), "global")

  n_valid <- rowSums(!is.na(diff_mat))
  exceed <- rowSums(diff_mat >= abs(obs_diff) - 1e-12, na.rm = TRUE)
  p_edge <- ifelse(is.na(obs_diff), NA_real_,
                   (exceed + 1) / (n_valid + 1))
  p_adj <- stats::p.adjust(p_edge, method = "BH")
  p_global <- (sum(global_null >= obs_global - 1e-12) + 1) / (n_perm + 1)

  edges <- data.frame(
    from  = pooled$from,
    to    = pooled$to,
    order = pooled$order,
    count = pooled$count,
    stringsAsFactors = FALSE
  )
  edges[[paste0("count_", names[1L])]] <- count_of(count_x_env)
  edges[[paste0("count_", names[2L])]] <- count_of(count_y_env)
  edges[[paste0("prob_", names[1L])]] <- obs_probs[, 1L]
  edges[[paste0("prob_", names[2L])]] <- obs_probs[, 2L]
  edges$diff <- obs_diff
  edges$p_value <- p_edge
  edges$p_adj <- p_adj
  edges$significant <- !is.na(p_adj) & p_adj < alpha
  edges$n_perm_used <- as.integer(n_valid)

  structure(
    list(
      edges = edges,
      global = list(statistic = obs_global, p_value = p_global),
      names = names,
      n_perm = n_perm,
      alpha = alpha,
      max_order = max_order,
      min_freq = min_freq,
      n_trajectories = stats::setNames(c(n_x, n_y), names),
      seed = seed
    ),
    class = "net_hon_compare"
  )
}

# ---------------------------------------------------------------------------
# S3 methods: net_hon_boot
# ---------------------------------------------------------------------------

#' Coerce a net_hon_boot to its tidy inference table
#'
#' @param x A `net_hon_boot` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param min_support Numeric in `[0, 1]` or NULL. Keep only rule edges
#'   with at least this bootstrap support.
#' @param order_min Integer or NULL. Keep only rule edges of at least this
#'   order (e.g. `2` for the genuinely higher-order rules).
#' @param sort_by `NULL` (order/from/to, default), `"count"`,
#'   `"probability"`, or `"support"` - sort the table by that column,
#'   largest first (ties broken by from/to so the order is deterministic).
#' @return A data.frame, one row per rule edge: `from`, `to`, `order`,
#'   `count`, `probability`, `ci_lower`, `ci_upper`, `support`,
#'   `n_boot_used`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_hon_boot <- function(x, row.names = NULL,
                                       optional = FALSE, ...,
                                       min_support = NULL,
                                       order_min = NULL,
                                       sort_by = NULL, top = NULL) {
  out <- x$edges
  if (!is.null(min_support)) {
    stopifnot("`min_support` must be a single number in [0, 1]" =
                is.numeric(min_support) && length(min_support) == 1L &&
                min_support >= 0 && min_support <= 1)
    out <- out[out$support >= min_support, , drop = FALSE]
  }
  if (!is.null(order_min)) {
    stopifnot("`order_min` must be a single integer >= 1" =
                is.numeric(order_min) && length(order_min) == 1L &&
                order_min >= 1)
    out <- out[out$order >= order_min, , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("count", "probability", "support"))
    out <- out[order(-out[[sort_by]], out$from, out$to), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Print method for net_hon_boot
#'
#' @param x A `net_hon_boot` object.
#' @param ... Additional arguments (ignored).
#' @return The input `x`, invisibly.
#' @export
print.net_hon_boot <- function(x, ...) {
  e <- x$edges
  ho <- e[e$order > 1L, , drop = FALSE]
  cat(sprintf("HON bootstrap: %d rule edges (%d higher-order) from %d sequences\n",
              nrow(e), nrow(ho), x$n_trajectories))
  cat(sprintf("  %d replicates, %.0f%% percentile CIs\n",
              x$n_boot, 100 * x$level))
  if (nrow(ho) > 0L) {
    cat(sprintf("  Higher-order rule support: min %.2f, median %.2f, max %.2f\n",
                min(ho$support), stats::median(ho$support), max(ho$support)))
  }
  cat("  Tidy table: as.data.frame(x); higher-order only:",
      "as.data.frame(x, order_min = 2)\n")
  invisible(x)
}

#' Summary method for net_hon_boot
#'
#' @param object A `net_hon_boot` object.
#' @param ... Additional arguments (ignored).
#' @return A data.frame, one row per rule order: `order`, `n_edges`,
#'   `mean_support`, `min_support`, `mean_ci_width`.
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @export
summary.net_hon_boot <- function(object, ...) {
  e <- object$edges
  by_ord <- split(e, e$order)
  out <- do.call(rbind, lapply(by_ord, function(d) {
    data.frame(
      order = d$order[1L],
      n_edges = nrow(d),
      mean_support = mean(d$support),
      min_support = min(d$support),
      mean_ci_width = mean(d$ci_upper - d$ci_lower),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  invisible(out)
}

#' Plot method for net_hon_boot
#'
#' Forest plot of the rule-edge probabilities with their bootstrap
#' percentile intervals, the most frequent edges first. Order is encoded
#' by both colour (Okabe-Ito) and point shape.
#'
#' @param x A `net_hon_boot` object.
#' @param top Integer. Number of edges to show (by descending count).
#'   Default `20`.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object, invisibly.
#' @export
plot.net_hon_boot <- function(x, top = 20L, ...) {
  e <- x$edges
  e <- e[order(-e$count, e$from, e$to), , drop = FALSE]
  e <- e[seq_len(min(top, nrow(e))), , drop = FALSE]
  e$label <- paste(e$from, "->", e$to)
  e$label <- factor(e$label, levels = rev(e$label))
  e$order_f <- factor(e$order)
  pal <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
           "#D55E00", "#CC79A7", "#999999", "#000000")
  p <- ggplot2::ggplot(e, ggplot2::aes(x = .data$probability,
                                       y = .data$label)) +
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$ci_lower, xmax = .data$ci_upper,
                   colour = .data$order_f),
      orientation = "y", width = 0.25, linewidth = 0.5) +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$order_f, shape = .data$order_f),
      size = 2.2) +
    ggplot2::scale_colour_manual(values = pal, name = "order") +
    ggplot2::scale_shape_manual(values = c(16, 17, 15, 18, 8),
                                name = "order") +
    ggplot2::labs(x = sprintf("conditional probability (%.0f%% CI)",
                              100 * x$level),
                  y = NULL) +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}

# ---------------------------------------------------------------------------
# S3 methods: net_hon_compare
# ---------------------------------------------------------------------------

#' Coerce a net_hon_compare to its tidy edge table
#'
#' @param x A `net_hon_compare` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param significant Logical. `TRUE` restricts to edges whose adjusted
#'   p-value falls below the object's `alpha`. Default `FALSE` (all
#'   edges).
#' @param sort_by `NULL` (order/from/to, default), `"abs_diff"`,
#'   `"count"`, or `"p_adj"` - sort by absolute difference or count
#'   (largest first) or adjusted p-value (smallest first), ties broken by
#'   from/to.
#' @return A data.frame, one row per pooled rule edge (see
#'   [compare_hon()] for the columns).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_hon_compare <- function(x, row.names = NULL,
                                          optional = FALSE, ...,
                                          significant = FALSE,
                                          sort_by = NULL, top = NULL) {
  out <- x$edges
  if (isTRUE(significant)) {
    out <- out[out$significant, , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("abs_diff", "count", "p_adj"))
    key <- switch(sort_by,
                  abs_diff = -abs(out$diff),
                  count    = -out$count,
                  p_adj    = out$p_adj)
    out <- out[order(key, out$from, out$to), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Print method for net_hon_compare
#'
#' @param x A `net_hon_compare` object.
#' @param ... Additional arguments (ignored).
#' @return The input `x`, invisibly.
#' @export
print.net_hon_compare <- function(x, ...) {
  e <- x$edges
  cat(sprintf("HON comparison: %s (%d sequences) vs %s (%d sequences)\n",
              x$names[1L], x$n_trajectories[1L],
              x$names[2L], x$n_trajectories[2L]))
  cat(sprintf("  %d pooled rule edges, %d permutations\n",
              nrow(e), x$n_perm))
  cat(sprintf("  Global weighted |diff|: %.4f, p = %.4g\n",
              x$global$statistic, x$global$p_value))
  cat(sprintf("  Significant edges (BH, alpha = %.2f): %d\n",
              x$alpha, sum(e$significant)))
  cat("  Tidy table: as.data.frame(x); significant only:",
      "as.data.frame(x, significant = TRUE)\n")
  invisible(x)
}

#' Summary method for net_hon_compare
#'
#' @param object A `net_hon_compare` object.
#' @param ... Additional arguments (ignored).
#' @return A data.frame, one row per rule order: `order`, `n_edges`,
#'   `n_significant`, `max_abs_diff`.
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @export
summary.net_hon_compare <- function(object, ...) {
  e <- object$edges
  by_ord <- split(e, e$order)
  out <- do.call(rbind, lapply(by_ord, function(d) {
    data.frame(
      order = d$order[1L],
      n_edges = nrow(d),
      n_significant = sum(d$significant),
      max_abs_diff = max(abs(d$diff), na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  invisible(out)
}

#' Plot method for net_hon_compare
#'
#' Difference plot of the per-edge probability differences, the largest
#' absolute differences first. Significance (BH-adjusted) is encoded by
#' both colour and point shape.
#'
#' @param x A `net_hon_compare` object.
#' @param top Integer. Number of edges to show (by descending absolute
#'   difference). Default `20`.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object, invisibly.
#' @export
plot.net_hon_compare <- function(x, top = 20L, ...) {
  e <- x$edges
  e <- e[!is.na(e$diff), , drop = FALSE]
  e <- e[order(-abs(e$diff), e$from, e$to), , drop = FALSE]
  e <- e[seq_len(min(top, nrow(e))), , drop = FALSE]
  e$label <- paste(e$from, "->", e$to)
  e$label <- factor(e$label, levels = rev(e$label))
  e$sig_f <- factor(ifelse(e$significant, "significant", "not significant"),
                    levels = c("significant", "not significant"))
  p <- ggplot2::ggplot(e, ggplot2::aes(x = .data$diff, y = .data$label)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "#999999") +
    ggplot2::geom_point(
      ggplot2::aes(colour = .data$sig_f, shape = .data$sig_f),
      size = 2.4) +
    ggplot2::scale_colour_manual(
      values = c(significant = "#D55E00", `not significant` = "#0072B2"),
      name = NULL, drop = FALSE) +
    ggplot2::scale_shape_manual(
      values = c(significant = 17, `not significant` = 16),
      name = NULL, drop = FALSE) +
    ggplot2::labs(
      x = sprintf("probability difference (%s - %s)",
                  x$names[1L], x$names[2L]),
      y = NULL) +
    ggplot2::theme_minimal(base_size = 12)
  print(p)
  invisible(p)
}
