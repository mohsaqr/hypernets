# ---- HON rule extraction core (private) + hg_get() for net_hon ----
#
# hon() wraps the Nestimate HON builder (see memory_wrappers.R).
# hypernets keeps a PRIVATE copy of the eager ("hon") rule-extraction core
# because its own inference verbs, hg_bootstrap() and hg_compare(),
# re-run extraction on reweighted counts (.hon_extract_rules_count, which
# Nestimate does not have) and `:::` is not allowed. Every helper below is
# a verbatim copy of Nestimate's same-named internal; identity is asserted
# in local_testing_and_equivalence/test-identity-nestimate-memory.R.
# The one addition is the long-format guard (.hon_guard_long_format()),
# which every sequence-taking verb runs through .ho_sequence_input().

# Separator for encoding tuple keys as strings (non-printable, avoids
# collision with any state label)
.HON_SEP <- "\x01"

# ---------------------------------------------------------------------------
# Encoding / decoding helpers
# ---------------------------------------------------------------------------

#' Encode a character vector (source tuple) into a single string key
#' @param x Character vector representing a source tuple.
#' @return Single string with elements joined by .HON_SEP.
#' @noRd
.hon_encode <- function(x) {
  paste(x, collapse = .HON_SEP)
}

#' Decode a string key back to its character vector tuple
#' @param key Single string created by .hon_encode().
#' @return Character vector.
#' @noRd
.hon_decode <- function(key) {
  strsplit(key, .HON_SEP, fixed = TRUE)[[1L]]
}

#' Return the tuple length for one or more encoded keys
#' @param keys Character vector of encoded keys.
#' @return Integer vector of lengths.
#' @noRd
.hon_key_len <- function(keys) {
  # Count separators + 1
  n_sep <- nchar(keys) - nchar(gsub(.HON_SEP, "", keys, fixed = TRUE))
  as.integer(n_sep + 1L)
}

# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

#' Parse and validate sequence input (private copy of Nestimate's)
#'
#' Converts data.frame or list input into a canonical list of character
#' vectors (one per trajectory). Strips trailing NAs from data.frame rows.
#' Optionally collapses adjacent duplicate states.
#'
#' @param data data.frame or list of character/numeric vectors.
#' @param collapse_repeats Logical. Remove adjacent duplicates.
#' @return List of character vectors.
#' @noRd
# One event per row (long) or one trajectory per row (wide)? Nothing in a
# data.frame's shape distinguishes them, and this parser reads every frame as
# wide. Handed a long table it would treat each ROW as a trajectory and
# silently promote actor ids and timestamps to states -- a two-actor,
# four-turn table became states "1", "2", "3", "4", "A", "B", "C", "s1", "s2"
# and eight trajectories instead of three states and two.
#
# Long-format column names are recognised by role (case-insensitive): a
# frame carrying columns of two or more roles -- a state column, an actor or
# session column, a time or position column -- passed where a wide frame is
# expected, is a long table whose `action`/`actor`/`time` arguments were
# forgotten. Refusing is deliberate -- routing it automatically would guess
# the caller's intent. A wide frame of states never names its columns this
# way (T1, V1, step_1, ...), so it passes. This guard is hypernets' addition;
# Nestimate's parser has none.
.HON_LONG_COLUMNS <- list(
  action = c("action", "code", "state", "event", "activity", "verb"),
  actor = c("actor", "user", "user_id", "userid", "student", "student_id",
            "participant", "person", "session", "session_id", "id"),
  time = c("time", "timestamp", "datetime", "date", "order", "position",
           "order_in_session", "seq_order")
)

.hon_guard_long_format <- function(data) {
  if (!is.data.frame(data)) return(invisible(NULL))
  cols <- names(data)
  hits <- lapply(.HON_LONG_COLUMNS, function(aliases) {
    cols[tolower(cols) %in% aliases]
  })
  roles <- names(hits)[lengths(hits) > 0L]
  if (length(roles) < 2L) return(invisible(NULL))
  found <- unlist(hits[roles], use.names = FALSE)
  suggestion <- paste(vapply(roles, function(r) {
    sprintf("%s = \"%s\"", r, hits[[r]][1L])
  }, character(1L)), collapse = ", ")
  stop(errorCondition(
    sprintf(paste0("`data` looks like a long event table (columns %s), but it ",
                   "would be read as wide -- one trajectory per row, so actor ",
                   "ids and times would become states. To use it, pass the ",
                   "long-format arguments, e.g. `%s`."),
            paste(sprintf("`%s`", found), collapse = ", "), suggestion),
    class = c("hypernets_long_format", "hypernets_bad_input"), call = NULL))
}

.hon_parse_input <- function(data, collapse_repeats = FALSE) {
  .hon_guard_long_format(data)
  if (is.matrix(data) && !is.numeric(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if (is.data.frame(data)) {
    stopifnot("data.frame must have at least one column" = ncol(data) >= 1L)
    stopifnot("data.frame must have at least one row" = nrow(data) >= 1L)
    trajectories <- lapply(seq_len(nrow(data)), function(i) {
      row_vals <- as.character(unlist(data[i, ], use.names = FALSE))
      # Strip trailing NAs
      non_na <- which(!is.na(row_vals))
      if (length(non_na) == 0L) return(character(0L))
      row_vals[seq_len(max(non_na))]
    })
    # Remove empty trajectories
    trajectories <- trajectories[vapply(trajectories, length, integer(1L)) > 0L]
  } else if (is.list(data)) {
    trajectories <- lapply(data, function(x) as.character(x))
  } else {
    stop("'data' must be a data.frame or list of character vectors") # nocov
  }

  if (collapse_repeats) {
    trajectories <- lapply(trajectories, function(traj) {
      if (length(traj) <= 1L) return(traj) # nocov
      keep <- c(TRUE, traj[-1L] != traj[-length(traj)])
      traj[keep]
    })
  }

  # Filter trajectories with < 2 states (no transitions possible)
  trajectories <- trajectories[vapply(trajectories, length, integer(1L)) >= 2L]
  trajectories
}

# ---------------------------------------------------------------------------
# Observation counting
# ---------------------------------------------------------------------------

#' Build observation counts from trajectories
#'
#' For each trajectory and each subsequence length from 2 to (max_order + 1),
#' counts how many times each (source_tuple -> target) transition occurs.
#'
#' @param trajectories List of character vectors.
#' @param max_order Integer. Maximum order to consider.
#' @return Environment mapping encoded source keys to named integer vectors
#'   of target counts.
#' @noRd
.hon_build_observations <- function(trajectories, max_order) {
  count <- new.env(hash = TRUE, parent = emptyenv())

  for (traj in trajectories) {
    n <- length(traj)
    # Subsequence lengths from 2 to min(max_order + 1, n)
    max_len <- min(max_order + 1L, n)
    for (order in seq.int(2L, max_len)) {
      n_subseqs <- n - order + 1L
      for (start in seq_len(n_subseqs)) {
        subseq <- traj[start:(start + order - 1L)]
        target <- subseq[order]
        source_key <- .hon_encode(subseq[-order])
        if (is.null(count[[source_key]])) {
          count[[source_key]] <- integer(0L)
        }
        existing <- count[[source_key]]
        if (is.na(existing[target])) {
          existing[target] <- 1L
        } else {
          existing[target] <- existing[target] + 1L
        }
        count[[source_key]] <- existing
      }
    }
  }

  count
}

# ---------------------------------------------------------------------------
# Distribution building
# ---------------------------------------------------------------------------

#' Build probability distributions from observation counts
#'
#' Zeros out counts below min_freq, then normalizes remaining counts to
#' probabilities.
#'
#' @param count Environment from .hon_build_observations().
#' @param min_freq Integer. Minimum count to keep a transition.
#' @return Environment mapping encoded source keys to named numeric vectors
#'   of probabilities.
#' @noRd
.hon_build_distributions <- function(count, min_freq) {
  distr <- new.env(hash = TRUE, parent = emptyenv())

  for (source_key in ls(count)) {
    counts <- count[[source_key]]
    # Zero out below min_freq - modify count in place so KLDThreshold
    # uses filtered totals (matching pyHON's BuildDistributions behavior)
    counts[counts < min_freq] <- 0L
    count[[source_key]] <- counts
    total <- sum(counts)
    if (total > 0L) {
      # Keep only positive entries
      pos <- counts[counts > 0L]
      distr[[source_key]] <- pos / total
    } else {
      distr[[source_key]] <- numeric(0L)
    }
  }

  distr
}

# ---------------------------------------------------------------------------
# KL-divergence
# ---------------------------------------------------------------------------

#' KL-divergence from distribution a to distribution b
#'
#' Computes sum over targets in a of P(a,t) * log2(P(a,t) / P(b,t)).
#' When P(b,t) = 0 for a target where P(a,t) > 0, contributes Inf.
#'
#' @param a Named numeric vector (probabilities).
#' @param b Named numeric vector (probabilities).
#' @return Numeric scalar.
#' @noRd
.hon_kld <- function(a, b) {
  if (length(a) == 0L) return(0.0)

  divergence <- 0.0
  for (target in names(a)) {
    p_a <- unname(a[target])
    p_b_raw <- b[target]
    p_b <- if (is.na(p_b_raw)) 0.0 else unname(p_b_raw)
    if (p_a > 0.0 && p_b > 0.0) {
      divergence <- divergence + p_a * log2(p_a / p_b)
    } else if (p_a > 0.0 && p_b == 0.0) {
      divergence <- Inf
    }
  }
  divergence
}

#' KL-divergence threshold for deciding whether to extend
#'
#' Threshold = new_order / log2(1 + sum(counts for ext_source)).
#' More data -> lower threshold -> easier to extend.
#' Higher order -> higher threshold -> harder to extend.
#'
#' @param new_order Integer. The proposed new order.
#' @param ext_source_key Encoded source key.
#' @param count Environment of observation counts.
#' @return Numeric scalar.
#' @noRd
.hon_kld_threshold <- function(new_order, ext_source_key, count) {
  total <- sum(count[[ext_source_key]])
  new_order / log2(1.0 + total)
}

# ---------------------------------------------------------------------------
# Source cache (suffix -> extended sources)
# ---------------------------------------------------------------------------

#' Build mapping from suffix tuples to their extended source tuples
#'
#' For each source tuple of length > 1 in the distribution, maps each suffix
#' of that tuple to the full source at the appropriate order.
#'
#' @param distr Environment of distributions.
#' @return Environment mapping encoded suffix keys to environments mapping
#'   order (as character) to character vectors of encoded extended source keys.
#' @noRd
.hon_build_source_cache <- function(distr) {
  cache <- new.env(hash = TRUE, parent = emptyenv())

  for (source_key in ls(distr)) {
    source <- .hon_decode(source_key)
    src_len <- length(source)
    if (src_len > 1L) {
      new_order <- src_len
      order_key <- as.character(new_order)
      # For each suffix (starting positions 2, 3, ..., src_len)
      for (start in seq.int(2L, src_len)) {
        curr <- source[start:src_len]
        curr_key <- .hon_encode(curr)
        if (is.null(cache[[curr_key]])) {
          cache[[curr_key]] <- new.env(hash = TRUE, parent = emptyenv())
        }
        order_env <- cache[[curr_key]]
        if (is.null(order_env[[order_key]])) {
          order_env[[order_key]] <- character(0L)
        }
        existing <- order_env[[order_key]]
        # Add if not already present (set behavior)
        if (!(source_key %in% existing)) {
          order_env[[order_key]] <- c(existing, source_key)
        }
      }
    }
  }

  cache
}

#' Get extensions for a source at a given order from the cache
#'
#' @param curr_key Encoded key for the current source.
#' @param new_order Integer order.
#' @param cache Environment from .hon_build_source_cache().
#' @return Character vector of encoded extended source keys, or character(0).
#' @noRd
.hon_get_extensions <- function(curr_key, new_order, cache) {
  if (is.null(cache[[curr_key]])) return(character(0L))
  order_env <- cache[[curr_key]]
  order_key <- as.character(new_order)
  if (is.null(order_env[[order_key]])) return(character(0L))
  order_env[[order_key]]
}

# ---------------------------------------------------------------------------
# Rule generation (recursive)
# ---------------------------------------------------------------------------

#' Add a source and all its prefixes to the rules
#'
#' @param source_key Encoded source key.
#' @param distr Environment of distributions.
#' @param rules Environment to store rules.
#' @noRd
.hon_add_to_rules <- function(source_key, distr, rules) {
  source <- .hon_decode(source_key)
  if (length(source) > 0L) {
    if (is.null(rules[[source_key]])) {
      rules[[source_key]] <- distr[[source_key]]
    }
    if (length(source) > 1L) {
      prev_key <- .hon_encode(source[-length(source)])
      .hon_add_to_rules(prev_key, distr, rules)
    }
  }
}

#' Recursively extend a rule to higher orders if KLD justifies it
#'
#' @param valid_key Encoded key of the currently valid (accepted) source.
#' @param curr_key Encoded key of the current source to try extending from.
#' @param order Current order level.
#' @param max_order Maximum order.
#' @param distr Environment of distributions.
#' @param count Environment of observation counts.
#' @param cache Source extension cache.
#' @param rules Environment to store rules.
#' @noRd
.hon_extend_rule <- function(valid_key, curr_key, order, max_order,
                             distr, count, cache, rules) {
  if (order >= max_order) {
    .hon_add_to_rules(valid_key, distr, rules)
  } else {
    valid_distr <- distr[[valid_key]]
    new_order <- order + 1L
    extended <- .hon_get_extensions(curr_key, new_order, cache)
    if (length(extended) == 0L) {
      .hon_add_to_rules(valid_key, distr, rules)
    } else {
      for (ext_key in extended) {
        ext_distr <- distr[[ext_key]]
        if (length(ext_distr) > 0L &&
            .hon_kld(ext_distr, valid_distr) >
            .hon_kld_threshold(new_order, ext_key, count)) {
          .hon_extend_rule(ext_key, ext_key, new_order, max_order,
                           distr, count, cache, rules)
        } else {
          .hon_extend_rule(valid_key, ext_key, new_order, max_order,
                           distr, count, cache, rules)
        }
      }
    }
  }
}

#' Extract all HON rules from trajectories
#'
#' Main rule extraction pipeline: builds observations, distributions,
#' source cache, and then generates rules via recursive KLD extension.
#'
#' @param trajectories List of character vectors.
#' @param max_order Integer.
#' @param min_freq Integer.
#' @return List with \code{rules} (environment) and \code{count} (environment).
#' @noRd
.hon_extract_rules <- function(trajectories, max_order, min_freq) {
  count <- .hon_build_observations(trajectories, max_order)
  .hon_extract_rules_count(count, max_order, min_freq)
}

#' Extract rules from a precomputed observation-count environment
#'
#' The counts-based core of .hon_extract_rules(), split out so the
#' inference verbs (hg_bootstrap, hg_compare) can re-run extraction on
#' reweighted counts without re-scanning trajectories.
#'
#' @param count Environment from .hon_build_observations() (or an
#'   equivalent weighted aggregation).
#' @param max_order,min_freq As in .hon_extract_rules().
#' @return list(rules, count).
#' @noRd
.hon_extract_rules_count <- function(count, max_order, min_freq) {
  distr <- .hon_build_distributions(count, min_freq)
  cache <- .hon_build_source_cache(distr)
  rules <- new.env(hash = TRUE, parent = emptyenv())

  # Start from each first-order source
  for (source_key in ls(distr)) {
    if (.hon_key_len(source_key) == 1L) {
      .hon_add_to_rules(source_key, distr, rules)
      .hon_extend_rule(source_key, source_key, 1L, max_order,
                       distr, count, cache, rules)
    }
  }

  list(rules = rules, count = count)
}

# ---------------------------------------------------------------------------
# Node naming
# ---------------------------------------------------------------------------

#' Convert a source/target tuple to readable arrow notation
#'
#' Examples:
#'   ("A",)  -> "A"
#'   ("A","B")  -> "A -> B"
#'   ("X","A","B") -> "X -> A -> B"
#'
#' @param seq Character vector (decoded tuple).
#' @return Character string in arrow notation.
#' @noRd
.hon_sequence_to_node <- function(seq) {
  if (length(seq) == 0L) return("")
  paste(seq, collapse = " -> ")
}

# ---------------------------------------------------------------------------
# hg_get() for net_hon
# ---------------------------------------------------------------------------

#' Tables of a higher-order network
#'
#' @param x A `net_hon` object from [hon()].
#' @param what `"rules"` (default) for the extracted higher-order rules,
#'   `"nodes"` for the node table of the higher-order topology, or
#'   `"pathways"` for the genuinely higher-order rules written as pathways
#'   (`"a b -> c"`: the context states, then the next state), most frequent
#'   first -- the form [simplicial()] turns into simplices.
#' @param ... For `what = "pathways"`: `min_count`, `min_prob` (drop rules
#'   below this count or probability) and `order` (keep one context order).
#'   Otherwise unused.
#' @param order_min Integer or `NULL`. Keep only rules whose source context
#'   is at least this order (e.g. `2` for the genuinely higher-order rules).
#'   Only for `what = "rules"`.
#' @param sort_by `NULL` (construction order, default), `"count"` or
#'   `"probability"` - sort largest first, ties broken by `from`/`to` so the
#'   order is deterministic. Only for `what = "rules"`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return A data.frame. For `what = "rules"`, one row per rule edge with
#'   columns `path`, `from`, `to`, `count`, `probability`, `from_order`,
#'   `to_order`. For `what = "nodes"`, one row per higher-order node with
#'   columns `id` and `node`. For `what = "pathways"`, one row per
#'   higher-order rule with the column `pathway`.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' net <- hon(seqs, max_order = 2)
#' hg_get(net, order_min = 2)
#' hg_get(net, what = "nodes")
#' hg_get(net, what = "pathways")
#' @export
hg_get.net_hon <- function(x, what = c("rules", "nodes", "pathways"), ...,
                           order_min = NULL, sort_by = NULL,
                           top = NULL) {
  what <- match.arg(what)
  if (what == "pathways") {
    return(.ho_top(.ho_pathway_table(Nestimate::pathways(x, ...)), top))
  }
  if (what == "nodes") {
    # the network's node table carries the label twice (label and name)
    out <- data.frame(id = x$nodes$id, node = x$nodes$label,
                      stringsAsFactors = FALSE)
    return(.ho_top(out, top))
  }
  out <- x$ho_edges
  if (!is.null(order_min)) {
    stopifnot("`order_min` must be a single integer >= 1" =
                is.numeric(order_min) && length(order_min) == 1L &&
                order_min >= 1)
    out <- out[out$from_order >= order_min, , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("count", "probability"))
    out <- out[order(-out[[sort_by]], out$from, out$to), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
