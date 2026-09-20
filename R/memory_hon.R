# ---- Higher-Order Networks (HON) ----

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

#' Parse and validate input for build_hon
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
# wide. Handed a long table it therefore treated each ROW as a trajectory and
# silently promoted actor ids and timestamps to states -- a two-actor,
# four-turn table became states "1", "2", "3", "4", "A", "B", "C", "s1", "s2"
# and eight trajectories instead of three states and two. No error, no
# warning, a confidently wrong model.
#
# The canonical long column names are the argument names the memory verbs
# already use, so their presence is a reliable signal: a frame carrying two or
# more of them, passed where a wide frame is expected, is a long table whose
# `action`/`actor`/`time` arguments were forgotten. Refusing is deliberate --
# routing it automatically would be a guess about the caller's intent, and the
# verbs that cannot take those arguments at all have no route to offer.
.HON_LONG_COLUMNS <- c("action", "actor", "time")

.hon_guard_long_format <- function(data, verb = NULL) {
  if (!is.data.frame(data)) return(invisible(NULL))
  found <- .HON_LONG_COLUMNS[.HON_LONG_COLUMNS %in% tolower(names(data))]
  if (length(found) < 2L) return(invisible(NULL))
  takes_long <- is.null(verb) || verb %in% c("build_hon", "bootstrap_hon",
                                             "compare_hon")
  remedy <- if (takes_long) {
    paste0("pass the long-format arguments, e.g. ",
           "`action = \"action\", actor = \"actor\", time = \"time\"`")
  } else {
    paste0("this verb reads wide sequences only -- split the table first, ",
           "e.g. `split(d$action, d$actor)`")
  }
  stop(errorCondition(
    sprintf(paste0("`data` looks like a long event table (columns %s), but it ",
                   "would be read as wide -- one trajectory per row, so actor ",
                   "ids and times would become states. To use it, %s."),
            paste(sprintf("`%s`", found), collapse = ", "), remedy),
    class = c("hypernets_long_format", "hypernets_bad_input"), call = NULL))
}

.hon_parse_input <- function(data, collapse_repeats = FALSE,
                             verb = NULL) {
  .hon_guard_long_format(data, verb)
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
# HON+ helpers (parameter-free, lazy observation building)
# ---------------------------------------------------------------------------

#' HON+: maximum-divergence singleton distribution
#'
#' Returns a named numeric(1) that places all probability mass on the least
#' probable target in \code{distr_vec}.  This is the distribution that would
#' maximally diverge from \code{distr_vec}, used as a pre-check in
#' \code{.honp_extend_rule()} to prune branches that cannot possibly exceed the
#' KLD threshold.
#'
#' @param distr_vec Named numeric vector of probabilities.
#' @return Named numeric(1) with value 1.0, named after the least probable
#'   target.
#' @noRd
.honp_max_divergence <- function(distr_vec) {
  min_target <- names(which.min(distr_vec))
  result <- 1.0
  names(result) <- min_target
  result
}

#' HON+: build order-1 counts, distributions, and starting points
#'
#' Scans all trajectories once to collect order-1 transition counts together
#' with the (trajectory-index, position) pairs needed for lazy higher-order
#' extension.
#'
#' @param trajectories List of character vectors.
#' @param min_freq Integer.  Minimum count to retain a transition.
#' @return List with three environments: \code{count}, \code{distr}, and
#'   \code{starting_points}.
#' @noRd
.honp_build_order1 <- function(trajectories, min_freq) {
  count <- new.env(hash = TRUE, parent = emptyenv())
  starting_points <- new.env(hash = TRUE, parent = emptyenv())

  for (ti in seq_along(trajectories)) {
    traj <- trajectories[[ti]]
    n <- length(traj)
    for (pos in seq_len(n - 1L)) {
      source_key <- .hon_encode(traj[pos])
      target <- traj[pos + 1L]

      # Count
      cur <- count[[source_key]]
      if (is.null(cur)) {
        cur <- integer(0L)
        names(cur) <- character(0L)
      }
      idx <- match(target, names(cur))
      if (is.na(idx)) {
        cur[target] <- 1L
      } else {
        cur[idx] <- cur[idx] + 1L
      }
      count[[source_key]] <- cur

      # StartingPoints: store (traj_index, position) pairs
      sp <- starting_points[[source_key]]
      if (is.null(sp)) sp <- list()
      sp[[length(sp) + 1L]] <- c(ti, pos)
      starting_points[[source_key]] <- sp
    }
  }

  # Build distributions (apply min_freq filtering then normalise)
  distr <- new.env(hash = TRUE, parent = emptyenv())
  for (source_key in ls(count)) {
    counts <- count[[source_key]]
    counts[counts < min_freq] <- 0L
    count[[source_key]] <- counts  # zero in place (for threshold calc)
    total <- sum(counts)
    if (total > 0L) {
      pos_counts <- counts[counts > 0L]
      distr[[source_key]] <- pos_counts / total
    } else {
      distr[[source_key]] <- numeric(0L)
    }
  }

  list(count = count, distr = distr, starting_points = starting_points)
}

#' HON+: lazily build higher-order counts for a source key
#'
#' For each starting point stored under \code{source_key}, looks one step back
#' in the trajectory to discover extended sources of order
#' \code{length(source) + 1}.  Populates \code{count}, \code{distr},
#' \code{starting_points}, and \code{source_to_ext} in place.
#'
#' @param source_key Encoded source key whose extensions are needed.
#' @param source_to_ext Environment mapping suffix keys to environments of
#'   extended source keys (set-like).
#' @param count Environment of observation counts (modified in place).
#' @param distr Environment of distributions (modified in place).
#' @param starting_points Environment of starting-point lists (modified in
#'   place).
#' @param trajectories List of character vectors.
#' @param min_freq Integer minimum frequency.
#' @return \code{NULL} invisibly.
#' @noRd
.honp_extend_observation <- function(source_key, source_to_ext,
                                     count, distr, starting_points,
                                     trajectories, min_freq) {
  source <- .hon_decode(source_key)
  order <- length(source)

  # Ensure the suffix order is already extended before going deeper
  if (order > 1L) {
    suffix_key <- .hon_encode(source[-1L])
    if (is.null(count[[suffix_key]]) || length(count[[suffix_key]]) == 0L) {
      .honp_extend_observation(suffix_key, source_to_ext, count, distr, # nocov start
                               starting_points, trajectories, min_freq) # nocov end
    }
  }

  sp <- starting_points[[source_key]]
  if (is.null(sp) || length(sp) == 0L) return(invisible(NULL)) # nocov

  # Accumulate counts for extended sources in a local environment
  local_count <- new.env(hash = TRUE, parent = emptyenv())

  for (point in sp) {
    ti <- point[1L]
    pos <- point[2L]
    traj <- trajectories[[ti]]

    # Look one step back and one step forward
    if (pos > 1L && (pos + order) <= length(traj)) {
      ext_source <- traj[(pos - 1L):(pos + order - 1L)]
      target <- traj[pos + order]
      ext_key <- .hon_encode(ext_source)

      # Accumulate count
      cur <- local_count[[ext_key]]
      if (is.null(cur)) {
        cur <- integer(0L)
        names(cur) <- character(0L)
      }
      idx <- match(target, names(cur))
      if (is.na(idx)) {
        cur[target] <- 1L
      } else {
        cur[idx] <- cur[idx] + 1L
      }
      local_count[[ext_key]] <- cur

      # Record starting point for the extended source
      esp <- starting_points[[ext_key]]
      if (is.null(esp)) esp <- list()
      esp[[length(esp) + 1L]] <- c(ti, pos - 1L)
      starting_points[[ext_key]] <- esp
    }
  }

  # Apply min_freq, build distributions, register in source_to_ext
  for (ext_key in ls(local_count)) {
    lc <- local_count[[ext_key]]
    lc[lc < min_freq] <- 0L
    count[[ext_key]] <- lc
    total <- sum(lc)
    if (total > 0L) {
      pos_counts <- lc[lc > 0L]
      distr[[ext_key]] <- pos_counts / total

      # Map suffix -> ext_source in source_to_ext
      suffix_key <- .hon_encode(.hon_decode(ext_key)[-1L])
      s2e <- source_to_ext[[suffix_key]]
      if (is.null(s2e)) {
        s2e <- new.env(hash = TRUE, parent = emptyenv())
        source_to_ext[[suffix_key]] <- s2e
      }
      s2e[[ext_key]] <- TRUE
    }
  }

  invisible(NULL)
}

#' HON+: lazy cache lookup for extended sources
#'
#' Returns the character vector of encoded extended-source keys for
#' \code{curr_key}.  If the cache has not been populated yet, calls
#' \code{.honp_extend_observation()} first.
#'
#' @param curr_key Encoded key of the current source.
#' @param source_to_ext Extension cache environment.
#' @param count Count environment.
#' @param distr Distribution environment.
#' @param starting_points Starting-points environment.
#' @param trajectories List of character vectors.
#' @param min_freq Integer minimum frequency.
#' @return Character vector of extended source keys that have non-empty
#'   distributions.
#' @noRd
.honp_extend_source_fast <- function(curr_key, source_to_ext,
                                      count, distr, starting_points,
                                      trajectories, min_freq) {
  s2e <- source_to_ext[[curr_key]]
  if (!is.null(s2e)) {
    ext_keys <- ls(s2e)
    has_distr <- vapply(ext_keys, function(k) {
      d <- distr[[k]]
      !is.null(d) && length(d) > 0L
    }, logical(1L))
    return(ext_keys[has_distr])
  }

  # Not cached yet - build lazily
  .honp_extend_observation(curr_key, source_to_ext, count, distr,
                           starting_points, trajectories, min_freq)

  s2e <- source_to_ext[[curr_key]]
  if (is.null(s2e)) return(character(0L))

  ext_keys <- ls(s2e)
  has_distr <- vapply(ext_keys, function(k) {
    d <- distr[[k]]
    !is.null(d) && length(d) > 0L
  }, logical(1L))
  ext_keys[has_distr]
}

#' HON+: add a source and all its prefixes to the rules (HON+ variant)
#'
#' Mirrors pyHON+'s \code{AddToRules()}: for each prefix of \code{source_key},
#' if the prefix is not yet in \code{distr} (or has an empty distribution),
#' triggers \code{.honp_extend_source_fast()} on the prefix's suffix to
#' populate it lazily before writing to \code{rules}.
#'
#' @param source_key Encoded source key.
#' @param distr Distribution environment.
#' @param rules Rules environment.
#' @param source_to_ext Extension cache environment.
#' @param count Count environment.
#' @param starting_points Starting-points environment.
#' @param trajectories List of character vectors.
#' @param min_freq Integer minimum frequency.
#' @noRd
.honp_add_to_rules <- function(source_key, distr, rules,
                                source_to_ext, count, starting_points,
                                trajectories, min_freq) {
  source <- .hon_decode(source_key)
  for (ord in seq_along(source)) {
    prefix <- source[seq_len(ord)]
    prefix_key <- .hon_encode(prefix)
    d <- distr[[prefix_key]]
    if (is.null(d) || length(d) == 0L) {
      if (ord > 1L) {
        suffix_key <- .hon_encode(prefix[-1L])
        .honp_extend_source_fast(suffix_key, source_to_ext, count, distr,
                                  starting_points, trajectories, min_freq)
      }
    }
    rules[[prefix_key]] <- distr[[prefix_key]]
  }
}

#' HON+: recursively extend a rule with MaxDivergence pre-check
#'
#' Extends rules to higher orders.  Before attempting to extend, checks whether
#' even the maximum-possible divergence from \code{curr_distr} would exceed the
#' KLD threshold.  If not, the current \code{valid_key} is accepted and no
#' further extension is tried (pruning).
#'
#' @param valid_key Encoded key of the currently accepted (valid) source.
#' @param curr_key Encoded key of the source being evaluated for extension.
#' @param order Current order level.
#' @param max_order Maximum allowed order.
#' @param distr Distribution environment.
#' @param count Count environment.
#' @param source_to_ext Extension cache environment.
#' @param starting_points Starting-points environment.
#' @param trajectories List of character vectors.
#' @param min_freq Integer minimum frequency.
#' @param rules Rules environment (modified in place).
#' @noRd
.honp_extend_rule <- function(valid_key, curr_key, order, max_order,
                              distr, count, source_to_ext, starting_points,
                              trajectories, min_freq, rules) {
  if (order >= max_order) {
    .honp_add_to_rules(valid_key, distr, rules, source_to_ext, count,
                        starting_points, trajectories, min_freq)
    return(invisible(NULL))
  }

  valid_distr <- distr[[valid_key]]
  curr_distr  <- distr[[curr_key]]

  # HON+ MaxDivergence pre-check
  if (!is.null(curr_distr) && length(curr_distr) > 0L) {
    max_div <- .honp_max_divergence(curr_distr)
    if (.hon_kld(max_div, valid_distr) <
        .hon_kld_threshold(order + 1L, curr_key, count)) {
      .honp_add_to_rules(valid_key, distr, rules, source_to_ext, count,
                          starting_points, trajectories, min_freq)
      return(invisible(NULL))
    }
  } else {
    .honp_add_to_rules(valid_key, distr, rules, source_to_ext, count,
                        starting_points, trajectories, min_freq)
    return(invisible(NULL))
  }

  new_order <- order + 1L
  extended <- .honp_extend_source_fast(curr_key, source_to_ext, count,
                                        distr, starting_points,
                                        trajectories, min_freq)

  if (length(extended) == 0L) {
    .honp_add_to_rules(valid_key, distr, rules, source_to_ext, count,
                        starting_points, trajectories, min_freq)
  } else {
    for (ext_key in extended) {
      ext_distr <- distr[[ext_key]]
      if (length(ext_distr) > 0L &&
          .hon_kld(ext_distr, valid_distr) >
          .hon_kld_threshold(new_order, ext_key, count)) {
        .honp_extend_rule(ext_key, ext_key, new_order, max_order,
                          distr, count, source_to_ext, starting_points,
                          trajectories, min_freq, rules)
      } else {
        .honp_extend_rule(valid_key, ext_key, new_order, max_order,
                          distr, count, source_to_ext, starting_points,
                          trajectories, min_freq, rules)
      }
    }
  }
}

#' HON+: main rule extraction pipeline
#'
#' Builds order-1 observations then iterates over all order-1 sources,
#' calling \code{.honp_extend_rule()} for each to discover higher-order rules
#' via the MaxDivergence pre-check.
#'
#' @param trajectories List of character vectors.
#' @param max_order Integer. Maximum order to consider.
#' @param min_freq Integer. Minimum transition count.
#' @return List with \code{rules} (environment) and \code{count} (environment).
#' @noRd
.honp_extract_rules <- function(trajectories, max_order, min_freq) {
  o1 <- .honp_build_order1(trajectories, min_freq)
  source_to_ext <- new.env(hash = TRUE, parent = emptyenv())
  rules <- new.env(hash = TRUE, parent = emptyenv())

  for (source_key in ls(o1$distr)) {
    if (.hon_key_len(source_key) == 1L) {
      .honp_add_to_rules(source_key, o1$distr, rules, source_to_ext,
                          o1$count, o1$starting_points, trajectories, min_freq)
      .honp_extend_rule(source_key, source_key, 1L, max_order,
                        o1$distr, o1$count, source_to_ext,
                        o1$starting_points, trajectories, min_freq, rules)
    }
  }

  list(rules = rules, count = o1$count)
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
#' inference verbs (bootstrap_hon, compare_hon) can re-run extraction on
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
# Network building
# ---------------------------------------------------------------------------

#' Build HON graph from rules
#'
#' Processes rules (sorted by source length), adds edges, and rewires
#' to connect higher-order nodes.
#'
#' @param rules Environment of rules.
#' @return Environment mapping encoded source keys to environments mapping
#'   encoded target keys to numeric weights.
#' @noRd
.hon_build_network <- function(rules) {
  graph <- new.env(hash = TRUE, parent = emptyenv())

  # Sort sources by tuple length (ascending)
  all_sources <- ls(rules)
  source_lens <- .hon_key_len(all_sources)
  sorted_sources <- all_sources[order(source_lens)]

  for (source_key in sorted_sources) {
    rule_distr <- rules[[source_key]]
    if (is.null(graph[[source_key]])) {
      graph[[source_key]] <- new.env(hash = TRUE, parent = emptyenv())
    }

    for (target in names(rule_distr)) {
      target_key <- .hon_encode(target)
      graph[[source_key]][[target_key]] <- rule_distr[target]
      source <- .hon_decode(source_key)
      if (length(source) > 1L) {
        .hon_rewire(graph, source_key, target_key)
      }
    }
  }

  .hon_rewire_tails(graph)
  graph
}

#' Rewire incoming edges for a higher-order source
#'
#' When a higher-order source (A,B) is added, redirect the edge from
#' its prefix (A) pointing to (B) so it now points to (A,B).
#'
#' @param graph Environment (the HON graph).
#' @param source_key Encoded source (length > 1).
#' @param target_key Encoded target.
#' @noRd
.hon_rewire <- function(graph, source_key, target_key) {
  source <- .hon_decode(source_key)
  prev_source_key <- .hon_encode(source[-length(source)])
  prev_target_key <- .hon_encode(source[length(source)])

  # Only rewire if the prev_source doesn't already point to this source
  prev_graph <- graph[[prev_source_key]]
  if (is.null(prev_graph) || !is.null(prev_graph[[source_key]])) {
    return(invisible(NULL))
  }

  # Check if prev_source has an edge to prev_target
  if (!is.null(prev_graph[[prev_target_key]])) {
    # Redirect: prev_source -> source (instead of -> prev_target)
    prev_graph[[source_key]] <- prev_graph[[prev_target_key]]
    rm(list = prev_target_key, envir = prev_graph)
  }

  invisible(NULL)
}

#' Rewire outgoing edges to point to highest-order targets
#'
#' For each edge (source -> target) where target is first-order, check if
#' a higher-order version of the target (source + target suffix) exists as
#' a node in the graph. If so, redirect.
#'
#' @param graph Environment (the HON graph).
#' @noRd
.hon_rewire_tails <- function(graph) {
  to_add <- list()
  to_remove <- list()

  for (source_key in ls(graph)) {
    source <- .hon_decode(source_key)
    source_graph <- graph[[source_key]]
    for (target_key in ls(source_graph)) {
      target <- .hon_decode(target_key)
      if (length(target) == 1L) {
        new_target <- c(source, target)
        found <- FALSE
        while (length(new_target) > 1L) {
          new_target_key <- .hon_encode(new_target)
          if (!is.null(graph[[new_target_key]])) {
            to_add[[length(to_add) + 1L]] <- list(
              source = source_key,
              target = new_target_key,
              weight = source_graph[[target_key]]
            )
            to_remove[[length(to_remove) + 1L]] <- list(
              source = source_key,
              target = target_key
            )
            found <- TRUE
            break
          }
          new_target <- new_target[-1L]
        }
      }
    }
  }

  # Apply additions
  for (item in to_add) {
    graph[[item$source]][[item$target]] <- item$weight
  }
  # Apply removals
  for (item in to_remove) {
    src_env <- graph[[item$source]]
    if (!is.null(src_env[[item$target]])) {
      rm(list = item$target, envir = src_env)
    }
  }

  invisible(NULL)
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
# Graph to edge list
# ---------------------------------------------------------------------------

#' Convert HON graph environment to edge list data.frame
#'
#' Returns a data.frame with unified column structure matching MOGen/HYPA:
#' \code{path} (full state sequence), \code{from} (context/conditioning states),
#' \code{to} (next state), \code{count}, \code{probability}, plus
#' \code{from_order} and \code{to_order}.
#'
#' @param graph Environment from .hon_build_network().
#' @param count_env Environment of observation counts (optional, for raw counts).
#' @return data.frame with columns: path, from, to, count, probability,
#'   from_order, to_order.
#' @noRd
.hon_graph_to_edgelist <- function(graph, count_env = NULL) {
  rows <- list()
  idx <- 0L

  for (source_key in ls(graph)) {
    source <- .hon_decode(source_key)
    from_node <- .hon_sequence_to_node(source)
    from_order <- length(source)
    source_graph <- graph[[source_key]]
    for (target_key in ls(source_graph)) {
      target <- .hon_decode(target_key)
      to_order <- length(target)
      prob <- source_graph[[target_key]]

      # Next state = last element of target tuple
      next_state <- target[length(target)]
      # Path = source context + next state
      path <- paste(c(source, next_state), collapse = " -> ")

      # Count lookup from the count environment
      edge_count <- NA_integer_
      if (!is.null(count_env) && !is.null(count_env[[source_key]])) {
        raw <- count_env[[source_key]][[next_state]]
        if (!is.null(raw) && !is.na(raw)) edge_count <- as.integer(raw)
      }

      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        path = path,
        from = from_node,
        to = next_state,
        count = edge_count,
        probability = as.numeric(prob),
        from_order = from_order,
        to_order = to_order,
        stringsAsFactors = FALSE
      )
    }
  }

  if (idx == 0L) {
    return(data.frame(
      path = character(0L),
      from = character(0L),
      to = character(0L),
      count = integer(0L),
      probability = numeric(0L),
      from_order = integer(0L),
      to_order = integer(0L),
      stringsAsFactors = FALSE
    ))
  }

  do.call(rbind, rows)
}

# ---------------------------------------------------------------------------
# Assembly
# ---------------------------------------------------------------------------

#' Assemble final net_hon output object
#'
#' @param graph Environment from .hon_build_network().
#' @param rules Environment from .hon_extract_rules().
#' @param trajectories List of character vectors (parsed input).
#' @param max_order Integer.
#' @param min_freq Integer.
#' @param count_env Environment of observation counts (for raw counts in edges).
#' @return net_hon S3 object.
#' @noRd
.hon_assemble_output <- function(graph, rules, trajectories, max_order,
                                 min_freq, count_env = NULL) {
  edges <- .hon_graph_to_edgelist(graph, count_env)

  # Collect all node names from graph (both sources and targets)
  all_nodes <- character(0L)
  for (source_key in ls(graph)) {
    source <- .hon_decode(source_key)
    all_nodes <- c(all_nodes, .hon_sequence_to_node(source))
    source_graph <- graph[[source_key]]
    for (target_key in ls(source_graph)) {
      target <- .hon_decode(target_key)
      all_nodes <- c(all_nodes, .hon_sequence_to_node(target))
    }
  }
  nodes <- sort(unique(all_nodes))

  # Compute max order observed
  rule_keys <- ls(rules)
  if (length(rule_keys) > 0L) {
    observed_max_order <- max(.hon_key_len(rule_keys))
  } else {
    observed_max_order <- 1L
  }

  # Build adjacency matrix from graph directly (using readable node names)
  n <- length(nodes)
  mat <- matrix(0.0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
  for (source_key in ls(graph)) {
    from_node <- .hon_sequence_to_node(.hon_decode(source_key))
    source_graph <- graph[[source_key]]
    for (target_key in ls(source_graph)) {
      to_node <- .hon_sequence_to_node(.hon_decode(target_key))
      mat[from_node, to_node] <- source_graph[[target_key]]
    }
  }

  # Unique first-order states
  first_order_states <- sort(unique(unlist(trajectories, use.names = FALSE)))

  # cograph-compatible fields
  cg <- .ho_cograph_fields(mat, nodes, method = "hon")

  result <- structure(
    c(list(
      weights = mat,
      matrix = mat,
      ho_edges = edges,
      edges = cg$edges,
      nodes = cg$nodes,
      n_nodes = n,
      n_edges = nrow(edges),
      first_order_states = first_order_states,
      max_order_requested = max_order,
      max_order_observed = observed_max_order,
      min_freq = min_freq,
      n_trajectories = length(trajectories),
      directed = TRUE
    ), cg[c("meta", "node_groups")]),
    class = c("net_hon", "cograph_network")
  )

  result
}

# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

#' Build a Higher-Order Network (HON)
#'
#' @description
#' Constructs a Higher-Order Network from sequential data, faithfully
#' implementing the BuildHON algorithm (Xu, Wickramarathne & Chawla, 2016).
#'
#' The algorithm detects when a first-order Markov model is insufficient
#' to capture sequential dependencies and automatically creates higher-order
#' nodes. Uses KL-divergence to determine whether extending a node's history
#' provides significantly different transition distributions.
#'
#' @param data One of:
#'   \itemize{
#'     \item \code{data.frame}: rows are trajectories, columns are time steps.
#'       Trailing \code{NA}s are stripped. All non-NA values are coerced to
#'       character.
#'     \item \code{list}: each element is a character (or coercible) vector
#'       representing one trajectory.
#'     \item \code{tna}: a tna object with sequence data. Numeric state
#'       IDs are automatically converted to label names.
#'     \item \code{netobject}: a netobject with sequence data.
#'     \item long \code{data.frame}: one event per row, together with
#'       \code{action} (and optionally \code{actor} and \code{time}).
#'   }
#' @param action,actor,time Long-format column names: \code{action} holds
#'   the categorical state of each event, \code{actor} groups events into
#'   trajectories (one per actor; \code{NULL} treats all rows as one), and
#'   \code{time} orders events within an actor (row order when
#'   \code{NULL}). Leave all three \code{NULL} for wide, list, or model
#'   input. Same interface as \code{\link{bootstrap_hon}()}.
#' @param max_order Integer. Maximum order of the HON. Default 5.
#'   The algorithm may produce lower-order nodes if the data do not
#'   justify higher orders.
#' @param min_freq Integer. Minimum frequency for a transition to be
#'   considered. Transitions observed fewer than \code{min_freq} times are
#'   treated as zero. Default 1.
#' @param collapse_repeats Logical. If \code{TRUE}, adjacent duplicate states
#'   within each trajectory are collapsed before analysis. Default \code{FALSE}.
#' @param method Character. Algorithm to use: \code{"hon+"} (default,
#'   parameter-free BuildHON+ with lazy observation building and MaxDivergence
#'   pruning) or \code{"hon"} (original BuildHON with eager observation
#'   building).
#'
#' @return An S3 object of class \code{"net_hon"} containing:
#' \describe{
#'   \item{matrix}{Weighted adjacency matrix (rows = from, cols = to).
#'     Rows and columns use readable arrow notation (e.g., \code{"A -> B"}).}
#'   \item{edges}{Data frame with columns: \code{path} (full state sequence,
#'     e.g., "A -> B -> C"), \code{from} (context/conditioning states),
#'     \code{to} (predicted next state), \code{count} (raw frequency),
#'     \code{probability} (transition probability), \code{from_order},
#'     \code{to_order}.}
#'   \item{nodes}{data.frame with columns \code{id}, \code{label},
#'     \code{name} (one row per HON node; \code{label}/\code{name} are the
#'     arrow-notation node names). Stored as a data.frame for
#'     \code{cograph_network} compatibility.}
#'   \item{n_nodes}{Number of HON nodes.}
#'   \item{n_edges}{Number of edges.}
#'   \item{first_order_states}{Character vector of unique original states.}
#'   \item{max_order_requested}{The \code{max_order} parameter used.}
#'   \item{max_order_observed}{Highest order actually present.}
#'   \item{min_freq}{The \code{min_freq} parameter used.}
#'   \item{n_trajectories}{Number of trajectories after parsing.}
#'   \item{directed}{Logical. Always \code{TRUE}.}
#' }
#'
#' @details
#' \strong{Node naming convention}: Higher-order nodes use readable arrow
#' notation. A first-order node is simply \code{"A"}. A second-order node
#' representing the context "came from A, now at B" is \code{"A -> B"}.
#' Third-order: \code{"A -> B -> C"}, etc.
#'
#' \strong{Algorithm overview}:
#' \enumerate{
#'   \item Count all subsequence transitions up to \code{max_order + 1}.
#'   \item Build probability distributions, filtering by \code{min_freq}.
#'   \item For each first-order source, recursively test whether extending
#'     the history (adding more context) produces a significantly different
#'     distribution (via KL-divergence vs. an adaptive threshold).
#'   \item Build the network from the accepted rules, rewiring edges so
#'     higher-order nodes are properly connected.
#' }
#'
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances},
#' 2(5), e1600028.
#'
#' Saebi, M., Xu, J., Kaplan, L. M., Ribeiro, B., & Chawla, N. V. (2020).
#' Efficient modeling of higher-order dependencies in networks: from algorithm
#' to application for anomaly detection. \emph{EPJ Data Science}, 9(1), 15.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon <- build_hon(seqs, max_order = 2)
#'
#' \donttest{
#' # From list of trajectories
#' trajs <- list(
#'   c("A", "B", "C", "D", "A"),
#'   c("A", "B", "D", "C", "A"),
#'   c("A", "B", "C", "D", "A")
#' )
#' hon <- build_hon(trajs, max_order = 3, min_freq = 1)
#' print(hon)
#' summary(hon)
#'
#' # From data.frame (rows = trajectories)
#' df <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"),
#'                  T3 = c("C", "D"), T4 = c("D", "C"))
#' hon <- build_hon(df, max_order = 2)
#' }
#'
#' @export
build_hon <- function(data, max_order = 5L, min_freq = 1L,
                      collapse_repeats = FALSE, method = "hon+",
                      action = NULL, actor = NULL, time = NULL) {
  # --- Long format: one event per row (same contract as bootstrap_hon) ---
  if (!is.null(action)) {
    method <- match.arg(method, c("hon+", "hon"))
    return(.hon_from_trajectories(
      .hi_parse(data, action = action, actor = actor, time = time,
                collapse_repeats = collapse_repeats),
      max_order = as.integer(max_order), min_freq = as.integer(min_freq),
      method = method))
  }
  stopifnot(
    "`actor` requires `action`" = is.null(actor),
    "`time` requires `action`" = is.null(time)
  )

  # --- Coerce tna/netobject to labeled data.frame ---
  data <- .coerce_sequence_input(data)

  # --- Input validation ---
  stopifnot(
    "'data' must be a data.frame or list" =
      is.data.frame(data) || is.list(data),
    "'data' must have at least one trajectory/row" =
      (is.data.frame(data) && nrow(data) >= 1L && ncol(data) >= 1L) ||
      (is.list(data) && !is.data.frame(data) && length(data) >= 1L),
    "'max_order' must be >= 1" = is.numeric(max_order) && max_order >= 1L,
    "'min_freq' must be >= 1" = is.numeric(min_freq) && min_freq >= 1L
  )

  method <- match.arg(method, c("hon+", "hon"))
  max_order <- as.integer(max_order)
  min_freq <- as.integer(min_freq)

  # --- Parse input ---
  trajectories <- .hon_parse_input(data, collapse_repeats = collapse_repeats,
                             verb = "build_hon")

  if (length(trajectories) == 0L) {
    stop("No valid trajectories (each must have at least 2 states)")
  }

  .hon_from_trajectories(trajectories, max_order, min_freq, method)
}

#' Build a HON from parsed trajectories
#'
#' The post-parsing core of build_hon(), shared by its wide/list path and
#' its long-format path so both produce byte-identical networks from the
#' same trajectories.
#'
#' @param trajectories List of character vectors.
#' @param max_order,min_freq Integers.
#' @param method "hon+" or "hon".
#' @return A `net_hon` object.
#' @noRd
.hon_from_trajectories <- function(trajectories, max_order, min_freq,
                                   method) {
  # --- Extract rules ---
  extraction <- if (method == "hon+") {
    .honp_extract_rules(trajectories, max_order, min_freq)
  } else {
    .hon_extract_rules(trajectories, max_order, min_freq)
  }
  rules <- extraction$rules
  count_env <- extraction$count

  # --- Build network ---
  graph <- .hon_build_network(rules)

  # --- Assemble output ---
  .hon_assemble_output(graph, rules, trajectories, max_order, min_freq,
                       count_env)
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_hon
#'
#' @param x A \code{net_hon} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon <- build_hon(seqs, max_order = 2)
#' print(hon)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' hon <- build_hon(seqs, max_order = 2L)
#' print(hon)
#' }
#'
#' @export
print.net_hon <- function(x, ...) {
  cat("Higher-Order Network (HON)\n")
  cat(sprintf("  Nodes:        %d (%d first-order states)\n",
              x$n_nodes, length(x$first_order_states)))
  cat(sprintf("  Edges:        %d\n", x$n_edges))
  cat(sprintf("  Max order:    %d (requested %d)\n",
              x$max_order_observed, x$max_order_requested))
  cat(sprintf("  Min freq:     %d\n", x$min_freq))
  cat(sprintf("  Trajectories: %d\n", x$n_trajectories))
  invisible(x)
}

#' Summary Method for net_hon
#'
#' @param object A \code{net_hon} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The edge data.frame \code{object$edges} (columns \code{path},
#'   \code{from}, \code{to}, \code{count}, \code{probability},
#'   \code{from_order}, \code{to_order}), returned visibly; the summary text
#'   is printed as a side effect.
#'
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon <- build_hon(seqs, max_order = 2)
#' summary(hon)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' hon <- build_hon(seqs, max_order = 2L)
#' summary(hon)
#' }
#'
#' @export
summary.net_hon <- function(object, ...) {
  cat("Higher-Order Network (HON) Summary\n")
  cat(sprintf("  Nodes: %d | Edges: %d | Trajectories: %d\n",
              object$n_nodes, object$n_edges, object$n_trajectories))
  cat(sprintf("  First-order states: %s\n",
              paste(object$first_order_states, collapse = ", ")))
  cat(sprintf("  Max order observed: %d (requested: %d)\n",
              object$max_order_observed, object$max_order_requested))
  cat(sprintf("  Min frequency: %d\n", object$min_freq))

  if (object$n_nodes > 0L) {
    node_orders <- vapply(object$nodes$label, function(nd) {
      length(strsplit(nd, " -> ", fixed = TRUE)[[1L]])
    }, integer(1L))
    order_tab <- table(node_orders)
    cat("  Node order distribution:\n")
    for (ord in names(order_tab)) {
      cat(sprintf("    Order %s: %d nodes\n", ord, order_tab[ord]))
    }
  }

  invisible(object$edges)
}


#' Coerce a net_hon to a tidy table
#'
#' @param x A `net_hon` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"rules"` (default) for the extracted higher-order rules, or
#'   `"nodes"` for the node table of the higher-order topology.
#' @param order_min Integer or `NULL`. Keep only rules whose source context
#'   is at least this order (e.g. `2` for the genuinely higher-order rules).
#'   Ignored when `what = "nodes"`.
#' @param sort_by `NULL` (construction order, default), `"count"` or
#'   `"probability"` - sort largest first, ties broken by `from`/`to` so the
#'   order is deterministic. Ignored when `what = "nodes"`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return A data.frame. For `what = "rules"`, one row per rule edge with
#'   columns `path`, `from`, `to`, `count`, `probability`, `from_order`,
#'   `to_order`. For `what = "nodes"`, one row per higher-order node with
#'   columns `id`, `label`, `name`.
#' @export
as.data.frame.net_hon <- function(x, row.names = NULL, optional = FALSE, ...,
                                  what = c("rules", "nodes"),
                                  order_min = NULL, sort_by = NULL,
                                  top = NULL) {
  what <- match.arg(what)
  if (what == "nodes") {
    out <- x$nodes
    rownames(out) <- NULL
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
