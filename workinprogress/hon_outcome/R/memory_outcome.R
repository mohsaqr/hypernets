# ---- hon_outcome: higher-order structure related to an actor outcome -----
#
# The package's first supervised verb. It is an ASSEMBLY AND INFERENCE layer:
# every predictor comes from a verb that already exists (hon_centrality(),
# build_honem()), and the only new mathematics is the regression itself,
# which is ordinary least squares / a canonical-link GLM with an optional
# cluster-robust covariance.
#
# The chain is:
#
#   1. hon_centrality(x, project = FALSE) and/or build_honem(x) attach a
#      value to every HIGHER-ORDER NODE of the fitted network.
#   2. each actor's own sequences are decoded against that same node set with
#      the longest-suffix rule (the rule BuildHON itself uses when it rewires
#      edges): at position i the actor occupies the longest-history node
#      "s[i-k+1] -> ... -> s[i]" that the network contains.
#   3. the actor's feature is the mean of the node value over that actor's
#      own visits -- a visit-frequency weighted mean, pooled over all of the
#      actor's sequences.
#   4. the actor-level table is regressed on the outcome.
#
# Step 2-3 is the whole aggregation. It is deliberately the only one on
# offer: an actor's exposure to the contexts the network singles out.
#
# What the intervals mean, stated plainly: the higher-order network is
# treated as FIXED. The features are estimated from the pooled data (which
# includes the actor's own sequences), so the reported standard errors
# condition on the estimated structure and do not propagate the uncertainty
# in it. They are the usual regression intervals for a fixed design, not
# intervals for "the effect of higher-order structure" in any wider sense.

# ---------------------------------------------------------------------------
# Input coercion
# ---------------------------------------------------------------------------

#' Per-actor sequences, keyed by actor
#'
#' @param sequences Named list, long data.frame, or wide data.frame with an
#'   id column.
#' @param by Name of the actor key.
#' @param action Long-format state column, or NULL.
#' @param time Long-format ordering column, or NULL.
#' @return Named list of character vectors, one per actor.
#' @noRd
.hoo_actor_sequences <- function(sequences, by, action, time) {
  if (!is.null(action)) {
    stopifnot(
      "`sequences` must be a data.frame when `action` is given" =
        is.data.frame(sequences),
      "`action` must name a column of `sequences`" =
        is.character(action) && length(action) == 1L &&
        action %in% names(sequences),
      "`by` must name a column of `sequences` when `action` is given" =
        by %in% names(sequences),
      "`time` must be NULL or name a column of `sequences`" =
        is.null(time) || (is.character(time) && length(time) == 1L &&
                            time %in% names(sequences))
    )
    state <- as.character(sequences[[action]])
    actor <- as.character(sequences[[by]])
    keep <- !is.na(actor) & !is.na(state)
    if (!any(keep)) {
      stop(errorCondition(
        "`sequences` has no rows with both an actor and an action.",
        class = "hypernets_bad_input", call = NULL))
    }
    state <- state[keep]
    actor <- actor[keep]
    # deterministic tie-breaking: row order is the explicit secondary key
    ord <- if (is.null(time)) {
      order(actor, seq_along(actor))
    } else {
      order(actor, sequences[[time]][keep], seq_along(actor))
    }
    return(split(state[ord], actor[ord]))
  }
  if (is.data.frame(sequences)) {
    stopifnot(
      "`by` must name a column of `sequences`, or `action` must be given" =
        by %in% names(sequences)
    )
    actor <- as.character(sequences[[by]])
    steps <- sequences[, setdiff(names(sequences), by), drop = FALSE]
    out <- lapply(seq_len(nrow(steps)), function(i) {
      v <- as.character(unlist(steps[i, ], use.names = FALSE))
      # wide format: NA is padding
      v[!is.na(v)]
    })
    names(out) <- actor
    return(out)
  }
  stopifnot(
    "`sequences` must be a named list, or a data.frame" =
      is.list(sequences) && !is.null(names(sequences)) &&
      !anyNA(names(sequences)) && all(nzchar(names(sequences)))
  )
  lapply(sequences, function(s) as.character(s[!is.na(s)]))
}

#' Outcome values keyed by actor
#'
#' @param outcome Named vector, or data.frame with the key and value columns.
#' @param by Name of the actor key.
#' @param outcome_col Name of the value column, or NULL.
#' @return Named numeric vector.
#' @noRd
.hoo_outcome_vector <- function(outcome, by, outcome_col) {
  if (is.data.frame(outcome)) {
    stopifnot(
      "`by` must name a column of `outcome`" = by %in% names(outcome),
      "`outcome_col` must be NULL or name a column of `outcome`" =
        is.null(outcome_col) ||
        (is.character(outcome_col) && length(outcome_col) == 1L &&
           outcome_col %in% names(outcome))
    )
    others <- setdiff(names(outcome), by)
    if (is.null(outcome_col)) {
      if (length(others) != 1L) {
        stop(errorCondition(
          paste0("`outcome` has ", length(others), " non-key columns (",
                 paste(others, collapse = ", "),
                 "); name the one to model with `outcome_col`."),
          class = "hypernets_bad_input", call = NULL))
      }
      outcome_col <- others
    }
    keys <- as.character(outcome[[by]])
    values <- outcome[[outcome_col]]
  } else {
    stopifnot(
      "`outcome` must be a data.frame or a named vector keyed by actor" =
        !is.null(names(outcome)) && length(outcome) >= 1L
    )
    keys <- names(outcome)
    values <- unname(outcome)
  }
  if (is.factor(values)) values <- as.character(values)
  if (is.logical(values)) values <- as.numeric(values)
  stopifnot(
    "the outcome column must be numeric (or logical)" = is.numeric(values),
    "`outcome` keys must not be missing" = !anyNA(keys)
  )
  if (anyDuplicated(keys) > 0L) {
    dup <- unique(keys[duplicated(keys)])
    stop(errorCondition(
      paste0("`outcome` has duplicate keys (", paste(utils::head(dup, 5L),
                                                     collapse = ", "),
             "); it must hold one value per actor."),
      class = "hypernets_bad_input", call = NULL))
  }
  stats::setNames(as.numeric(values), keys)
}

#' Cluster labels keyed by actor
#'
#' @param nested_in Column name of `outcome`, or a named vector.
#' @param outcome The original outcome argument.
#' @param by Name of the actor key.
#' @param actors Character vector of actors in the analysis sample.
#' @return Character vector, one per actor in `actors`.
#' @noRd
.hoo_cluster_vector <- function(nested_in, outcome, by, actors) {
  if (is.character(nested_in) && length(nested_in) == 1L &&
      is.data.frame(outcome) && nested_in %in% names(outcome)) {
    lookup <- stats::setNames(as.character(outcome[[nested_in]]),
                              as.character(outcome[[by]]))
  } else {
    stopifnot(
      "`nested_in` must name a column of `outcome`, or be a named vector" =
        !is.null(names(nested_in)) && length(nested_in) >= 1L
    )
    lookup <- stats::setNames(as.character(nested_in), names(nested_in))
  }
  out <- unname(lookup[actors])
  if (anyNA(out)) {
    stop(errorCondition(
      paste0("`nested_in` is missing for ", sum(is.na(out)),
             " of the ", length(actors),
             " actors in the analysis sample; a cluster-robust covariance ",
             "needs every observation assigned to a cluster."),
      class = "hypernets_bad_input", call = NULL))
  }
  out
}

# ---------------------------------------------------------------------------
# Actor visits over the higher-order node set
# ---------------------------------------------------------------------------

#' Longest-suffix decoding of one sequence onto a higher-order node set
#'
#' At position `i` the actor occupies the longest-history node
#' `s[i-k+1] -> ... -> s[i]` present in `index`, exactly the rule BuildHON
#' uses to rewire edges onto higher-order nodes.
#'
#' @param s Character vector (one sequence).
#' @param index Named integer vector mapping node name to column.
#' @param max_order Integer, longest history to try.
#' @return Integer vector of node indices, `NA` where nothing matched.
#' @noRd
.hoo_decode <- function(s, index, max_order) {
  n <- length(s)
  if (n == 0L) return(integer(0L))
  k_max <- min(as.integer(max_order), n)
  rolls <- if (k_max <= 1L) {
    list(s)
  } else {
    Reduce(function(prev, k) {
      cur <- rep(NA_character_, n)
      pos <- seq.int(k, n)
      cur[pos] <- paste(prev[pos - 1L], s[pos], sep = " -> ")
      cur
    }, seq.int(2L, k_max), accumulate = TRUE, init = s)
  }
  hits <- lapply(rolls, function(nm) unname(index[nm]))
  # longer history wins: later (higher-order) hits overwrite earlier ones
  Reduce(function(acc, new) {
    ok <- !is.na(new)
    acc[ok] <- new[ok]
    acc
  }, hits)
}

#' Visit counts of every actor over the higher-order node set
#'
#' @param seqs Named list of character vectors.
#' @param nodes Character vector of higher-order node names.
#' @param max_order Integer.
#' @return List with `counts` (actors x nodes integer matrix), `n_events`
#'   and `n_unmatched` (named integer vectors).
#' @noRd
.hoo_visit_counts <- function(seqs, nodes, max_order) {
  index <- stats::setNames(seq_along(nodes), nodes)
  per <- lapply(seqs, function(s) {
    idx <- .hoo_decode(s, index, max_order)
    matched <- idx[!is.na(idx)]
    list(counts = tabulate(matched, nbins = length(nodes)),
         n_events = length(s),
         n_unmatched = length(s) - length(matched))
  })
  counts <- do.call(rbind, lapply(per, `[[`, "counts"))
  dimnames(counts) <- list(names(seqs), nodes)
  list(
    counts = counts,
    n_events = vapply(per, `[[`, integer(1L), "n_events"),
    n_unmatched = vapply(per, `[[`, integer(1L), "n_unmatched")
  )
}

#' Node-level predictor values from the verbs that already exist
#'
#' @param x A `net_hon`.
#' @param features Character vector, a subset of "centrality" / "honem".
#' @param centrality_type Character vector passed to `hon_centrality()`.
#' @param honem_dim Integer embedding dimension.
#' @param weighted Logical, passed to `hon_centrality()`.
#' @param max_paths Passed to `hon_centrality()`.
#' @return Numeric matrix, one row per higher-order node, one column per
#'   predictor, with the node names as row names.
#' @noRd
.hoo_node_values <- function(x, features, centrality_type, honem_dim,
                             weighted, max_paths) {
  nodes <- rownames(x$matrix)
  parts <- list()
  if ("centrality" %in% features) {
    cen <- hon_centrality(x, type = centrality_type, project = FALSE,
                          weighted = weighted, max_paths = max_paths)
    m <- as.matrix(cen[, centrality_type, drop = FALSE])
    rownames(m) <- cen$node
    parts$centrality <- m[nodes, , drop = FALSE]
  }
  if ("honem" %in% features) {
    emb <- build_honem(x, dim = honem_dim)
    m <- emb$embeddings
    colnames(m) <- paste0("honem_", seq_len(ncol(m)))
    parts$honem <- m[nodes, , drop = FALSE]
  }
  do.call(cbind, unname(parts))
}

# ---------------------------------------------------------------------------
# Least squares / GLM with an optional cluster-robust covariance
# ---------------------------------------------------------------------------

#' Pivot-aware inverse of the cross-product implied by a QR decomposition
#'
#' @param qr_obj A `qr` object of X (or of sqrt(W) X).
#' @param k Integer number of columns.
#' @return The k x k matrix `(X'X)^-1` in the ORIGINAL column order.
#' @noRd
.hoo_xtx_inv <- function(qr_obj, k) {
  r <- qr.R(qr_obj)[seq_len(k), seq_len(k), drop = FALSE]
  inv <- chol2inv(r)
  piv <- qr_obj$pivot[seq_len(k)]
  out <- matrix(0, k, k)
  out[piv, piv] <- inv
  out
}

#' Fit the actor-level model and build its covariance
#'
#' Gaussian: ordinary least squares by QR. Binomial / Poisson: IRLS through
#' `stats::glm.fit()` with the canonical link.
#'
#' Cluster-robust covariance (when `cluster` is given):
#'   V = c * B (sum_g u_g u_g') B,  u_g = sum_{i in g} x_i e_i,
#' with bread B = (X'X)^-1 (gaussian) or (X'WX)^-1 (GLM, W the IRLS working
#' weights), e = y - mu, and the CR1S ("Stata") small-sample factor
#'   c = G/(G-1) * (n-1)/(n-k).
#' This is `sandwich::vcovCL(type = "HC1", cadjust = TRUE)`.
#'
#' @param X Design matrix including the intercept column.
#' @param y Numeric response.
#' @param family "gaussian", "binomial" or "poisson".
#' @param cluster Character vector of cluster labels, or NULL.
#' @return List with coefficients, vcov, df, distribution and fit statistics.
#' @noRd
.hoo_fit <- function(X, y, family, cluster, vcov = "CR3") {
  n <- nrow(X)
  k <- ncol(X)
  converged <- TRUE
  if (identical(family, "gaussian")) {
    qr_obj <- qr(X)
    beta <- qr.coef(qr_obj, y)
    mu <- as.numeric(X %*% beta)
    resid <- y - mu
    bread <- .hoo_xtx_inv(qr_obj, k)
    work_w <- rep(1, n)
    dispersion <- sum(resid^2) / (n - k)
    vcov_model <- dispersion * bread
    rss <- sum(resid^2)
    tss <- sum((y - mean(y))^2)
    r_squared <- if (tss > 0) 1 - rss / tss else NA_real_
    fit_stats <- list(
      sigma = sqrt(dispersion),
      r_squared = r_squared,
      adj_r_squared = if (is.na(r_squared)) NA_real_ else
        1 - (1 - r_squared) * (n - 1) / (n - k),
      deviance = rss,
      null_deviance = tss,
      aic = n * log(rss / n) + 2 * (k + 1) + n + n * log(2 * pi)
    )
  } else {
    fam <- switch(family,
                  binomial = stats::binomial(),
                  poisson = stats::poisson())
    fit <- stats::glm.fit(x = X, y = y, family = fam,
                          control = stats::glm.control(epsilon = 1e-12,
                                                        maxit = 100L))
    # epsilon is tighter than glm()'s 1e-8: IRLS stops on the change in
    # deviance, which is quadratic at the optimum, so 1e-8 pins the
    # coefficients only to ~1e-4 and the CR3 leverage inherits that error.
    converged <- isTRUE(fit$converged)
    beta <- fit$coefficients
    mu <- fit$fitted.values
    resid <- y - mu
    bread <- .hoo_xtx_inv(fit$qr, k)
    work_w <- fit$weights
    vcov_model <- bread
    fit_stats <- list(
      sigma = NA_real_,
      r_squared = NA_real_,
      adj_r_squared = NA_real_,
      deviance = fit$deviance,
      null_deviance = fit$null.deviance,
      aic = fit$aic
    )
  }
  names(beta) <- colnames(X)

  n_clusters <- NA_integer_
  if (is.null(cluster)) {
    vcov_use <- vcov_model
    vcov_type <- "model"
    df <- if (identical(family, "gaussian")) n - k else Inf
  } else {
    groups <- split(seq_len(n), cluster)
    n_clusters <- length(groups)
    if (identical(vcov, "CR3")) {
      # Cluster jackknife. Each cluster's score is inflated by
      # (I - H_gg)^-1, H_gg = X_g B X_g', which undoes the leverage that
      # shrinks a cluster's own residuals toward zero. No G/(G-1) factor:
      # the leverage term is the small-sample correction, and this is
      # exactly sandwich::vcovCL(type = "HC3").
      # For a GLM the bread is (X'WX)^-1, so the leverage must carry the
      # working weights too: H_gg = X_g B X_g' W_g. Dropping W_g is invisible
      # for gaussian (W = I) but mis-scales binomial/Poisson CR3.
      u <- t(vapply(groups, function(idx) {
        xg <- X[idx, , drop = FALSE]
        wg <- work_w[idx]
        leverage <- diag(length(idx)) - xg %*% bread %*% t(xg * wg)
        # I - H_gg is similar to the symmetric PSD matrix
        # I - W^1/2 X B X' W^1/2, so that matrix's smallest eigenvalue is the
        # right singularity test. rcond() is not: for a singleton cluster
        # the matrix is 1x1 and rcond() returns 1 however close to zero the
        # entry is, so a leverage-1 cluster slips through.
        xs <- xg * sqrt(wg)
        smallest <- min(eigen(diag(length(idx)) - xs %*% bread %*% t(xs),
                              symmetric = TRUE, only.values = TRUE)$values)
        if (smallest <= sqrt(.Machine$double.eps)) {
          stop(errorCondition(
            paste0("cluster '", cluster[idx[1L]], "' has ", length(idx),
                   " observation(s) and leverage 1: the jackknife is ",
                   "undefined there. Use `vcov = \"CR1\"`."),
            class = c("hypernets_rank_deficient", "hypernets_bad_input"),
            call = NULL))
        }
        crossprod(xg, qr.solve(leverage, resid[idx]))[, 1L]
      }, numeric(k)))
      vcov_use <- bread %*% crossprod(u) %*% bread
      vcov_type <- "CR3"
    } else {
      u <- rowsum(X * resid, group = cluster, reorder = TRUE)
      adj <- (n_clusters / (n_clusters - 1)) * ((n - 1) / (n - k))
      vcov_use <- adj * (bread %*% crossprod(u) %*% bread)
      vcov_type <- "CR1S"
    }
    df <- n_clusters - 1
  }
  dimnames(vcov_use) <- list(colnames(X), colnames(X))

  se <- sqrt(diag(vcov_use))
  statistic <- beta / se
  distribution <- if (is.finite(df)) "t" else "normal"
  p <- if (is.finite(df)) {
    2 * stats::pt(-abs(statistic), df = df)
  } else {
    2 * stats::pnorm(-abs(statistic))
  }
  crit <- if (is.finite(df)) stats::qt(0.975, df = df) else stats::qnorm(0.975)

  c(list(
    coefficients = beta,
    std_error = se,
    statistic = statistic,
    p = p,
    conf_low = beta - crit * se,
    conf_high = beta + crit * se,
    vcov = vcov_use,
    vcov_model = vcov_model,
    vcov_type = vcov_type,
    df = df,
    distribution = distribution,
    n_clusters = n_clusters,
    fitted = mu,
    residuals = resid,
    converged = converged
  ), fit_stats)
}

# ---------------------------------------------------------------------------
# Interpretation layer: per-SD effects, percentile contrasts, labels
# ---------------------------------------------------------------------------
#
# Both tables are exact functions of the stored fit -- no refit. The per-SD
# effect rescales each coefficient by sd(x_j) / s_j, where s_j is the scale the
# feature entered the design on (sd(x_j) under `standardize = TRUE`, 1 on the
# raw scale); a positive rescaling of one column maps its coefficient, its
# standard error and both interval ends by the same factor, so the interval is
# exactly the one a refit on z-scored features would report. The contrast is
# a linear combination c' beta of the coefficients, with variance c' V c.

#' A short plain-language name for an argument, for captions
#'
#' @param expr The unevaluated argument (from `substitute()`).
#' @param value The evaluated argument.
#' @param by Key column to ignore when `value` is a data.frame, or NULL.
#' @param col An explicit column name, or NULL.
#' @param fallback Label when nothing better is available.
#' @return A single string.
#' @noRd
.hoo_arg_label <- function(expr, value, by, col, fallback) {
  if (is.character(col) && length(col) == 1L && !is.na(col)) return(col)
  if (is.data.frame(value)) {
    others <- setdiff(names(value), by)
    if (length(others) == 1L) return(others)
  }
  if (is.symbol(expr)) return(deparse(expr))
  fallback
}

#' Human-readable names for the model's features
#'
#' @param features Character vector of feature (column) names.
#' @return Character vector of the same length.
#' @noRd
.hoo_feature_label <- function(features) {
  known <- c(pagerank = "PageRank exposure",
             betweenness = "Betweenness exposure",
             closeness = "Closeness exposure")
  out <- unname(known[features])
  is_honem <- grepl("^honem_[0-9]+$", features)
  out[is_honem] <- paste("HONEM dimension", sub("^honem_", "",
                                                features[is_honem]))
  out[is.na(out)] <- features[is.na(out)]
  out
}

#' What one unit of the per-SD effect is, per family
#'
#' @param family "gaussian", "binomial" or "poisson".
#' @return A single string.
#' @noRd
.hoo_effect_scale <- function(family) {
  switch(family, gaussian = "outcome units", binomial = "odds ratio",
         poisson = "rate ratio")
}

#' Per-SD effect columns for the terms table
#'
#' @param x A `net_outcome`.
#' @param terms The terms table (with or without the intercept).
#' @return `terms` with `feature_sd`, `effect_per_sd`, `per_sd_low`,
#'   `per_sd_high` and `scale` appended (`NA` on the intercept row).
#' @noRd
.hoo_add_per_sd <- function(x, terms) {
  feature_sd <- unname(x$feature_sd[terms$feature])
  # coefficient per design unit -> per SD of the raw feature
  factor <- feature_sd / unname(x$scale[terms$feature])
  back <- if (identical(x$family, "gaussian")) identity else exp
  terms$feature_sd <- feature_sd
  terms$effect_per_sd <- back(terms$estimate * factor)
  terms$per_sd_low <- back(terms$conf_low * factor)
  terms$per_sd_high <- back(terms$conf_high * factor)
  terms$scale <- ifelse(is.na(feature_sd), NA_character_,
                        .hoo_effect_scale(x$family))
  terms
}

#' Predicted outcome at the 25th vs 75th percentile of each feature
#'
#' Every other feature is held at its sample mean. The linear predictor at a
#' design row c is eta = c' beta with Var(eta) = c' V c (V the stored,
#' possibly cluster-robust, covariance). Predicted values and their intervals
#' are formed on the link scale and mapped through the inverse link. The
#' difference mu(c_high) - mu(c_low) is on the outcome scale; its standard
#' error is the delta method with gradient
#' g = mu'(eta_high) c_high - mu'(eta_low) c_low, which for the identity link
#' reduces to the exact |beta_j| (q75 - q25) / s_j interval. For a GLM the
#' link-scale contrast is also reported exponentiated (odds / rate ratio).
#'
#' @param x A `net_outcome`.
#' @return One row per feature (see `as.data.frame.net_outcome()`).
#' @noRd
.hoo_contrasts <- function(x) {
  fam <- switch(x$family, gaussian = stats::gaussian(),
                binomial = stats::binomial(), poisson = stats::poisson())
  beta <- x$terms$estimate
  names(beta) <- x$terms$feature
  V <- x$vcov[names(beta), names(beta)]
  raw <- as.matrix(x$features[, x$feature_names, drop = FALSE])
  # the design holds (raw - centre) / scale; others sit at their sample mean
  to_design <- function(v, j) (v - x$centre[[j]]) / x$scale[[j]]
  base_row <- c(1, vapply(x$feature_names, function(j)
    to_design(mean(raw[, j]), j), numeric(1L)))
  rows <- lapply(x$feature_names, function(j) {
    q <- unname(stats::quantile(raw[, j], c(0.25, 0.75), type = 7))
    c_low <- base_row
    c_high <- base_row
    c_low[[j]] <- to_design(q[1L], j)
    c_high[[j]] <- to_design(q[2L], j)
    eta <- c(sum(c_low * beta), sum(c_high * beta))
    se_eta <- sqrt(c(drop(c_low %*% V %*% c_low),
                     drop(c_high %*% V %*% c_high)))
    mu <- fam$linkinv(eta)
    grad <- fam$mu.eta(eta[2L]) * c_high - fam$mu.eta(eta[1L]) * c_low
    se_diff <- sqrt(drop(grad %*% V %*% grad))
    difference <- mu[2L] - mu[1L]
    # link-scale contrast (c_high - c_low)' beta, exponentiated for a GLM:
    # the odds / rate ratio of the two actors, whose interval stays inside
    # the parameter space where the Wald difference interval need not
    dc <- c_high - c_low
    link_diff <- sum(dc * beta)
    se_link <- sqrt(drop(dc %*% V %*% dc))
    ratio <- if (identical(x$family, "gaussian")) rep(NA_real_, 3L) else
      exp(link_diff + c(0, -1, 1) * x$crit * se_link)
    data.frame(
      feature = j,
      low_value = q[1L],
      high_value = q[2L],
      predicted_low = mu[1L],
      predicted_low_ci_low = fam$linkinv(eta[1L] - x$crit * se_eta[1L]),
      predicted_low_ci_high = fam$linkinv(eta[1L] + x$crit * se_eta[1L]),
      predicted_high = mu[2L],
      predicted_high_ci_low = fam$linkinv(eta[2L] - x$crit * se_eta[2L]),
      predicted_high_ci_high = fam$linkinv(eta[2L] + x$crit * se_eta[2L]),
      difference = difference,
      std_error = se_diff,
      conf_low = difference - x$crit * se_diff,
      conf_high = difference + x$crit * se_diff,
      ratio = ratio[1L],
      ratio_low = ratio[2L],
      ratio_high = ratio[3L],
      scale = switch(x$family, gaussian = "outcome units",
                     binomial = "probability", poisson = "expected count"),
      stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' One line of model context for a plot subtitle
#'
#' @param x A `net_outcome`.
#' @return A single string.
#' @noRd
.hoo_model_context <- function(x) {
  who <- if (identical(x$vcov_type, "model")) {
    sprintf("%d actors, model-based covariance", x$n)
  } else {
    sprintf("%d actors in %d clusters (%s), %s cluster-robust covariance",
            x$n, x$cluster_n, x$cluster_label, x$vcov_type)
  }
  fit_line <- if (identical(x$family, "gaussian")) {
    sprintf("R-squared %.2f", x$r_squared)
  } else {
    sprintf("deviance explained %.1f%%",
            100 * (1 - x$deviance / x$null_deviance))
  }
  sprintf("%s\n%s; %s over %d features", who, fit_line, x$p_adjust,
          length(x$feature_names))
}

#' How the intervals were computed, in words
#'
#' @param x A `net_outcome`.
#' @return A single string.
#' @noRd
.hoo_interval_words <- function(x) {
  if (identical(x$vcov_type, "model")) {
    "95% CIs, model-based"
  } else {
    sprintf("95%% CIs, cluster-robust %s over %d clusters (%s)",
            x$vcov_type, x$cluster_n, x$cluster_label)
  }
}

# ---------------------------------------------------------------------------
# Main verb
# ---------------------------------------------------------------------------

#' Relate higher-order structure to an actor-level outcome
#'
#' @description
#' Regresses an outcome measured once per actor on per-actor features built
#' from a fitted higher-order network. It is an assembly and inference layer:
#' the predictors come from \code{\link{hon_centrality}()} and
#' \code{\link{build_honem}()}, and the only new step is the regression.
#'
#' @details
#' \strong{How a per-actor feature is built.} Node-level values are attached
#' to every higher-order node of \code{x} -- a centrality from
#' \code{hon_centrality(x, project = FALSE)}, a HONEM coordinate from
#' \code{build_honem(x)}, or both. Each actor's own sequences are then
#' decoded against that same node set with the \emph{longest-suffix} rule: at
#' position \eqn{i} the actor occupies the longest-history node
#' \code{"s[i-k+1] -> ... -> s[i]"} that the network contains, which is the
#' rule BuildHON itself uses when it rewires edges onto higher-order nodes.
#' The actor's feature is the mean of the node value over that actor's own
#' visits, pooling every sequence that actor contributed:
#' \deqn{f_a = \sum_v c_{av} \, \mathrm{value}(v) / \sum_v c_{av}}
#' with \eqn{c_{av}} the number of times actor \eqn{a} occupies node
#' \eqn{v}. Positions whose state never entered the network (pruned by
#' \code{min_freq}, or unseen) match no node, are excluded from the mean, and
#' are counted in \code{n_unmatched} in the features table.
#'
#' \strong{Standardization.} With \code{standardize = TRUE} (the default)
#' every feature is z-scored over the analysis sample, so an estimate is the
#' change in the outcome per one standard deviation of that feature, and
#' features on wildly different scales (a PageRank near \eqn{10^{-3}}, a HONEM
#' coordinate near 1) are comparable. The outcome is never transformed.
#'
#' \strong{Reading the result.} Whatever \code{standardize} is, the accessor
#' reports the effect of a one-SD increase in each feature,
#' \eqn{\hat\beta_j \, \mathrm{sd}(x_j) / s_j} with \eqn{s_j} the scale
#' the feature entered the design on, and scales the interval by the same
#' positive factor -- exactly the estimate and interval a refit on z-scored
#' features gives, because the (cluster-robust) covariance is equivariant to
#' rescaling a column. For binomial and Poisson fits it is exponentiated to an
#' odds or rate ratio per SD. \code{as.data.frame(fit, what = "contrasts")}
#' turns the same fit into two actors a reader can picture: one at the 25th
#' and one at the 75th percentile of a feature, every other feature at its
#' mean, with the predicted outcome for each and the delta-method interval
#' \eqn{\sqrt{g' V g}} for their difference. The comparison is always
#' across actors: a feature is the average value of the states an actor
#' visited, so a one-SD contrast compares actors whose trajectories run
#' through more central states with actors whose trajectories do not.
#'
#' \strong{Inference.} Gaussian outcomes are fitted by ordinary least squares
#' (QR); \code{family = "binomial"} and \code{"poisson"} by IRLS with the
#' canonical link. Without \code{nested_in} the covariance is the usual
#' model-based one, with \eqn{t_{n-k}} intervals for the gaussian family and
#' normal intervals otherwise. With \code{nested_in} the covariance is the
#' cluster-robust sandwich
#' \deqn{V = c \, B \left(\sum_g u_g u_g'\right) B, \quad
#'       u_g = \sum_{i \in g} x_i e_i}
#' with bread \eqn{B = (X'X)^{-1}} (gaussian) or \eqn{(X'WX)^{-1}} (GLM,
#' \eqn{W} the IRLS working weights), \eqn{e = y - \hat\mu}, and the CR1S
#' ("Stata") small-sample correction \eqn{c = \frac{G}{G-1}\cdot
#' \frac{n-1}{n-k}}; intervals then use \eqn{t_{G-1}}. This is exactly
#' \code{sandwich::vcovCL(type = "HC1", cadjust = TRUE)}, and the test suite
#' checks it against that package where it is installed. Cluster-robust
#' standard errors are badly biased with few clusters, so fewer than
#' \code{min_clusters} clusters raises a \code{hypernets_few_clusters}
#' warning.
#'
#' \strong{Measured coverage.} In the package's own simulation
#' (\code{tests/testthat/test-memory_outcome.R}, 2000 replicates with the
#' design held fixed and only the outcome redrawn) the model-based 95%
#' interval covered a known effect 95.3% of the time, and the default
#' cluster-robust interval (\code{vcov = "CR3"}) covered 95.3% with 40
#' clusters \emph{when the predictor is itself cluster-correlated} -- the
#' hardest case, and the one where ignoring the clustering is catastrophic
#' (the model-based interval covered only 72.0% there).
#'
#' \code{vcov = "CR1"} covers 93.3% on the same design. That shortfall is a
#' property of the CR1 estimator at a moderate number of clusters, not of
#' this implementation: the identical simulation run through
#' \code{sandwich::vcovCL()} gives the same 93.3%, with 94.1% for CR2 and
#' 95.3% for CR3. CR3 is therefore the default, and CR1 is kept because it is
#' what most published results report -- read a CR1 interval from a few dozen
#' clusters as slightly optimistic.
#'
#' \strong{What the intervals do and do not cover.} The higher-order network
#' is treated as fixed. The features are estimated from the pooled data,
#' which includes each actor's own sequences, so the intervals are the usual
#' regression intervals \emph{conditional on the estimated structure}; they
#' do not propagate the uncertainty in the network itself. Use
#' \code{\link{bootstrap_hon}()} for that uncertainty.
#'
#' \strong{Collinearity.} A rank-deficient design raises
#' \code{hypernets_rank_deficient} rather than silently dropping terms, and a
#' scaled condition index above 30 (Belsley, Kuh and Welsch 1980) warns with
#' the same class.
#'
#' \strong{Conditions raised.} \code{hypernets_bad_input} (an input that
#' cannot be modelled at all: a malformed outcome, duplicate actor keys, a
#' missing cluster label, fewer actors than terms, fewer clusters than terms);
#' \code{hypernets_rank_deficient} (error on an exactly collinear or constant
#' design, warning on a near-collinear one);
#' \code{hypernets_dropped_actors} (warning: events that matched no node, or
#' actors excluded from the fit -- both are listed in the result);
#' \code{hypernets_few_clusters} (warning); \code{hypernets_no_converge}
#' (warning: IRLS hit its iteration bound); \code{hypernets_degenerate_fit}
#' (warning: separation in a binomial fit, whose standard errors are then
#' meaningless). Nothing is dropped or rescued silently.
#'
#' @param x A \code{net_hon} from \code{\link{build_hon}()}.
#' @param outcome One value per actor: a named numeric vector keyed by actor,
#'   or a data.frame holding the actor key in column \code{by} and the
#'   outcome in \code{outcome_col} (or in its only other column).
#' @param sequences The sequences the network was built from, keyed by actor:
#'   a named list of character vectors, a long data.frame together with
#'   \code{action} (and optionally \code{time}), or a wide data.frame whose
#'   \code{by} column holds the actor and whose remaining columns are time
#'   steps. Required, because a \code{net_hon} does not carry the actor
#'   identity of the trajectories it was built from.
#' @param by Name of the actor key, used for the column in \code{outcome},
#'   \code{sequences} and the returned features table. Default \code{"actor"}.
#' @param features Which per-actor predictors to build: any of
#'   \code{"centrality"} and \code{"honem"}, or a data.frame of ready-made
#'   per-actor features holding the actor key in column \code{by} and one
#'   numeric column per feature.
#' @param nested_in Optional clustering variable for cluster-robust standard
#'   errors: the name of a column of \code{outcome}, or a vector named by
#'   actor. \code{NULL} (default) uses model-based standard errors.
#' @param family \code{"gaussian"} (default), \code{"binomial"} or
#'   \code{"poisson"}; the last two use the canonical link.
#' @param outcome_col Name of the outcome column when \code{outcome} is a
#'   data.frame with more than one non-key column.
#' @param action,time Long-format column names of \code{sequences}:
#'   \code{action} holds the state of each event and \code{time} orders
#'   events within an actor (row order when \code{NULL}). Leave both
#'   \code{NULL} for a named list or a wide data.frame.
#' @param centrality_type Which centralities to attach to the higher-order
#'   nodes, passed to \code{\link{hon_centrality}()}. Betweenness and
#'   closeness enumerate shortest paths and are the expensive ones.
#' @param honem_dim Integer. Number of HONEM dimensions to use as features
#'   when \code{"honem"} is requested. Default 2.
#' @param weighted Logical. Passed to \code{\link{hon_centrality}()}: use
#'   the higher-order edge weights instead of the binary pattern. Default
#'   \code{FALSE}, matching \code{hon_centrality()}.
#' @param standardize Logical. Z-score each feature over the analysis sample.
#'   Default \code{TRUE}.
#' @param vcov Cluster covariance estimator, used only when `nested_in` is
#'   given. `"CR3"` (default) is the cluster jackknife, matching
#'   `sandwich::vcovCL(type = "HC3")`; `"CR1"` is the Stata estimator,
#'   `sandwich::vcovCL(type = "HC1", cadjust = TRUE)`. CR1 under-covers at
#'   moderate cluster counts -- measured 0.933 against a nominal 0.95 at 40
#'   clusters with a cluster-correlated regressor, where CR3 reached 0.953 --
#'   which is why CR3 is the default. CR1 is kept because it is what most
#'   published results report. CR3 is undefined for a cluster whose leverage
#'   is 1 (typically a singleton cluster) and raises there.
#' @param p_adjust Multiplicity correction applied across the feature
#'   coefficients, passed to \code{\link[stats]{p.adjust}}. Default
#'   \code{"BH"}; \code{"none"} switches it off and is recorded as such.
#' @param min_clusters Integer. Warn below this many clusters. Default 20.
#' @param max_paths Passed to \code{\link{hon_centrality}()}.
#' @param ... Ignored; present so the verb can grow arguments without
#'   breaking positional calls.
#'
#' @return An object of class \code{net_outcome}. Reach into it with
#'   \code{\link[=as.data.frame.net_outcome]{as.data.frame}()}, which returns
#'   one row per feature with \code{feature}, \code{estimate},
#'   \code{std_error}, \code{conf_low}, \code{conf_high}, \code{statistic},
#'   \code{p}, \code{p_adj} and \code{n}, followed by the interpretable
#'   per-SD columns \code{feature_sd}, \code{effect_per_sd},
#'   \code{per_sd_low}, \code{per_sd_high} and \code{scale} (the effect of a
#'   one-SD increase in the feature, in outcome units for gaussian, as an
#'   odds or rate ratio otherwise); \code{what = "terms"} adds the
#'   intercept, \code{what = "contrasts"} gives the predicted outcome at the
#'   25th against the 75th percentile of each feature,
#'   \code{what = "features"} returns the per-actor design table, and
#'   \code{what = "dropped"} the actors that were excluded and why.
#'   \code{summary()} reports the model, the covariance and the correction;
#'   \code{plot()} draws the per-SD effects (default), the percentile
#'   contrasts (\code{what = "contrasts"}) or the raw coefficients
#'   (\code{what = "raw"}).
#'
#' @references
#' Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). \emph{Regression
#' Diagnostics: Identifying Influential Data and Sources of Collinearity}.
#' Wiley.
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate. \emph{Journal of the Royal Statistical Society B}, 57(1), 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#'
#' Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to
#' cluster-robust inference. \emph{Journal of Human Resources}, 50(2),
#' 317-372. \doi{10.3368/jhr.50.2.317}
#'
#' Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
#' generalized linear models. \emph{Biometrika}, 73(1), 13-22.
#' \doi{10.1093/biomet/73.1.13}
#'
#' Scholtes, I., Wider, N., & Garas, A. (2016). Higher-order aggregate
#' networks in the analysis of temporal networks. \emph{European Physical
#' Journal B}, 89, 61. \doi{10.1140/epjb/e2016-60663-0}
#'
#' @seealso \code{\link{hon_centrality}()}, \code{\link{build_honem}()},
#'   \code{\link{bootstrap_hon}()}
#'
#' @examples
#' # 60 actors, each a different mixture of four work motifs; the motifs
#' # carry the second-order dependence the network will pick up.
#' set.seed(42)
#' motifs <- list(c("plan", "code", "test"), c("debug", "code", "debug"),
#'                c("read", "plan", "code"), c("test", "read", "plan"))
#' actors <- sprintf("s%02d", seq_len(60))
#' mix <- matrix(stats::runif(60 * 4), nrow = 60)
#' seqs <- stats::setNames(
#'   lapply(seq_len(60), function(i)
#'     unlist(motifs[sample(4, 10, replace = TRUE, prob = mix[i, ])],
#'            use.names = FALSE)),
#'   actors)
#' hon <- build_hon(seqs, max_order = 2)
#'
#' scores <- stats::setNames(60 + 8 * mix[, 1] + stats::rnorm(60, sd = 2),
#'                           actors)
#' fit <- hon_outcome(hon, outcome = scores, sequences = seqs)
#' fit
#' as.data.frame(fit, sort_by = "p_adj")
#' as.data.frame(fit, what = "contrasts")
#'
#' \donttest{
#' # actors nested in classes: cluster-robust standard errors
#' classes <- stats::setNames(rep(sprintf("c%02d", seq_len(20)), each = 3),
#'                            actors)
#' nested <- hon_outcome(hon, outcome = scores, sequences = seqs,
#'                       nested_in = classes)
#' summary(nested)
#' plot(nested)
#' plot(nested, what = "contrasts")
#' }
#'
#' @export
hon_outcome <- function(x, outcome, sequences, by = "actor",
                        features = "centrality", nested_in = NULL,
                        family = c("gaussian", "binomial", "poisson"),
                        outcome_col = NULL, action = NULL, time = NULL,
                        centrality_type = c("pagerank", "betweenness",
                                            "closeness"),
                        honem_dim = 2L, weighted = FALSE,
                        standardize = TRUE, vcov = c("CR3", "CR1"),
                        p_adjust = "BH", min_clusters = 20L,
                        max_paths = 1e6, ...) {
  stopifnot(
    "`x` must be a net_hon object from build_hon()" = inherits(x, "net_hon"),
    "`by` must be a single non-empty string" =
      is.character(by) && length(by) == 1L && !is.na(by) && nzchar(by),
    "`features` must be a character vector or a per-actor data.frame" =
      is.data.frame(features) || is.character(features),
    "`honem_dim` must be a single whole number >= 1" =
      is.numeric(honem_dim) && length(honem_dim) == 1L &&
      is.finite(honem_dim) && honem_dim == round(honem_dim) && honem_dim >= 1,
    "`weighted` must be TRUE or FALSE" =
      is.logical(weighted) && length(weighted) == 1L && !is.na(weighted),
    "`standardize` must be TRUE or FALSE" =
      is.logical(standardize) && length(standardize) == 1L &&
      !is.na(standardize),
    "`p_adjust` must be a single string" =
      is.character(p_adjust) && length(p_adjust) == 1L && !is.na(p_adjust),
    "`min_clusters` must be a single whole number >= 2" =
      is.numeric(min_clusters) && length(min_clusters) == 1L &&
      is.finite(min_clusters) && min_clusters == round(min_clusters) &&
      min_clusters >= 2
  )
  family <- match.arg(family)
  vcov <- match.arg(vcov)
  # plain-language names for plot captions, read before `outcome` and
  # `nested_in` are evaluated
  outcome_label <- if (missing(outcome)) "outcome" else
    .hoo_arg_label(substitute(outcome), outcome, by, outcome_col,
                   fallback = "outcome")
  cluster_label <- .hoo_arg_label(substitute(nested_in), nested_in, NULL,
                                  if (is.character(nested_in) &&
                                      length(nested_in) == 1L) nested_in,
                                  fallback = "clusters")
  p_adjust <- match.arg(p_adjust, stats::p.adjust.methods)
  honem_dim <- as.integer(honem_dim)

  # ---- per-actor features ------------------------------------------------
  if (is.data.frame(features)) {
    feature_spec <- "supplied"
    stopifnot(
      "`by` must name a column of the `features` data.frame" =
        by %in% names(features)
    )
    feature_names <- setdiff(names(features), by)
    if (length(feature_names) == 0L) {
      stop(errorCondition(
        "The `features` data.frame has no feature columns beside the key.",
        class = "hypernets_bad_input", call = NULL))
    }
    is_num <- vapply(features[feature_names], is.numeric, logical(1L))
    if (!all(is_num)) {
      stop(errorCondition(
        paste0("Every feature column must be numeric; these are not: ",
               paste(feature_names[!is_num], collapse = ", "), "."),
        class = "hypernets_bad_input", call = NULL))
    }
    actors <- as.character(features[[by]])
    if (anyDuplicated(actors) > 0L) {
      stop(errorCondition(
        "The `features` data.frame must hold one row per actor.",
        class = "hypernets_bad_input", call = NULL))
    }
    feature_matrix <- as.matrix(features[, feature_names, drop = FALSE])
    rownames(feature_matrix) <- actors
    support <- data.frame(
      key = actors, n_events = NA_integer_, n_unmatched = NA_integer_,
      stringsAsFactors = FALSE)
  } else {
    features <- unique(features)
    unknown <- setdiff(features, c("centrality", "honem"))
    if (length(unknown) > 0L || length(features) == 0L) {
      stop(errorCondition(
        paste0("`features` must be a subset of \"centrality\", \"honem\", ",
               "or a data.frame; got: ",
               paste(unknown, collapse = ", "), "."),
        class = "hypernets_bad_input", call = NULL))
    }
    if (missing(sequences)) {
      stop(errorCondition(
        paste0("`sequences` is required: a net_hon does not carry the actor ",
               "identity of the trajectories it was built from, so the ",
               "per-actor features cannot be derived from `x` alone."),
        class = "hypernets_bad_input", call = NULL))
    }
    allowed_types <- eval(formals(hon_outcome)$centrality_type)
    unknown_types <- setdiff(centrality_type, allowed_types)
    if (length(unknown_types) > 0L) {
      stop(errorCondition(
        paste0("unknown `centrality_type`: ",
               paste(unknown_types, collapse = ", "), ". Choose from: ",
               paste(allowed_types, collapse = ", "), "."),
        class = "hypernets_bad_input", call = NULL))
    }
    centrality_type <- match.arg(centrality_type, several.ok = TRUE)
    feature_spec <- paste(features, collapse = " + ")
    seqs <- .hoo_actor_sequences(sequences, by = by, action = action,
                                 time = time)
    if (length(seqs) == 0L) {
      stop(errorCondition("`sequences` holds no actors.",
                          class = "hypernets_bad_input", call = NULL))
    }
    if (anyDuplicated(names(seqs)) > 0L) {
      # several sequences per actor are pooled into one visit profile
      seqs <- lapply(split(seqs, names(seqs)), function(g)
        unlist(g, use.names = FALSE))
    }
    node_values <- .hoo_node_values(x, features, centrality_type, honem_dim,
                                    weighted = weighted,
                                    max_paths = max_paths)
    visits <- .hoo_visit_counts(seqs, rownames(x$matrix),
                                x$max_order_observed)
    total <- rowSums(visits$counts)
    feature_matrix <- matrix(NA_real_, nrow = nrow(visits$counts),
                             ncol = ncol(node_values),
                             dimnames = list(rownames(visits$counts),
                                             colnames(node_values)))
    ok <- total > 0
    if (any(ok)) {
      feature_matrix[ok, ] <- (visits$counts[ok, , drop = FALSE] %*%
                                 node_values) / total[ok]
    }
    feature_names <- colnames(feature_matrix)
    actors <- rownames(feature_matrix)
    support <- data.frame(
      key = actors,
      n_events = as.integer(visits$n_events),
      n_unmatched = as.integer(visits$n_unmatched),
      stringsAsFactors = FALSE)
    if (sum(visits$n_unmatched) > 0L) {
      warning(warningCondition(
        sprintf(paste0("%d of %d events matched no node of the network and ",
                       "were excluded from the feature means (see ",
                       "as.data.frame(x, what = \"features\")$n_unmatched)."),
                sum(visits$n_unmatched), sum(visits$n_events)),
        class = "hypernets_dropped_actors"))
    }
  }

  # ---- align the outcome, and record what is dropped ----------------------
  y_all <- .hoo_outcome_vector(outcome, by = by, outcome_col = outcome_col)
  reason <- rep(NA_character_, length(actors))
  reason[!(actors %in% names(y_all))] <- "no outcome value"
  y <- unname(y_all[actors])
  reason[is.na(reason) & is.na(y)] <- "outcome is NA"
  bad_feature <- !stats::complete.cases(feature_matrix)
  reason[is.na(reason) & bad_feature] <- "no feature value (no matched visit)"
  extra <- setdiff(names(y_all), actors)

  dropped <- rbind(
    data.frame(key = actors[!is.na(reason)],
               reason = reason[!is.na(reason)], stringsAsFactors = FALSE),
    data.frame(key = extra,
               reason = rep("no sequences", length(extra)),
               stringsAsFactors = FALSE))
  names(dropped)[1L] <- by
  rownames(dropped) <- NULL
  if (nrow(dropped) > 0L) {
    warning(warningCondition(
      sprintf(paste0("%d actor(s) dropped from the fit (see ",
                     "as.data.frame(x, what = \"dropped\")): %s."),
              nrow(dropped),
              paste(utils::head(unique(dropped$reason), 3L),
                    collapse = "; ")),
      class = "hypernets_dropped_actors"))
  }

  keep <- is.na(reason)
  actors_fit <- actors[keep]
  y_fit <- y[keep]
  fmat <- feature_matrix[keep, , drop = FALSE]
  n <- length(y_fit)
  k <- ncol(fmat) + 1L
  if (n <= k) {
    stop(errorCondition(
      sprintf(paste0("Only %d actor(s) remain for %d model terms; an ",
                     "outcome model needs more actors than terms."), n, k),
      class = "hypernets_bad_input", call = NULL))
  }
  # exact equality is the contract here: a Bernoulli outcome is 0 or 1, and a
  # count is a whole number -- neither is the result of a floating-point
  # computation, so a tolerance would only let bad data through
  if (identical(family, "binomial") && !all(y_fit %in% c(0, 1))) {
    stop(errorCondition(
      "`family = \"binomial\"` needs an outcome of 0s and 1s.",
      class = "hypernets_bad_input", call = NULL))
  }
  if (identical(family, "poisson") &&
      (any(y_fit < 0) || any(abs(y_fit - round(y_fit)) > 0))) {
    stop(errorCondition(
      "`family = \"poisson\"` needs a non-negative integer outcome.",
      class = "hypernets_bad_input", call = NULL))
  }

  # ---- standardize, then guard the design --------------------------------
  centre <- rep(0, ncol(fmat))
  scale_by <- rep(1, ncol(fmat))
  sds <- apply(fmat, 2L, stats::sd)
  constant <- !(sds > sqrt(.Machine$double.eps) * max(1, max(abs(fmat))))
  if (any(constant)) {
    stop(errorCondition(
      paste0("These features are constant over the analysis sample and ",
             "cannot enter the model: ",
             paste(colnames(fmat)[constant], collapse = ", "),
             ". Drop them, or widen the sample."),
      class = "hypernets_rank_deficient", call = NULL))
  }
  if (standardize) {
    centre <- colMeans(fmat)
    scale_by <- sds
    fmat <- sweep(sweep(fmat, 2L, centre, "-"), 2L, scale_by, "/")
  }
  X <- cbind(`(Intercept)` = 1, fmat)
  colnames(X) <- c("(Intercept)", feature_names)
  qr_rank <- qr(X)$rank
  if (qr_rank < ncol(X)) {
    stop(errorCondition(
      sprintf(paste0("The design matrix is rank deficient (rank %d of %d ",
                     "columns): the features are exactly collinear. No term ",
                     "is dropped silently -- remove a feature, or supply ",
                     "your own `features` data.frame."),
              qr_rank, ncol(X)),
      class = "hypernets_rank_deficient", call = NULL))
  }
  # collinearity AMONG THE FEATURES: centred, then scaled to unit length, so
  # the diagnostic is invariant to `standardize` and to the units a feature
  # happens to carry (Belsley, Kuh & Welsch 1980).
  centred <- sweep(fmat, 2L, colMeans(fmat), "-")
  norms <- sqrt(colSums(centred^2))
  sv <- svd(sweep(centred, 2L, norms, "/"))$d
  condition_index <- max(sv) / min(sv)
  if (condition_index > 30) {
    warning(warningCondition(
      sprintf(paste0("Scaled condition index %.1f (> 30): the features are ",
                     "near-collinear, so individual estimates are poorly ",
                     "identified even though the fit succeeds."),
              condition_index),
      class = "hypernets_rank_deficient"))
  }

  # ---- cluster labels ----------------------------------------------------
  cluster <- NULL
  if (!is.null(nested_in)) {
    cluster <- .hoo_cluster_vector(nested_in, outcome, by, actors_fit)
    n_cl <- length(unique(cluster))
    if (n_cl < 2L) {
      stop(errorCondition(
        paste0("`nested_in` gives a single cluster; a cluster-robust ",
               "covariance needs at least two."),
        class = "hypernets_bad_input", call = NULL))
    }
    if (n_cl <= ncol(X)) {
      stop(errorCondition(
        sprintf(paste0("`nested_in` gives %d clusters for %d model terms; ",
                       "the cluster-robust meat is singular at G <= k."),
                n_cl, ncol(X)),
        class = "hypernets_bad_input", call = NULL))
    }
    if (n_cl < min_clusters) {
      warning(warningCondition(
        sprintf(paste0("Only %d clusters: cluster-robust standard errors are ",
                       "downward biased with few clusters and the t(%d) ",
                       "intervals may under-cover (Cameron & Miller 2015)."),
                n_cl, n_cl - 1L),
        class = "hypernets_few_clusters"))
    }
  }

  fit <- .hoo_fit(X, y_fit, family = family, cluster = cluster,
                  vcov = vcov)
  if (!fit$converged) {
    warning(warningCondition(
      paste0("The IRLS fit did not converge in 100 iterations; the ",
             "estimates below are not a converged maximum likelihood fit."),
      class = "hypernets_no_converge"))
  }
  if (identical(family, "binomial")) {
    if (any(fit$fitted < 1e-8 | fit$fitted > 1 - 1e-8)) {
      warning(warningCondition(
        paste0("Some fitted probabilities are numerically 0 or 1 (complete ",
               "or quasi-complete separation); the standard errors below are ",
               "not trustworthy."),
        class = "hypernets_degenerate_fit"))
    }
  }

  terms_tab <- data.frame(
    feature = names(fit$coefficients),
    estimate = unname(fit$coefficients),
    std_error = unname(fit$std_error),
    conf_low = unname(fit$conf_low),
    conf_high = unname(fit$conf_high),
    statistic = unname(fit$statistic),
    p = unname(fit$p),
    p_adj = NA_real_,
    n = rep(n, length(fit$coefficients)),
    stringsAsFactors = FALSE)
  is_feature <- terms_tab$feature != "(Intercept)"
  terms_tab$p_adj[is_feature] <-
    stats::p.adjust(terms_tab$p[is_feature], method = p_adjust)
  rownames(terms_tab) <- NULL

  features_tab <- data.frame(key = actors_fit, stringsAsFactors = FALSE)
  names(features_tab) <- by
  features_tab <- cbind(features_tab,
                        as.data.frame(feature_matrix[keep, , drop = FALSE],
                                      stringsAsFactors = FALSE))
  features_tab$outcome <- y_fit
  if (!is.null(cluster)) features_tab$cluster <- cluster
  sup <- support[match(actors_fit, support$key), , drop = FALSE]
  features_tab$n_events <- sup$n_events
  features_tab$n_unmatched <- sup$n_unmatched
  rownames(features_tab) <- NULL

  structure(list(
    terms = terms_tab,
    features = features_tab,
    dropped = dropped,
    feature_names = feature_names,
    key = by,
    n = n,
    n_terms = ncol(X),
    family = family,
    link = switch(family, gaussian = "identity", binomial = "logit",
                  poisson = "log"),
    standardize = standardize,
    centre = stats::setNames(as.numeric(centre), feature_names),
    scale = stats::setNames(as.numeric(scale_by), feature_names),
    feature_sd = stats::setNames(as.numeric(sds), feature_names),
    crit = if (is.finite(fit$df)) stats::qt(0.975, df = fit$df) else
      stats::qnorm(0.975),
    outcome_label = outcome_label,
    cluster_label = cluster_label,
    p_adjust = p_adjust,
    vcov_type = fit$vcov_type,
    vcov = fit$vcov,
    small_sample = switch(fit$vcov_type,
      CR1S = "CR1S: G/(G-1) * (n-1)/(n-k)",
      CR3  = "CR3: cluster jackknife, (I - H_gg)^-1 per cluster",
      "none"),
    nested_in = if (is.null(nested_in)) NULL else by,
    cluster_n = fit$n_clusters,
    df = fit$df,
    distribution = fit$distribution,
    condition_index = condition_index,
    r_squared = fit$r_squared,
    adj_r_squared = fit$adj_r_squared,
    sigma = fit$sigma,
    deviance = fit$deviance,
    null_deviance = fit$null_deviance,
    aic = fit$aic,
    converged = fit$converged,
    feature_spec = feature_spec,
    centrality_type = if (is.data.frame(features)) NULL else centrality_type,
    call = match.call()
  ), class = "net_outcome")
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print method for net_outcome
#'
#' @param x A \code{net_outcome} from \code{\link{hon_outcome}()}.
#' @param digits Number of significant digits in the printed table.
#' @param ... Ignored.
#' @return \code{x}, invisibly.
#' @inherit hon_outcome examples
#' @export
print.net_outcome <- function(x, digits = 3L, ...) {
  cat(sprintf("Higher-order structure -> outcome  [%s, %s covariance]\n",
              x$family,
              if (identical(x$vcov_type, "model")) "model-based" else
                "cluster-robust"))
  cat(sprintf("  %d actors / %d terms / features: %s\n",
              x$n, x$n_terms, x$feature_spec))
  if (!identical(x$vcov_type, "model")) {
    cat(sprintf("  %d clusters, %s, %s(%g) intervals\n", x$cluster_n,
                x$small_sample, x$distribution, x$df))
  }
  tab <- x$terms
  num <- setdiff(names(tab), c("feature", "n"))
  tab[num] <- lapply(tab[num], function(v) signif(v, digits))
  print(tab, row.names = FALSE)
  cat(sprintf("\n  p_adj: %s across %d features\n", x$p_adjust,
              length(x$feature_names)))
  invisible(x)
}

#' Summary method for net_outcome
#'
#' @description
#' Prints the model, what the covariance and the multiplicity correction
#' were, the collinearity diagnostic, and the actors that were dropped.
#'
#' @param object A \code{net_outcome} from \code{\link{hon_outcome}()}.
#' @param ... Ignored.
#' @return The tidy coefficient table (one row per feature), invisibly.
#' @inherit hon_outcome examples
#' @export
summary.net_outcome <- function(object, ...) {
  cat("Higher-order outcome model\n")
  cat(sprintf("  family              : %s (link: %s)\n", object$family,
              object$link))
  cat(sprintf("  actors in the fit   : %d (%d dropped)\n", object$n,
              nrow(object$dropped)))
  cat(sprintf("  features            : %s [%s]\n", object$feature_spec,
              paste(object$feature_names, collapse = ", ")))
  cat(sprintf("  predictors          : %s\n",
              if (object$standardize)
                "z-scored (estimate = change per 1 SD)" else "raw scale"))
  cat(sprintf("  covariance          : %s\n",
              if (identical(object$vcov_type, "model")) "model-based"
              else sprintf("cluster-robust (%s) over %d clusters",
                           object$vcov_type, object$cluster_n)))
  cat(sprintf("  small-sample factor : %s\n", object$small_sample))
  cat(sprintf("  intervals           : 95%%, %s(%g)\n", object$distribution,
              object$df))
  cat(sprintf("  multiplicity        : %s across %d features\n",
              object$p_adjust, length(object$feature_names)))
  cat(sprintf("  scaled condition    : %.1f%s\n", object$condition_index,
              if (object$condition_index > 30) "  (> 30: near-collinear)"
              else ""))
  if (identical(object$family, "gaussian")) {
    cat(sprintf("  R-squared           : %.4f (adjusted %.4f)\n",
                object$r_squared, object$adj_r_squared))
    cat(sprintf("  residual SD         : %.4f\n", object$sigma))
  } else {
    cat(sprintf("  deviance            : %.2f of null %.2f (AIC %.2f)\n",
                object$deviance, object$null_deviance, object$aic))
    cat(sprintf("  converged           : %s\n", object$converged))
  }
  cat("\n")
  out <- as.data.frame(object)
  print(out, row.names = FALSE)
  cat(paste0("\n  The network is treated as fixed: these intervals condition",
             " on the\n  estimated higher-order structure and do not",
             " propagate its uncertainty.\n"))
  invisible(out)
}

#' Coerce a net_outcome to a tidy table
#'
#' @param x A \code{net_outcome} from \code{\link{hon_outcome}()}.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Ignored.
#' @param what \code{"coefficients"} (default) for one row per feature,
#'   \code{"terms"} for the same table with the intercept,
#'   \code{"contrasts"} for the predicted outcome at the 25th vs the 75th
#'   percentile of each feature, \code{"features"} for the per-actor feature
#'   table on its raw scale, or \code{"dropped"} for the actors excluded from
#'   the fit and why.
#' @param significant Logical. For \code{"coefficients"} and \code{"terms"},
#'   keep only the rows with \code{p_adj < alpha}. Default \code{FALSE}.
#' @param alpha Threshold used by \code{significant}. Default 0.05.
#' @param sort_by \code{NULL} (default, model order), or one of
#'   \code{"estimate"} (largest absolute estimate first), \code{"p"} or
#'   \code{"p_adj"} (smallest first).
#' @param top Integer or \code{NULL}. Return only the first \code{top} rows,
#'   applied after every filter and after \code{sort_by}.
#' @return A base data.frame.
#'
#'   For \code{"coefficients"} and \code{"terms"}, one row per model term
#'   with \code{feature}, \code{estimate}, \code{std_error}, \code{conf_low},
#'   \code{conf_high}, \code{statistic}, \code{p}, \code{p_adj} (the
#'   multiplicity-corrected p-value; \code{NA} for the intercept, which is not
#'   one of the tested features) and \code{n} (actors in the fit) -- all on
#'   the scale the feature entered the model -- followed by the
#'   interpretable, per-SD columns: \code{feature_sd} (the standard deviation
#'   of the raw per-actor feature across the analysed actors),
#'   \code{effect_per_sd} (the effect of a one-SD increase in that feature,
#'   the others held fixed), \code{per_sd_low} / \code{per_sd_high} (its 95%
#'   interval) and \code{scale}, which names the unit: \code{"outcome units"}
#'   for the gaussian family (the coefficient times \code{feature_sd}),
#'   \code{"odds ratio"} for binomial and \code{"rate ratio"} for Poisson
#'   (the exponentiated coefficient times \code{feature_sd}, interval ends
#'   exponentiated too). The per-SD columns are \code{NA} on the intercept
#'   row.
#'
#'   For \code{"contrasts"}, one row per feature: \code{feature},
#'   \code{low_value} / \code{high_value} (the 25th and 75th percentiles of
#'   the raw per-actor feature), \code{predicted_low} / \code{predicted_high}
#'   (the predicted outcome at those two values with every other feature at
#'   its sample mean), each with a 95% interval
#'   (\code{predicted_low_ci_low}, \code{predicted_low_ci_high},
#'   \code{predicted_high_ci_low}, \code{predicted_high_ci_high}; formed on
#'   the link scale and mapped through the inverse link), \code{difference}
#'   (\code{predicted_high - predicted_low}, on the outcome scale),
#'   \code{std_error} (its delta-method standard error from the stored,
#'   possibly cluster-robust, covariance), \code{conf_low} /
#'   \code{conf_high} (its 95% Wald interval), \code{ratio} /
#'   \code{ratio_low} / \code{ratio_high} (binomial and Poisson only,
#'   \code{NA} for gaussian: the odds or rate ratio of the 75th- against the
#'   25th-percentile actor, \eqn{\exp\{\beta_j (q_{75} - q_{25}) / s_j\}}, with
#'   its interval exponentiated from the link scale) and \code{scale}
#'   (\code{"outcome units"}, \code{"probability"} or
#'   \code{"expected count"}). A difference of probabilities or counts is a
#'   Wald interval on the outcome scale and, when the uncertainty is large,
#'   can extend past the possible range; the ratio interval cannot.
#'
#'   For \code{"features"}, one row per actor with the key column, one
#'   column per feature on its RAW (un-standardized) scale, \code{outcome},
#'   the \code{cluster} when one was given, and the \code{n_events} /
#'   \code{n_unmatched} support counts. For \code{"dropped"}, one row per
#'   excluded actor with the key column and \code{reason}.
#' @inherit hon_outcome examples
#' @export
as.data.frame.net_outcome <- function(x, row.names = NULL, optional = FALSE,
                                      ..., what = c("coefficients", "terms",
                                                    "contrasts", "features",
                                                    "dropped"),
                                      significant = FALSE, alpha = 0.05,
                                      sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  stopifnot(
    "`significant` must be TRUE or FALSE" =
      is.logical(significant) && length(significant) == 1L &&
      !is.na(significant),
    "`alpha` must be a single number strictly between 0 and 1" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
      alpha > 0 && alpha < 1
  )
  if (what %in% c("features", "dropped", "contrasts")) {
    if (!is.null(sort_by) || isTRUE(significant)) {
      stop(errorCondition(
        sprintf("`sort_by` and `significant` do not apply to what = \"%s\".",
                what),
        class = "hypernets_bad_input", call = NULL))
    }
    out <- switch(what, features = x$features, dropped = x$dropped,
                  contrasts = .hoo_contrasts(x))
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  out <- .hoo_add_per_sd(x, x$terms)
  if (identical(what, "coefficients")) {
    out <- out[out$feature != "(Intercept)", , drop = FALSE]
  }
  if (isTRUE(significant)) {
    out <- out[!is.na(out$p_adj) & out$p_adj < alpha, , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("estimate", "p", "p_adj"))
    ord <- switch(sort_by,
                  estimate = order(-abs(out$estimate), out$feature),
                  p = order(out$p, out$feature),
                  p_adj = order(out$p_adj, out$feature))
    out <- out[ord, , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Plot method for net_outcome
#'
#' @description
#' Three views of the same fit. Every feature is a per-actor average: the
#' mean centrality (or HONEM coordinate) of the higher-order states an actor
#' visited, weighted by how often they visited each. The comparison is
#' therefore ACROSS actors -- actors whose trajectories run through more
#' central states against actors whose trajectories do not.
#'
#' \describe{
#'   \item{\code{"per_sd"} (default)}{The effect of a one-standard-deviation
#'     increase in each feature, the other features held fixed, on a common
#'     axis: outcome units for a gaussian fit (reference line at 0), odds or
#'     rate ratios on a log axis for binomial and Poisson fits (reference
#'     line at 1). The SD is printed in the feature label.}
#'   \item{\code{"contrasts"}}{A dumbbell per feature: the predicted outcome
#'     for an actor at the 25th percentile of that feature against one at the
#'     75th, the other features at their means, each with its 95% interval.}
#'   \item{\code{"raw"}}{The coefficient plot on the scale the features
#'     entered the model (\code{standardize}).}
#' }
#'
#' Features whose BH-corrected p-value clears \code{alpha} are marked by both
#' shape and fill, never colour alone; the subtitle gives the number of
#' actors, the clusters and covariance, and the fit (R-squared or deviance
#' explained); the caption says how to read the plot. Okabe-Ito palette.
#'
#' @param x A \code{net_outcome} from \code{\link{hon_outcome}()}.
#' @param what \code{"per_sd"} (default), \code{"contrasts"} or
#'   \code{"raw"}; see Description.
#' @param alpha Threshold on \code{p_adj} used for the marking. Default 0.05.
#' @param sort_by \code{"estimate"} (default) orders the features by the size
#'   of their effect; \code{NULL} keeps model order.
#' @param ... Ignored.
#' @return A ggplot object. An unknown \code{what} raises
#'   \code{hypernets_bad_input}.
#' @inherit hon_outcome examples
#' @export
plot.net_outcome <- function(x, what = c("per_sd", "contrasts", "raw"),
                             alpha = 0.05, sort_by = "estimate", ...) {
  stopifnot(
    "`alpha` must be a single number strictly between 0 and 1" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
      alpha > 0 && alpha < 1
  )
  allowed <- eval(formals(plot.net_outcome)$what)
  if (!is.character(what) || length(what) < 1L ||
      !all(what %in% allowed)) {
    stop(errorCondition(
      sprintf("`what` must be one of %s.",
              paste0("\"", allowed, "\"", collapse = ", ")),
      class = "hypernets_bad_input", call = NULL))
  }
  what <- what[[1L]]
  if (!is.null(sort_by)) sort_by <- match.arg(sort_by, "estimate")
  switch(what,
         per_sd = .hoo_plot_per_sd(x, alpha, sort_by),
         contrasts = .hoo_plot_contrasts(x, sort_by),
         raw = .hoo_plot_raw(x, alpha, sort_by))
}

#' Format numbers one at a time, so a vector is not padded to a common width
#'
#' @param v Numeric vector.
#' @param digits Significant digits.
#' @return Character vector.
#' @noRd
.hoo_num <- function(v, digits) {
  vapply(v, function(z) format(signif(z, digits)), character(1L))
}

#' Significance marking shared by the coefficient-style plots
#'
#' @param p_adj Numeric vector of corrected p-values.
#' @param alpha Threshold.
#' @return A factor with two fixed levels.
#' @noRd
.hoo_mark <- function(p_adj, alpha) {
  lev <- c(sprintf("BH-corrected p < %g", alpha), "not significant")
  factor(ifelse(!is.na(p_adj) & p_adj < alpha, lev[1L], lev[2L]),
         levels = lev)
}

#' Filled circle for significant, open triangle for not: shape AND fill
#'
#' @param lev The two marking levels.
#' @return A list of ggplot2 scales.
#' @noRd
.hoo_mark_scales <- function(lev) {
  list(
    ggplot2::scale_shape_manual(values = stats::setNames(c(21L, 24L), lev),
                                drop = FALSE),
    ggplot2::scale_fill_manual(values = stats::setNames(
      c("#0072B2", "white"), lev), drop = FALSE))
}

#' Per-SD effect plot
#'
#' @param x A `net_outcome`.
#' @return A ggplot object.
#' @noRd
.hoo_plot_per_sd <- function(x, alpha, sort_by) {
  d <- as.data.frame(x)
  ratio <- !identical(x$family, "gaussian")
  size <- if (ratio) abs(log(d$effect_per_sd)) else abs(d$effect_per_sd)
  if (!is.null(sort_by)) d <- d[order(-size, d$feature), , drop = FALSE]
  d$marked <- .hoo_mark(d$p_adj, alpha)
  d$label <- sprintf("%s\n(1 SD = %s)", .hoo_feature_label(d$feature),
                     .hoo_num(d$feature_sd, 2L))
  d$label <- factor(d$label, levels = rev(d$label))
  lev <- levels(d$marked)
  x_lab <- switch(x$family,
    gaussian = sprintf("Change in %s per 1-SD increase", x$outcome_label),
    binomial = sprintf("Odds ratio of %s per 1-SD increase (log scale)",
                       x$outcome_label),
    poisson = sprintf("Rate ratio of %s per 1-SD increase (log scale)",
                      x$outcome_label))
  caption <- sprintf(paste0(
    "Each actor's feature is the average value of the higher-order states ",
    "they visited.\nEach point compares two actors one SD apart on that ",
    "feature, the other features held fixed (%s).\nBars: %s, ",
    "unadjusted; the marker shows the BH-corrected test."),
    if (ratio) "1 = no difference" else "0 = no difference",
    .hoo_interval_words(x))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$effect_per_sd,
                                       y = .data$label)) +
    ggplot2::geom_vline(xintercept = if (ratio) 1 else 0,
                        linetype = "dashed", colour = "#999999") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$per_sd_low, xmax = .data$per_sd_high),
      colour = "#000000", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$marked,
                                     fill = .data$marked),
                        colour = "#000000", size = 3, stroke = 0.7,
                        show.legend = TRUE) +
    .hoo_mark_scales(lev) +
    ggplot2::labs(x = x_lab, y = NULL, shape = NULL, fill = NULL,
                  title = sprintf("How %s differs with higher-order exposure",
                                  x$outcome_label),
                  subtitle = .hoo_model_context(x),
                  caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.title.position = "plot",
                   plot.caption.position = "plot")
  if (ratio) p <- p + ggplot2::scale_x_log10()
  p
}

#' Percentile-contrast dumbbell plot
#'
#' @param x A `net_outcome`.
#' @return A ggplot object.
#' @noRd
.hoo_plot_contrasts <- function(x, sort_by) {
  d <- as.data.frame(x, what = "contrasts")
  if (!is.null(sort_by)) {
    d <- d[order(-abs(d$difference), d$feature), , drop = FALSE]
  }
  labels <- sprintf("%s\n(%s vs %s)", .hoo_feature_label(d$feature),
                    .hoo_num(d$low_value, 3L), .hoo_num(d$high_value, 3L))
  pos <- rev(seq_len(nrow(d)))
  lev <- c("25th percentile", "75th percentile")
  long <- data.frame(
    y = c(pos + 0.12, pos - 0.12),
    y_mid = c(pos, pos),
    at = factor(rep(lev, each = nrow(d)), levels = lev),
    predicted = c(d$predicted_low, d$predicted_high),
    ci_low = c(d$predicted_low_ci_low, d$predicted_high_ci_low),
    ci_high = c(d$predicted_low_ci_high, d$predicted_high_ci_high))
  segs <- data.frame(y = pos, from = d$predicted_low, to = d$predicted_high)
  unit <- switch(x$family, gaussian = "", binomial = "probability of ",
                 poisson = "expected ")
  caption <- sprintf(paste0(
    "Two hypothetical actors per row: one at the 25th and one at the 75th ",
    "percentile of that feature,\nevery other feature at its mean. ",
    "Bars: %s."), .hoo_interval_words(x))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$predicted, y = .data$y)) +
    ggplot2::geom_segment(data = segs,
                          ggplot2::aes(x = .data$from, xend = .data$to,
                                       y = .data$y, yend = .data$y),
                          colour = "#999999", linewidth = 1.2,
                          inherit.aes = FALSE) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = .data$ci_low,
                                         xmax = .data$ci_high,
                                         colour = .data$at),
                            linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(colour = .data$at, shape = .data$at),
                        size = 3) +
    ggplot2::scale_colour_manual(values = stats::setNames(
      c("#E69F00", "#0072B2"), lev)) +
    ggplot2::scale_shape_manual(values = stats::setNames(c(16L, 17L), lev)) +
    ggplot2::scale_y_continuous(breaks = pos, labels = labels,
                                expand = ggplot2::expansion(add = 0.5)) +
    ggplot2::labs(
      x = sprintf("Predicted %s%s", unit, x$outcome_label), y = NULL,
      colour = "Actor at the feature's", shape = "Actor at the feature's",
      title = sprintf("Predicted %s: low vs high exposure", x$outcome_label),
      subtitle = .hoo_model_context(x), caption = caption) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid.minor.y = ggplot2::element_blank(),
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.title.position = "plot",
                   plot.caption.position = "plot")
}

#' Coefficient plot on the model scale
#'
#' @param x A `net_outcome`.
#' @return A ggplot object.
#' @noRd
.hoo_plot_raw <- function(x, alpha, sort_by) {
  d <- as.data.frame(x, sort_by = sort_by)
  d$marked <- .hoo_mark(d$p_adj, alpha)
  lev <- levels(d$marked)
  d$feature <- factor(d$feature, levels = rev(d$feature))
  ggplot2::ggplot(d, ggplot2::aes(x = .data$estimate, y = .data$feature)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "#999999") +
    ggplot2::geom_linerange(
      ggplot2::aes(xmin = .data$conf_low, xmax = .data$conf_high),
      colour = "#000000", linewidth = 0.6) +
    ggplot2::geom_point(ggplot2::aes(shape = .data$marked,
                                     fill = .data$marked),
                        colour = "#000000", size = 3, stroke = 0.7,
                        show.legend = TRUE) +
    .hoo_mark_scales(lev) +
    ggplot2::labs(
      x = if (x$standardize) "Coefficient (per 1 SD of the feature, link scale)"
      else "Coefficient (per raw unit of the feature, link scale)",
      y = NULL, shape = NULL, fill = NULL,
      title = "Model coefficients",
      subtitle = .hoo_model_context(x),
      caption = sprintf("Bars: %s.", .hoo_interval_words(x))) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   plot.caption = ggplot2::element_text(hjust = 0),
                   plot.title.position = "plot",
                   plot.caption.position = "plot")
}
