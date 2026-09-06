# ---- Markov Order Test (Exact within-context permutation) ----
#
# Principled test of Markov order for categorical sequence data.
#
# At each order k, test H0: "process is order (k-1)" against H1: "process
# is order k" by reframing as a conditional-independence test on
# (k+1)-grams (x, w, s) where:
#   - w is the (k-1)-gram context
#   - x is the additional (k-th-back) state
#   - s is the next state
# Under H0, s is independent of x given w.
#
# Test statistic: classical LRT for conditional independence
#   G^2 = 2 sum_{x,w,s} N(x,w,s) log[ N(x,w,s) / E(x,w,s) ]
# with E(x,w,s) = N(x,w,.) * N(.,w,s) / N(.,w,.),
# df = sum_w (r_w - 1)(c_w - 1).
#
# Null distribution: EXACT within-w permutation. For each context w,
# the successor labels are exchangeable under H0; shuffling s within
# each w-group gives an exact reference distribution. No plug-in MLE
# bias (which inflated the parametric-bootstrap null on sparse k+1-gram
# tables). No refitting per replicate.


# ---------------------------------------------------------------------------
# Internal: Extract (x, w, s) tuples from trajectories
# ---------------------------------------------------------------------------

#' @noRd
.mot_extract_tuples <- function(trajectories, k) {
  stopifnot(k >= 1L)
  win <- k + 1L

  tuples <- lapply(trajectories, function(traj) {
    n <- length(traj)
    if (n < win) return(NULL)
    pos <- seq_len(n - k)
    x <- traj[pos]
    s <- traj[pos + k]
    w <- if (k == 1L) {
      rep("", length(pos))
    } else {
      vapply(pos, function(p) paste(traj[(p + 1L):(p + k - 1L)],
                                    collapse = .HON_SEP),
             character(1L))
    }
    data.frame(x = x, w = w, s = s, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, tuples[!vapply(tuples, is.null, logical(1L))])
  if (is.null(out)) {
    data.frame(x = character(0L), w = character(0L), s = character(0L),
               stringsAsFactors = FALSE)
  } else out
}


# ---------------------------------------------------------------------------
# Internal: G^2 statistic and df for one w-table
# ---------------------------------------------------------------------------

#' @noRd
.mot_g2_one_w <- function(x, s) {
  tab <- table(x, s)
  if (nrow(tab) < 2L || ncol(tab) < 2L) return(c(stat = 0, df = 0))
  rs <- rowSums(tab); cs <- colSums(tab); total <- sum(tab)
  if (total == 0) return(c(stat = 0, df = 0))
  expt <- outer(rs, cs) / total
  nz <- tab > 0 & expt > 0
  g2 <- 2 * sum(tab[nz] * log(tab[nz] / expt[nz]))
  df <- (nrow(tab) - 1L) * (ncol(tab) - 1L)
  c(stat = g2, df = df)
}


#' @noRd
.mot_g2 <- function(tuples) {
  if (nrow(tuples) == 0L) return(list(stat = 0, df = 0L))
  by_w <- split(seq_len(nrow(tuples)), tuples$w)
  parts <- vapply(by_w, function(idx) {
    .mot_g2_one_w(tuples$x[idx], tuples$s[idx])
  }, numeric(2L))
  list(stat = sum(parts["stat", ]),
       df   = as.integer(sum(parts["df", ])))
}


# ---------------------------------------------------------------------------
# Internal: Within-w permutation null
# ---------------------------------------------------------------------------

#' @noRd
.mot_permute_null <- function(tuples, n_perm, parallel = FALSE, n_cores = 2L) {
  if (nrow(tuples) == 0L) return(numeric(0L))
  by_w <- split(seq_len(nrow(tuples)), tuples$w)
  s_orig <- tuples$s
  x <- tuples$x

  one_perm <- function(.) {
    s <- s_orig
    shuffled <- lapply(by_w, function(idx) {
      if (length(idx) > 1L) s[idx] <<- sample(s[idx])
      NULL
    })
    # Recompute G^2 directly on shuffled s (x, w unchanged)
    parts <- vapply(by_w, function(idx) {
      .mot_g2_one_w(x[idx], s[idx])
    }, numeric(2L))
    sum(parts["stat", ])
  }

  map_fn <- if (parallel && .Platform$OS.type != "windows") {
    function(X, FUN) parallel::mclapply(X, FUN, mc.cores = n_cores)
  } else lapply

  unlist(map_fn(seq_len(n_perm), one_perm), use.names = FALSE)
}


# ---------------------------------------------------------------------------
# Internal: Fit models + compute per-order log-likelihood (for IC panel)
# ---------------------------------------------------------------------------

#' @noRd
.mot_fit_models <- function(trajectories, max_order) {
  marginal <- .mogen_marginal(trajectories)
  trans_mats <- vector("list", max_order + 1L)
  trans_mats[[1L]] <- marginal
  layer_dofs <- integer(max_order + 1L)
  layer_dofs[1L] <- length(marginal) - 1L

  lapply(seq_len(max_order), function(k) {
    kg <- .mogen_count_kgrams(trajectories, k)
    tm <- .mogen_transition_matrix(kg$nodes, kg$edges)
    trans_mats[[k + 1L]] <<- tm
    layer_dofs[k + 1L] <<- .mogen_layer_dof(tm)
    NULL
  })

  logliks <- vapply(0L:max_order, function(k) {
    .mogen_log_likelihood(trajectories, k, trans_mats[seq_len(k + 1L)])
  }, numeric(1L))

  list(trans_mats = trans_mats, marginal = marginal,
       layer_dofs = layer_dofs, logliks = logliks)
}


# ---------------------------------------------------------------------------
# Main: markov_order_test()
# ---------------------------------------------------------------------------

#' Test the Markov order of a sequential process
#'
#' @description
#' Principled test of whether a categorical sequence is best described
#' as a \eqn{k}-th order Markov chain. At each order
#' \eqn{k = 1, \ldots,} \code{max_order}, the function computes the
#' classical likelihood-ratio statistic (\eqn{G^2}) for the conditional
#' independence \eqn{s \perp x \mid w}, where \eqn{w} is the
#' \eqn{(k-1)}-gram context, \eqn{x} is the extra (k-th-back) state
#' added at order \eqn{k}, and \eqn{s} is the next state. Under
#' \eqn{H_0} (order-\eqn{(k-1)} is correct), \eqn{s} is independent of
#' \eqn{x} given \eqn{w}.
#'
#' The null distribution is obtained by an \strong{exact within-\eqn{w}
#' permutation test}: for each context \eqn{w} the successor labels are
#' exchangeable under \eqn{H_0}, so shuffling \eqn{s} within each
#' \eqn{w}-group yields an exact reference distribution for \eqn{G^2}.
#' No plug-in MLE bias and no refitting per replicate. An asymptotic
#' \eqn{\chi^2} p-value is reported alongside for reference.
#'
#' The optimal order is the smallest \eqn{k} that is \strong{not}
#' significantly better than \eqn{k - 1} at level \code{alpha}: we
#' keep raising the order while the test rejects, and stop at the first
#' non-rejection.
#'
#' @param data A data.frame (wide format, one sequence per row) or list
#'   of character vectors (one per trajectory). NAs are treated as end
#'   of sequence.
#' @param max_order Integer. Highest Markov order to test. Default 3.
#' @param n_perm Integer. Number of within-\eqn{w} permutations per order.
#'   Default 500.
#' @param alpha Numeric. Significance level for order selection. Default 0.05.
#' @param parallel Logical. Use \code{parallel::mclapply} for permutations.
#'   Default \code{FALSE} (set \code{TRUE} only on Unix-like systems).
#' @param n_cores Integer. Cores for parallel execution. Default 2.
#' @param seed Optional integer seed for reproducibility.
#' @return An object of class \code{net_markov_order} with elements:
#' \describe{
#'   \item{optimal_order}{Integer. Selected order via sequential permutation test.}
#'   \item{bic_order}{Integer. Order minimising BIC (reported by
#'     \code{print}/\code{plot}; not used in the permutation selection).}
#'   \item{aic_order}{Integer. Order minimising AIC (reported by
#'     \code{print}/\code{plot}; not used in the permutation selection).}
#'   \item{test_table}{Tidy data.frame, one row per order tested with
#'     columns \code{order}, \code{loglik}, \code{AIC}, \code{BIC},
#'     \code{df}, \code{g2}, \code{p_permutation}, \code{p_asymptotic},
#'     \code{significant}. \code{AIC}/\code{BIC} are the
#'     information-criterion values used for the \code{ic} plot panel and
#'     to derive \code{aic_order}/\code{bic_order}; the order-0 row has
#'     \code{NA} for the test columns (\code{df}, \code{g2},
#'     \code{p_permutation}, \code{p_asymptotic}, \code{significant}).}
#'   \item{permutation_null}{List of numeric vectors (length \code{max_order}),
#'     one empirical null \eqn{G^2} distribution per order.}
#'   \item{logliks}{Named numeric vector of log-likelihoods per order
#'     (for AIC / BIC panel only, not used in the test).}
#'   \item{layer_dofs}{Named integer vector of model degrees of freedom per
#'     order (free parameters added at each layer), used to compute the
#'     AIC / BIC columns.}
#'   \item{transition_matrices}{List of fitted transition matrices.}
#'   \item{states}{Character vector of observed state labels.}
#'   \item{n_sequences, n_observations}{Data summary.}
#'   \item{n_perm, alpha, max_order}{Call settings.}
#' }
#'
#' @examples
#' \donttest{
#' # First-order Markov data: test should select order 1
#' set.seed(1)
#' states <- letters[1:4]
#' tm <- matrix(runif(16), 4, 4, dimnames = list(states, states))
#' tm <- tm / rowSums(tm)
#' seqs <- lapply(1:30, function(.) {
#'   s <- character(50); s[1] <- sample(states, 1)
#'   for (i in 2:50) s[i] <- sample(states, 1, prob = tm[s[i - 1], ])
#'   s
#' })
#' res <- markov_order_test(seqs, max_order = 3, n_perm = 300, seed = 1)
#' as.data.frame(res)
#' summary(res)
#' plot(res)
#' }
#' @export
markov_order_test <- function(data, max_order = 3L, n_perm = 500L, alpha = 0.05,
                               parallel = FALSE, n_cores = 2L, seed = NULL) {
  if (inherits(data, "netobject_group")) {
    out <- lapply(data, function(net) {
      seq_data <- net$data
      if (is.null(seq_data)) {
        stop("netobject_group member has no $data; rebuild with build_network().",
             call. = FALSE)
      }
      markov_order_test(seq_data, max_order = max_order, n_perm = n_perm,
                        alpha = alpha, parallel = parallel,
                        n_cores = n_cores, seed = seed)
    })
    class(out) <- c("net_markov_order_group", "list")
    return(out)
  }
  if (inherits(data, "netobject")) {
    if (is.null(data$data)) {
      stop("netobject has no $data; rebuild with build_network().",
           call. = FALSE)
    }
    data <- data$data
  }
  max_order <- as.integer(max_order)
  n_perm    <- as.integer(n_perm)
  stopifnot(
    "'max_order' must be >= 1" = max_order >= 1L,
    "'n_perm' must be >= 1"    = n_perm >= 1L,
    "'alpha' must be in (0, 1)" = alpha > 0 && alpha < 1,
    is.logical(parallel)
  )
  if (!is.null(seed)) set.seed(as.integer(seed))

  # Accept a prepare() result directly: use its $sequence_data (wide df)
  if (is.list(data) && !is.data.frame(data) && "sequence_data" %in% names(data)) {
    data <- data$sequence_data
  }

  trajectories <- .hon_parse_input(data, collapse_repeats = FALSE)
  n_seqs <- length(trajectories)
  stopifnot("need >= 1 sequence after parsing" = n_seqs >= 1L)
  seq_lens <- vapply(trajectories, length, integer(1L))

  longest <- max(seq_lens)
  if (max_order >= longest) {
    max_order <- longest - 1L
    message(sprintf("max_order capped at %d (longest sequence = %d)",
                    max_order, longest))
  }

  fits <- .mot_fit_models(trajectories, max_order)

  # --- Per-order test ---
  per_order <- lapply(seq_len(max_order), function(k) {
    tuples <- .mot_extract_tuples(trajectories, k)
    g2 <- .mot_g2(tuples)
    if (nrow(tuples) > 0L && g2$df > 0L) {
      null_stats <- .mot_permute_null(tuples, n_perm, parallel, n_cores)
      p_perm  <- (sum(null_stats >= g2$stat) + 1) / (length(null_stats) + 1)
      p_asymp <- stats::pchisq(g2$stat, g2$df, lower.tail = FALSE)
    } else {
      null_stats <- numeric(0L)
      p_perm <- NA_real_; p_asymp <- NA_real_
    }
    list(stat = g2$stat, df = g2$df, null = null_stats,
         p_perm = p_perm, p_asymp = p_asymp, n_tuples = nrow(tuples))
  })

  # Sequential selection
  stat_vec   <- vapply(per_order, `[[`, numeric(1L), "stat")
  df_vec     <- vapply(per_order, `[[`, numeric(1L), "df")
  p_perm_vec <- vapply(per_order, `[[`, numeric(1L), "p_perm")
  p_asy_vec  <- vapply(per_order, `[[`, numeric(1L), "p_asymp")

  sig_vec <- !is.na(p_perm_vec) & p_perm_vec < alpha
  optimal_order <- 0L
  cont <- TRUE
  lapply(seq_len(max_order), function(k) {
    if (cont && sig_vec[k]) optimal_order <<- k
    else cont <<- FALSE
    NULL
  })

  orders_all <- 0L:max_order
  cum_dof    <- cumsum(fits$layer_dofs)
  aic_vec    <- 2 * cum_dof - 2 * fits$logliks
  bic_vec    <- log(sum(seq_lens)) * cum_dof - 2 * fits$logliks

  test_table <- data.frame(
    order         = orders_all,
    loglik        = as.numeric(fits$logliks),
    AIC           = as.numeric(aic_vec),
    BIC           = as.numeric(bic_vec),
    df            = c(NA_integer_, as.integer(df_vec)),
    g2            = c(NA_real_,    as.numeric(stat_vec)),
    p_permutation = c(NA_real_,    as.numeric(p_perm_vec)),
    p_asymptotic  = c(NA_real_,    as.numeric(p_asy_vec)),
    significant   = c(NA,          as.logical(sig_vec)),
    stringsAsFactors = FALSE,
    row.names     = NULL
  )

  null_list <- lapply(per_order, `[[`, "null")
  names(null_list) <- paste0("order_", seq_len(max_order))
  names(fits$logliks) <- paste0("order_", orders_all)
  names(fits$layer_dofs) <- paste0("order_", orders_all)

  bic_order <- as.integer(test_table$order[which.min(test_table$BIC)])
  aic_order <- as.integer(test_table$order[which.min(test_table$AIC)])

  structure(
    list(
      optimal_order       = as.integer(optimal_order),
      bic_order           = bic_order,
      aic_order           = aic_order,
      test_table          = test_table,
      permutation_null    = null_list,
      logliks             = fits$logliks,
      layer_dofs          = fits$layer_dofs,
      transition_matrices = fits$trans_mats,
      states              = names(fits$marginal),
      n_sequences         = n_seqs,
      n_observations      = sum(seq_lens),
      n_perm              = n_perm,
      alpha               = alpha,
      max_order           = max_order
    ),
    class = "net_markov_order"
  )
}


# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_markov_order
#'
#' @param x A \code{net_markov_order} object.
#' @param ... Ignored.
#' @return The input object, invisibly.
#' @inherit markov_order_test examples
#' @export
print.net_markov_order <- function(x, ...) {
  cat(sprintf("Markov Order Test  [within-w permutation, n_perm = %d, alpha = %.3f]\n",
              x$n_perm, x$alpha))
  cat(sprintf("  %d sequences / %d observations / %d states\n\n",
              x$n_sequences, x$n_observations, length(x$states)))
  cat(sprintf("  Selected order  BIC: %d   AIC: %d   permutation-LRT: %d\n\n",
              x$bic_order, x$aic_order, x$optimal_order))
  tt <- x$test_table
  tt$loglik <- round(tt$loglik, 2)
  tt$AIC    <- round(tt$AIC,    2)
  tt$BIC    <- round(tt$BIC,    2)
  tt$g2     <- round(tt$g2,     2)
  print(tt, row.names = FALSE)
  invisible(x)
}

#' Print method for `net_markov_order_group`
#'
#' @param x A `net_markov_order_group` (named list of `net_markov_order`
#'   results, one per group).
#' @param ... Forwarded to `print.net_markov_order` for each element.
#' @return `x` invisibly.
#' @export
print.net_markov_order_group <- function(x, ...) {
  cat(sprintf("Markov Order Test -- %d groups: %s\n\n",
              length(x), paste(names(x), collapse = ", ")))
  for (nm in names(x)) {
    cat(sprintf("--- %s ---\n", nm))
    print(x[[nm]], ...)
    cat("\n")
  }
  invisible(x)
}


#' Summary Method for net_markov_order
#'
#' @param object A \code{net_markov_order} object.
#' @param ... Ignored.
#' @return The tidy \code{test_table} data.frame, with the selected
#'   \code{optimal_order} attached as an attribute.
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @inherit markov_order_test examples
#' @export
summary.net_markov_order <- function(object, ...) {
  out <- object$test_table
  attr(out, "optimal_order") <- object$optimal_order
  attr(out, "bic_order")     <- object$bic_order
  attr(out, "aic_order")     <- object$aic_order
  attr(out, "alpha")         <- object$alpha
  attr(out, "n_perm")        <- object$n_perm
  invisible(out)
}


#' Plot Method for net_markov_order
#'
#' @description
#' Two-panel professional visualization:
#' \itemize{
#'   \item Panel A: log-likelihood, AIC, BIC across tested orders with
#'     the selected order highlighted (both the permutation-selected
#'     order and the BIC-minimizing order are marked).
#'   \item Panel B: permutation null density per order with the observed
#'     \eqn{G^2} as a vertical marker; colored by rejection at
#'     \code{alpha}.
#' }
#' Uses the Okabe-Ito colorblind-safe palette.
#'
#' @param x A \code{net_markov_order} object.
#' @param panel Which panel(s) to render: \code{"both"}, \code{"ic"},
#'   or \code{"permutation"}. Default \code{"both"}.
#' @param combined When \code{panel = "both"} and \code{combined = TRUE}
#'   (default), the two panels are drawn side-by-side. If \pkg{gridExtra}
#'   is installed they are arranged into a single drawable/saveable
#'   gtable (returned); otherwise base \code{grid} viewports draw both
#'   panels and a named list of the two ggplots is returned invisibly.
#'   When \code{FALSE}, returns that named list (\code{ic},
#'   \code{permutation}) without drawing. Ignored when
#'   \code{panel != "both"}.
#' @param ... Ignored.
#' @return A ggplot (single panel); for \code{panel = "both"}, either a
#'   \code{gridExtra} gtable (when \pkg{gridExtra} is installed) or a
#'   named list of two ggplots (\code{ic}, \code{permutation}) drawn
#'   side-by-side and returned invisibly.
#' @inherit markov_order_test examples
#' @export
plot.net_markov_order <- function(x,
                                   panel = c("both", "ic", "permutation"),
                                   combined = TRUE,
                                   ...) {
  panel <- match.arg(panel)
  stopifnot(is.logical(combined), length(combined) == 1L)
  pal <- c(reject = "#D55E00", accept = "#0072B2",
           ic_ll = "#009E73", ic_aic = "#E69F00", ic_bic = "#56B4E9",
           selected = "#D55E00")

  tt <- x$test_table
  ic_df <- rbind(
    data.frame(order = tt$order, metric = "log-likelihood", value = tt$loglik),
    data.frame(order = tt$order, metric = "AIC",            value = tt$AIC),
    data.frame(order = tt$order, metric = "BIC",            value = tt$BIC)
  )
  ic_df$metric <- factor(ic_df$metric,
                         levels = c("log-likelihood", "AIC", "BIC"))

  p_ic <- ggplot2::ggplot(ic_df,
      ggplot2::aes(x = .data$order, y = .data$value,
                   color = .data$metric)) +
    ggplot2::geom_vline(xintercept = x$optimal_order,
                        linetype = "dashed", color = pal["selected"],
                        linewidth = 0.6) +
    ggplot2::geom_vline(xintercept = x$bic_order,
                        linetype = "dotted", color = pal["ic_bic"],
                        linewidth = 0.6) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(size = 2.4) +
    ggplot2::scale_color_manual(values = unname(pal[c("ic_ll",
                                                       "ic_aic",
                                                       "ic_bic")]),
                                name = NULL) +
    ggplot2::scale_x_continuous(breaks = tt$order) +
    ggplot2::facet_wrap(~ .data$metric, scales = "free_y", ncol = 1) +
    ggplot2::labs(x = "Markov order",
                  y = NULL,
                  title = "Model fit across orders",
                  subtitle = sprintf(
                    "Permutation LRT selects order %d (alpha = %.2f); BIC minimized at %d",
                    x$optimal_order, x$alpha, x$bic_order)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  # --- Permutation null ---
  null_df <- do.call(rbind, lapply(seq_along(x$permutation_null), function(k) {
    nn <- x$permutation_null[[k]]
    if (length(nn) == 0L) return(NULL)
    data.frame(order = k, stat = nn, stringsAsFactors = FALSE)
  }))
  obs_df <- data.frame(
    order = seq_len(x$max_order),
    g2    = tt$g2[-1L],
    p     = tt$p_permutation[-1L],
    sig   = tt$significant[-1L]
  )
  obs_df$decision <- ifelse(obs_df$sig, "reject H0", "accept H0")
  if (!is.null(null_df)) {
    null_df$order_fac <- factor(sprintf("order %d vs %d", null_df$order,
                                         null_df$order - 1L),
                                 levels = sprintf("order %d vs %d",
                                                   seq_len(x$max_order),
                                                   seq_len(x$max_order) - 1L))
  }
  obs_df$order_fac <- factor(sprintf("order %d vs %d", obs_df$order,
                                      obs_df$order - 1L),
                              levels = sprintf("order %d vs %d",
                                                seq_len(x$max_order),
                                                seq_len(x$max_order) - 1L))

  p_perm <- ggplot2::ggplot(null_df,
      ggplot2::aes(x = .data$stat)) +
    ggplot2::geom_density(fill = "#CCCCCC", color = "#888888",
                          alpha = 0.6) +
    ggplot2::geom_vline(data = obs_df,
                        ggplot2::aes(xintercept = .data$g2,
                                     color = .data$decision),
                        linewidth = 0.9) +
    ggplot2::geom_text(data = obs_df,
                       ggplot2::aes(x = .data$g2, y = 0,
                                    label = sprintf("p = %.3f", .data$p),
                                    color = .data$decision),
                       hjust = -0.1, vjust = -0.8, size = 3.2,
                       show.legend = FALSE) +
    ggplot2::scale_color_manual(values = c("reject H0" = unname(pal["reject"]),
                                            "accept H0" = unname(pal["accept"])),
                                name = NULL) +
    ggplot2::facet_wrap(~ .data$order_fac, scales = "free", ncol = 1) +
    ggplot2::labs(x = expression(G^2 ~ "statistic"),
                  y = "density",
                  title = "Permutation null distributions",
                  subtitle = sprintf(
                    "Within-context permutation, n_perm = %d; observed G^2 as vertical line",
                    x$n_perm)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   strip.text = ggplot2::element_text(face = "bold"),
                   plot.title = ggplot2::element_text(face = "bold"),
                   panel.grid.minor = ggplot2::element_blank())

  if (panel == "ic") return(p_ic)
  if (panel == "permutation") return(p_perm)

  if (!combined) return(invisible(list(ic = p_ic, permutation = p_perm)))

  # Side-by-side. When gridExtra is available, return its arranged gtable (a
  # drawable/saveable grob -- `ggsave(plot(fit))` keeps working). Otherwise
  # fall back to base `grid` viewports (always available): both panels still
  # render, but only a named list of the two ggplots can be returned.
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(p_ic, p_perm, ncol = 2))
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1L, 2L)))
  print(p_ic,   vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 1L))
  print(p_perm, vp = grid::viewport(layout.pos.row = 1L, layout.pos.col = 2L))
  grid::popViewport()
  invisible(list(ic = p_ic, permutation = p_perm))
}

#' Coerce a net_markov_order to a tidy table
#'
#' @param x A `net_markov_order` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"orders"` (default) for the per-order test table, or
#'   `"null"` for the permutation null distribution behind it.
#' @return A data.frame. For `what = "orders"`, one row per candidate order
#'   with the log-likelihood, information criteria, likelihood-ratio
#'   statistic and permutation p-value. For `what = "null"`, one row per
#'   permutation replicate: `order`, `replicate`, `g2` (the replicate's
#'   likelihood-ratio statistic under the null).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_markov_order <- function(x, row.names = NULL,
                                           optional = FALSE, ...,
                                           what = c("orders", "null"),
                                           top = NULL) {
  what <- match.arg(what)
  if (what == "orders") {
    out <- x$test_table
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  nulls <- x$permutation_null
  if (!length(nulls)) {
    return(data.frame(order = integer(0), replicate = integer(0),
                      g2 = numeric(0), stringsAsFactors = FALSE))
  }
  ord <- suppressWarnings(as.integer(names(nulls)))
  if (anyNA(ord)) ord <- seq_along(nulls)
  out <- do.call(rbind, Map(function(v, k) {
    data.frame(order = rep.int(k, length(v)),
               replicate = seq_along(v),
               g2 = as.numeric(v),
               stringsAsFactors = FALSE)
  }, nulls, ord))
  rownames(out) <- NULL
  .ho_top(out, top)
}
