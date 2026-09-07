# ---- Path Dependence: Per-Context Order-k vs Order-1 Diagnostic ----

#' @noRd
.pd_extract_data <- function(x) {
  if (is.matrix(x))     return(as.data.frame(x, stringsAsFactors = FALSE))
  if (is.data.frame(x)) return(x)
  if (inherits(x, "netobject") || inherits(x, "cograph_network")) {
    if (!is.null(x$data) && (is.data.frame(x$data) || is.matrix(x$data)))
      return(as.data.frame(x$data, stringsAsFactors = FALSE))
    stop(
      "netobject does not carry the source sequence data. Pass the wide ",
      "sequence data.frame directly to path_dependence().",
      call. = FALSE
    )
  }
  stop(
    "'x' must be a wide sequence data.frame, matrix, or netobject built ",
    "with raw data preserved.",
    call. = FALSE
  )
}

#' @noRd
.pd_kgram_counts <- function(data, k) {
  ## Wide-format only. Vectorised k-gram extraction across columns,
  ## results pooled across rows. NAs in any of the k positions zero
  ## out the contribution.
  m <- as.matrix(data)
  storage.mode(m) <- "character"
  ncols <- ncol(m)
  if (ncols < k) {
    stop(sprintf("Sequences too short for order %d (need >=%d columns, got %d).",
                 k - 1L, k, ncols), call. = FALSE)
  }
  windows <- vapply(seq_len(ncols - k + 1L), function(t) {
    do.call(paste, c(lapply(seq_len(k), function(j) m[, t + j - 1L]),
                     list(sep = "\x01")))
  }, character(nrow(m)))
  pooled <- as.vector(windows)
  pooled <- pooled[!grepl("NA", pooled, fixed = TRUE)]
  table(pooled)
}

#' @noRd
.pd_split_kgram <- function(kgram_names, k) {
  parts  <- strsplit(kgram_names, "\x01", fixed = TRUE)
  ctx    <- vapply(parts,
                   function(p) paste(p[seq_len(k - 1L)], collapse = " -> "),
                   character(1))
  nxt    <- vapply(parts, function(p) p[k], character(1))
  list(context = ctx, next_state = nxt)
}

#' @noRd
.pd_kl <- function(p, q, base = 2) {
  ## Returns KL(p || q) in `base` units. Defines 0 log 0 := 0; if p[i] > 0
  ## but q[i] == 0, returns Inf.
  msk <- p > 0
  if (any(p[msk] > 0 & q[msk] == 0)) return(Inf)
  sum(p[msk] * log(p[msk] / q[msk], base = base))
}

#' @noRd
.pd_entropy <- function(p, base = 2) {
  msk <- p > 0
  -sum(p[msk] * log(p[msk], base = base))
}


# ---- path_dependence() ----

#' Per-Context Path Dependence at Order k
#'
#' @description
#' Diagnoses where a chain's order-1 Markov assumption fails by comparing,
#' for each order-k context \eqn{(s_1, \ldots, s_{k-1})}, the empirical
#' next-state distribution \eqn{P(s_k \mid s_1, \ldots, s_{k-1})} against
#' the order-1 prediction \eqn{P(s_k \mid s_{k-1})} that uses only the most
#' recent state. Returns a tidy per-context table sorted by Kullback-Leibler
#' divergence so the analyst can see exactly *which* histories carry extra
#' predictive information.
#'
#' @param x A wide sequence data.frame / matrix (rows = actors, columns =
#'   time-steps), or a \code{netobject} that carries the source data.
#' @param order Integer. Order of the conditioning context. \code{order = 2}
#'   (default) compares 2-step memory against 1-step; \code{order = 3}
#'   compares 3-step memory; etc. Must be a whole number; a non-integer
#'   value is an error rather than being silently truncated.
#' @param min_count Integer. Drop contexts seen fewer than this many times.
#'   Default 5. Very small samples produce noisy KL estimates.
#' @param base Numeric. Logarithm base for entropy and KL. Default 2 (bits).
#'
#' @return An object of class \code{"net_path_dependence"} with
#' \describe{
#'   \item{contexts}{tidy data.frame, one row per order-k context, sorted
#'     by KL descending. Columns: \code{context} (e.g. "A -> B"), \code{n}
#'     (count), \code{H_order1} (entropy of \eqn{P(\cdot \mid s_{k-1})}),
#'     \code{H_orderk} (entropy of \eqn{P(\cdot \mid \mathrm{context})}),
#'     \code{H_drop} (= \code{H_order1} - \code{H_orderk}), \code{KL}
#'     (= \eqn{D_{KL}(P_k \| P_1)}), \code{top_o1} (most likely next state
#'     under order-1), \code{top_ok} (most likely next state under order-k),
#'     \code{flips} (logical: did the most likely next state change?).}
#'   \item{chain}{list with chain-level summaries: \code{KL_weighted}
#'     (count-weighted mean KL across contexts), \code{H_drop_weighted}
#'     (count-weighted mean entropy drop), \code{n_contexts}, \code{n_flips}
#'     (contexts where the most-likely next state changed).}
#'   \item{order}{integer}
#'   \item{base}{numeric}
#'   \item{min_count}{integer}
#'   \item{states}{character vector}
#' }
#'
#' @details
#' For each context \eqn{c = (s_1, \ldots, s_{k-1})} occurring at least
#' \code{min_count} times, the function computes:
#' \itemize{
#'   \item the empirical conditional \eqn{P_k(\cdot \mid c)} from k-gram counts;
#'   \item the order-1 prediction \eqn{P_1(\cdot \mid s_{k-1})} from the
#'     most recent state alone (the bigram-marginal estimator);
#'   \item the entropy drop \eqn{H(P_1) - H(P_k)} - bits of uncertainty
#'     removed by extending memory by one step in this specific context;
#'   \item the Kullback-Leibler divergence \eqn{D_{KL}(P_k \,\|\, P_1)} -
#'     bits of "surprise" if you used the order-1 model when the order-k
#'     model is true.
#' }
#'
#' \code{KL = 0} means longer history adds no information for that context.
#' \code{H_drop > 0} means longer history sharpens the prediction;
#' \code{H_drop < 0} indicates the order-k context happens to spread
#' probability across more outcomes than order-1 alone (small-sample noise
#' or genuine context-induced uncertainty - inspect \code{n}).
#'
#' Contexts where \code{flips = TRUE} are the substantively interesting
#' ones: the longer history changes the *modal* prediction, not just its
#' confidence.
#'
#' Pair this with \code{\link{markov_order_test}} (which decides whether
#' order-k is needed *globally*) to see the chain-level decision broken
#' down per context.
#'
#' @examples
#' \donttest{
#' # Sequences with a genuine order-2 rule: A then B is always followed by C
#' set.seed(7)
#' seqs <- t(replicate(60, {
#'   s <- character(12)
#'   s[1:2] <- sample(c("A", "B", "C"), 2, replace = TRUE)
#'   for (j in 3:12) {
#'     s[j] <- if (s[j - 2] == "A" && s[j - 1] == "B") "C"
#'             else sample(c("A", "B", "C"), 1)
#'   }
#'   s
#' }))
#' seqs <- as.data.frame(seqs)
#' pd <- path_dependence(seqs, order = 2, min_count = 3)
#' print(pd)
#' summary(pd)
#' plot(pd)
#' }
#'
#' @seealso \code{\link{markov_order_test}}, \code{\link{build_mogen}}
#'
#' @references
#' Cover, T.M. & Thomas, J.A. (2006). \emph{Elements of Information
#' Theory}, 2nd ed., chapters 2 and 4. Wiley. (KL divergence and
#' conditional entropy.)
#'
#' @export
path_dependence <- function(x, order = 2L, min_count = 5L, base = 2) {
  stopifnot(
    "'order' must be a whole number" =
      is.numeric(order) && length(order) == 1L &&
      order == round(order)
  )
  order     <- as.integer(order)
  min_count <- as.integer(min_count)
  stopifnot(
    "'order' must be >= 2" = order >= 2L,
    "'min_count' must be >= 1" = min_count >= 1L,
    is.numeric(base), length(base) == 1L, base > 0, base != 1
  )

  data <- .pd_extract_data(x)

  ## Counts at order k (full context + next) and at order 1 (last state + next)
  c_k <- .pd_kgram_counts(data, k = order + 1L)
  c_1 <- .pd_kgram_counts(data, k = 2L)

  if (length(c_k) == 0L)
    stop("No complete order-", order, " contexts found in data.",
         call. = FALSE)

  ## Order-1 conditional table P_1(next | last)
  parts1 <- .pd_split_kgram(names(c_1), k = 2L)
  o1_df  <- data.frame(
    last = parts1$context, next_state = parts1$next_state,
    count = as.integer(c_1), stringsAsFactors = FALSE
  )
  o1_totals <- tapply(o1_df$count, o1_df$last, sum)
  o1_df$p   <- o1_df$count / o1_totals[o1_df$last]

  states <- sort(unique(c(parts1$context, parts1$next_state)))

  build_p1 <- function(last_state) {
    out <- setNames(numeric(length(states)), states)
    sub <- o1_df[o1_df$last == last_state, ]
    out[sub$next_state] <- sub$p
    out
  }
  p1_cache <- vapply(states, build_p1, numeric(length(states)))
  ## p1_cache is states x states matrix: rows = next_state, cols = last_state

  ## Order-k conditional table
  parts_k <- .pd_split_kgram(names(c_k), k = order + 1L)
  ok_df   <- data.frame(
    context = parts_k$context, next_state = parts_k$next_state,
    count   = as.integer(c_k), stringsAsFactors = FALSE
  )
  ok_df$last_state <- vapply(strsplit(ok_df$context, " -> ", fixed = TRUE),
                             function(p) p[length(p)], character(1))
  ok_totals <- tapply(ok_df$count, ok_df$context, sum)
  ok_df$p_k <- ok_df$count / ok_totals[ok_df$context]

  ## Aggregate per context, dropping low-count contexts
  ctxs <- names(ok_totals)[ok_totals >= min_count]

  ctx_rows <- do.call(rbind, lapply(ctxs, function(c) {
    sub  <- ok_df[ok_df$context == c, ]
    last <- sub$last_state[1L]
    p_k  <- setNames(numeric(length(states)), states)
    p_k[sub$next_state] <- sub$p_k
    p_1  <- p1_cache[, last]

    H1 <- .pd_entropy(p_1, base = base)
    Hk <- .pd_entropy(p_k, base = base)
    KL <- .pd_kl(p_k, p_1, base = base)

    top_o1 <- states[which.max(p_1)]
    top_ok <- states[which.max(p_k)]

    data.frame(
      context  = c,
      n        = as.integer(ok_totals[c]),
      H_order1 = H1,
      H_orderk = Hk,
      H_drop   = H1 - Hk,
      KL       = KL,
      top_o1   = top_o1,
      top_ok   = top_ok,
      flips    = top_o1 != top_ok,
      stringsAsFactors = FALSE
    )
  }))
  ctx_rows <- ctx_rows[order(-ctx_rows$KL), ]
  rownames(ctx_rows) <- NULL

  finite_KL <- is.finite(ctx_rows$KL)
  total_n   <- sum(ctx_rows$n[finite_KL])
  KL_weighted     <- if (total_n > 0)
    sum(ctx_rows$KL[finite_KL]     * ctx_rows$n[finite_KL]) / total_n else 0
  H_drop_weighted <- if (total_n > 0)
    sum(ctx_rows$H_drop[finite_KL] * ctx_rows$n[finite_KL]) / total_n else 0

  structure(
    list(
      contexts   = ctx_rows,
      chain      = list(
        KL_weighted     = KL_weighted,
        H_drop_weighted = H_drop_weighted,
        n_contexts      = nrow(ctx_rows),
        n_flips         = sum(ctx_rows$flips)
      ),
      order      = order,
      base       = base,
      min_count  = min_count,
      states     = states
    ),
    class = "net_path_dependence"
  )
}


#' Print method for `net_path_dependence`
#'
#' @param x A `net_path_dependence` object.
#' @param top Integer. Number of top contexts to show. Default 10.
#' @param digits Integer. Digits to round numeric output. Default 3.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.net_path_dependence <- function(x, top = 10L, digits = 3L, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf("Path Dependence (order %d vs order 1, %s)\n\n",
              x$order, unit))
  cat(sprintf("Contexts: %d (min_count = %d).  Modal-prediction flips: %d.\n",
              x$chain$n_contexts, x$min_count, x$chain$n_flips))
  cat(sprintf("Chain-level KL_weighted      = %.*f %s\n",
              digits, x$chain$KL_weighted, unit))
  cat(sprintf("Chain-level H_drop_weighted  = %.*f %s\n\n",
              digits, x$chain$H_drop_weighted, unit))

  n_top <- min(top, nrow(x$contexts))
  cat(sprintf("Top %d contexts by KL:\n", n_top))
  tbl <- x$contexts[seq_len(n_top), ]
  tbl$H_order1 <- round(tbl$H_order1, digits)
  tbl$H_orderk <- round(tbl$H_orderk, digits)
  tbl$H_drop   <- round(tbl$H_drop,   digits)
  tbl$KL       <- round(tbl$KL,       digits)
  print(tbl, row.names = FALSE)
  invisible(x)
}


#' Summary method for `net_path_dependence`
#'
#' @param object A `net_path_dependence` object.
#' @param ... Ignored.
#' @return A `summary.net_path_dependence` with the full sorted table and
#'   chain-level summaries.
#' @export
summary.net_path_dependence <- function(object, ...) {
  structure(
    list(
      contexts   = object$contexts,
      chain      = object$chain,
      order      = object$order,
      base       = object$base,
      min_count  = object$min_count,
      n_states   = length(object$states)
    ),
    class = "summary.net_path_dependence"
  )
}

#' Print method for `summary.net_path_dependence`
#'
#' @param x A `summary.net_path_dependence` object.
#' @param digits Integer. Digits to round numeric output. Default 3.
#' @param ... Ignored.
#' @return `x` invisibly.
#' @export
print.summary.net_path_dependence <- function(x, digits = 3L, ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"

  cat(sprintf("Path Dependence Summary (order %d vs order 1, %s)\n\n",
              x$order, unit))
  cat(sprintf("All %d contexts (sorted by KL):\n", nrow(x$contexts)))
  tbl <- x$contexts
  tbl$H_order1 <- round(tbl$H_order1, digits)
  tbl$H_orderk <- round(tbl$H_orderk, digits)
  tbl$H_drop   <- round(tbl$H_drop,   digits)
  tbl$KL       <- round(tbl$KL,       digits)
  print(tbl, row.names = FALSE)

  cat("\nChain-level:\n")
  ch <- data.frame(
    KL_weighted     = round(x$chain$KL_weighted,     digits),
    H_drop_weighted = round(x$chain$H_drop_weighted, digits),
    n_contexts      = x$chain$n_contexts,
    n_flips         = x$chain$n_flips
  )
  print(ch, row.names = FALSE)
  invisible(x)
}


#' Plot method for `net_path_dependence`
#'
#' @description
#' Lollipop chart of per-context KL divergence, sorted descending. Point
#' size is proportional to context count; points where the modal next
#' state flips between orders are marked with an X to highlight
#' substantively meaningful order-2 effects.
#'
#' @param x A `net_path_dependence` object.
#' @param top Integer. Number of contexts to show (top by KL). Default 15.
#' @param title Character. Plot title.
#' @param ... Ignored.
#' @return A ggplot object.
#' @export
plot.net_path_dependence <- function(x,
                                     top = 15L,
                                     title = NULL,
                                     ...) {
  unit <- switch(as.character(x$base),
                 "2"  = "bits",
                 "10" = "hartleys",
                 sprintf("log_%g", x$base))
  if (isTRUE(all.equal(x$base, exp(1)))) unit <- "nats"
  if (is.null(title))
    title <- sprintf("Path dependence: order %d vs order 1", x$order)

  df <- x$contexts
  df <- df[seq_len(min(top, nrow(df))), ]
  df$context <- factor(df$context, levels = rev(df$context))

  ggplot2::ggplot(df, ggplot2::aes(x = .data$KL, y = .data$context)) +
    ggplot2::geom_segment(ggplot2::aes(xend = 0, yend = .data$context),
                          colour = "grey70") +
    ggplot2::geom_point(ggplot2::aes(size = .data$n,
                                     colour = .data$flips)) +
    ggplot2::geom_text(
      data = df[df$flips, ],
      ggplot2::aes(label = sprintf("%s -> %s", .data$top_o1, .data$top_ok)),
      hjust = -0.15, size = 3.2, colour = "#D55E00") +
    ggplot2::scale_colour_manual(
      values = c(`FALSE` = "#0072B2", `TRUE` = "#D55E00"),
      labels = c("modal next-state same",
                 "modal next-state changes"),
      name = NULL) +
    ggplot2::scale_size_continuous(name = "context count") +
    # The flip labels are drawn to the RIGHT of each point (hjust < 0), so the
    # top-ranked context's label runs past the panel and ggplot clips it --
    # the highest-KL row is exactly the one a reader most wants to read.
    # Widening the figure cannot fix this; the panel has to make room.
    ggplot2::scale_x_continuous(
      expand = ggplot2::expansion(mult = c(0.02, 0.35))) +
    ggplot2::labs(
      x = sprintf("KL(P_order%d || P_order1)  (%s)", x$order, unit),
      y = NULL,
      title = title,
      subtitle = sprintf(
        "Chain-level weighted KL = %.3f %s; %d / %d contexts flip the modal next state",
        x$chain$KL_weighted, unit, x$chain$n_flips, x$chain$n_contexts)) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      plot.title         = ggplot2::element_text(face = "bold"),
      legend.position    = "bottom"
    )
}

#' Coerce a net_path_dependence to its tidy per-context table
#'
#' @param x A `net_path_dependence` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param min_count Integer or `NULL`. Keep only contexts observed at least
#'   this many times.
#' @param sort_by `NULL` (construction order, default) or `"KL"` - sort by
#'   the order-k vs order-1 divergence, largest first.
#' @return A data.frame, one row per context, carrying the context label,
#'   its observed count, and the Kullback-Leibler divergence between its
#'   order-k and order-1 continuation distributions.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_path_dependence <- function(x, row.names = NULL,
                                              optional = FALSE, ...,
                                              min_count = NULL,
                                              sort_by = NULL, top = NULL) {
  out <- x$contexts
  if (!is.null(min_count)) {
    stopifnot("`min_count` must be a single integer >= 1" =
                is.numeric(min_count) && length(min_count) == 1L &&
                min_count >= 1)
    cnt <- intersect(c("count", "n"), names(out))
    if (!length(cnt)) {
      stop(errorCondition("no count column found on this result",
                          class = "hypernets_bad_input", call = NULL))
    }
    out <- out[out[[cnt[1L]]] >= min_count, , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, "KL")
    out <- out[order(-out[[sort_by]], out$context), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
