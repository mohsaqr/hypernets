# Hypergeometric anomaly detection for hypergraph co-occurrence.
#
# The memory family's build_hypa() scores transitions of a De Bruijn graph
# against a hypergeometric null whose propensity is the product of the
# endpoints' strengths (LaRock et al. 2020). The same null applies to a
# hypergraph: a pair of nodes co-occurs in some number of hyperedges, and the
# question is whether that count exceeds what the two hyperdegrees predict.
#
# Two things keep this tractable where build_hypa() is not. Only pairs that
# actually co-occur are scored, so the propensity is held sparse rather than
# as the dense n x n outer product build_hypa() materialises. And `min_count`
# drops cells that cannot reach significance before they enter the
# multiplicity correction.

#' Hypergeometric anomaly detection for co-occurring node pairs
#'
#' Scores every co-occurring pair of nodes against a hypergeometric null in
#' which a pair's propensity to share hyperedges is the product of the two
#' nodes' hyperdegrees. A pair that co-occurs far more often than that
#' product predicts is over-represented; far less often, under-represented.
#' This is the hypergraph counterpart of [build_hypa()], which applies the
#' same null to the transitions of a higher-order network, and it returns the
#' same columns so the two read alike.
#'
#' The test is analytic, so it needs no resampling. That is what makes it
#' usable where [hg_null_test()] is not: [hg_null_test()] answers whether the
#' hypergraph as a whole carries more repetition than chance, by Monte Carlo,
#' and returns one row; this returns one row per pair.
#'
#' @details
#' For a pair \eqn{(u, v)} the propensity is \eqn{\Xi_{uv} = d_u d_v}, the
#' product of the hyperdegrees. With \eqn{N = \sum \Xi} over scored pairs,
#' \eqn{K = \Xi_{uv}} and \eqn{n} the total co-occurrence mass, the observed
#' count is referred to a \eqn{\mathrm{Hypergeometric}(N, K, n)} law:
#' `p_under` is \eqn{P(X \le f)} and `p_over` is \eqn{P(X \ge f)}.
#'
#' `min_count` is not a convenience. A pair seen once or twice cannot reach
#' significance after correction, so scoring it spends multiplicity budget
#' that the pairs which can reach significance then have to pay for. Raising
#' the floor raises the power of the whole test. The floor selects the
#' hypotheses by count, never by p-value, so restricting the correction to
#' them is legitimate. \eqn{N} and the draw count are still taken over every
#' co-occurring pair, so the null describes the whole hypergraph.
#'
#' Time is handled by scoring a cumulative snapshot rather than by slicing
#' into periods. Disjoint periods split a pair's evidence across cells while
#' multiplying the number of tests, so both terms move the wrong way; a
#' sequence of `hypergraph_snapshot(hg, at = ..., mode = "cumulative")`
#' accumulates evidence instead and shows when a pair becomes anomalous.
#'
#' @param hg A `net_hypergraph`, or a snapshot of a temporal one.
#' @param min_count Minimum observed co-occurrence count for a pair to be
#'   scored (default `5L`). Pairs below the floor are not tested and do not
#'   enter the multiplicity correction.
#' @param alpha Significance threshold for the `anomaly` column
#'   (default `0.05`), applied to the adjusted p-values.
#' @param p_adjust Multiplicity correction, any [stats::p.adjust] method or
#'   `"none"` (default `"BH"`).
#' @param top Optional number of rows to keep, the most significantly
#'   over-represented first. `NULL` (default) keeps every scored pair.
#' @return A base `data.frame`, one row per scored pair, most significantly
#'   over-represented first (by `p_adjusted_over`, ties broken by `ratio`),
#'   with columns `from`, `to`, `observed` (co-occurrence count),
#'   `expected`, `ratio` (observed / expected), `p_under`, `p_over`,
#'   `p_adjusted_under`, `p_adjusted_over` and `anomaly` (`"over"`,
#'   `"under"` or `"normal"`).
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies in
#' time series data on networks. *Proceedings of the 2020 SIAM International
#' Conference on Data Mining*, 460-468.
#' @seealso [build_hypa()] for the same null on a higher-order network,
#'   [hg_null_test()] for the Monte Carlo whole-hypergraph tests.
#' @examples
#' memberships <- data.frame(
#'   actor = c("a","b","c", "a","b","d", "a","b","e", "c","d","e", "a","b","f"),
#'   group = rep(c("g1","g2","g3","g4","g5"), each = 3)
#' )
#' hg <- group_hypergraph(memberships, actor = "actor", group = "group")
#' hg_hypa(hg, min_count = 2L)
#' @export
hg_hypa <- function(hg, min_count = 5L, alpha = 0.05, p_adjust = "BH",
                    top = NULL) {
  .thg_check_hg(hg)
  stopifnot(
    "`min_count` must be a single count >= 1" =
      length(min_count) == 1L && is.finite(min_count) && min_count >= 1,
    "`alpha` must be a single number in (0, 1)" =
      length(alpha) == 1L && is.finite(alpha) && alpha > 0 && alpha < 1,
    "`top` must be NULL or a single count >= 1" =
      is.null(top) || (length(top) == 1L && is.finite(top) && top >= 1)
  )
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% c(stats::p.adjust.methods, "none")) {
    stop(errorCondition(
      "`p_adjust` must be one of stats::p.adjust.methods or \"none\"",
      class = "hypernets_bad_input", call = NULL
    ))
  }

  b <- .thg_binary(hg$incidence)
  co <- Matrix::tcrossprod(b)              # node x node co-occurrence counts
  co <- Matrix::triu(co, k = 1L)           # each unordered pair once
  all_pairs <- Matrix::summary(methods::as(co, "TsparseMatrix"))
  hits <- all_pairs[all_pairs$x >= min_count, , drop = FALSE]
  if (nrow(hits) == 0L) {
    stop(errorCondition(
      sprintf("no pair co-occurs at least %d times", as.integer(min_count)),
      class = "hypernets_empty_result", call = NULL
    ))
  }

  # The null is the whole co-occurrence ensemble, as in build_hypa(): N and
  # the draw count come from every co-occurring pair, not only the pairs above
  # the floor. Estimating them from the survivors alone would make the null
  # self-fulfilling -- with one survivor, expected would always equal observed.
  degree <- as.numeric(Matrix::rowSums(b))
  total <- sum(degree[all_pairs$i] * degree[all_pairs$j])
  draws <- sum(all_pairs$x)
  xi <- degree[hits$i] * degree[hits$j]
  observed <- hits$x

  # Same inclusive tails and clamping as .hypa_compute_scores(), so the two
  # families report the same quantity.
  n_total <- round(total)
  k_white <- pmax(pmin(round(xi), n_total), 0)
  n_draw <- min(draws, n_total)
  p_under <- stats::phyper(observed, k_white, n_total - k_white, n_draw)
  p_over <- stats::phyper(observed - 1, k_white, n_total - k_white, n_draw,
                          lower.tail = FALSE)
  expected <- n_draw * k_white / n_total

  nodes <- hg$nodes
  out <- data.frame(
    from = nodes[hits$i],
    to = nodes[hits$j],
    observed = as.integer(observed),
    expected = expected,
    ratio = ifelse(expected > 0, observed / expected, NA_real_),
    p_under = p_under,
    p_over = p_over,
    stringsAsFactors = FALSE
  )
  if (identical(p_adjust, "none")) {
    out$p_adjusted_under <- out$p_under
    out$p_adjusted_over <- out$p_over
  } else {
    out$p_adjusted_under <- stats::p.adjust(out$p_under, method = p_adjust)
    out$p_adjusted_over <- stats::p.adjust(out$p_over, method = p_adjust)
  }
  out$anomaly <- ifelse(
    out$p_adjusted_over < alpha, "over",
    ifelse(out$p_adjusted_under < alpha, "under", "normal")
  )
  # Ordered by significance, not by `ratio`. Ratio is maximised by the rarest
  # pairs sitting just above `min_count` -- a pair seen 5 times against an
  # expectation of 0.005 has a ratio near 1000 and tells you little -- so
  # ranking by it buries the heavily-cited pairs that actually depart from the
  # null. Sorting by the adjusted p-value puts weight and departure together.
  out <- out[order(out$p_adjusted_over, -out$ratio, out$from, out$to), ,
             drop = FALSE]
  rownames(out) <- NULL
  if (!is.null(top)) out <- utils::head(out, as.integer(top))
  out
}
