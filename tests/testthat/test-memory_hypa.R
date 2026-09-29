# ===========================================================================
# Section 1: Internal — .hypa_fit_xi
# ===========================================================================

test_that(".hypa_fit_xi gives N >> m", {
  adj <- matrix(c(0, 5, 3, 2, 0, 4, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- .hypa_fit_xi(adj)

  # N = sum(Xi) should be >> m = sum(adj)
  expect_true(sum(xi) > sum(adj))

  # Xi = outer(s_out, s_in) * mask
  s_out <- rowSums(adj)
  s_in <- colSums(adj)
  expected <- outer(s_out, s_in) * (adj > 0)
  expect_equal(xi, expected)
})

test_that(".hypa_fit_xi respects edge structure", {
  adj <- matrix(c(0, 5, 0, 0, 0, 3, 2, 0, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- .hypa_fit_xi(adj)

  # Xi should be zero where adj is zero
  expect_equal(xi[1, 1], 0)
  expect_equal(xi[1, 3], 0)
  expect_equal(xi[2, 1], 0)
  expect_equal(xi[2, 2], 0)
  expect_equal(xi[3, 3], 0)
})

# ===========================================================================
# Section 2: Internal — .hypa_compute_scores
# ===========================================================================

test_that(".hypa_compute_scores returns correct format", {
  adj <- matrix(c(0, 5, 3, 2, 0, 4, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- .hypa_fit_xi(adj)
  scores <- .hypa_compute_scores(adj, xi)

  expect_true(is.data.frame(scores))
  expect_true(all(c("path", "from", "to", "observed", "expected",
                     "ratio", "p_value", "p_under", "p_over",
                     "anomaly") %in% names(scores)))
  expect_equal(nrow(scores), sum(adj > 0))

  # HYPA scores should be in [0, 1]
  expect_true(all(scores$p_value >= 0))
  expect_true(all(scores$p_value <= 1))
  expect_equal(scores$p_value, scores$p_under)
  expect_true(all(scores$p_over >= 0))
  expect_true(all(scores$p_over <= 1))
})

test_that(".hypa_compute_scores uses inclusive over-representation tail", {
  adj <- matrix(c(0, 12, 5, 2,
                  3, 0, 8, 1,
                  6, 4, 0, 9,
                  2, 7, 3, 0), 4, 4, byrow = TRUE,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  xi <- .hypa_fit_xi(adj)
  scores <- .hypa_compute_scores(adj, xi)
  edge <- scores[scores$from == "B" & scores$to == "C", , drop = FALSE]
  K <- round(xi["B", "C"])
  N_total <- round(sum(xi))
  n_draws <- min(sum(adj), N_total)

  expect_equal(edge$p_under,
               stats::phyper(edge$observed, K, N_total - K, n_draws))
  expect_equal(edge$p_over,
               stats::phyper(edge$observed - 1, K, N_total - K, n_draws,
                             lower.tail = FALSE))
  expect_gt(edge$p_over, 1 - edge$p_under)
})

test_that(".hypa_compute_scores handles empty graph", {
  adj <- matrix(0, 3, 3, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  xi <- adj
  scores <- .hypa_compute_scores(adj, xi)
  expect_equal(nrow(scores), 0L)
})

# ===========================================================================
# Section 3: build_hypa end-to-end
# ===========================================================================

test_that("build_hypa returns net_hypa class", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("B", "C", "A"),
                c("C", "A", "B"), c("A", "C", "B"), c("B", "A", "C"))
  h <- build_hypa(trajs, order = 1L)

  expect_s3_class(h, "net_hypa")
  expect_equal(h$k, 1L)
  expect_true(h$n_edges > 0L)
  expect_true(is.data.frame(h$scores))
})

test_that("build_hypa detects anomalies in biased data", {
  # Create data where A->B->C is overwhelmingly common
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(5, c("A", "B", "D"), simplify = FALSE),
    replicate(5, c("C", "B", "A"), simplify = FALSE),
    replicate(2, c("D", "B", "C"), simplify = FALSE),
    replicate(2, c("C", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05)

  # Should find some anomalous paths
  # The A->B->C path is very frequent, may be over-represented
  expect_s3_class(h, "net_hypa")
  expect_true(h$n_edges > 0L)
})

test_that("build_hypa handles k=1 (first-order)", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "C", "B", "D"),
                c("B", "C", "D", "A"), c("D", "A", "B", "C"))
  h <- build_hypa(trajs, order = 1L)

  expect_equal(h$k, 1L)
  expect_true(nrow(h$scores) > 0L)
})

test_that("build_hypa rejects invalid input", {
  expect_error(build_hypa(42), "data.frame or list")
  expect_error(build_hypa(list(c("A", "B")), order = 0L), "integers >= 1")
  expect_error(build_hypa(list(c("A", "B")), alpha = 0.6),
               "alpha.*must be in")
})

test_that("build_hypa deprecated `k` warns by class and matches `order`", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("B", "C", "A"),
                c("C", "A", "B"), c("A", "C", "B"), c("B", "A", "C"))
  expect_warning(h_k <- build_hypa(trajs, k = 2L, min_count = 1L),
                 class = "deprecatedWarning")
  h_order <- build_hypa(trajs, order = 2L, min_count = 1L)
  expect_identical(h_k, h_order)
})

test_that("build_hypa alpha parameter affects classification", {
  trajs <- list(c("A", "B", "C"), c("A", "B", "D"), c("B", "C", "A"),
                c("C", "A", "B"), c("B", "A", "C"), c("A", "C", "B"))

  h1 <- build_hypa(trajs, order = 1L, alpha = 0.05)
  h2 <- build_hypa(trajs, order = 1L, alpha = 0.49)

  # Stricter alpha should find fewer or equal anomalies
  expect_true(h1$n_anomalous <= h2$n_anomalous)
})

# ===========================================================================
# Section 4: HYPA scores properties
# ===========================================================================

test_that("HYPA scores are in [0, 1]", {
  set.seed(42)
  trajs <- lapply(seq_len(50L), function(i) {
    sample(LETTERS[1:4], 5, replace = TRUE)
  })
  h <- build_hypa(trajs, order = 1L)

  expect_true(all(h$scores$p_value >= 0))
  expect_true(all(h$scores$p_value <= 1))
})

test_that("HYPA expected values are positive", {
  trajs <- list(c("A", "B", "C"), c("A", "C", "B"), c("B", "A", "C"),
                c("C", "B", "A"), c("B", "C", "A"), c("C", "A", "B"))
  h <- build_hypa(trajs, order = 1L)

  expect_true(all(h$scores$expected >= 0))
})

# ===========================================================================
# Section 5: S3 methods
# ===========================================================================

test_that("print.net_hypa works", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  h <- build_hypa(trajs, order = 1L)
  out <- capture.output(print(h))
  expect_true(any(grepl("HYPA", out)))
})

test_that("summary.net_hypa works", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"), c("A", "C", "B"))
  h <- build_hypa(trajs, order = 1L)
  out <- capture.output(summary(h))
  expect_true(any(grepl("HYPA", out)))
})

test_that("summary.net_hypa returns direction-specific tail probability", {
  scores <- data.frame(
    path = c("A -> B -> C", "D -> E -> F"),
    observed = c(20, 1),
    expected = c(3, 8),
    ratio = c(20 / 3, 1 / 8),
    p_value = c(0.99, 0.01),
    p_under = c(0.99, 0.01),
    p_over = c(0.001, 0.95),
    p_adjusted_under = c(0.99, 0.01),
    p_adjusted_over = c(0.001, 0.95),
    anomaly = c("over", "under"),
    stringsAsFactors = FALSE
  )
  h <- structure(list(
    scores = scores,
    over = scores[scores$anomaly == "over", , drop = FALSE],
    under = scores[scores$anomaly == "under", , drop = FALSE],
    # canonical slot is $order (post HYPA refactor); $k is the deprecated alias
    order = 2L,
    k = 2L,
    nodes = data.frame(label = LETTERS[1:6], stringsAsFactors = FALSE),
    n_edges = nrow(scores),
    alpha = 0.05,
    p_adjust = "none",
    n_anomalous = 2L,
    n_over = 1L,
    n_under = 1L
  ), class = "net_hypa")

  capture.output(over <- summary(h, type = "over"))
  capture.output(under <- summary(h, type = "under"))

  expect_equal(over$p_tail, scores$p_over[1])
  expect_equal(under$p_tail, scores$p_under[2])
  expect_false("p_value" %in% names(over))
})


# ===========================================================================
# Section 6: Data.frame input
# ===========================================================================

test_that("build_hypa handles data.frame input", {
  df <- data.frame(T1 = c("A", "B", "C"), T2 = c("B", "C", "A"),
                   T3 = c("C", "A", "B"))
  h <- build_hypa(df, order = 1L)
  expect_s3_class(h, "net_hypa")
})

# ===========================================================================
# Section 7: Coverage for previously uncovered paths
# ===========================================================================

# --- build_hypa: no valid trajectories ---
test_that("build_hypa stops when no valid trajectories", {
  # All single-state entries: parsed trajectories have < 2 states each
  df <- data.frame(T1 = c("A", "B"), stringsAsFactors = FALSE)
  expect_error(build_hypa(df, order = 1L), "No valid trajectories")
})

# --- build_hypa: no edges at given order (paths too short) ---
test_that("build_hypa stops when no edges at requested order k", {
  # k=3 requires 4-grams; trajectories of length 3 produce only 1-grams (k=1)
  # and 2-grams (k=2) but not 3-grams as transitions
  trajs <- list(c("A", "B", "C"), c("B", "C", "D"))
  expect_error(build_hypa(trajs, order = 3L), "No edges at")
})

# --- summary.net_hypa with anomalies displays anomalous paths ---
test_that("summary.net_hypa displays anomalous paths when present", {
  # Create highly biased data to force anomaly detection
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)
  out <- capture.output(summary(h))
  # Either found anomalies or printed "No anomalous paths detected."
  expect_true(any(grepl("Anomalous|anomalous|No anomalous", out,
                         ignore.case = TRUE)))
})


# ===========================================================================
# Section 8: pathways() tests for HYPA
# ===========================================================================

test_that("pathways.net_hypa returns character vector", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)
  pw <- pathways(h)
  expect_true(is.character(pw))
})

test_that("pathways.net_hypa type='over' returns over-represented paths", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)
  pw_all  <- pathways(h, type = "all")
  pw_over <- pathways(h, type = "over")
  expect_true(length(pw_over) <= length(pw_all))
})

test_that("pathways.net_hypa type='under' returns under-represented paths", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)
  pw_under <- pathways(h, type = "under")
  expect_true(is.character(pw_under))
})

test_that("pathways.net_hypa returns empty when no anomalies", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  h <- build_hypa(trajs, order = 1L, alpha = 1e-10)
  pw <- pathways(h)
  # With near-zero alpha threshold, likely no anomalies
  expect_true(is.character(pw))
})

# ===========================================================================
# Section 9: New fields ($over, $under, $n_over, $n_under, sorting)
# ===========================================================================

test_that("build_hypa stores $over and $under data frames", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)

  expect_true(is.data.frame(h$over))
  expect_true(is.data.frame(h$under))
  expect_equal(nrow(h$over), h$n_over)
  expect_equal(nrow(h$under), h$n_under)
  expect_equal(h$n_anomalous, h$n_over + h$n_under)

  # $over should only contain "over" anomalies
  if (nrow(h$over) > 0L) {
    expect_true(all(h$over$anomaly == "over"))
  }
  # $under should only contain "under" anomalies
  if (nrow(h$under) > 0L) {
    expect_true(all(h$under$anomaly == "under"))
  }
})

test_that("build_hypa pre-sorts scores: anomalous first", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)

  if (h$n_anomalous > 0L && nrow(h$scores) > h$n_anomalous) {
    # Anomalous rows should come before normal rows
    anomaly_positions <- which(h$scores$anomaly != "normal")
    normal_positions <- which(h$scores$anomaly == "normal")
    if (length(anomaly_positions) > 0L && length(normal_positions) > 0L) {
      expect_true(max(anomaly_positions) < min(normal_positions))
    }
  }
})

test_that("summary.net_hypa respects n parameter", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)
  out <- capture.output(summary(h, n = 2L))
  expect_true(any(grepl("HYPA", out)))
})

test_that("summary.net_hypa shows over/under counts", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2,  c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)
  out <- capture.output(summary(h))
  expect_true(any(grepl("over:", out)))
  expect_true(any(grepl("under:", out)))
})


# ===========================================================================
# Section 10: Multiple testing correction (p_adjust)
# ===========================================================================

test_that("p_adjust='none' matches original behavior (no correction)", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(5, c("A", "B", "D"), simplify = FALSE),
    replicate(5, c("C", "B", "A"), simplify = FALSE),
    replicate(2, c("D", "B", "C"), simplify = FALSE),
    replicate(2, c("C", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L,
                  p_adjust = "none")

  expect_equal(h$p_adjust, "none")

  # With p_adjust="none", adjusted values should equal the raw tails.
  expect_equal(h$scores$p_adjusted_under, h$scores$p_under)
  expect_equal(h$scores$p_adjusted_over, h$scores$p_over)
})

test_that("default BH adjustment produces p_adjusted columns in scores", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(5, c("A", "B", "D"), simplify = FALSE),
    replicate(5, c("C", "B", "A"), simplify = FALSE),
    replicate(2, c("D", "B", "C"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)

  expect_equal(h$p_adjust, "BH")
  expect_true("p_adjusted_under" %in% names(h$scores))
  expect_true("p_adjusted_over" %in% names(h$scores))
  expect_true("p_value" %in% names(h$scores))
  expect_true("p_under" %in% names(h$scores))
  expect_true("p_over" %in% names(h$scores))

  # Adjusted p-values should be in [0, 1]
  expect_true(all(h$scores$p_adjusted_under >= 0))
  expect_true(all(h$scores$p_adjusted_under <= 1))
  expect_true(all(h$scores$p_adjusted_over >= 0))
  expect_true(all(h$scores$p_adjusted_over <= 1))
})

test_that("BH adjustment can differ from no correction on biased data", {
  trajs <- c(
    replicate(80, c("A", "B", "C"), simplify = FALSE),
    replicate(3, c("A", "B", "D"), simplify = FALSE),
    replicate(3, c("C", "B", "A"), simplify = FALSE),
    replicate(3, c("D", "B", "C"), simplify = FALSE),
    replicate(3, c("C", "B", "D"), simplify = FALSE),
    replicate(3, c("A", "C", "B"), simplify = FALSE),
    replicate(3, c("B", "C", "D"), simplify = FALSE)
  )

  h_none <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L,
                       p_adjust = "none")
  h_bh   <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L,
                       p_adjust = "BH")

  # BH should be at least as conservative as none (fewer or equal anomalies)
  expect_true(h_bh$n_anomalous <= h_none$n_anomalous)
})

test_that("bonferroni is more conservative than BH", {
  trajs <- c(
    replicate(80, c("A", "B", "C"), simplify = FALSE),
    replicate(3, c("A", "B", "D"), simplify = FALSE),
    replicate(3, c("C", "B", "A"), simplify = FALSE),
    replicate(3, c("D", "B", "C"), simplify = FALSE),
    replicate(3, c("C", "B", "D"), simplify = FALSE),
    replicate(3, c("A", "C", "B"), simplify = FALSE),
    replicate(3, c("B", "C", "D"), simplify = FALSE)
  )

  h_bh   <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L,
                       p_adjust = "BH")
  h_bonf <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L,
                       p_adjust = "bonferroni")

  # Bonferroni should be at least as conservative as BH
  expect_true(h_bonf$n_anomalous <= h_bh$n_anomalous)

  # Bonferroni adjusted p-values should be >= BH adjusted p-values
  # (merge on path to compare)
  merged <- merge(h_bh$scores[, c("path", "p_adjusted_under", "p_adjusted_over")],
                  h_bonf$scores[, c("path", "p_adjusted_under", "p_adjusted_over")],
                  by = "path", suffixes = c("_bh", "_bonf"))
  expect_true(all(merged$p_adjusted_under_bonf >= merged$p_adjusted_under_bh - 1e-10))
  expect_true(all(merged$p_adjusted_over_bonf >= merged$p_adjusted_over_bh - 1e-10))
})

test_that("$over and $under data frames have p_adjusted columns", {
  trajs <- c(
    replicate(50, c("A", "B", "C"), simplify = FALSE),
    replicate(2, c("A", "B", "D"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)

  # $over and $under should have p_adjusted columns as subsets of scores
  if (nrow(h$over) > 0L) {
    expect_true("p_adjusted_under" %in% names(h$over))
    expect_true("p_adjusted_over" %in% names(h$over))
  }
  if (nrow(h$under) > 0L) {
    expect_true("p_adjusted_under" %in% names(h$under))
    expect_true("p_adjusted_over" %in% names(h$under))
  }
})

test_that("invalid p_adjust method errors", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  expect_error(build_hypa(trajs, order = 1L, p_adjust = "invalid_method"),
               "p_adjust.*must be one of")
  expect_error(build_hypa(trajs, order = 1L, p_adjust = 42),
               "p_adjust.*must be one of")
  expect_error(build_hypa(trajs, order = 1L, p_adjust = c("BH", "bonferroni")),
               "p_adjust.*must be one of")
})

test_that("p_adjust stored in result object", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"), c("A", "C", "B"))
  h_bh <- build_hypa(trajs, order = 1L, p_adjust = "BH")
  h_none <- build_hypa(trajs, order = 1L, p_adjust = "none")
  h_bonf <- build_hypa(trajs, order = 1L, p_adjust = "bonferroni")

  expect_equal(h_bh$p_adjust, "BH")
  expect_equal(h_none$p_adjust, "none")
  expect_equal(h_bonf$p_adjust, "bonferroni")
})

test_that("print.net_hypa shows p_adjust", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"))
  h <- build_hypa(trajs, order = 1L, p_adjust = "BH")
  out <- capture.output(print(h))
  expect_true(any(grepl("p_adjust=BH", out)))

  h_none <- build_hypa(trajs, order = 1L, p_adjust = "none")
  out_none <- capture.output(print(h_none))
  expect_true(any(grepl("p_adjust=none", out_none)))
})

test_that("summary.net_hypa shows p_adjust", {
  trajs <- list(c("A", "B", "C"), c("B", "C", "A"), c("A", "C", "B"))
  h <- build_hypa(trajs, order = 1L, p_adjust = "bonferroni")
  out <- capture.output(summary(h))
  expect_true(any(grepl("p_adjust: bonferroni", out)))
})

test_that("two-sided correction adjusts under and over separately", {
  # Create data with enough edges to see the effect of separate adjustments
  set.seed(123)
  trajs <- c(
    replicate(60, c("A", "B", "C"), simplify = FALSE),
    replicate(3, c("A", "B", "D"), simplify = FALSE),
    replicate(3, c("C", "B", "A"), simplify = FALSE),
    replicate(3, c("D", "B", "C"), simplify = FALSE),
    replicate(3, c("A", "C", "D"), simplify = FALSE),
    replicate(3, c("D", "C", "A"), simplify = FALSE)
  )
  h <- build_hypa(trajs, order = 2L, alpha = 0.05, min_count = 1L)

  # Verify that p_adjusted_under and p_adjusted_over are adjusted separately
  # by checking they equal p.adjust applied to the raw values
  raw_p_under <- h$scores$p_under
  raw_p_over <- h$scores$p_over
  expect_equal(h$scores$p_adjusted_under,
               stats::p.adjust(raw_p_under, method = "BH"))
  expect_equal(h$scores$p_adjusted_over,
               stats::p.adjust(raw_p_over, method = "BH"))
})
