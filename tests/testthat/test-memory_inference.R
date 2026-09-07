# ---- bootstrap_hon() / compare_hon() tests --------------------------------

.hi_det_seqs <- function(n = 6L) {
  replicate(n, c("a", "b", "c", "a", "b", "c"), simplify = FALSE)
}

.hi_planted <- function(n_each = 8L) {
  c(replicate(n_each, rep(c("a", "b", "c"), 5), simplify = FALSE),
    replicate(n_each, rep(c("x", "b", "d"), 5), simplify = FALSE))
}

test_that("weighted count aggregation equals re-counting the multiset", {
  set.seed(21)
  trajectories <- replicate(6, sample(letters[1:4], 10, replace = TRUE),
                            simplify = FALSE)
  sc <- hypernets:::.hi_seq_counts(trajectories, max_order = 3L)
  results <- lapply(1:10, function(i) {
    w <- sample(0:3, 6, replace = TRUE)
    if (sum(w) == 0L) w[1L] <- 1L
    env_fast <- hypernets:::.hi_count_env(sc, w)
    # genuine exception to the no-loop rule is not needed: rep() expands
    # the multiset that the slow reference re-counts
    multiset <- rep(trajectories, times = w)
    env_slow <- hypernets:::.hon_build_observations(multiset, 3L)
    keys <- sort(ls(env_slow))
    expect_identical(sort(ls(env_fast)), keys)
    for (k in keys) {
      # loop over keys: comparing two environments key by key
      fast <- env_fast[[k]]
      slow <- env_slow[[k]]
      expect_identical(fast[sort(names(fast))], slow[sort(names(slow))])
    }
    w
  })
  expect_length(results, 10L)
})

test_that("deterministic sequences give degenerate CIs and full support", {
  bs <- bootstrap_hon(.hi_det_seqs(), n_boot = 100, max_order = 2, seed = 1)
  df <- as.data.frame(bs)
  expect_true(all(df$probability == 1))
  expect_true(all(df$ci_lower == 1) && all(df$ci_upper == 1))
  expect_true(all(df$support == 1))
  expect_identical(unique(df$n_boot_used), 100L)
})

test_that("planted second-order rules have high bootstrap support", {
  bs <- bootstrap_hon(.hi_planted(), n_boot = 200, max_order = 3, seed = 2)
  ho <- as.data.frame(bs, order_min = 2)
  expect_identical(sort(ho$from), c("a -> b", "x -> b"))
  expect_true(all(ho$support > 0.9))
  expect_true(all(ho$probability == 1))
})

test_that("bootstrap CI covers a known conditional probability", {
  # iid states: P(next = 'b' | current = 'a') equals the marginal P('b')
  p_b <- 0.6
  covered <- vapply(1:10, function(s) {
    set.seed(1000 + s)
    seqs <- replicate(12, sample(c("a", "b"), 40, replace = TRUE,
                                 prob = c(1 - p_b, p_b)),
                      simplify = FALSE)
    bs <- bootstrap_hon(seqs, n_boot = 200, max_order = 1, seed = s)
    df <- subset(as.data.frame(bs), from == "a" & to == "b")
    df$ci_lower <= p_b && p_b <= df$ci_upper
  }, logical(1L))
  # 95% nominal coverage: allow at most 2 misses in 10 seeded runs
  expect_gte(sum(covered), 8L)
})

test_that("observed probabilities match a hand count", {
  seqs <- list(c("a", "b", "a", "c"), c("a", "b", "a", "b"))
  bs <- bootstrap_hon(seqs, n_boot = 20, max_order = 1, seed = 3)
  bs_table <- as.data.frame(bs)
  df <- subset(bs_table, from == "a")
  # transitions from 'a': a->b (3), a->c (1)
  expect_identical(df$count, c(3L, 1L))
  expect_equal(df$probability, c(0.75, 0.25), tolerance = 1e-12)
})

test_that("parallel and serial bootstrap are identical under a seed", {
  set.seed(10)
  seqs <- replicate(8, sample(c("a", "b", "c"), 15, replace = TRUE),
                    simplify = FALSE)
  b_ser <- bootstrap_hon(seqs, n_boot = 60, max_order = 2, seed = 9)
  b_par <- bootstrap_hon(seqs, n_boot = 60, max_order = 2, seed = 9,
                         parallel = TRUE)
  expect_identical(b_ser$edges, b_par$edges)
  b_rep <- bootstrap_hon(seqs, n_boot = 60, max_order = 2, seed = 9)
  expect_identical(b_ser$edges, b_rep$edges)
})

test_that("long and list input give the same bootstrap", {
  seqs <- list(s1 = rep(c("a", "b", "c"), 4), s2 = rep(c("a", "b", "c"), 4),
               s3 = rep(c("c", "b", "a"), 4), s4 = rep(c("c", "b", "a"), 4))
  long <- data.frame(
    code = unlist(seqs, use.names = FALSE),
    id   = rep(names(seqs), times = lengths(seqs)),
    t    = unlist(lapply(lengths(seqs), seq_len), use.names = FALSE),
    stringsAsFactors = FALSE
  )
  b_list <- bootstrap_hon(seqs, n_boot = 50, max_order = 2, seed = 4)
  b_long <- bootstrap_hon(long, action = "code", actor = "id", time = "t",
                          n_boot = 50, max_order = 2, seed = 4)
  expect_identical(b_list$edges, b_long$edges)
})

test_that("compare_hon detects a planted rule difference", {
  x <- replicate(10, rep(c("a", "b", "c"), 5), simplify = FALSE)
  y <- replicate(10, rep(c("a", "b", "d"), 5), simplify = FALSE)
  cmp <- compare_hon(x, y, n_perm = 199, max_order = 2, seed = 4,
                     names = c("early", "late"))
  expect_lt(cmp$global$p_value, 0.05)
  sig <- as.data.frame(cmp, significant = TRUE)
  expect_true(all(c("c", "d") %in% sig$to))
  expect_equal(subset(sig, to == "c")$diff, 1, tolerance = 1e-12)
  expect_identical(names(cmp$n_trajectories), c("early", "late"))
})

test_that("compare_hon is calibrated under the null", {
  rejections <- vapply(1:10, function(s) {
    set.seed(2000 + s)
    gen <- function() replicate(8, sample(c("a", "b", "c"), 15,
                                          replace = TRUE),
                                simplify = FALSE)
    cmp <- compare_hon(gen(), gen(), n_perm = 99, max_order = 2, seed = s)
    cmp$global$p_value < 0.05
  }, logical(1L))
  # nominal 5% level: 10 null runs should almost never reject 4+ times
  expect_lte(sum(rejections), 3L)
})

test_that("permutation p-values are valid (add-one) and BH-adjusted", {
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  cmp <- compare_hon(x, y, n_perm = 49, max_order = 2, seed = 7)
  df <- as.data.frame(cmp)
  expect_true(all(df$p_value > 0 & df$p_value <= 1, na.rm = TRUE))
  expect_true(all(df$p_adj >= df$p_value, na.rm = TRUE))
  # identical cohorts: no differences at all
  expect_true(all(abs(df$diff) < 1e-12, na.rm = TRUE))
  expect_false(any(df$significant))
})

test_that("methods: print, summary, plot, accessor filters", {
  bs <- bootstrap_hon(.hi_planted(), n_boot = 50, max_order = 3, seed = 5)
  expect_invisible(print(bs))
  expect_output(print(bs), "HON bootstrap")
  s <- summary(bs)
  expect_s3_class(s, "data.frame")
  expect_identical(names(s), c("order", "n_edges", "mean_support",
                               "min_support", "mean_ci_width"))
  high_support <- as.data.frame(bs, min_support = 0.95)
  expect_true(all(high_support$support >= 0.95))
  grDevices::pdf(NULL)
  p <- plot(bs, top = 5)
  expect_s3_class(p, "ggplot")
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "d"), 4), simplify = FALSE)
  cmp <- compare_hon(x, y, n_perm = 49, max_order = 2, seed = 6)
  expect_invisible(print(cmp))
  expect_output(print(cmp), "HON comparison")
  cmp_summary <- summary(cmp)
  expect_s3_class(cmp_summary, "data.frame")
  p2 <- plot(cmp, top = 5)
  expect_s3_class(p2, "ggplot")
  grDevices::dev.off()
})

test_that("error paths: invalid arguments", {
  seqs <- .hi_det_seqs(4L)
  expect_error(bootstrap_hon(seqs, n_boot = 1), "n_boot")
  expect_error(bootstrap_hon(seqs, level = 1.5), "level")
  expect_error(bootstrap_hon(list(c("a", "b"))), "at least 2 sequences")
  expect_error(bootstrap_hon(seqs, actor = "id"), "actor")
  expect_error(compare_hon(seqs, list(c("a", "b"))), "at least 2 sequences")
  expect_error(compare_hon(seqs, seqs, names = c("g", "g")), "names")
  expect_error(compare_hon(seqs, seqs, n_perm = 1), "n_perm")
  expect_error(
    bootstrap_hon(data.frame(a = "x"), action = "missing"), "action")
})

test_that("accessor sort_by orders deterministically", {
  bs <- bootstrap_hon(.hi_planted(), n_boot = 30, max_order = 3, seed = 8)
  df <- as.data.frame(bs, sort_by = "count")
  expect_true(all(diff(df$count) <= 0))
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "d"), 4), simplify = FALSE)
  cmp <- compare_hon(x, y, n_perm = 49, max_order = 2, seed = 6)
  dfc <- as.data.frame(cmp, sort_by = "abs_diff")
  d_sorted <- abs(dfc$diff)
  d_sorted <- d_sorted[!is.na(d_sorted)]  # NA diffs sort last by order()
  expect_true(all(diff(d_sorted) <= 1e-12))
  expect_error(as.data.frame(bs, sort_by = "nope"), "arg")
})
