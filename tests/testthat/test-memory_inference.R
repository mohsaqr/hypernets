# ---- hg_bootstrap() / hg_compare() tests --------------------------------

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
  bs <- hg_bootstrap(.hi_det_seqs(), n_boot = 100, max_order = 2, seed = 1)
  df <- hg_get(bs)
  expect_true(all(df$probability == 1))
  expect_true(all(df$ci_lower == 1) && all(df$ci_upper == 1))
  expect_true(all(df$support == 1))
  expect_identical(unique(df$n_boot_used), 100L)
})

test_that("planted second-order rules have high bootstrap support", {
  bs <- hg_bootstrap(.hi_planted(), n_boot = 200, max_order = 3, seed = 2)
  ho <- hg_get(bs, order_min = 2)
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
    bs <- hg_bootstrap(seqs, n_boot = 200, max_order = 1, seed = s)
    df <- subset(hg_get(bs), from == "a" & to == "b")
    df$ci_lower <= p_b && p_b <= df$ci_upper
  }, logical(1L))
  # 95% nominal coverage: allow at most 2 misses in 10 seeded runs
  expect_gte(sum(covered), 8L)
})

test_that("observed probabilities match a hand count", {
  seqs <- list(c("a", "b", "a", "c"), c("a", "b", "a", "b"))
  bs <- hg_bootstrap(seqs, n_boot = 20, max_order = 1, seed = 3)
  bs_table <- hg_get(bs)
  df <- subset(bs_table, from == "a")
  # transitions from 'a': a->b (3), a->c (1)
  expect_identical(df$count, c(3L, 1L))
  expect_equal(df$probability, c(0.75, 0.25), tolerance = 1e-12)
})

test_that("parallel and serial bootstrap are identical under a seed", {
  set.seed(10)
  seqs <- replicate(8, sample(c("a", "b", "c"), 15, replace = TRUE),
                    simplify = FALSE)
  b_ser <- hg_bootstrap(seqs, n_boot = 60, max_order = 2, seed = 9)
  b_par <- hg_bootstrap(seqs, n_boot = 60, max_order = 2, seed = 9,
                         parallel = TRUE)
  expect_identical(b_ser$edges, b_par$edges)
  b_rep <- hg_bootstrap(seqs, n_boot = 60, max_order = 2, seed = 9)
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
  b_list <- hg_bootstrap(seqs, n_boot = 50, max_order = 2, seed = 4)
  b_long <- hg_bootstrap(long, action = "code", actor = "id", time = "t",
                          n_boot = 50, max_order = 2, seed = 4)
  expect_identical(b_list$edges, b_long$edges)
})


# hg_compare() takes a group model; build one from two lists of sequences.
.pair <- function(x, y, max_order, labels = c("x", "y")) {
  hon(c(x, y), group = rep(labels, c(length(x), length(y))),
      max_order = max_order)
}

test_that("hg_compare detects a planted rule difference", {
  x <- replicate(10, rep(c("a", "b", "c"), 5), simplify = FALSE)
  y <- replicate(10, rep(c("a", "b", "d"), 5), simplify = FALSE)
  cmp <- hg_compare(.pair(x, y, 2, c("early", "late")), n_perm = 199,
                    seed = 4)
  expect_lt(cmp$global$p_value, 0.05)
  sig <- hg_get(cmp, significant = TRUE)
  expect_true(all(c("c", "d") %in% sig$to))
  expect_equal(subset(sig, to == "c")$diff, 1, tolerance = 1e-12)
  expect_identical(names(cmp$n_trajectories), c("early", "late"))
})

test_that("hg_compare is calibrated under the null", {
  rejections <- vapply(1:10, function(s) {
    set.seed(2000 + s)
    gen <- function() replicate(8, sample(c("a", "b", "c"), 15,
                                          replace = TRUE),
                                simplify = FALSE)
    cmp <- hg_compare(.pair(gen(), gen(), 2), n_perm = 99, seed = s)
    cmp$global$p_value < 0.05
  }, logical(1L))
  # nominal 5% level: 10 null runs should almost never reject 4+ times
  expect_lte(sum(rejections), 3L)
})

test_that("permutation p-values are valid (add-one) and BH-adjusted", {
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  cmp <- hg_compare(.pair(x, y, 2), n_perm = 49, seed = 7)
  df <- hg_get(cmp)
  expect_true(all(df$p_value > 0 & df$p_value <= 1, na.rm = TRUE))
  expect_true(all(df$p_adj >= df$p_value, na.rm = TRUE))
  # identical cohorts: no differences at all
  expect_true(all(abs(df$diff) < 1e-12, na.rm = TRUE))
  expect_false(any(df$significant))
})

test_that("methods: print, summary, plot, accessor filters", {
  bs <- hg_bootstrap(.hi_planted(), n_boot = 50, max_order = 3, seed = 5)
  expect_invisible(print(bs))
  expect_output(print(bs), "HON bootstrap")
  s <- summary(bs)
  expect_s3_class(s, "hypernets_summary")
  expect_identical(s$edges, hg_get(bs))
  expect_identical(names(s$by_order), c("order", "n_edges", "mean_support",
                                        "min_support", "mean_ci_width"))
  high_support <- hg_get(bs, min_support = 0.95)
  expect_true(all(high_support$support >= 0.95))
  grDevices::pdf(NULL)
  p <- plot(bs, top = 5)
  expect_s3_class(p, "ggplot")
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "d"), 4), simplify = FALSE)
  cmp <- hg_compare(.pair(x, y, 2), n_perm = 49, seed = 6)
  expect_invisible(print(cmp))
  expect_output(print(cmp), "HON comparison")
  cmp_summary <- summary(cmp)
  expect_s3_class(cmp_summary$by_order, "data.frame")
  expect_named(cmp_summary$overall, c("statistic", "p_value", "n_permutations"))
  p2 <- plot(cmp, top = 5)
  expect_s3_class(p2, "ggplot")
  grDevices::dev.off()
})

test_that("error paths: invalid arguments", {
  seqs <- .hi_det_seqs(4L)
  expect_error(hg_bootstrap(seqs, n_boot = 1), "n_boot")
  expect_error(hg_bootstrap(seqs, level = 1.5), "level")
  expect_error(hg_bootstrap(list(c("a", "b"))), "at least 2 sequences")
  expect_error(hg_bootstrap(seqs, actor = "id"), "actor")
  expect_error(hg_compare(.pair(seqs, list(c("a", "b")), 2)),
               "at least 2 sequences")
  expect_error(hg_compare(.pair(seqs, seqs, 2), groups = c("x", "x")),
               class = "hypernets_bad_input")
  expect_error(hg_compare(.pair(seqs, seqs, 2), n_perm = 1), "n_perm")
  expect_error(hg_compare(seqs), class = "hypernets_bad_input")
  expect_error(
    hg_bootstrap(data.frame(a = "x"), action = "missing"), "action")
})

test_that("accessor sort_by orders deterministically", {
  bs <- hg_bootstrap(.hi_planted(), n_boot = 30, max_order = 3, seed = 8)
  df <- hg_get(bs, sort_by = "count")
  expect_true(all(diff(df$count) <= 0))
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "d"), 4), simplify = FALSE)
  cmp <- hg_compare(.pair(x, y, 2), n_perm = 49, seed = 6)
  dfc <- hg_get(cmp, sort_by = "abs_diff")
  d_sorted <- abs(dfc$diff)
  d_sorted <- d_sorted[!is.na(d_sorted)]  # NA diffs sort last by order()
  expect_true(all(diff(d_sorted) <= 1e-12))
  expect_error(hg_get(bs, sort_by = "nope"), "arg")
})


test_that("hon(group =) builds one network per group; hg_compare reads it", {
  x <- replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE)
  y <- replicate(6, rep(c("a", "b", "d"), 4), simplify = FALSE)
  z <- replicate(6, rep(c("c", "b", "a"), 4), simplify = FALSE)
  model <- hon(c(x, y, z), group = rep(c("x", "y", "z"), each = 6),
               max_order = 2)
  expect_s3_class(model, "net_hon_group")
  expect_identical(names(model), c("x", "y", "z"))
  expect_identical(model[["y"]], hon(y, max_order = 2))
  stacked <- hg_get(model)
  expect_identical(unique(stacked$group), c("x", "y", "z"))
  expect_output(print(model), "3 groups")
  # three groups need an explicit pair; the pair equals the engine
  expect_error(hg_compare(model, n_perm = 9), class = "hypernets_bad_input")
  cmp <- hg_compare(model, groups = c("x", "z"), n_perm = 19, seed = 3)
  ref <- .hg_compare_pair(x, z, n_perm = 19, max_order = 2, min_freq = 1,
                          names = c("x", "z"), seed = 3)
  expect_identical(cmp, ref)
  # long data with a group column
  long <- data.frame(id = rep(paste0("s", 1:12), each = 12),
                     t = rep(1:12, 12), code = unlist(c(x, y)),
                     period = rep(c("early", "late"), each = 72))
  by_period <- hon(long, actor = "id", action = "code", time = "t",
                   group = "period", max_order = 2)
  expect_identical(names(by_period), c("early", "late"))
  expect_identical(by_period[["late"]]$matrix, hon(y, max_order = 2)$matrix)
  # bootstrap per group
  bg <- hg_bootstrap(by_period, n_boot = 10, seed = 2)
  expect_s3_class(bg, "net_hon_boot_group")
  expect_identical(unique(summary(bg)$edges$group), c("early", "late"))
  expect_error(hon(long, actor = "id", action = "code", time = "t",
                   group = c("a", "b")), class = "hypernets_bad_input")
})
