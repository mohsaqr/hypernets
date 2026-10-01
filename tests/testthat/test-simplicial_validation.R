# Statistically validated window complexes (Musciotto, Battiston & Mantegna
# 2021): the exact intersection tail, the validation table and the complex.

test_that("the intersection tail of two sets is the hypergeometric tail", {
  cases <- expand.grid(n_obs = c(0L, 3L, 9L), d1 = c(10L, 25L),
                       d2 = c(12L, 30L), total = 60L)
  ours <- vapply(seq_len(nrow(cases)), \(i) {
    .svh_upper_tail(cases$n_obs[i], c(cases$d1[i], cases$d2[i]),
                    cases$total[i])
  }, numeric(1L))
  exact <- stats::phyper(cases$n_obs - 1L, cases$d1,
                         cases$total - cases$d1, cases$d2,
                         lower.tail = FALSE)
  expect_equal(ours, pmin(1, exact), tolerance = 1e-12)
})

test_that("the intersection tail does not depend on the order of the sets", {
  degrees <- c(30L, 12L, 25L, 18L)
  orders <- list(1:4, 4:1, c(2L, 4L, 1L, 3L), c(3L, 1L, 4L, 2L))
  tails <- vapply(orders, \(o) .svh_upper_tail(4L, degrees[o], 60L),
                  numeric(1L))
  expect_equal(tails, rep(tails[1L], length(tails)), tolerance = 1e-12)
})

test_that("the intersection distribution matches a simulation of the null", {
  set.seed(20260930)
  degrees <- c(30L, 25L, 20L)
  total <- 60L
  sizes <- replicate(20000L, {
    picks <- lapply(degrees, \(d) sample.int(total, d))
    length(Reduce(intersect, picks))
  })
  thresholds <- 2:9
  simulated <- vapply(thresholds, \(t) mean(sizes >= t), numeric(1L))
  exact <- vapply(thresholds, \(t) .svh_upper_tail(t, degrees, total),
                  numeric(1L))
  # binomial Monte Carlo error of a proportion from 20000 draws
  mc_se <- sqrt(exact * (1 - exact) / 20000)
  expect_true(all(abs(simulated - exact) <= 4 * mc_se + 1e-4))
})

test_that("the swap chain keeps every window size and action frequency", {
  set.seed(3)
  sets <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L, 4L), c(4L, 5L), c(1L, 5L))
  m <- .svh_incidence(sets, c(4L, 2L, 3L, 1L, 2L), 5L)
  shuffled <- .thg_swap_chain(m, attempts = 500L)
  expect_identical(rowSums(shuffled), rowSums(m))
  expect_identical(colSums(shuffled), colSums(m))
  expect_false(identical(shuffled, m))
})

# Every binary matrix with the margins of `m`, by enumerating the rows'
# possible action sets; the swap chain samples these uniformly.
.all_matrices <- function(m) {
  n_col <- ncol(m)
  row_options <- lapply(rowSums(m), \(k) utils::combn(n_col, k, simplify = FALSE))
  grid <- expand.grid(lapply(row_options, seq_along))
  mats <- lapply(seq_len(nrow(grid)), \(g) {
    rows <- Map(\(opts, j) tabulate(opts[[j]], nbins = n_col), row_options,
                unlist(grid[g, ]))
    do.call(rbind, rows)
  })
  Filter(\(x) identical(as.numeric(colSums(x)), as.numeric(colSums(m))), mats)
}

test_that("the swap null matches the exact null by enumeration", {
  set.seed(11)
  sets <- list(c(1L, 2L), c(1L, 3L), c(2L, 3L), c(3L, 4L), c(1L, 2L, 3L))
  counts <- c(2L, 1L, 1L, 1L, 1L)
  m <- .svh_incidence(sets, counts, 4L)
  targets <- .svh_row_keys(.svh_incidence(sets[1:4], rep(1L, 4L), 4L))
  exact_counts <- vapply(.all_matrices(m), \(x) {
    tabulate(match(.svh_row_keys(x), targets), nbins = length(targets))
  }, numeric(length(targets)))
  exact_mean <- rowMeans(exact_counts)
  exact_tail <- rowMeans(exact_counts >= counts[1:4])
  n_draws <- 4000L
  draws <- .svh_swap_counts(m, targets, n_draws)
  # Monte Carlo error of a mean and a proportion from n_draws
  # (autocorrelated) draws, with a factor of 3 for the autocorrelation
  se_mean <- sqrt(apply(exact_counts, 1L, stats::var) / n_draws) * 3
  se_tail <- sqrt(exact_tail * (1 - exact_tail) / n_draws) * 3
  expect_true(all(abs(rowMeans(draws) - exact_mean) <= 4 * se_mean + 1e-3))
  expect_true(all(abs(rowMeans(draws >= counts[1:4]) - exact_tail) <=
                    4 * se_tail + 1e-3))
})

# 30 actions; the triple a1-a2-a3 is planted in many sessions, the rest of
# every session is uniform noise over the alphabet
.planted_sequences <- function() {
  set.seed(7)
  alphabet <- sprintf("a%02d", seq_len(30L))
  lapply(seq_len(200L), \(i) {
    noise <- sample(alphabet, 12L, replace = TRUE)
    if (i <= 120L) c(noise[1:6], "a01", "a02", "a03", noise[7:12]) else noise
  })
}

.validation_columns <- c("members", "size", "count", "n_windows", "expected",
                         "z", "p_value", "p_adj", "n_tests", "significant")

# 6 actions, so 999 shuffles resolve every size (choose(6, 3) = 20 sets);
# the triple a1-a2-a3 is planted in half of the sessions
.small_planted <- function() {
  set.seed(4)
  alphabet <- sprintf("a%d", seq_len(6L))
  lapply(seq_len(80L), \(i) {
    noise <- sample(alphabet, 10L, replace = TRUE)
    if (i <= 40L) c(noise[1:5], "a1", "a2", "a3", noise[6:10]) else noise
  })
}

test_that("swap validation keeps the planted set and builds the complex", {
  sc <- simplicial(.small_planted(), type = "window", window = 3L,
                   validate = TRUE, seed = 1)
  v <- hg_get(sc, what = "validation")
  expect_s3_class(v, "data.frame")
  expect_named(v, .validation_columns)
  expect_true(subset(v, members == "a1, a2, a3")$significant)
  expect_gt(subset(v, members == "a1, a2, a3")$z, 3)
  # the multiple-testing base is every possible set of that size
  expect_identical(unique(subset(v, size == 3L)$n_tests), choose(6, 3))
  expect_true(all(v$p_adj >= v$p_value))
  expect_true(all(v$p_value >= 1 / 1000))
  # dim = 2 keeps the sets of three actions
  triangles <- hg_get(sc, what = "validation", dimension = 2L)
  expect_true(all(triangles$size == 3L))
  expect_identical(nrow(triangles), sum(v$size == 3L))
  # the complex holds exactly the validated sets and their faces
  simplices <- hg_get(sc)
  validated <- subset(v, significant)$members
  expect_true(all(validated %in% simplices$members))
  top_sets <- subset(simplices, size >= 2L)$members
  covered <- vapply(top_sets, \(m) {
    parts <- strsplit(m, ", ", fixed = TRUE)[[1L]]
    any(vapply(strsplit(validated, ", ", fixed = TRUE),
               \(s) all(parts %in% s), logical(1L)))
  }, logical(1L))
  expect_true(all(covered))
  # every action stays a vertex
  expect_identical(nrow(hg_get(sc, dimension = 0L)), 6L)
})

test_that("too few shuffles for the number of possible sets is reported", {
  expect_warning(simplicial(.planted_sequences(), type = "window",
                            validate = TRUE, n_null = 99L, seed = 1),
                 class = "hypernets_low_resolution")
})

test_that("the swap validation is reproducible and leaves the RNG alone", {
  seqs <- .small_planted()
  set.seed(99)
  before <- .Random.seed
  a <- simplicial(seqs, type = "window", validate = TRUE, n_null = 499L,
                  seed = 5)
  expect_identical(.Random.seed, before)
  b <- simplicial(seqs, type = "window", validate = TRUE, n_null = 499L,
                  seed = 5)
  expect_identical(hg_get(a, what = "validation"),
                   hg_get(b, what = "validation"))
})

test_that("the hypergeometric null gives the exact tail of Musciotto et al.", {
  seqs <- .planted_sequences()
  # 30 actions: the triples are dense (warned), the pairs sparse
  expect_warning(
    sc <- simplicial(seqs, type = "window", window = 3L, validate = TRUE,
                     null = "hypergeometric"),
    class = "hypernets_dense_cooccurrence")
  v <- hg_get(sc, what = "validation")
  expect_named(v, .validation_columns)
  expect_true(subset(v, members == "a01, a02, a03")$significant)
  expect_output(print(sc), "hypergeometric null")
  # a pair's p-value is the hypergeometric tail of the two actions' numbers
  # of two-action windows, counted here independently of the package code
  wh <- window_hypergraph(seqs, window = 3L)
  pairs <- lengths(wh$hyperedges) == 2L
  pair_sets <- lapply(wh$hyperedges[pairs], \(s) wh$nodes[s])
  pair_counts <- wh$window_counts[pairs]
  degree <- \(a) sum(pair_counts[vapply(pair_sets, \(s) a %in% s, logical(1L))])
  n_pairs <- sum(pair_counts)
  observed <- pair_counts[vapply(pair_sets, \(s) setequal(s, c("a01", "a02")),
                                 logical(1L))]
  row <- subset(v, members == "a01, a02")
  expect_identical(row$count, as.integer(observed))
  expect_identical(row$n_windows, as.integer(n_pairs))
  expect_equal(row$p_value,
               stats::phyper(observed - 1L, degree("a01"),
                             n_pairs - degree("a01"), degree("a02"),
                             lower.tail = FALSE), tolerance = 1e-12)
  expect_equal(row$expected, degree("a01") * degree("a02") / n_pairs,
               tolerance = 1e-12)
})

test_that("the hypergeometric null warns on dense co-occurrence only", {
  expect_warning(
    simplicial(human_long, type = "window", window = 3L, validate = TRUE,
               null = "hypergeometric", actor = "session_id",
               action = "code", time = "order_in_session"),
    class = "hypernets_dense_cooccurrence")
  # 400 actors meeting in random triples, three triples repeated: sparse,
  # so no warning, and the repeated triples are found
  set.seed(1)
  actors <- sprintf("n%03d", seq_len(400L))
  meetings <- c(replicate(3000L, sample(actors, 3L), simplify = FALSE),
                rep(list(actors[1:3], actors[4:6], actors[7:9]), each = 6L))
  # each meeting is one sequence of its members; window = 3 covers it
  expect_no_warning(
    sc <- simplicial(meetings, type = "window", window = 3L, validate = TRUE,
                     null = "hypergeometric"))
  v <- hg_get(sc, what = "validation")
  expect_setequal(subset(v, significant)$members,
                  c("n001, n002, n003", "n004, n005, n006",
                    "n007, n008, n009"))
})

test_that("a validated complex satisfies the Euler-Poincare identity", {
  sc <- simplicial(.small_planted(), type = "window", window = 3L,
                   validate = TRUE, n_null = 499L, seed = 2)
  betti <- hg_betti(sc)
  expect_identical(as.integer(hg_euler(sc)),
                   as.integer(sum((-1L)^(seq_along(betti) - 1L) * betti)))
})

test_that("validation prints its summary and persists over its counts", {
  sc <- simplicial(.small_planted(), type = "window", window = 3L,
                   validate = TRUE, seed = 3)
  expect_output(print(sc), "swap null, 999 shuffles")
  ph <- hg_homology(sc, max_dim = 2L)
  expect_s3_class(ph, "persistent_homology")
})

test_that("validation rejects a count threshold and a missing table", {
  seqs <- .planted_sequences()
  expect_error(simplicial(seqs, type = "window", validate = TRUE,
                          min_count = 5L),
               class = "hypernets_bad_input")
  plain <- simplicial(seqs, type = "window", window = 3L)
  expect_error(hg_get(plain, what = "validation"),
               class = "hypernets_bad_input")
  expect_error(simplicial(seqs, type = "window", validate = TRUE, alpha = 2))
  expect_error(simplicial(seqs, type = "window", validate = TRUE,
                          n_null = 5L))
  expect_error(simplicial(seqs, type = "window", validate = TRUE,
                          null = "binomial"))
})
