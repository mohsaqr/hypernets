testthat::skip_on_cran()

# ---- markov_stability(): dwell, recurrence and first passage -------------
#
# The calibration fixture is a three-state ergodic chain whose stationary
# distribution and whole mean-first-passage matrix are solvable in closed
# form, so the tests below compare against hand-derived rationals rather
# than against the function's own output.
#
#        A    B    C
#   A  0.7  0.2  0.1
#   B  0.3  0.5  0.2
#   C  0.2  0.3  0.5
#
# Stationary distribution, from pi P = pi with sum(pi) = 1:
#   eq(A): 0.3 piA = 0.3 piB + 0.2 piC  =>  piA = piB + (2/3) piC
#   eq(B): 0.5 piB = 0.2 piA + 0.3 piC  =>  piB = 0.4 piA + 0.6 piC
#   solving:  piA = (19/9) piC,  piB = (13/9) piC,  sum = (41/9) piC = 1
#   =>  pi = (19, 13, 9) / 41
#
# Mean first passage times, from the first-step equations
#   m_ij = 1 + sum_{k != j} P_ik m_kj  (i != j),  m_jj = 1 / pi_j:
#   m_BA = 1.4 / 0.38 = 70/19      m_CA = 2 + 0.6 * 70/19 = 80/19
#   m_AB = 1.2 / 0.26 = 60/13      m_CB = 2 + 0.4 * 60/13 = 50/13
#   m_AC = 1.4 / 0.18 = 70/9       m_BC = 2 + 0.6 * 70/9  = 20/3

.ms_P <- matrix(
  c(0.7, 0.2, 0.1,
    0.3, 0.5, 0.2,
    0.2, 0.3, 0.5),
  nrow = 3L, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
)

# Hand-derived stationary distribution and MFPT matrix (see the header).
.ms_pi <- c(A = 19, B = 13, C = 9) / 41
.ms_M <- matrix(
  c(41 / 19, 60 / 13, 70 / 9,
    70 / 19, 41 / 13, 20 / 3,
    80 / 19, 50 / 13, 41 / 9),
  nrow = 3L, byrow = TRUE,
  dimnames = list(c("A", "B", "C"), c("A", "B", "C"))
)

.ms_seqs <- function() split(human_long$code, human_long$session_id)


# ---- Calibration against the closed-form solution ------------------------

test_that("stationary distribution matches the hand-solved chain", {
  ms <- markov_stability(.ms_P)
  d <- as.data.frame(ms)
  expect_equal(d$stationary_prob, unname(.ms_pi), tolerance = 1e-12)
  # mean recurrence time is 1 / pi, exactly 41/19, 41/13, 41/9
  expect_equal(d$return_time, 41 / c(19, 13, 9), tolerance = 1e-12)
})

test_that("mean first passage times match the first-step equations", {
  ms <- markov_stability(.ms_P)
  pt <- as.data.frame(ms, what = "passage_time")
  got <- matrix(pt$steps, nrow = 3L,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_equal(got, .ms_M, tolerance = 1e-12)
})

test_that("persistence, sojourn and the passage means follow from those", {
  ms <- markov_stability(.ms_P)
  d <- as.data.frame(ms)
  expect_equal(d$persistence, c(0.7, 0.5, 0.5), tolerance = 1e-12)
  # expected consecutive dwell = 1 / (1 - P_ii): 10/3, 2, 2
  expect_equal(d$sojourn_time, c(10 / 3, 2, 2), tolerance = 1e-12)
  # off-diagonal row and column means of the hand-derived MFPT matrix
  expect_equal(d$avg_time_to_others,
               c((60 / 13 + 70 / 9) / 2,
                 (70 / 19 + 20 / 3) / 2,
                 (80 / 19 + 50 / 13) / 2),
               tolerance = 1e-12)
  expect_equal(d$avg_time_from_others,
               c((70 / 19 + 80 / 19) / 2,
                 (60 / 13 + 50 / 13) / 2,
                 (70 / 9 + 20 / 3) / 2),
               tolerance = 1e-12)
})


# ---- Invariants ----------------------------------------------------------

test_that("pi is a stationary distribution of P", {
  ms <- markov_stability(.ms_P)
  p <- as.data.frame(ms)$stationary_prob
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_true(all(p > 0))
  # the defining property, checked directly rather than assumed
  expect_equal(as.vector(p %*% .ms_P), p, tolerance = 1e-10)
})

test_that("the passage-time diagonal is the mean recurrence time", {
  ms <- markov_stability(.ms_P)
  pt <- as.data.frame(ms, what = "passage_time")
  diag_rows <- pt[pt$from == pt$to, , drop = FALSE]
  st <- as.data.frame(ms, what = "stationary")
  expect_equal(diag_rows$steps[match(st$state, diag_rows$from)],
               st$return_time, tolerance = 1e-12)
})

test_that("the Kemeny constant is the same from every starting state", {
  # sum_j pi_j m_ij over j != i does not depend on i in an ergodic chain
  ms <- markov_stability(.ms_P)
  st <- as.data.frame(ms, what = "stationary")$stationary_prob
  pt <- as.data.frame(ms, what = "passage_time")
  M <- matrix(pt$steps, nrow = 3L)
  k <- vapply(seq_len(3L),
              function(i) sum(M[i, -i] * st[-i]), numeric(1L))
  expect_equal(k, rep(k[1L], 3L), tolerance = 1e-10)
})

test_that("relabelling the states permutes the result and changes nothing", {
  ord <- c(3L, 1L, 2L)
  ms_ref <- as.data.frame(markov_stability(.ms_P))
  ms_perm <- as.data.frame(markov_stability(.ms_P[ord, ord, drop = FALSE]))
  expect_identical(ms_perm$state, ms_ref$state[ord])
  numeric_cols <- setdiff(names(ms_ref), "state")
  # eigen() on a permuted matrix is not bit-identical, so a tolerance is
  # required here; it is the only place in this file that needs one.
  expect_equal(ms_perm[numeric_cols], ms_ref[ord, numeric_cols, drop = FALSE],
               tolerance = 1e-10, ignore_attr = "row.names")
})

test_that("an unnamed matrix gets synthetic state names, same numbers", {
  P <- unname(.ms_P)
  d <- as.data.frame(markov_stability(P))
  expect_identical(d$state, c("S1", "S2", "S3"))
  expect_equal(d$stationary_prob, unname(.ms_pi), tolerance = 1e-12)
})


# ---- Determinism ---------------------------------------------------------

test_that("markov_stability is deterministic and leaves the RNG alone", {
  set.seed(20260920L)
  before <- get(".Random.seed", envir = globalenv())
  a <- markov_stability(.ms_P)
  b <- markov_stability(.ms_P)
  after <- get(".Random.seed", envir = globalenv())
  expect_identical(a, b)
  # the verb is pure linear algebra: it must not consume the caller's stream
  expect_identical(after, before)
})


# ---- Error paths, asserted by condition class ----------------------------

test_that("a state with no outgoing transition is a classed error", {
  P_dead <- .ms_P
  P_dead["B", ] <- 0
  expect_error(markov_stability(P_dead), class = "hypernets_bad_input")
  # the message names the offending state so the user can act on it
  expect_error(markov_stability(P_dead), "B")
})

test_that("unsupported input, non-square and single-state are classed errors", {
  expect_error(markov_stability("not a matrix"),
               class = "hypernets_bad_input")
  expect_error(markov_stability(matrix(0.5, nrow = 2L, ncol = 3L)),
               class = "hypernets_bad_input")
  expect_error(markov_stability(matrix(1, 1L, 1L)),
               class = "hypernets_bad_input")
  expect_error(markov_stability(matrix(c(0.5, NA, 0.5, 0.5), 2L, 2L)),
               class = "hypernets_bad_input")
})

test_that("normalize = FALSE rejects a non-stochastic matrix", {
  expect_error(markov_stability(.ms_P * 2, normalize = FALSE),
               class = "hypernets_bad_input")
  expect_error(markov_stability(.ms_P, normalize = NA),
               "`normalize` must be")
})

test_that("normalize = TRUE rescales rows with a classed warning", {
  expect_warning(ms <- markov_stability(.ms_P * 2),
                 class = "hypernets_renormalized")
  expect_equal(as.data.frame(ms)$stationary_prob, unname(.ms_pi),
               tolerance = 1e-12)
})


# ---- Higher-order input: the reason the verb lives in hypernets ----------

test_that("a net_hon gives one row per higher-order state", {
  hon <- build_hon(.ms_seqs(), max_order = 2L, min_freq = 50L)
  ms <- markov_stability(hon)
  d <- as.data.frame(ms)
  expect_identical(d$state, colnames(hon$weights))
  expect_equal(sum(d$stationary_prob), 1, tolerance = 1e-10)
  # at least one state carries memory, i.e. this is not a first-order result
  expect_true(any(grepl(" -> ", d$state, fixed = TRUE)))
})

test_that("sequence input is the first-order chain build_hon() builds", {
  seqs <- .ms_seqs()
  from_seqs <- markov_stability(seqs)
  from_hon <- markov_stability(build_hon(seqs, max_order = 1L, min_freq = 1L))
  expect_identical(from_seqs, from_hon)
})


# ---- Taxonomy: accessor, print, summary, plot ----------------------------

test_that("as.data.frame returns tidy tables for every `what`", {
  ms <- markov_stability(.ms_P)
  states <- as.data.frame(ms)
  expect_identical(class(states), "data.frame")
  expect_named(states, c("state", "persistence", "stationary_prob",
                         "return_time", "sojourn_time",
                         "avg_time_to_others", "avg_time_from_others"))
  expect_identical(attr(states, "row.names"), seq_len(3L))

  pt <- as.data.frame(ms, what = "passage_time")
  expect_named(pt, c("from", "to", "steps"))
  expect_identical(nrow(pt), 9L)

  st <- as.data.frame(ms, what = "stationary")
  expect_named(st, c("state", "stationary_prob", "return_time"))
})

test_that("sort_by and top compose, and top applies last", {
  ms <- markov_stability(.ms_P)
  sorted <- as.data.frame(ms, sort_by = "stationary_prob")
  expect_identical(sorted$stationary_prob,
                   sort(sorted$stationary_prob, decreasing = TRUE))
  topped <- as.data.frame(ms, sort_by = "stationary_prob", top = 2L)
  expect_identical(topped, `rownames<-`(utils::head(sorted, 2L), NULL))
  expect_error(as.data.frame(ms, sort_by = "nope"), "`sort_by` must be")
  expect_error(as.data.frame(ms, top = 0L), "`top` must be")
})

test_that("print, summary and plot are all defined", {
  ms <- markov_stability(.ms_P)
  expect_output(print(ms), "Markov Stability")
  expect_output(s <- summary(ms), "attractor")
  expect_s3_class(s, "data.frame")
  expect_identical(attr(s, "attractor"), "A")
  expect_s3_class(plot(ms), "ggplot")
  expect_s3_class(plot(ms, metrics = "sojourn_time", top = 2L), "ggplot")
})

test_that("an absorbing state makes the chain reducible, and that is said", {
  P <- matrix(c(0.5, 0.5,
                0.0, 1.0), nrow = 2L, byrow = TRUE,
              dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(markov_stability(P), class = "hypernets_not_ergodic")
  # the class composes, so a caller guarding on bad input still catches it
  expect_error(markov_stability(P), class = "hypernets_bad_input")
})

test_that("a pruned higher-order network can be reducible, and is diagnosed", {
  # The regression this exists for: at min_freq = 20 this corpus leaves an
  # absorbing higher-order state, (I - P + Pi) is singular, and solve()
  # would fail with a LAPACK message that says nothing about the data.
  hon <- build_hon(.ms_seqs(), max_order = 2L, min_freq = 20L)
  expect_error(markov_stability(hon), class = "hypernets_not_ergodic")
  expect_error(markov_stability(hon), "reducible")
})

test_that("a diagonal that is 1 in double precision gives Inf sojourn", {
  # P_AA is 1 to the last bit but the A -> B edge exists, so the chain is
  # irreducible while 1 - P_AA underflows to 0.
  P <- matrix(c(1 - 1e-17, 1e-17,
                0.5,       0.5), nrow = 2L, byrow = TRUE,
              dimnames = list(c("A", "B"), c("A", "B")))
  skip_if_not(identical(P[1L, 1L], 1))
  # pi underflows to (1, 0) here, which the defensive guard reports
  expect_warning(ms <- markov_stability(P), class = "hypernets_not_ergodic")
  expect_identical(as.data.frame(ms)$sojourn_time[1L], Inf)
  expect_s3_class(plot(ms), "ggplot")
})
