testthat::skip_on_cran()

# ---- hg_markov_stability(): dwell, recurrence and first passage -------------
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
  ms <- hg_markov_stability(.ms_P)
  # the stationary table is read from Nestimate's unrounded passage-time
  # component; the "states" table is Nestimate's, rounded (checked below)
  d <- hg_get(ms, what = "stationary")
  expect_equal(d$stationary, unname(.ms_pi), tolerance = 1e-12)
  # mean recurrence time is 1 / pi, exactly 41/19, 41/13, 41/9
  expect_equal(d$return_time, 41 / c(19, 13, 9), tolerance = 1e-12)
})

test_that("mean first passage times match the first-step equations", {
  ms <- hg_markov_stability(.ms_P)
  pt <- hg_get(ms, what = "passage_time")
  got <- matrix(pt$steps, nrow = 3L,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_equal(got, .ms_M, tolerance = 1e-12)
})

test_that("persistence, sojourn and the passage means follow from those", {
  ms <- hg_markov_stability(.ms_P)
  # Nestimate stores this table rounded: probabilities to 4 decimals, times
  # to 2. The hand-derived values are rounded the same way, then compared.
  d <- hg_get(ms)
  expect_equal(d$persistence, c(0.7, 0.5, 0.5), tolerance = 1e-12)
  expect_equal(d$stationary, round(unname(.ms_pi), 4), tolerance = 1e-12)
  # expected consecutive dwell = 1 / (1 - P_ii): 10/3, 2, 2
  expect_equal(d$sojourn_time, round(c(10 / 3, 2, 2), 2), tolerance = 1e-12)
  # off-diagonal row and column means of the hand-derived MFPT matrix
  expect_equal(d$avg_time_to_others,
               round(c((60 / 13 + 70 / 9) / 2,
                       (70 / 19 + 20 / 3) / 2,
                       (80 / 19 + 50 / 13) / 2), 2),
               tolerance = 1e-12)
  expect_equal(d$avg_time_from_others,
               round(c((70 / 19 + 80 / 19) / 2,
                       (60 / 13 + 50 / 13) / 2,
                       (70 / 9 + 20 / 3) / 2), 2),
               tolerance = 1e-12)
})


# ---- Invariants ----------------------------------------------------------

test_that("pi is a stationary distribution of P", {
  ms <- hg_markov_stability(.ms_P)
  p <- hg_get(ms, what = "stationary")$stationary
  expect_equal(sum(p), 1, tolerance = 1e-12)
  expect_true(all(p > 0))
  # the defining property, checked directly rather than assumed
  expect_equal(as.vector(p %*% .ms_P), p, tolerance = 1e-10)
})

test_that("the passage-time diagonal is the mean recurrence time", {
  ms <- hg_markov_stability(.ms_P)
  pt <- hg_get(ms, what = "passage_time")
  diag_rows <- pt[pt$from == pt$to, , drop = FALSE]
  st <- hg_get(ms, what = "stationary")
  expect_equal(diag_rows$steps[match(st$state, diag_rows$from)],
               st$return_time, tolerance = 1e-12)
})

test_that("the Kemeny constant is the same from every starting state", {
  # sum_j pi_j m_ij over j != i does not depend on i in an ergodic chain
  ms <- hg_markov_stability(.ms_P)
  st <- hg_get(ms, what = "stationary")$stationary
  pt <- hg_get(ms, what = "passage_time")
  M <- matrix(pt$steps, nrow = 3L)
  k <- vapply(seq_len(3L),
              function(i) sum(M[i, -i] * st[-i]), numeric(1L))
  expect_equal(k, rep(k[1L], 3L), tolerance = 1e-10)
})

test_that("relabelling the states permutes the result and changes nothing", {
  ord <- c(3L, 1L, 2L)
  ms_ref <- hg_get(hg_markov_stability(.ms_P))
  ms_perm <- hg_get(hg_markov_stability(.ms_P[ord, ord, drop = FALSE]))
  expect_identical(ms_perm$state, ms_ref$state[ord])
  numeric_cols <- setdiff(names(ms_ref), "state")
  # eigen() on a permuted matrix is not bit-identical, so a tolerance is
  # required here; it is the only place in this file that needs one.
  expect_equal(ms_perm[numeric_cols], ms_ref[ord, numeric_cols, drop = FALSE],
               tolerance = 1e-10, ignore_attr = "row.names")
})

test_that("an unnamed matrix gets synthetic state names, same numbers", {
  P <- unname(.ms_P)
  d <- hg_get(hg_markov_stability(P), what = "stationary")
  expect_identical(d$state, c("S1", "S2", "S3"))
  expect_equal(d$stationary, unname(.ms_pi), tolerance = 1e-12)
})


# ---- Determinism ---------------------------------------------------------

test_that("markov_stability is deterministic and leaves the RNG alone", {
  set.seed(20260920L)
  before <- get(".Random.seed", envir = globalenv())
  a <- hg_markov_stability(.ms_P)
  b <- hg_markov_stability(.ms_P)
  after <- get(".Random.seed", envir = globalenv())
  expect_identical(a, b)
  # the verb is pure linear algebra: it must not consume the caller's stream
  expect_identical(after, before)
})


# ---- Error paths ----------------------------------------------------------
#
# hg_markov_stability() is Nestimate's; it raises unclassed conditions, so these
# assert that the condition is raised (and, where useful, what it names).

test_that("a state with no outgoing transition is an error naming it", {
  P_dead <- .ms_P
  P_dead["B", ] <- 0
  expect_error(hg_markov_stability(P_dead), "B")
})

test_that("unsupported input, non-square and NA matrices are errors", {
  expect_error(hg_markov_stability("not a matrix"))
  expect_error(hg_markov_stability(matrix(0.5, nrow = 2L, ncol = 3L)))
  expect_error(hg_markov_stability(matrix(c(0.5, NA, 0.5, 0.5), 2L, 2L)))
})

test_that("normalize = FALSE rejects a non-stochastic matrix", {
  expect_error(hg_markov_stability(.ms_P * 2, normalize = FALSE))
})

test_that("normalize = TRUE rescales rows with a warning", {
  expect_warning(ms <- hg_markov_stability(.ms_P * 2), "normaliz")
  expect_equal(hg_get(ms, what = "stationary")$stationary,
               unname(.ms_pi), tolerance = 1e-12)
})


# ---- Higher-order input: the reason the verb lives in hypernets ----------

test_that("a net_hon gives one row per higher-order state", {
  hon <- hon(.ms_seqs(), max_order = 2L, min_freq = 50L)
  ms <- hg_markov_stability(hon)
  d <- hg_get(ms, what = "stationary")
  expect_identical(d$state, colnames(hon$weights))
  expect_equal(sum(d$stationary), 1, tolerance = 1e-10)
  # at least one state carries memory, i.e. this is not a first-order result
  expect_true(any(grepl(" -> ", d$state, fixed = TRUE)))
})


# ---- Taxonomy: accessor, print, summary, plot ----------------------------

test_that("as.data.frame returns tidy tables for every `what`", {
  ms <- hg_markov_stability(.ms_P)
  states <- hg_get(ms)
  expect_identical(class(states), "data.frame")
  expect_named(states, c("state", "persistence", "stationary",
                         "return_time", "sojourn_time",
                         "avg_time_to_others", "avg_time_from_others"))
  expect_identical(attr(states, "row.names"), seq_len(3L))

  pt <- hg_get(ms, what = "passage_time")
  expect_named(pt, c("from", "to", "steps"))
  expect_identical(nrow(pt), 9L)

  st <- hg_get(ms, what = "stationary")
  expect_named(st, c("state", "stationary", "return_time"))
})

test_that("sort_by and top compose, and top applies last", {
  ms <- hg_markov_stability(.ms_P)
  sorted <- hg_get(ms, sort_by = "stationary")
  expect_identical(sorted$stationary,
                   sort(sorted$stationary, decreasing = TRUE))
  topped <- hg_get(ms, sort_by = "stationary", top = 2L)
  expect_identical(topped, `rownames<-`(utils::head(sorted, 2L), NULL))
  expect_error(hg_get(ms, sort_by = "nope"), "`sort_by` must be")
  expect_error(hg_get(ms, top = 0L), "`top` must be")
})

test_that("a diagonal that is 1 in double precision is an absorbing state", {
  # P_AA is 1 to the last bit, so in double precision A cannot be left and
  # the stationary distribution underflows to (1, 0): the chain is refused
  # instead of reporting an infinite sojourn and meaningless passage times.
  P <- matrix(c(1 - 1e-17, 1e-17,
                0.5,       0.5), nrow = 2L, byrow = TRUE,
              dimnames = list(c("A", "B"), c("A", "B")))
  skip_if_not(identical(P[1L, 1L], 1))
  expect_error(hg_markov_stability(P), class = "hypernets_not_ergodic")
})


# ---- Accessor: ascending order and the from/to passage filter ------------

# The accessor exactly as it stood before `decreasing`, `from` and `to` were
# added (0.5.0), reading Nestimate's object layout (the passage matrix and
# stationary vector live in its `mpt` component). The default and descending
# calls must still reproduce it.
.ms_accessor_050 <- function(x, what = "states", sort_by = NULL, top = NULL) {
  out <- switch(
    what,
    states = .ho_rename(x$stability, c(stationary_prob = "stationary")),
    stationary = data.frame(state = x$mpt$states,
                            stationary = unname(x$mpt$stationary),
                            return_time = unname(x$mpt$return_times),
                            stringsAsFactors = FALSE),
    passage_time = {
      M <- x$mpt$matrix
      data.frame(from = rep(rownames(M), times = ncol(M)),
                 to = rep(colnames(M), each = nrow(M)),
                 steps = as.vector(M), stringsAsFactors = FALSE)
    }
  )
  if (!is.null(sort_by)) {
    out <- out[order(out[[sort_by]], decreasing = TRUE), , drop = FALSE]
  }
  rownames(out) <- NULL
  if (!is.null(top)) out <- utils::head(out, top)
  rownames(out) <- NULL
  out
}

test_that("default and descending calls are unchanged from 0.5.0", {
  hon <- hon(.ms_seqs(), max_order = 2L, min_freq = 50L)
  ms <- hg_markov_stability(hon)
  expect_identical(hg_get(ms), .ms_accessor_050(ms))
  expect_identical(hg_get(ms, what = "passage_time"),
                   .ms_accessor_050(ms, what = "passage_time"))
  expect_identical(hg_get(ms, what = "stationary"),
                   .ms_accessor_050(ms, what = "stationary"))
  expect_identical(hg_get(ms, sort_by = "stationary", top = 5L),
                   .ms_accessor_050(ms, sort_by = "stationary", top = 5L))
  expect_identical(
    hg_get(ms, what = "passage_time", sort_by = "steps", top = 5L),
    .ms_accessor_050(ms, what = "passage_time", sort_by = "steps", top = 5L))
})

test_that("decreasing = FALSE sorts ascending, the reverse of descending", {
  ms <- hg_markov_stability(.ms_P)
  # every column below is untied on the calibration chain
  asc <- hg_get(ms, sort_by = "stationary", decreasing = FALSE)
  expect_identical(asc$state, c("C", "B", "A"))           # 9 < 13 < 19 / 41
  desc <- hg_get(ms, sort_by = "stationary")
  expect_identical(asc, `rownames<-`(desc[rev(seq_len(nrow(desc))), ], NULL))

  pt_asc <- hg_get(ms, what = "passage_time", sort_by = "steps",
                          decreasing = FALSE)
  expect_equal(pt_asc$steps, sort(as.vector(.ms_M)), tolerance = 1e-12)
  pt_desc <- hg_get(ms, what = "passage_time", sort_by = "steps")
  expect_identical(pt_asc,
                   `rownames<-`(pt_desc[rev(seq_len(nrow(pt_desc))), ], NULL))
  # the shortest passage is the recurrence of A, 41/19
  expect_identical(c(pt_asc$from[1L], pt_asc$to[1L]), c("A", "A"))

  # top applies after the ascending sort: the two shortest sojourns
  low <- hg_get(ms, sort_by = "sojourn_time", decreasing = FALSE,
                       top = 2L)
  expect_identical(low$state, c("B", "C"))                # both 1 / 0.5 = 2
})

test_that("ties are broken by the row key in both directions", {
  P <- matrix(c(0.5, 0.25, 0.25,
                0.25, 0.5, 0.25,
                0.25, 0.25, 0.5), nrow = 3L, byrow = TRUE,
              dimnames = list(c("Z", "Y", "X"), c("Z", "Y", "X")))
  ms <- hg_markov_stability(P)
  # all three persistences are 0.5: the key (state, ascending) decides
  expect_identical(hg_get(ms, sort_by = "persistence")$state,
                   c("X", "Y", "Z"))
  expect_identical(
    hg_get(ms, sort_by = "persistence", decreasing = FALSE)$state,
    c("X", "Y", "Z"))
})

test_that("from and to restrict the passage table to hand-derived rows", {
  ms <- hg_markov_stability(.ms_P)
  a_out <- hg_get(ms, what = "passage_time", from = "A")
  expect_identical(a_out$to, c("A", "B", "C"))
  expect_equal(a_out$steps, c(41 / 19, 60 / 13, 70 / 9), tolerance = 1e-12)

  into_c <- hg_get(ms, what = "passage_time", to = "C")
  expect_identical(into_c$from, c("A", "B", "C"))
  expect_equal(into_c$steps, c(70 / 9, 20 / 3, 41 / 9), tolerance = 1e-12)

  ab_bc <- hg_get(ms, what = "passage_time", from = c("A", "B"),
                         to = c("B", "C"), sort_by = "steps",
                         decreasing = FALSE)
  expect_identical(paste(ab_bc$from, ab_bc$to),
                   c("B B", "A B", "B C", "A C"))
  expect_equal(ab_bc$steps, c(41 / 13, 60 / 13, 20 / 3, 70 / 9),
               tolerance = 1e-12)
})

test_that("from/to on the wrong table or an unknown state is a classed error", {
  ms <- hg_markov_stability(.ms_P)
  expect_error(hg_get(ms, from = "A"), class = "hypernets_bad_input")
  expect_error(hg_get(ms, what = "stationary", to = "A"),
               class = "hypernets_bad_input")
  expect_error(hg_get(ms, what = "passage_time", from = "Q"),
               class = "hypernets_bad_input")
  expect_error(hg_get(ms, what = "passage_time", to = 1),
               class = "hypernets_bad_input")
  expect_error(hg_get(ms, decreasing = NA), "`decreasing` must be")
})


test_that("plot() draws every state for each measure, and the passage heatmap", {
  ms <- hg_markov_stability(human_long, actor = "session_id", action = "code",
                            time = "order_in_session")
  p <- plot(ms)
  expect_s3_class(p, "ggplot")
  # every state has a bar in every one of the six panels
  expect_identical(nrow(p$data), 6L * nrow(hg_get(ms)))
  expect_setequal(levels(p$data$state), hg_get(ms)$state)
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_no_error(plot(ms, what = "passage_time"))
  expect_error(plot(ms, what = "nope"))
})
