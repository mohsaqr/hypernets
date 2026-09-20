testthat::skip_on_cran()

# The memory verbs read a data.frame as WIDE: one trajectory per row. Handed a
# LONG event table they used to promote actor ids and timestamps to states and
# return a confidently wrong model with no error and no warning. These tests
# pin the guard that now refuses, and -- just as importantly -- pin that the
# guard does not fire on the wide frames the verbs are supposed to accept.

.lf_long <- function() {
  data.frame(
    actor  = rep(c("s1", "s2"), each = 4L),
    time   = rep(1:4, 2L),
    action = c("A", "B", "A", "C", "B", "B", "C", "A"),
    stringsAsFactors = FALSE
  )
}

.lf_wide <- function() {
  data.frame(
    t1 = c("A", "B", "A"), t2 = c("B", "B", "C"),
    t3 = c("A", "C", "A"), t4 = c("C", "A", "B"),
    stringsAsFactors = FALSE
  )
}

test_that("a long table is refused instead of silently read as wide", {
  long <- .lf_long()
  expect_error(build_hon(long), class = "hypernets_long_format")
  expect_error(build_hon(long), class = "hypernets_bad_input")
  # every memory verb that parses a wide frame is covered, including the ones
  # that cannot take the long-format arguments at all
  expect_error(build_mogen(long), class = "hypernets_long_format")
  expect_error(build_hypa(long), class = "hypernets_long_format")
  expect_error(markov_order_test(long), class = "hypernets_long_format")
  expect_error(bootstrap_hon(long, n_boot = 2L), class = "hypernets_long_format")
})

test_that("the message names the columns and the remedy for that verb", {
  long <- .lf_long()
  long_msg <- conditionMessage(tryCatch(build_hon(long), error = function(e) e))
  expect_match(long_msg, "`action`")
  expect_match(long_msg, "`actor`")
  expect_match(long_msg, "`time`")
  # a verb that accepts the long arguments is told to pass them
  expect_match(long_msg, "action = ", fixed = TRUE)
  # one that does not is told to split instead, not sent to a dead end
  wide_only <- conditionMessage(
    tryCatch(markov_order_test(long), error = function(e) e)
  )
  expect_match(wide_only, "wide sequences only")
  expect_false(grepl("pass the long-format arguments", wide_only))
})

test_that("the long-format arguments still work and give the right model", {
  long <- .lf_long()
  fit <- build_hon(long, action = "action", actor = "actor", time = "time")
  expect_s3_class(fit, "net_hon")
  # three states and two trajectories -- not nine states and eight, which is
  # what the silent wide reading produced
  expect_identical(sort(rownames(fit$weights)), c("A", "B", "C"))
  expect_identical(fit$n_trajectories, 2L)
  # identical to building from the hand-split list
  split_fit <- build_hon(split(long$action, long$actor))
  expect_identical(fit$weights, split_fit$weights)
})

test_that("the guard does not fire on the wide frames the verbs accept", {
  wide <- .lf_wide()
  expect_s3_class(build_hon(wide), "net_hon")
  expect_s3_class(build_mogen(wide), "net_mogen")
  # a list of sequences is untouched by the guard
  expect_s3_class(build_hon(list(c("A", "B", "C"), c("B", "C", "A"))), "net_hon")
  # one canonical name alone is not enough: a wide frame may legitimately
  # have a column called "time" without being a long table
  one_name <- .lf_wide()
  names(one_name)[1L] <- "time"
  expect_s3_class(build_hon(one_name), "net_hon")
})

test_that("detection is on names, case-insensitively, and needs two of three", {
  expect_null(hypernets:::.hon_guard_long_format(list(c("A", "B"))))
  expect_null(hypernets:::.hon_guard_long_format(.lf_wide()))
  shouty <- .lf_long()
  names(shouty) <- toupper(names(shouty))
  expect_error(hypernets:::.hon_guard_long_format(shouty),
               class = "hypernets_long_format")
  two_of_three <- .lf_long()[, c("actor", "action")]
  expect_error(hypernets:::.hon_guard_long_format(two_of_three),
               class = "hypernets_long_format")
})

test_that("hg_sequences output is refused bare and accepted with arguments", {
  # the verb that makes long tables easy to obtain is the reason the guard
  # matters, so the two are tested together
  corpus <- data.frame(
    student = rep(c("s1", "s2"), each = 4L),
    turn = rep(1:4, 2L),
    # "graph" is a deliberate bridge word: doc-mode corpora are connected only
    # if shared words chain every document, and hg_cluster() refuses otherwise
    text = c("slope rate graph", "the data graph", "slope rate graph",
             "model residual graph", "the data graph", "the data graph",
             "model residual graph", "slope rate graph"),
    stringsAsFactors = FALSE
  )
  hg <- text_hypergraph(corpus, column = "text", nodes = "doc", min_count = 1L)
  topics <- hg_cluster(hg, k = 2L)
  seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
  expect_identical(names(seqs), c("actor", "time", "action"))
  expect_error(build_hon(seqs), class = "hypernets_long_format")
  fit <- build_hon(seqs, action = "action", actor = "actor", time = "time")
  expect_identical(fit$n_trajectories, 2L)
})
