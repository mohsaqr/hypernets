testthat::skip_on_cran()

# A data.frame is read as a long event table when `action` is given or a
# column is named `action` (with `time` and `session` / `session_id` found by
# name too, as build_network() finds them). A long table whose columns carry
# other names (code, user, timestamp, ...) and that is passed without
# `action =` would be read as wide, promoting actor ids and timestamps to
# states. These tests pin the guard that refuses it, and pin that it does
# not fire on the wide frames the verbs accept.

.lf_long <- function() {
  data.frame(
    user      = rep(c("s1", "s2"), each = 4L),
    timestamp = rep(1:4, 2L),
    code      = c("A", "B", "A", "C", "B", "B", "C", "A"),
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
  expect_error(hg_bootstrap(long, n_boot = 2L), class = "hypernets_long_format")
  expect_error(hg_bootstrap(long, n_boot = 2L), class = "hypernets_bad_input")
  expect_error(hon(long), class = "hypernets_long_format")
})

test_that("the message names the columns and the remedy", {
  long <- .lf_long()
  long_msg <- conditionMessage(tryCatch(hg_bootstrap(long, n_boot = 2L),
                                        error = function(e) e))
  expect_match(long_msg, "`code`")
  expect_match(long_msg, "`user`")
  expect_match(long_msg, "`timestamp`")
  expect_match(long_msg, "action = ", fixed = TRUE)
})

test_that("the long-format arguments still work and give the right model", {
  long <- .lf_long()
  fit <- hg_bootstrap(long, action = "code", actor = "user",
                       time = "timestamp", n_boot = 5L, max_order = 2L,
                       seed = 1L)
  expect_s3_class(fit, "net_hon_boot")
  # two trajectories -- not eight, which is what a wide reading produced
  expect_identical(fit$n_trajectories, 2L)
  # identical to the hand-split list
  split_fit <- hg_bootstrap(split(long$code, long$user), n_boot = 5L,
                             max_order = 2L, seed = 1L)
  expect_identical(hg_get(fit), hg_get(split_fit))
})

test_that("the guard does not fire on the wide frames the verbs accept", {
  wide <- .lf_wide()
  expect_s3_class(hg_bootstrap(wide, n_boot = 2L, seed = 1L), "net_hon_boot")
  # a list of sequences is untouched by the guard
  expect_s3_class(hg_bootstrap(list(c("A", "B", "C"), c("B", "C", "A")),
                                n_boot = 2L, seed = 1L), "net_hon_boot")
  # one canonical name alone is not enough: a wide frame may legitimately
  # have a column called "time" without being a long table
  one_name <- .lf_wide()
  names(one_name)[1L] <- "time"
  expect_s3_class(hg_bootstrap(one_name, n_boot = 2L, seed = 1L),
                  "net_hon_boot")
})

test_that("detection is on names, case-insensitively, and needs two of three", {
  expect_null(hypernets:::.hon_guard_long_format(list(c("A", "B"))))
  expect_null(hypernets:::.hon_guard_long_format(.lf_wide()))
  shouty <- .lf_long()
  names(shouty) <- toupper(names(shouty))
  expect_error(hypernets:::.hon_guard_long_format(shouty),
               class = "hypernets_long_format")
  two_of_three <- .lf_long()[, c("user", "code")]
  expect_error(hypernets:::.hon_guard_long_format(two_of_three),
               class = "hypernets_long_format")
})

test_that("hg_sequences output is one sequence bare and per actor with arguments", {
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
  # its columns are named action and time, so they are found by name; with
  # no actor the events form one sequence, which bootstrap cannot resample
  expect_message(expect_error(hg_bootstrap(seqs, n_boot = 2L)),
                 class = "hypernets_single_sequence")
  fit <- hg_bootstrap(seqs, action = "action", actor = "actor", time = "time",
                       n_boot = 2L, seed = 1L)
  expect_identical(fit$n_trajectories, 2L)
})
