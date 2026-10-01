testthat::skip_on_cran()

# ---- Memory-network wrappers: identity with the Nestimate estimators ----
#
# hon(), mogen(), hypa() on sequences, markov_order(),
# memory() and hg_markov_stability() compute nothing: they coerce
# the input and return the estimator's object. Every test below asserts
# identical() against the direct Nestimate call, on every input form the
# wrappers accept: a list, a wide frame, a long event table, a netobject and
# a tna model.

.mw_fixture <- function() {
  set.seed(3)
  seqs <- lapply(seq_len(12L), function(i) {
    sample(c("a", "b", "c", "d"), sample(6:10, 1L), replace = TRUE)
  })
  names(seqs) <- sprintf("u%02d", seq_len(12L))
  long <- do.call(rbind, lapply(names(seqs), function(n) {
    data.frame(who = n, when = seq_along(seqs[[n]]) * 10, what = seqs[[n]],
               stringsAsFactors = FALSE)
  }))
  long <- long[sample(nrow(long)), ]          # row order must not matter
  rownames(long) <- NULL
  wide <- .ho_wide_sequences(seqs)
  net <- Nestimate::build_network(long, method = "relative", actor = "who",
                                  action = "what", time = "when",
                                  time_threshold = Inf)
  list(seqs = seqs, long = long, wide = wide, net = net)
}

test_that("the long route builds the sequences build_network() builds", {
  # the same call, argument for argument, in both packages
  same <- function(data, ...) {
    net <- suppressMessages(
      Nestimate::build_network(data, method = "relative", ...))
    ours <- .ho_sequence_input(data, ...)
    expect_identical(ours, net$data)
    nrow(ours)
  }
  sessions <- human_long
  n <- c(
    position = same(sessions, action = "code", actor = "session_id",
                    time = "order_in_session"),
    timestamp = same(sessions, action = "code", actor = "session_id",
                     time = "timestamp"),
    no_gaps = same(sessions, action = "code", actor = "session_id",
                   time = "timestamp", time_threshold = FALSE),
    project = same(sessions, action = "code", actor = "project",
                   time = "order_in_session"),
    project_session = same(sessions, action = "code", actor = "project",
                           session = "session_id",
                           time = "order_in_session"),
    two_actors = same(sessions, action = "code",
                      actor = c("project", "session_id"),
                      time = "order_in_session")
  )
  # a gap of more than 900 seconds starts a new sequence; the session column
  # is found by its name when `session` is NULL
  expect_identical(unname(n), c(429L, 526L, 429L, 429L, 429L, 429L))
  # an `action` column is found by its name too
  renamed <- human_long
  names(renamed)[names(renamed) == "code"] <- "action"
  expect_identical(.ho_sequence_input(renamed, actor = "session_id",
                                      time = "order_in_session"),
                   .ho_sequence_input(human_long, action = "code",
                                      actor = "session_id",
                                      time = "order_in_session"))
  fx <- .mw_fixture()
  expect_identical(same(fx$long, action = "what", actor = "who",
                        time = "when"), 12L)
})

test_that("session = FALSE switches the session detection off", {
  by_project <- .ho_sequence_input(human_long, action = "code",
                                   actor = "project",
                                   time = "order_in_session",
                                   session = FALSE)
  # without the session column there is nothing to detect: the same call
  without <- human_long[setdiff(names(human_long), "session_id")]
  net <- Nestimate::build_network(without, method = "relative",
                                  action = "code", actor = "project",
                                  time = "order_in_session")
  expect_identical(by_project, net$data)
  expect_identical(nrow(by_project), length(unique(human_long$project)))
})

test_that("missing actors and sessions are refused, by class", {
  fx <- .mw_fixture()
  long <- fx$long
  long$who[3L] <- NA
  expect_error(hon(long, action = "what", actor = "who", time = "when"),
               class = "hypernets_bad_input")
  expect_error(hon(fx$long, action = "what", actor = "who", time = "when",
                   time_threshold = -1))
})

test_that("hon() is the HON builder on every input form", {
  fx <- .mw_fixture()
  ref <- Nestimate::build_hon(fx$seqs, max_order = 2L)
  expect_identical(.ho_unresult(hon(fx$seqs, max_order = 2L)), ref)
  expect_identical(.ho_unresult(hon(fx$wide, max_order = 2L)),
                   Nestimate::build_hon(fx$wide, max_order = 2L))
  expect_identical(.ho_unresult(hon(fx$long, action = "what", actor = "who", time = "when",
                       max_order = 2L)),
                   Nestimate::build_hon(fx$net, max_order = 2L))
  expect_identical(.ho_unresult(hon(fx$long, action = "what", actor = "who", time = "when",
                       max_order = 2L)), ref)
  expect_identical(.ho_unresult(hon(fx$net, max_order = 2L)),
                   Nestimate::build_hon(fx$net, max_order = 2L))
  expect_identical(.ho_unresult(hon(fx$seqs, max_order = 3L, min_freq = 2L,
                       collapse_repeats = TRUE, method = "hon")),
                   Nestimate::build_hon(fx$seqs, max_order = 3L, min_freq = 2L,
                                        collapse_repeats = TRUE,
                                        method = "hon"))
})

# A tna model as the tna package stores it: integer-coded sequences with
# their labels (built here so the test needs no extra package).
.mw_tna <- function(wide) {
  labels <- sort(unique(stats::na.omit(unlist(wide, use.names = FALSE))))
  codes <- vapply(wide, function(col) match(col, labels), integer(nrow(wide)))
  codes <- matrix(codes, nrow(wide), dimnames = list(NULL, names(wide)))
  attr(codes, "labels") <- labels
  structure(list(data = codes, labels = labels), class = "tna")
}

test_that("hon() reads a tna model", {
  fx <- .mw_fixture()
  tm <- .mw_tna(fx$wide)
  expect_identical(.ho_unresult(hon(tm, max_order = 2L)),
                   Nestimate::build_hon(tm, max_order = 2L))
  expect_identical(hon(tm, max_order = 2L), hon(fx$seqs, max_order = 2L))
})

test_that("honem() is the HONEM builder", {
  fx <- .mw_fixture()
  h <- hon(fx$seqs, max_order = 2L)
  expect_identical(.ho_unresult(honem(h, dim = 3L, max_power = 5L)),
                   Nestimate::build_honem(h, dim = 3L, max_power = 5L))
})

test_that("mogen() is the MOGen builder on every input form", {
  fx <- .mw_fixture()
  expect_identical(.ho_unresult(mogen(fx$seqs, max_order = 2L)),
                   Nestimate::build_mogen(fx$seqs, max_order = 2L))
  expect_identical(.ho_unresult(mogen(fx$wide, max_order = 2L, criterion = "bic")),
                   Nestimate::build_mogen(fx$wide, max_order = 2L,
                                          criterion = "bic"))
  expect_identical(.ho_unresult(mogen(fx$long, action = "what", actor = "who",
                         time = "when", max_order = 2L)),
                   Nestimate::build_mogen(fx$net, max_order = 2L))
  expect_identical(.ho_unresult(mogen(fx$net, max_order = 2L)),
                   Nestimate::build_mogen(fx$net, max_order = 2L))
})

test_that("hypa() on sequences is the HYPA builder on every input form", {
  fx <- .mw_fixture()
  expect_identical(.unview(hypa(fx$seqs, min_count = 1L)),
                   Nestimate::build_hypa(fx$seqs, min_count = 1L))
  expect_identical(.unview(hypa(fx$wide, order = c(1L, 2L), min_count = 2L,
                           p_adjust = "holm")),
                   Nestimate::build_hypa(fx$wide, order = c(1L, 2L),
                                         min_count = 2L, p_adjust = "holm"))
  expect_identical(.unview(hypa(fx$long, action = "what", actor = "who",
                           time = "when", min_count = 1L)),
                   Nestimate::build_hypa(fx$net, min_count = 1L))
  expect_identical(.unview(hypa(fx$net, min_count = 1L)),
                   Nestimate::build_hypa(fx$net, min_count = 1L))
})

test_that("markov_order() is the Markov order test on every input form", {
  fx <- .mw_fixture()
  ref <- Nestimate::markov_order_test(fx$seqs, max_order = 2L, n_perm = 20L,
                                      seed = 1L)
  expect_identical(.ho_unresult(markov_order(fx$seqs, max_order = 2L, n_perm = 20L,
                                   seed = 1L)), ref)
  expect_identical(.ho_unresult(markov_order(fx$wide, max_order = 2L, n_perm = 20L,
                                   seed = 1L)),
                   Nestimate::markov_order_test(fx$wide, max_order = 2L,
                                                n_perm = 20L, seed = 1L))
  expect_identical(.ho_unresult(markov_order(fx$long, action = "what", actor = "who",
                                   time = "when", max_order = 2L,
                                   n_perm = 20L, seed = 1L)),
                   Nestimate::markov_order_test(fx$net, max_order = 2L,
                                                n_perm = 20L, seed = 1L))
  expect_identical(.ho_unresult(markov_order(fx$net, max_order = 2L, n_perm = 20L,
                                   seed = 1L)),
                   Nestimate::markov_order_test(fx$net, max_order = 2L,
                                                n_perm = 20L, seed = 1L))
})

test_that("markov_order() decodes a tna model to its sequences", {
  fx <- .mw_fixture()
  expect_identical(markov_order(.mw_tna(fx$wide), max_order = 2L,
                                n_perm = 20L, seed = 1L),
                   markov_order(fx$seqs, max_order = 2L, n_perm = 20L,
                                seed = 1L))
})

test_that("memory() is the path-dependence estimator on every input form", {
  fx <- .mw_fixture()
  ref <- Nestimate::path_dependence(fx$wide, order = 2L, min_count = 1L)
  expect_identical(.ho_unresult(memory(fx$wide, order = 2L, min_count = 1L)),
                   ref)
  # the estimator reads only frames: a list is padded into the same frame
  expect_identical(.ho_unresult(memory(fx$seqs, order = 2L, min_count = 1L)),
                   ref)
  expect_identical(.ho_unresult(memory(fx$long, action = "what", actor = "who",
                                      time = "when", min_count = 1L)),
                   Nestimate::path_dependence(fx$net, min_count = 1L))
  expect_identical(.ho_unresult(memory(fx$net, min_count = 1L, base = exp(1))),
                   Nestimate::path_dependence(fx$net, min_count = 1L,
                                              base = exp(1)))
})

test_that("hg_markov_stability() is the stability estimator on every input form", {
  fx <- .mw_fixture()
  P <- matrix(c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), 3, 3,
              byrow = TRUE, dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_identical(.ho_unresult(hg_markov_stability(P)), Nestimate::markov_stability(P))
  expect_identical(.ho_unresult(hg_markov_stability(fx$wide)),
                   Nestimate::markov_stability(fx$wide))
  expect_identical(.ho_unresult(hg_markov_stability(fx$seqs)),
                   Nestimate::markov_stability(fx$wide))
  expect_identical(.ho_unresult(hg_markov_stability(fx$long, action = "what", actor = "who",
                                       time = "when")),
                   Nestimate::markov_stability(fx$net))
  expect_identical(.ho_unresult(hg_markov_stability(fx$net)),
                   Nestimate::markov_stability(fx$net))
})

test_that("session splits an actor's events into separate sequences", {
  fx <- .mw_fixture()
  long <- fx$long
  long$sess <- ifelse(long$when > 50, "s2", "s1")
  net <- Nestimate::build_network(long, method = "relative", actor = "who",
                                  session = "sess", action = "what",
                                  time = "when", time_threshold = Inf)
  expect_identical(.ho_unresult(hon(long, action = "what", actor = "who", session = "sess",
                       time = "when", max_order = 2L)),
                   Nestimate::build_hon(net, max_order = 2L))
  seqs <- .ho_sequence_input(long, action = "what", actor = "who",
                             session = "sess", time = "when")
  # every actor with events on both sides of the split gives two sequences
  expect_identical(nrow(seqs), 24L)
  expect_identical(sum(!is.na(as.matrix(seqs))), nrow(long))
})

test_that("a long table passed without action is refused, by class", {
  fx <- .mw_fixture()
  canonical <- data.frame(who = c("p", "p", "q"), state = c("a", "b", "a"),
                          timestamp = 1:3)
  lapply(list(
    function(d) hon(d),
    function(d) mogen(d),
    function(d) hypa(d),
    function(d) markov_order(d),
    function(d) memory(d),
    function(d) hg_markov_stability(d),
    function(d) hg_bootstrap(d, n_boot = 2L)
  ), function(f) {
    expect_error(f(canonical), class = "hypernets_long_format")
    expect_error(f(human_long), class = "hypernets_bad_input")
  })
  # a wide frame of states passes the guard
  expect_s3_class(hon(fx$wide, max_order = 2L), "net_hon")
})

test_that("long-format column arguments are validated", {
  fx <- .mw_fixture()
  expect_error(hon(fx$seqs, actor = "who"), class = "hypernets_bad_input")
  expect_error(hon(fx$wide, time = "when"), class = "hypernets_bad_input")
  expect_error(hon(fx$long, action = "nope"), class = "hypernets_bad_input")
  expect_error(hon(fx$long, action = "what", actor = "nope"),
               class = "hypernets_bad_input")
  expect_error(hon(fx$seqs, action = "what"), class = "hypernets_bad_input")
})

test_that("the wrappers are invariant to the row order of a long table", {
  fx <- .mw_fixture()
  shuffled <- fx$long[rev(seq_len(nrow(fx$long))), ]
  expect_identical(
    hon(shuffled, action = "what", actor = "who", time = "when",
        max_order = 2L),
    hon(fx$long, action = "what", actor = "who", time = "when",
        max_order = 2L))
})

test_that("hg_markov_stability() refuses a chain that is not irreducible", {
  absorbing <- matrix(c(0.5, 0.4, 0.1,
                        0.3, 0.6, 0.1,
                        0.0, 0.0, 1.0), 3, byrow = TRUE,
                      dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_error(hg_markov_stability(absorbing),
               class = "hypernets_not_ergodic")
  err <- tryCatch(hg_markov_stability(absorbing), error = identity)
  expect_match(conditionMessage(err), "transient: A, B")
  expect_match(conditionMessage(err), "absorbing: C")
  # an irreducible chain is unchanged: the estimator's own result
  regular <- matrix(c(0.5, 0.3, 0.2, 0.2, 0.5, 0.3, 0.1, 0.4, 0.5), 3,
                    byrow = TRUE,
                    dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  expect_identical(.ho_unresult(hg_markov_stability(regular)),
                   Nestimate::markov_stability(regular))
})

test_that("a state with no outgoing transition is refused by class", {
  dead_end <- matrix(c(0.5, 0.5, 0,
                       0.3, 0.3, 0.4,
                       0, 0, 0), 3, byrow = TRUE,
                     dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  err <- tryCatch(hg_markov_stability(dead_end), error = identity)
  expect_s3_class(err, "hypernets_not_ergodic")
  expect_match(conditionMessage(err), "no outgoing transition from: C")
})
