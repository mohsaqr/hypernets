testthat::skip_on_cran()

# ---- build_hon(by =): one BuildHON model per unit (net_hon_group) --------

.hg_long <- function() {
  subset(human_long, project %in% c("Project_10", "Project_11", "Project_19",
                                    "Project_12", "Project_15"))
}

.hg_fit <- function(...) {
  withCallingHandlers(
    build_hon(.hg_long(), action = "code", actor = "session_id",
              time = "timestamp", by = "project", max_order = 2L, ...),
    hypernets_insufficient_support = function(w) invokeRestart("muffleWarning"))
}

test_that("build_hon default output is unchanged (frozen fixture)", {
  fx <- readRDS(test_path("fixtures", "build_hon_frozen.rds"))
  seqs <- list(c("A", "B", "C", "D", "A"), c("A", "B", "D", "C", "A"),
               c("B", "C", "D", "A", "B", "C"))
  df <- data.frame(T1 = c("A", "A", "B"), T2 = c("B", "B", "C"),
                   T3 = c("C", "D", "A"), T4 = c("D", "C", NA))
  long <- subset(human_long,
                 project %in% c("Project_10", "Project_11", "Project_19"))
  ho <- c(rep(list(c("A", "B", "C")), 6), rep(list(c("D", "B", "E")), 6))
  expect_identical(build_hon(seqs), fx$list_default)
  expect_identical(build_hon(seqs, max_order = 2L), fx$list_order2)
  expect_identical(build_hon(seqs, method = "hon", max_order = 3L),
                   fx$list_hon)
  expect_identical(build_hon(df, max_order = 2L), fx$wide_default)
  expect_identical(build_hon(long, action = "code", actor = "session_id",
                             time = "timestamp", max_order = 2L),
                   fx$long_default)
  expect_identical(build_hon(long, action = "code", actor = "session_id",
                             time = "timestamp", max_order = 2L,
                             collapse_repeats = TRUE, min_freq = 2L),
                   fx$long_collapse)
  expect_identical(build_hon(ho, max_order = 3L), fx$list_higher)
})

test_that("every per-unit model is identical() to build_hon on that unit alone", {
  long <- .hg_long()
  fits <- .hg_fit()
  expect_s3_class(fits, "net_hon_group")
  fitted_units <- as.data.frame(fits, what = "actors")
  fitted_units <- fitted_units$actor[fitted_units$fitted]
  expect_gte(length(fitted_units), 3L)
  lapply(fitted_units, function(u) {
    alone <- build_hon(long[long$project == u, , drop = FALSE],
                       action = "code", actor = "session_id",
                       time = "timestamp", max_order = 2L)
    expect_identical(fits$models[[u]], alone, info = u)
  })
})

test_that("per-unit identity holds for wide and list input and for hon/min_freq", {
  set.seed(11)
  states <- c("a", "b", "c")
  seqs <- replicate(24, sample(states, 12, replace = TRUE), simplify = FALSE)
  units <- rep(c("u1", "u2", "u3"), each = 8)
  fits <- build_hon(seqs, by = units, max_order = 3L, min_freq = 2L,
                    method = "hon")
  lapply(c("u1", "u2", "u3"), function(u) {
    expect_identical(fits$models[[u]],
                     build_hon(seqs[units == u], max_order = 3L, min_freq = 2L,
                               method = "hon"), info = u)
  })

  wide <- as.data.frame(do.call(rbind, seqs), stringsAsFactors = FALSE)
  wide$student <- units
  fits_w <- build_hon(wide, by = "student", max_order = 2L,
                      collapse_repeats = TRUE)
  lapply(c("u1", "u2", "u3"), function(u) {
    rows <- wide[wide$student == u, setdiff(names(wide), "student")]
    expect_identical(fits_w$models[[u]],
                     build_hon(rows, max_order = 2L, collapse_repeats = TRUE),
                     info = u)
  })
})

test_that("support guard fires by class, names the units and excludes them", {
  w <- NULL
  fits <- withCallingHandlers(
    build_hon(.hg_long(), action = "code", actor = "session_id",
              time = "timestamp", by = "project", max_order = 2L),
    hypernets_insufficient_support = function(cnd) {
      w <<- cnd
      invokeRestart("muffleWarning")
    })
  expect_s3_class(w, "hypernets_insufficient_support")
  actors <- as.data.frame(fits, what = "actors")
  below <- actors$actor[!actors$sufficient]
  expect_setequal(w$actors, below)
  expect_true(all(c("Project_12", "Project_15") %in% below))
  expect_false(any(below %in% names(fits$models)))
  expect_false(any(actors$fitted[!actors$sufficient]))
  # the default floor is S(S-1) on the pooled state space
  n_states <- length(unique(.hg_long()$code))
  expect_identical(unique(actors$min_transitions),
                   as.integer(n_states * (n_states - 1L)))
  expect_warning(
    build_hon(.hg_long(), action = "code", actor = "session_id",
              time = "timestamp", by = "project", max_order = 2L),
    class = "hypernets_insufficient_support")
})

test_that("keep_sparse fits units below the floor and flags them", {
  expect_warning(
    fits <- build_hon(.hg_long(), action = "code", actor = "session_id",
                      time = "timestamp", by = "project", max_order = 2L,
                      keep_sparse = TRUE),
    class = "hypernets_insufficient_support")
  actors <- as.data.frame(fits, what = "actors")
  sparse <- actors[!actors$sufficient, ]
  expect_true(all(sparse$fitted[sparse$n_transitions >= 1L]))
  expect_true("Project_12" %in% names(fits$models))
})

test_that("a stated floor is honoured and silences the guard when met", {
  expect_no_warning(
    fits <- build_hon(.hg_long(), action = "code", actor = "session_id",
                      time = "timestamp", by = "project", max_order = 2L,
                      min_transitions = 1L))
  actors <- as.data.frame(fits, what = "actors")
  expect_true(all(actors$min_transitions == 1L))
  expect_true(all(actors$fitted))
})

test_that("no fittable unit is an error of the support class", {
  seqs <- list(c("a", "b"), c("b", "c"), c("c", "a"))
  expect_error(build_hon(seqs, by = c("x", "y", "z")),
               class = "hypernets_insufficient_support")
})

test_that("as.data.frame edges: one row per unit per rule; probabilities sum to 1", {
  fits <- .hg_fit()
  edges <- as.data.frame(fits)
  expect_identical(class(edges), "data.frame")
  expect_identical(attr(edges, "row.names"), seq_len(nrow(edges)))
  expect_named(edges, c("actor", "path", "from", "to", "count",
                        "probability", "from_order", "to_order"))
  expect_false(anyDuplicated(edges[c("actor", "path")]) > 0L)
  n_rules <- vapply(fits$models, \(m) nrow(m$ho_edges), integer(1L))
  expect_identical(as.vector(table(edges$actor)[names(n_rules)]),
                   unname(n_rules))
  sums <- tapply(edges$probability, paste(edges$actor, edges$from), sum)
  expect_true(all(abs(sums - 1) < sqrt(.Machine$double.eps)))

  ho <- as.data.frame(fits, order_min = 2L)
  expect_true(all(ho$from_order >= 2L))
  one <- as.data.frame(fits, actor = "Project_19")
  expect_identical(unique(one$actor), "Project_19")
  top <- as.data.frame(fits, sort_by = "count", top = 5L)
  expect_identical(nrow(top), 5L)
  expect_false(is.unsorted(rev(top$count)))
})

test_that("as.data.frame actors: one row per unit with the support columns", {
  fits <- .hg_fit()
  actors <- as.data.frame(fits, what = "actors")
  expect_identical(class(actors), "data.frame")
  expect_identical(actors$actor, sort(unique(.hg_long()$project)))
  expect_true(all(c("actor", "n_sequences", "n_transitions", "n_states",
                    "min_transitions", "sufficient", "fitted", "n_nodes",
                    "n_edges", "max_order_observed") %in% names(actors)))
  expect_identical(is.na(actors$n_nodes), !actors$fitted)
  # support counts match an independent count on the raw events
  long <- .hg_long()
  p10 <- long[long$project == "Project_10", ]
  per_session <- table(p10$session_id)
  expect_identical(actors$n_transitions[actors$actor == "Project_10"],
                   as.integer(sum(per_session[per_session >= 2] - 1L)))
  expect_identical(summary(fits), actors)
})

test_that("print and plot views run and delegate to cograph", {
  fits <- .hg_fit()
  expect_output(print(fits), "one BuildHON model per unit")
  p <- plot(fits, top = 10L)
  expect_s3_class(p, "ggplot")
  expect_s3_class(plot(fits, what = "support"), "ggplot")
  pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_invisible(plot(fits, what = "network", actor = "Project_19"))
})

test_that("bad input is refused by class", {
  long <- .hg_long()
  expect_error(build_hon(long, action = "code", actor = "session_id",
                         by = "nope"), class = "hypernets_bad_input")
  expect_error(build_hon(list(c("a", "b")), by = c("x", "y")),
               class = "hypernets_bad_input")
  expect_error(build_hon(list(c("a", "b"), c("b", "a")), by = c("x", NA)),
               class = "hypernets_bad_input")
  expect_error(build_hon(long, action = "code", actor = "session_id",
                         by = "project", min_transitions = 0),
               class = "hypernets_bad_input")
  expect_error(build_hon(list(c("a", "b", "a")), min_transitions = 5L),
               class = "hypernets_bad_input")
  expect_error(build_hon(list(c("a", "b", "a")), keep_sparse = TRUE),
               class = "hypernets_bad_input")
  fits <- .hg_fit()
  expect_error(as.data.frame(fits, actor = "Project_99"),
               class = "hypernets_bad_input")
  expect_error(plot(fits, what = "network"), class = "hypernets_bad_input")
  expect_error(plot(fits, what = "bogus"), class = "hypernets_bad_input")
})

test_that("unit order is invariant to row permutation of the input", {
  long <- .hg_long()
  set.seed(3)
  shuffled_units <- sample(unique(long$project))
  reordered <- do.call(rbind, split(long, long$project)[shuffled_units])
  a <- .hg_fit()
  b <- withCallingHandlers(
    build_hon(reordered, action = "code", actor = "session_id",
              time = "timestamp", by = "project", max_order = 2L),
    hypernets_insufficient_support = function(w) invokeRestart("muffleWarning"))
  expect_identical(names(a$models), names(b$models))
  expect_identical(as.data.frame(a), as.data.frame(b))
})
