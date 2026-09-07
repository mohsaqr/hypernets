# temporal_hypergraph() on sequence tables: a session is a hyperedge, a
# state's membership is a contact at its position, and a hyperedge spans
# the hull of its memberships. Constant-time data must be untouched.

wide <- data.frame(
  T1 = c("a", "b", "c"), T2 = c("b", "b", "c"),
  T3 = c("a", NA, "d"), T4 = c("c", NA, NA),
  stringsAsFactors = FALSE
)

test_that("a wide sequence table is read as one hyperedge per session", {
  th <- temporal_hypergraph(wide)
  expect_s3_class(th, "net_temporal_hypergraph")
  expect_true(th$params$membership_times)
  expect_identical(th$format, "interval")
  expect_identical(th$time_unit, "step")
  expect_identical(th$times, as.numeric(1:4))
  edges <- as.data.frame(th, what = "edges")
  expect_identical(edges$edge, paste0("sequence_", 1:3))
  expect_identical(edges$start, c(1, 1, 1))
  expect_identical(edges$end, c(4, 2, 3))
  memberships <- as.data.frame(th)
  expect_identical(nrow(memberships), 9L)
  expect_identical(memberships$start, memberships$end)
  expect_identical(memberships$start[memberships$edge == "sequence_1"], as.numeric(1:4))
})

test_that("a list of vectors and a long log with per-row times give the same memberships", {
  th <- temporal_hypergraph(wide)
  as_list <- temporal_hypergraph(list(
    sequence_1 = c("a", "b", "a", "c"), sequence_2 = c("b", "b"),
    sequence_3 = c("c", "c", "d")
  ))
  expect_identical(as.data.frame(as_list), as.data.frame(th))
  long <- data.frame(
    code = c("a", "b", "a", "c", "b", "b", "c", "c", "d"),
    session = rep(paste0("sequence_", 1:3), c(4, 2, 3)),
    step = c(1:4, 1:2, 1:3),
    stringsAsFactors = FALSE
  )
  as_long <- temporal_hypergraph(long, actor = "code", group = "session", time = "step")
  expect_identical(as.data.frame(as_long), as.data.frame(th))
  expect_true(as_long$params$membership_times)
})

test_that("snapshots keep the memberships present in the window", {
  th <- temporal_hypergraph(wide)
  seen_by_two <- as.data.frame(hypergraph_snapshot(th, at = 2, mode = "cumulative"))
  expect_identical(seen_by_two$states, c("a, b", "b", "c"))
  steps_two_three <- as.data.frame(hypergraph_snapshot(th, at = 2, window = 2))
  expect_identical(steps_two_three$states, c("a, b", "b", "c, d"))
  at_three <- as.data.frame(hypergraph_snapshot(th, at = 3))
  expect_identical(at_three$hyperedge, c("sequence_1", "sequence_3"))
  expect_identical(at_three$states, c("a", "d"))
})

test_that("a window snapshot of every session equals the sequence windows", {
  set.seed(7)
  sequences <- replicate(20, sample(letters[1:5], sample(3:9, 1), replace = TRUE),
                         simplify = FALSE)
  th <- temporal_hypergraph(sequences)
  window <- 3
  snaps <- hypergraph_snapshots(th, start = 1, end = 9, step = 1, window = window)
  # every full window of every sequence, as the set of its states
  expected <- unlist(lapply(sequences, function(s) {
    starts <- seq_len(max(length(s) - window + 1L, 0L))
    vapply(starts, function(i) paste(sort(unique(s[i:(i + window - 1L)])), collapse = ", "),
           character(1L))
  }))
  observed <- unlist(lapply(seq_len(9 - window + 1L), function(t) {
    hg <- snaps[[as.character(t)]]
    tab <- as.data.frame(hg)
    # sessions whose spell covers the whole window are the full windows
    full <- vapply(tab$hyperedge, function(e) th$edge_data$end[th$edge_data$edge == e] >= t + window - 1,
                   logical(1L))
    tab$states[full]
  }))
  expect_identical(sort(observed), sort(expected))
})

test_that("hg_growth counts memberships present, not whole sessions", {
  th <- temporal_hypergraph(wide)
  active <- hg_growth(th)
  expect_identical(active$n_memberships, c(3L, 3L, 2L, 1L))
  expect_identical(active$n_nodes, c(3L, 2L, 2L, 1L))
  expect_identical(active$n_edges, c(3L, 3L, 2L, 1L))
  expect_identical(active$n_edges_distinct, c(3L, 2L, 2L, 1L))
  cumulative <- hg_growth(th, mode = "cumulative")
  expect_identical(cumulative$n_memberships, c(3L, 6L, 8L, 9L))
  expect_identical(cumulative$n_nodes, c(3L, 3L, 4L, 4L))
})

test_that("constant-time data are untouched and snapshots do not move", {
  seats <- data.frame(
    case = rep(c("A", "B"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
    start = rep(c(1, 2), each = 3), end = rep(c(3, 2), each = 3)
  )
  th <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                            start = "start", end = "end")
  expect_false(th$params$membership_times)
  expect_identical(th$format, "interval")
  expect_identical(th$times, c(1, 2, 3))
  expect_identical(hypergraph_snapshot(th, 2)$n_hyperedges, 2L)
  expect_identical(hypergraph_snapshot(th, 3)$n_hyperedges, 1L)
})

test_that("a relational table without ends is not mistaken for a sequence table", {
  no_ends <- data.frame(person = c("a", "b"), when = c(1, 2))
  expect_error(temporal_hypergraph(no_ends), class = "hypernets_bad_input")
  categorical_with_clock <- data.frame(x = c("a", "b"), time = c("1", "2"))
  expect_error(temporal_hypergraph(categorical_with_clock), class = "hypernets_bad_input")
  expect_error(temporal_hypergraph(list(character(0))), class = "hypernets_bad_input")
})

test_that("print names the membership times", {
  th <- temporal_hypergraph(wide)
  expect_snapshot(print(th))
})
