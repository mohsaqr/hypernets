test_that("growing temporal hypergraphs produce cumulative active snapshots", {
  dat <- data.frame(
    member = c("a", "b", "b", "c", "c", "d"),
    event = rep(c("e1", "e2", "e3"), each = 2),
    year = rep(1:3, each = 2),
    source = rep(c("s1", "s2", "s3"), each = 2)
  )
  thg <- temporal_hypergraph(dat, "member", "event", time = "year",
                             source = "source")
  expect_s3_class(thg, "net_temporal_hypergraph")
  expect_equal(thg$times, 1:3)
  expect_equal(hypergraph_snapshot(thg, 2)$n_hyperedges, 2)
  expect_equal(length(hypergraph_snapshots(thg)), 3)
  expect_equal(as.data.frame(thg), thg$memberships)
  sm <- summary(thg)
  expect_equal(sm$n_nodes, 4)
  expect_equal(sm$n_hyperedges, 3)
  expect_equal(sm$n_memberships, 6)
  expect_equal(sm$mean_edge_size, 2)
  expect_equal(sm$evolution, "growing")
})

test_that("interval snapshots use closed intervals and wide member columns", {
  dat <- data.frame(
    case = c("A", "B"), p = c("p1", "p2"), a1 = c("a1", "a1"),
    a2 = c("a2", "a3"), start = c(1, 2), end = c(3, 2)
  )
  thg <- temporal_hypergraph(dat, c("p", "a1", "a2"), "case",
                             start = "start", end = "end",
                             evolution = "interval")
  expect_equal(hypergraph_snapshot(thg, 2)$n_hyperedges, 2)
  expect_equal(hypergraph_snapshot(thg, 3)$n_hyperedges, 1)
  expect_equal(hypergraph_snapshot(thg, 3, mode = "cumulative")$n_hyperedges, 2)
})

test_that("snapshot duplicate handling records multiplicity", {
  dat <- data.frame(
    member = rep(c("a", "b", "c"), 2),
    event = rep(c("e1", "e2"), each = 3), time = 1
  )
  thg <- temporal_hypergraph(dat, "member", "event", time = "time")
  multi <- hypergraph_snapshot(thg, 1, multiedges = TRUE)
  simple <- hypergraph_snapshot(thg, 1, multiedges = FALSE)
  expect_equal(multi$n_hyperedges, 2)
  expect_equal(simple$n_hyperedges, 1)
  expect_equal(simple$edge_multiplicity, 2L)
})

test_that("temporal constructors reject inconsistent or invalid spells", {
  bad <- data.frame(member = c("a", "b"), event = "e1", time = c(1, 2))
  expect_error(temporal_hypergraph(bad, "member", "event", time = "time"),
               class = "honets_bad_input")
  interval <- data.frame(member = "a", event = "e1", start = 2, end = 1)
  expect_error(temporal_hypergraph(interval, "member", "event",
                                   start = "start", end = "end",
                                   evolution = "interval"),
               class = "honets_bad_input")
})
