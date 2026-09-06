test_that("growing temporal hypergraphs produce cumulative active snapshots", {
  dat <- data.frame(
    member = c("a", "b", "b", "c", "c", "d"),
    event = rep(c("e1", "e2", "e3"), each = 2),
    year = rep(1:3, each = 2),
    source = rep(c("s1", "s2", "s3"), each = 2)
  )
  thg <- temporal_hypergraph(dat, actor = "member", cooccur_by = "event",
                             time = "year")
  expect_s3_class(thg, "net_temporal_hypergraph")
  expect_equal(thg$times, 1:3)
  expect_equal(thg$evolution, "growing")
  expect_equal(hypergraph_snapshot(thg, 2)$n_hyperedges, 2)
  expect_equal(length(hypergraph_snapshots(thg)), 3)
  expect_equal(as.data.frame(thg), thg$memberships)
  # a column constant within each hyperedge is a hyperedge attribute
  expect_identical(thg$edge_data$source, c("s1", "s2", "s3"))
  sm <- summary(thg)
  expect_equal(sm$n_nodes, 4)
  expect_equal(sm$n_hyperedges, 3)
  expect_equal(sm$n_memberships, 6)
  expect_equal(sm$mean_edge_size, 2)
  expect_equal(sm$evolution, "growing")
})

test_that("interval snapshots use closed intervals", {
  dat <- data.frame(
    case = rep(c("A", "B"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
    start = rep(c(1, 2), each = 3), end = rep(c(3, 2), each = 3)
  )
  thg <- temporal_hypergraph(dat, actor = "arbitrator", cooccur_by = "case",
                             start = "start", end = "end")
  expect_equal(thg$evolution, "interval")
  expect_equal(hypergraph_snapshot(thg, 2)$n_hyperedges, 2)
  expect_equal(hypergraph_snapshot(thg, 3)$n_hyperedges, 1)
  expect_equal(hypergraph_snapshot(thg, 3, mode = "cumulative")$n_hyperedges, 2)
})

test_that("an edge list gives hyperedges of size two", {
  contacts <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
                         time = 1:3, kind = c("call", "mail", "call"))
  thg <- temporal_hypergraph(contacts, from = "from", to = "to", time = "time")
  expect_equal(thg$edges, c("e1", "e2", "e3"))
  expect_equal(unname(thg$edge_data$kind), c("call", "mail", "call"))
  snap <- hypergraph_snapshot(thg, 2)
  expect_equal(snap$n_hyperedges, 2)
  expect_equal(snap$size_distribution, c(size_2 = 2L))
  static <- group_hypergraph(contacts, from = "from", to = "to")
  expect_equal(static$n_hyperedges, 3)
  expect_equal(static$nodes, c("a", "b", "c"))
  expect_error(temporal_hypergraph(contacts, from = "from", to = "to",
                                   actor = "from", time = "time"),
               class = "honets_bad_input")
})

test_that("snapshot duplicate handling records multiplicity", {
  dat <- data.frame(
    member = rep(c("a", "b", "c"), 2),
    event = rep(c("e1", "e2"), each = 3), time = 1
  )
  thg <- temporal_hypergraph(dat, actor = "member", cooccur_by = "event",
                             time = "time")
  multi <- hypergraph_snapshot(thg, 1, multiedges = TRUE)
  simple <- hypergraph_snapshot(thg, 1, multiedges = FALSE)
  expect_equal(multi$n_hyperedges, 2)
  expect_equal(simple$n_hyperedges, 1)
  expect_equal(simple$edge_multiplicity, 2L)
})

test_that("temporal constructors reject inconsistent or invalid spells", {
  bad <- data.frame(member = c("a", "b"), event = "e1", time = c(1, 2))
  expect_error(temporal_hypergraph(bad, actor = "member", cooccur_by = "event",
                                   time = "time"),
               class = "honets_bad_input")
  interval <- data.frame(member = "a", event = "e1", start = 2, end = 1)
  expect_error(temporal_hypergraph(interval, actor = "member", cooccur_by = "event",
                                   start = "start", end = "end"),
               class = "honets_bad_input")
  expect_error(temporal_hypergraph(bad, actor = "member", cooccur_by = "event"),
               class = "honets_bad_input")
  expect_error(temporal_hypergraph(bad, actor = "member", cooccur_by = "event",
                                   time = "time", start = "time"),
               class = "honets_bad_input")
})
