test_that("s-edge centralities equal graph centrality on the s-line graph", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c", "c", "d"),
               event = rep(c("e1", "e2", "e3"), each = 2)),
    "member", "event"
  )
  out <- hg_edge_centrality(h, s = 1)
  between <- out$value[out$measure == "betweenness"]
  names(between) <- out$edge[out$measure == "betweenness"]
  close <- out$value[out$measure == "closeness"]
  names(close) <- out$edge[out$measure == "closeness"]
  expect_equal(unname(between[c("e1", "e2", "e3")]), c(0, 1, 0))
  expect_equal(unname(close[c("e1", "e2", "e3")]), c(2 / 3, 1, 2 / 3))
})

test_that("edge centrality applies over temporal snapshots", {
  dat <- data.frame(
    member = c("a", "b", "b", "c", "c", "d"),
    event = rep(c("e1", "e2", "e3"), each = 2),
    time = rep(1:3, each = 2)
  )
  thg <- temporal_hypergraph(dat, actor = "member", group = "event", time = "time")
  out <- hg_edge_centrality(thg, measure = "closeness", snapshot_mode = "cumulative")
  expect_equal(unique(out$time), as.character(1:3))
  expect_equal(nrow(out), 1 + 2 + 3)
})

test_that("edge centrality uses the paper's disconnected-graph normalization", {
  h <- group_hypergraph(
    data.frame(
      member = c("a", "b", "b", "c", "x", "y"),
      event = rep(c("e1", "e2", "isolated"), each = 2)
    ),
    "member", "event"
  )
  out <- hg_edge_centrality(h, s = 1, normalized = TRUE)
  close <- out[out$measure == "closeness", c("edge", "value")]
  close <- stats::setNames(close$value, close$edge)

  # NetworkX/HypergraphX use the Wasserman--Faust correction: the two
  # connected vertices have raw closeness 1, multiplied by (2 - 1)/(3 - 1).
  expect_equal(unname(close[c("e1", "e2", "isolated")]), c(0.5, 0.5, 0))
})

test_that("subhypergraph centrality has the closed form for one dyad", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b"), event = "e1"), "member", "event"
  )
  out <- hypergraph_centrality(h, type = "subhypergraph")
  expect_equal(out$subhypergraph, rep(log(cosh(1)), 2), tolerance = 1e-12)
})
