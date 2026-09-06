testthat::skip_on_cran()

# Three tribunals: A {p1, a1, a2} on [1, 4], B {p2, a1, a3} on [2, 3],
# C {p1, a4, a5} on [4, 6].
.tribunals <- function() {
  data.frame(
    case = rep(c("A", "B", "C"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
    constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3)
  )
}

.interval_thg <- function() {
  temporal_hypergraph(.tribunals(), actor = "arbitrator", cooccur_by = "case",
                      start = "constituted", end = "concluded")
}

# Growing: e1 {a, b} at 1, e2 {b, c} at 2, e3 {c, d} at 3, e4 {a, b} at 3.
.growing_thg <- function(nodes = NULL) {
  dat <- data.frame(
    member = c("a", "b", "b", "c", "c", "d", "a", "b"),
    event = rep(c("e1", "e2", "e3", "e4"), each = 2),
    time = rep(c(1, 2, 3, 3), each = 2)
  )
  temporal_hypergraph(dat, actor = "member", cooccur_by = "event", time = "time",
                      nodes = nodes)
}

test_that("hg_growth counts active and cumulative structure of interval data", {
  g <- hg_growth(.interval_thg())
  expect_s3_class(g, "honets_series")
  expect_s3_class(g, "data.frame")
  expect_identical(g$time, c(1, 2, 3, 4, 6))
  expect_identical(g$n_nodes, c(3L, 5L, 5L, 5L, 3L))
  expect_identical(g$n_edges, c(1L, 2L, 2L, 2L, 1L))
  expect_identical(g$n_edges_distinct, g$n_edges)
  expect_identical(g$n_memberships, c(3L, 6L, 6L, 6L, 3L))
  expect_identical(g$n_nodes_cumulative, c(3L, 5L, 5L, 7L, 7L))
  expect_identical(g$n_edges_cumulative, c(1L, 2L, 2L, 3L, 3L))
})

test_that("hg_growth on growing data is cumulative and counts distinct sets", {
  g <- hg_growth(.growing_thg())
  expect_identical(g$time, c(1, 2, 3))
  expect_identical(g$n_nodes, c(2L, 3L, 4L))
  expect_identical(g$n_edges, c(1L, 2L, 4L))
  expect_identical(g$n_edges_distinct, c(1L, 2L, 3L))
  expect_identical(g$n_memberships, c(2L, 4L, 8L))
  expect_false("n_nodes_cumulative" %in% names(g))
})

test_that("a node universe with entry times drives node counts and snapshots", {
  universe <- data.frame(node = c("a", "b", "c", "d", "z"),
                         start = c(1, 1, 2, 3, 0))
  thg <- .growing_thg(nodes = universe)
  expect_identical(thg$nodes, c("a", "b", "c", "d", "z"))
  g <- hg_growth(thg)
  expect_identical(g$n_nodes, c(3L, 4L, 5L))
  snap <- hypergraph_snapshot(thg, at = 1)
  expect_identical(snap$nodes, c("a", "b", "z"))
  expect_identical(snap$n_hyperedges, 1L)
  expect_identical(hypergraph_snapshot(thg, mode = "all")$n_nodes, 5L)
  # the node table's first column is the name; `time` also gives the entry
  named <- data.frame(key = c("a", "b", "c", "d", "z"), time = c(1, 1, 2, 3, 0))
  dat <- data.frame(member = c("a", "b", "b", "c", "c", "d", "a", "b"),
                    event = rep(c("e1", "e2", "e3", "e4"), each = 2),
                    time = rep(c(1, 2, 3, 3), each = 2))
  renamed <- temporal_hypergraph(dat, actor = "member", cooccur_by = "event",
                                 time = "time", nodes = named)
  expect_identical(renamed$node_data, thg$node_data)
  expect_error(temporal_hypergraph(dat, actor = "member", cooccur_by = "event",
                                   time = "time", nodes = list(1)),
               class = "honets_bad_input")
  # a bare universe keeps every node in every snapshot
  bare <- .growing_thg(nodes = c("a", "b", "c", "d", "z"))
  expect_identical(hypergraph_snapshot(bare, at = 1)$n_nodes, 5L)
  expect_identical(hg_growth(bare)$n_nodes, c(2L, 3L, 4L))
  expect_error(.growing_thg(nodes = c("a", "b")), class = "honets_bad_input")
})

test_that("component statistics follow the shared-hyperedge connectivity", {
  # A and B share a1 at times 2-3; A and C share p1 at time 4: always one
  # component, whose diameter is 2 whenever two tribunals are active.
  g <- hg_growth(.interval_thg(), components = TRUE)
  expect_identical(g$n_components, rep(1L, 5))
  expect_equal(g$largest_component, rep(1, 5))
  expect_identical(g$diameter, c(1L, 2L, 2L, 2L, 1L))

  snap <- hypergraph_snapshot(.interval_thg(), at = 4)
  comp <- hg_measures(snap, what = "components")
  expect_identical(names(comp), c("component", "n_nodes", "n_edges", "share", "diameter"))
  expect_identical(comp$n_nodes, 5L)
  expect_identical(comp$n_edges, 2L)
  expect_identical(comp$diameter, 2L)

  # two tribunals with no member in common are two components
  apart <- temporal_hypergraph(
    data.frame(case = rep(c("A", "B"), each = 3),
               arbitrator = c("p1", "a1", "a2", "p2", "a3", "a4"),
               start = 1, end = 2),
    actor = "arbitrator", cooccur_by = "case", start = "start", end = "end"
  )
  two <- hg_measures(hypergraph_snapshot(apart, at = 1), what = "components")
  expect_identical(two$component, 1:2)
  expect_identical(two$n_nodes, c(3L, 3L))
  expect_identical(two$n_edges, c(1L, 1L))
  expect_equal(two$share, c(0.5, 0.5))
  expect_identical(hg_growth(apart, components = TRUE)$n_components, c(2L, 2L))
})

test_that("diameter matches a hand-computed path", {
  path <- matrix(0, 4, 4)
  path[cbind(1:3, 2:4)] <- 1
  path <- path + t(path)
  expect_identical(honets:::.thg_diameter(path > 0), 3L)
  expect_identical(honets:::.thg_diameter(matrix(TRUE, 1, 1)), 0L)
  two <- diag(2) > 0
  expect_true(is.na(honets:::.thg_diameter(two)))
})

test_that("series and distribution tables plot", {
  g <- hg_growth(.interval_thg())
  expect_s3_class(plot(g), "ggplot")
  expect_s3_class(plot(g, columns = c("n_nodes", "n_edges"), facets = FALSE), "ggplot")
  expect_error(plot(g, columns = "nope"), class = "honets_bad_input")
  snap <- hypergraph_snapshot(.interval_thg(), at = 2)
  d <- hg_measures(snap, what = "distribution", measure = "n_neighbors")
  expect_s3_class(d, "honets_distribution")
  expect_identical(attr(d, "measure"), "n_neighbors")
  expect_equal(d$ccdf[[1L]], 1)
  expect_s3_class(plot(d), "ggplot")
  expect_s3_class(plot(d, log = FALSE), "ggplot")
})

test_that("temporal summary reports the observation window", {
  s <- summary(.interval_thg())
  expect_identical(s$first_time, 1)
  expect_identical(s$last_time, 6)
  expect_equal(s$mean_duration, mean(c(3, 1, 2)))
  expect_true(is.na(summary(.growing_thg())$mean_duration))
})

test_that("hg_growth rejects non-temporal input", {
  expect_error(hg_growth(42), class = "honets_bad_input")
  expect_error(hg_growth(.interval_thg(), at = NA), class = "honets_bad_input")
})
