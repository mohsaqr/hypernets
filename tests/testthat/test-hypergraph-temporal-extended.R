testthat::skip_on_cran()

# Three tribunals as seats: A {p1, a1, a2} on [1, 4], B {p2, a1, a3} on
# [2, 3], C {p1, a4, a5} on [4, 6]; sector and pending are per case, the
# seat varies within a case.
.tribunal_seats <- function() {
  data.frame(
    case = rep(c("A", "B", "C"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
    seat = rep(c("president", "arbitrator_1", "arbitrator_2"), 3),
    constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3),
    sector = rep(c("oil", "gas", "oil"), each = 3),
    pending = rep(c(FALSE, FALSE, TRUE), each = 3)
  )
}

test_that("constant columns become hyperedge attributes and travel into subsets", {
  thg <- temporal_hypergraph(
    .tribunal_seats(), actor = "arbitrator", group = "case",
    start = "constituted", end = "concluded"
  )
  # `seat` varies within a case and is not a hyperedge attribute
  expect_identical(names(thg$edge_data),
                   c("edge", "start", "end", "sector", "pending"))
  expect_identical(thg$params$attributes, c("sector", "pending"))
  snap <- hypergraph_snapshot(thg, at = 4)
  expect_identical(snap$edge_data$edge, colnames(snap$incidence))
  expect_identical(snap$edge_data$sector, c("oil", "oil"))
  sub <- hg_subset(snap, edges = "C")
  expect_identical(sub$edge_data$pending, TRUE)
  sub_plot <- plot(sub, color_by = "sector", linetype_by = "pending")
  expect_s3_class(sub_plot, "ggplot")
})

test_that("sparse temporal hypergraphs give sparse snapshots identical to dense", {
  dat <- .tribunal_seats()
  universe <- c("p1", "p2", "a1", "a2", "a3", "a4", "a5", "zz")
  build <- function(sparse) {
    temporal_hypergraph(dat, actor = "arbitrator", group = "case",
                        start = "constituted", end = "concluded",
                        sparse = sparse, nodes = universe)
  }
  dense <- hypergraph_snapshot(build(FALSE), at = 2)
  sparse <- hypergraph_snapshot(build(TRUE), at = 2)
  expect_true(methods::is(sparse$incidence, "sparseMatrix"))
  expect_false(methods::is(dense$incidence, "sparseMatrix"))
  expect_identical(as.matrix(sparse$incidence) * 1, dense$incidence * 1)
  expect_identical(sparse$hyperedges, dense$hyperedges)
  expect_identical(sparse$nodes, dense$nodes)
  expect_true("zz" %in% sparse$nodes)
  # empty snapshot keeps the universe, in both storages
  none <- hypergraph_snapshot(build(TRUE), at = 0.5)
  expect_identical(none$n_hyperedges, 0L)
  expect_identical(none$n_nodes, 8L)
  expect_true(methods::is(none$incidence, "sparseMatrix"))
})

test_that("duplicate collapse works on sparse incidence and keeps metadata", {
  dat <- data.frame(
    member = rep(c("a", "b", "c"), 3),
    event = rep(c("e1", "e2", "e3"), each = 3), time = 1,
    label = rep(c("x", "y", "z"), each = 3)
  )
  thg <- temporal_hypergraph(dat, actor = "member", group = "event",
                             time = "time", sparse = TRUE)
  simple <- hypergraph_snapshot(thg, 1, multiedges = FALSE)
  expect_identical(simple$n_hyperedges, 1L)
  expect_identical(simple$edge_multiplicity, 3L)
  expect_identical(simple$edge_data$label, "x")
  expect_identical(simple$hyperedges, list(1:3))
})

test_that("hg_subset(size = ) keeps hyperedges by member count", {
  hg <- group_hypergraph(data.frame(
    member = c("a", "b", "c", "a", "b", "c", "d", "e"),
    group = c("e1", "e1", "e1", "e2", "e2", "e3", "e3", "e3")
  ), actor = "member", group = "group")
  three <- hg_subset(hg, size = 3)
  expect_identical(three$n_hyperedges, 2L)
  expect_identical(sort(colnames(three$incidence)), c("e1", "e3"))
  two_or_three <- hg_subset(hg, size = 2:3)
  expect_identical(two_or_three$n_hyperedges, 3L)
  expect_error(hg_subset(hg, size = 0), class = "hypernets_bad_input")
  expect_error(hg_subset(hg), class = "hypernets_bad_input")
})
