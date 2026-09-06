testthat::skip_on_cran()

# e1 = {a, b, c}, e2 = {b, c, d}, e3 = {d, e}; sources s1, s1, s2.
.subset_fixture <- function(sparse = FALSE) {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3")
  )
  hg <- group_hypergraph(dat, "member", "event", sparse = sparse)
  hg$edge_data <- data.frame(edge = c("e1", "e2", "e3"),
                             source = c("s1", "s1", "s2"),
                             stringsAsFactors = FALSE)
  hg$edge_multiplicity <- c(1L, 2L, 1L)
  hg
}

test_that("hg_subset keeps named hyperedges and their nodes", {
  out <- hg_subset(.subset_fixture(), edges = c("e1", "e2"))
  expect_s3_class(out, "net_hypergraph")
  expect_identical(colnames(out$incidence), c("e1", "e2"))
  expect_identical(out$nodes, c("a", "b", "c", "d"))
  expect_identical(out$n_hyperedges, 2L)
  expect_identical(out$n_nodes, 4L)
  expect_identical(out$hyperedges, list(1:3, 2:4))
  expect_identical(out$size_distribution, c(size_3 = 2L))
  expect_identical(out$edge_data$edge, c("e1", "e2"))
  expect_identical(out$edge_multiplicity, c(1L, 2L))
})

test_that("hg_subset by nodes keeps the induced sub-hypergraph", {
  out <- hg_subset(.subset_fixture(), nodes = c("b", "c", "d", "e"))
  expect_identical(colnames(out$incidence), c("e2", "e3"))
  expect_identical(out$nodes, c("b", "c", "d", "e"))
  # a node set with no hyperedge inside gives an edgeless hypergraph
  empty <- hg_subset(.subset_fixture(), nodes = c("a", "e"))
  expect_identical(empty$n_hyperedges, 0L)
  expect_identical(empty$n_nodes, 2L)
})

test_that("hg_subset by hyperedge attribute uses the edge metadata", {
  out <- hg_subset(.subset_fixture(), where = c(source = "s2"))
  expect_identical(colnames(out$incidence), "e3")
  expect_identical(out$nodes, c("d", "e"))
  kept <- hg_subset(.subset_fixture(), where = list(source = "s2"),
                    drop_isolated = FALSE)
  expect_identical(kept$nodes, c("a", "b", "c", "d", "e"))
  both <- hg_subset(.subset_fixture(), where = list(source = c("s1", "s2")))
  expect_identical(both$n_hyperedges, 3L)
  no_source <- group_hypergraph(
    data.frame(member = c("a", "b"), event = "e1"), "member", "event"
  )
  expect_error(hg_subset(no_source, where = c(source = "s1")), class = "honets_bad_input")
  expect_error(hg_subset(.subset_fixture(), where = "s1"), class = "honets_bad_input")
})

test_that("INVARIANT: sparse and dense subsets agree", {
  dense <- hg_subset(.subset_fixture(), edges = c("e2", "e3"))
  sparse <- hg_subset(.subset_fixture(sparse = TRUE), edges = c("e2", "e3"))
  expect_true(methods::is(sparse$incidence, "sparseMatrix"))
  expect_identical(as.matrix(sparse$incidence) * 1, dense$incidence * 1)
  expect_identical(sparse$hyperedges, dense$hyperedges)
  expect_identical(sparse$nodes, dense$nodes)
})

test_that("hg_subset takes a ranking table through its edge column", {
  hg <- .subset_fixture()
  ranking <- hg_edge_centrality(hg, s = 1, measure = "betweenness", top = 2)
  from_table <- hg_subset(hg, edges = ranking)
  from_names <- hg_subset(hg, edges = ranking$edge)
  expect_identical(from_table, from_names)
  expect_error(hg_subset(hg, edges = data.frame(x = 1)), class = "honets_bad_input")
})

test_that("hg_subset rejects unknown names and empty selectors", {
  hg <- .subset_fixture()
  expect_error(hg_subset(hg), class = "honets_bad_input")
  expect_error(hg_subset(hg, edges = "nope"), class = "honets_bad_input")
  expect_error(hg_subset(hg, nodes = "nope"), class = "honets_bad_input")
  expect_error(hg_subset(42, edges = "e1"), class = "honets_bad_input")
})
