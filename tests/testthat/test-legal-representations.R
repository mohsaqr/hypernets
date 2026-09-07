testthat::skip_on_cran()

# A {p1, a1, a2}, B {p2, a1, a3}, C {p1, a1, a2} (a repeat of A).
.rep_tribunals <- function(sparse = FALSE) {
  group_hypergraph(
    data.frame(
      member = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a1", "a2"),
      case = rep(c("A", "B", "C"), each = 3)
    ),
    "member", "case", sparse = sparse
  )
}

# e1 {b, c} cited by a, e2 {b} cited by a, e3 {c} cited by b.
.rep_citations <- function(sparse = FALSE) {
  hg <- group_hypergraph(
    data.frame(member = c("b", "c", "b", "c"),
               block = c("e1", "e1", "e2", "e3")),
    "member", "block", sparse = sparse
  )
  hg$edge_data <- data.frame(edge = c("e1", "e2", "e3"),
                             source = c("a", "a", "b"),
                             stringsAsFactors = FALSE)
  hg
}

test_that("clique representations match hand-computed Table 2 statistics", {
  out <- hg_representations(.rep_tribunals())
  expect_identical(out$representation, c("bg", "mg", "bh", "mh"))
  expect_identical(out$type, rep(c("graph", "hypergraph"), each = 2))
  expect_identical(out$n_nodes, rep(5L, 4))
  expect_identical(out$n_edges, c(6L, 9L, 2L, 3L))
  expect_equal(out$mean_degree, c(12 / 5, 18 / 5, 6 / 5, 9 / 5))
  # multigraph degrees: p1 4, a1 6, a2 4, p2 2, a3 2
  expect_equal(out$median_degree, c(2, 4, 1, 2))
})

test_that("citation representations count directed edges and degrees", {
  out <- hg_representations(.rep_citations(), graph = "citation")
  expect_identical(out$n_edges, c(3L, 4L, 3L, 3L))
  expect_identical(out$n_nodes[1:2], c(3L, 3L))
  # bg out-degrees a 2, b 1, c 0 and in-degrees a 0, b 1, c 2;
  # mg out-degrees a 3, b 1, c 0 and in-degrees a 0, b 2, c 2
  expect_equal(out$mean_degree[1:2], c(1, 4 / 3))
  expect_equal(out$median_out_degree[1:2], c(1, 1))
  expect_equal(out$median_in_degree[1:2], c(1, 2))
  expect_true(all(is.na(out$median_out_degree[3:4])))
})

test_that("INVARIANT: sparse and dense representations agree", {
  sparse_tribunals <- hg_representations(.rep_tribunals(sparse = TRUE))
  dense_tribunals <- hg_representations(.rep_tribunals())
  expect_equal(sparse_tribunals, dense_tribunals)
  sparse_citations <- hg_representations(.rep_citations(sparse = TRUE), graph = "citation")
  dense_citations <- hg_representations(.rep_citations(), graph = "citation")
  expect_equal(sparse_citations, dense_citations)
})

test_that("citation projection is the source-to-member graph", {
  hg <- .rep_citations()
  directed <- hg_project(hg, method = "citation", what = "matrix", directed = TRUE)
  expect_identical(rownames(directed), c("a", "b", "c"))
  expect_equal(directed["a", "b"], 2)
  expect_equal(directed["a", "c"], 1)
  expect_equal(directed["b", "c"], 1)
  expect_equal(sum(directed), 4)
  binary <- hg_project(hg, method = "citation", what = "matrix",
                       directed = TRUE, duplicate_edges = "collapse")
  expect_equal(sum(binary), 3)
  undirected <- hg_project(hg, method = "citation", what = "matrix")
  expect_true(isSymmetric(unname(as.matrix(undirected))))
  edges <- hg_project(hg, method = "citation")
  expect_identical(names(edges), c("from", "to", "weight"))
  expect_identical(edges$from, c("a", "a", "b"))
  expect_equal(edges$weight, c(2, 1, 1))
  directed_edges <- hg_project(hg, method = "citation", directed = TRUE)
  expect_identical(directed_edges$to, c("b", "c", "c"))
  # explicit sources override, and missing sources are an error
  bare <- .rep_citations()
  bare$edge_data <- NULL
  expect_error(hg_project(bare, method = "citation"), class = "hypernets_bad_input")
  explicit <- hg_project(bare, method = "citation", what = "matrix",
                         edge_source = c(e1 = "a", e2 = "a", e3 = "b"))
  expect_equal(as.matrix(explicit), as.matrix(undirected))
  # the source may be named as an attribute column of the hyperedges
  named <- .rep_citations()
  names(named$edge_data)[2L] <- "citing"
  named_projection <- hg_project(named, method = "citation", what = "matrix",
                                 edge_source = "citing")
  expect_equal(as.matrix(named_projection), as.matrix(undirected))
  expect_error(hg_project(named, method = "citation"), class = "hypernets_bad_input")
  expect_error(hg_project(hg, method = "clique", directed = TRUE),
               class = "hypernets_bad_input")
  expect_error(hg_project(hg, method = "citation", weighted = FALSE),
               class = "hypernets_bad_input")
})

test_that("hg_measures reports neighbours, distributions and components", {
  hg <- .rep_tribunals()
  nodes <- hg_measures(hg, what = "nodes")
  expect_identical(names(nodes), c("node", "hyperdegree", "strength",
                                   "max_edge_size", "n_neighbors"))
  expect_identical(nodes$n_neighbors[nodes$node == "a1"], 4L)
  expect_identical(nodes$n_neighbors[nodes$node == "a3"], 2L)
  sparse_nodes <- hg_measures(.rep_tribunals(sparse = TRUE), what = "nodes")
  expect_equal(sparse_nodes, nodes)
  d <- hg_measures(hg, what = "distribution", measure = "size")
  expect_identical(d$value, 3)
  expect_equal(d$ccdf, 1)
  comp <- hg_measures(hg, what = "components")
  expect_identical(comp$n_nodes, 5L)
  expect_identical(comp$n_edges, 3L)
  expect_identical(comp$diameter, 2L)
})
