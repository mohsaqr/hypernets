# Same fixture as test-project.R:
#   e1 = {a, b, c}   e2 = {c, d}   e3 = {a, d}   e4 = {a, b, c, d}
# Every vertex pair shares some hyperedge, so the vertex adjacency is K4.
#   size             3, 2, 2, 4
#   n_incident_edges 3, 3, 3, 3   (every pair of hyperedges overlaps)
#   n_neighbors      1, 2, 2, 0   (K4 minus the hyperedge's own members)
edge_fixture <- function(weight = 1) {
  long <- data.frame(
    vertex = c("a", "b", "c", "c", "d", "a", "d", "a", "b", "c", "d"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e3", "e3",
             "e4", "e4", "e4", "e4"),
    w = weight
  )
  group_hypergraph(long, actor = "vertex", group = "edge",
                              weight = "w")
}

test_that("hyperedge measures match the hand-computed values", {
  got <- hg_edges(edge_fixture())
  expect_identical(names(got), c("edge", "size", "weight",
                                 "n_incident_edges", "n_neighbors"))
  expect_identical(got$edge, c("e1", "e2", "e3", "e4"))
  expect_identical(got$size, c(3L, 2L, 2L, 4L))
  expect_identical(got$n_incident_edges, c(3L, 3L, 3L, 3L))
  expect_identical(got$n_neighbors, c(1L, 2L, 2L, 0L))
  expect_equal(got$weight, c(3, 2, 2, 4))
})

test_that("weight sums the incidence, size counts membership", {
  got <- hg_edges(edge_fixture(weight = c(2, 1, 3, 1, 4, 2, 1, 1, 5, 1, 2)))
  expect_identical(got$size, c(3L, 2L, 2L, 4L))
  expect_equal(got$weight, c(2 + 1 + 3, 1 + 4, 2 + 1, 1 + 5 + 1 + 2))
})

test_that("INVARIANT: size agrees with hg_measures, n_incident with the line graph", {
  hg <- edge_fixture()
  edge_table <- hg_edges(hg)
  edge_measures <- hg_measures(hg, what = "edges")
  expect_identical(edge_table$size, edge_measures$size)
  # n_incident_edges is the s = 1 line-graph degree.
  line <- hg_line_graph(hg, s = 1)
  degree <- table(c(line$from, line$to))
  expect_identical(as.integer(degree[edge_table$edge]),
                   edge_table$n_incident_edges)
})

test_that("INVARIANT: isolated and singleton hyperedges are handled", {
  long <- data.frame(vertex = c("a", "b", "c", "z"),
                     edge = c("e1", "e1", "e1", "solo"), w = 1)
  hg <- group_hypergraph(long, actor = "vertex", group = "edge",
                                    weight = "w")
  got <- hg_edges(hg)
  # "solo" holds one vertex, shares it with nothing, and reaches nobody.
  expect_identical(got$size, c(3L, 1L))
  expect_identical(got$n_incident_edges, c(0L, 0L))
  expect_identical(got$n_neighbors, c(0L, 0L))
})

test_that("INVARIANT: measures are invariant to input row order", {
  long <- data.frame(
    vertex = c("d", "a", "c", "b", "a", "c", "d", "b", "a", "c", "d"),
    edge = c("e4", "e3", "e1", "e4", "e1", "e2", "e2", "e1",
             "e4", "e4", "e3"),
    w = 1
  )
  shuffled <- group_hypergraph(long, actor = "vertex",
                                          group = "edge", weight = "w")
  ordered_edges <- hg_edges(edge_fixture())
  shuffled_edges <- hg_edges(shuffled)
  expect_equal(ordered_edges, shuffled_edges)
})

test_that("INVARIANT: sparse and dense incidences agree", {
  docs <- c(a = "salt and soup and night", b = "soup and stars",
            c = "stars and night", d = "night and pepper and salt")
  dense_hg <- text_hypergraph(docs)
  dense_edges <- hg_edges(dense_hg)
  sparse_hg <- text_hypergraph(docs, sparse = TRUE)
  sparse_edges <- hg_edges(sparse_hg)
  expect_equal(dense_edges, sparse_edges)
})

test_that("the distribution is a proper CCDF", {
  got <- hg_edges(edge_fixture(), what = "distribution")
  expect_identical(names(got), c("value", "n", "proportion", "ccdf"))
  expect_equal(got$value, c(2, 3, 4))
  expect_identical(got$n, c(2L, 1L, 1L))
  expect_equal(got$proportion, c(0.5, 0.25, 0.25))
  expect_equal(got$ccdf, c(1, 0.5, 0.25))
})

test_that("INVARIANT: CCDF starts at 1, is non-increasing, proportions sum to 1", {
  hg <- text_hypergraph(covid_abstracts[seq_len(20), ], column = "abstract",
                        id = "doc")
  got <- hg_edges(hg, what = "distribution", measure = "n_neighbors")
  expect_equal(got$ccdf[1], 1)
  expect_true(all(diff(got$ccdf) < 0))
  expect_equal(sum(got$proportion), 1)
  edge_table <- hg_edges(hg)
  expect_equal(sum(got$n), nrow(edge_table))
})

test_that("bad input raises classed conditions", {
  expect_error(hg_edges(42), class = "hypernets_bad_input")
  expect_error(hg_edges(edge_fixture(), what = "nope"), "arg")
  expect_error(hg_edges(edge_fixture(), measure = "nope"), "arg")
})
