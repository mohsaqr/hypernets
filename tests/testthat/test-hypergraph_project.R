# Fixture with hand-computable projections.
#
#   e1 = {a, b, c}   e2 = {c, d}   e3 = {a, d}   e4 = {a, b, c, d}
#
# association weights, w({u,v}) = sum_e 1/(|e| - 1):
#   ab = 1/2 + 1/3 = 5/6      ac = 1/2 + 1/3 = 5/6      ad = 1 + 1/3 = 4/3
#   bc = 1/2 + 1/3 = 5/6      bd =       1/3            cd = 1 + 1/3 = 4/3
# clique weights (all incidence weights are 1, so weighted == binary):
#   ab = 2  ac = 2  ad = 2  bc = 2  bd = 1  cd = 2
# pairwise hyperedge overlaps:
#   e1e2 = 1  e1e3 = 1  e1e4 = 3  e2e3 = 1  e2e4 = 2  e3e4 = 2
fixture <- function(weight = 1) {
  long <- data.frame(
    vertex = c("a", "b", "c", "c", "d", "a", "d", "a", "b", "c", "d"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e3", "e3",
             "e4", "e4", "e4", "e4"),
    w = weight
  )
  group_hypergraph(long, actor = "vertex", group = "edge",
                              weight = "w")
}

test_that("association weights match the hand-computed values", {
  got <- hg_project(fixture(), method = "association")
  expect_identical(names(got), c("from", "to", "weight"))
  expect_identical(got$from, c("a", "a", "a", "b", "b", "c"))
  expect_identical(got$to, c("b", "c", "d", "c", "d", "d"))
  expect_equal(got$weight, c(5 / 6, 5 / 6, 4 / 3, 5 / 6, 1 / 3, 4 / 3))
})

test_that("clique weights match the hand-computed values", {
  got <- hg_project(fixture(), method = "clique")
  expect_identical(got$from, c("a", "a", "a", "b", "b", "c"))
  expect_identical(got$to, c("b", "c", "d", "c", "d", "d"))
  expect_equal(got$weight, c(2, 2, 2, 2, 1, 2))
})

test_that("INVARIANT: association degree equals hyperdegree", {
  # Each hyperedge adds exactly weight 1 around each of its members, so the
  # weighted degree in the projection is the number of incident hyperedges.
  # This is the random-walk property the 1/(|e| - 1) normalisation exists for.
  hg <- fixture()
  w <- hg_project(hg, method = "association", what = "matrix")
  expect_equal(as.numeric(rowSums(w)),
               as.numeric(rowSums(hg$incidence != 0)))
  expect_equal(as.numeric(rowSums(w)), c(3, 2, 3, 3))
})

test_that("INVARIANT: association degree counts non-singleton hyperedges", {
  # "pepper" occurs in one document only, so it is a singleton hyperedge with
  # no pairs to carry weight: the degree identity counts hyperedges of size
  # >= 2, not every incident hyperedge.
  hg <- text_hypergraph(c(
    a = "salt and soup and night", b = "soup and stars",
    c = "stars and night", d = "night and pepper and salt"
  ))
  membership <- hg$incidence != 0
  expect_true(any(colSums(membership) == 1))
  shared <- membership[, colSums(membership) >= 2, drop = FALSE]
  w <- hg_project(hg, method = "association", what = "matrix")
  expect_equal(as.numeric(rowSums(w)), as.numeric(rowSums(shared)))
})

test_that("INVARIANT: projections are invariant to input row order", {
  hg <- fixture()
  long <- data.frame(
    vertex = c("d", "a", "c", "b", "a", "c", "d", "b", "a", "c", "d"),
    edge = c("e4", "e3", "e1", "e4", "e1", "e2", "e2", "e1",
             "e4", "e4", "e3"),
    w = 1
  )
  shuffled <- group_hypergraph(long, actor = "vertex",
                                          group = "edge", weight = "w")
  association <- hg_project(hg, method = "association")
  association_shuffled <- hg_project(shuffled, method = "association")
  expect_equal(association, association_shuffled)
  line <- hg_line_graph(hg)
  line_shuffled <- hg_line_graph(shuffled)
  expect_equal(line, line_shuffled)
})

test_that("clique projection matches clique_expansion()", {
  hg <- fixture(weight = c(2, 1, 3, 1, 4, 2, 1, 1, 5, 1, 2))
  ours <- hg_project(hg, method = "clique", what = "matrix")
  expansion <- clique_expansion(hg)
  theirs <- expansion$weights
  expect_equal(unname(as.matrix(ours)), unname(theirs))
})

test_that("weighted = FALSE drops the incidence weights", {
  weighted_hg <- fixture(weight = c(2, 1, 3, 1, 4, 2, 1, 1, 5, 1, 2))
  binary <- hg_project(weighted_hg, method = "clique", weighted = FALSE)
  unweighted <- hg_project(fixture(), method = "clique")
  expect_equal(binary, unweighted)
})

test_that("s-line graph overlaps match, and s filters them", {
  hg <- fixture()
  s1 <- hg_line_graph(hg)
  expect_identical(s1$from, c("e1", "e1", "e1", "e2", "e2", "e3"))
  expect_identical(s1$to, c("e2", "e3", "e4", "e3", "e4", "e4"))
  expect_equal(s1$weight, c(1, 1, 3, 1, 2, 2))

  s2 <- hg_line_graph(hg, s = 2)
  expect_identical(s2$from, c("e1", "e2", "e3"))
  expect_identical(s2$to, c("e4", "e4", "e4"))
  expect_equal(s2$weight, c(3, 2, 2))

  s3 <- hg_line_graph(hg, s = 3)
  expect_equal(s3$weight, 3)
  s4 <- hg_line_graph(hg, s = 4)
  expect_identical(nrow(s4), 0L)
  expect_identical(names(s4),
                   c("from", "to", "weight"))
})

test_that("INVARIANT: the s = 1 line graph is the dual's clique projection", {
  hg <- fixture()
  line <- hg_line_graph(hg, s = 1)
  dual <- dual_hypergraph(hg)
  dual_projection <- hg_project(dual, weighted = FALSE)
  expect_equal(line, dual_projection)
})

test_that("INVARIANT: sparse and dense incidences agree", {
  docs <- c(a = "salt and soup and night", b = "soup and stars",
            c = "stars and night", d = "night and pepper and salt")
  dense <- text_hypergraph(docs)
  sparse <- text_hypergraph(docs, sparse = TRUE)
  expect_true(methods::is(sparse$incidence, "sparseMatrix"))
  association_dense <- hg_project(dense, method = "association")
  association_sparse <- hg_project(sparse, method = "association")
  expect_equal(association_dense, association_sparse)
  clique_dense <- hg_project(dense, method = "clique")
  clique_sparse <- hg_project(sparse, method = "clique")
  expect_equal(clique_dense, clique_sparse)
  line_dense <- hg_line_graph(dense, s = 2)
  line_sparse <- hg_line_graph(sparse, s = 2)
  expect_equal(line_dense, line_sparse)
})

test_that("singleton hyperedges contribute nothing instead of dividing by zero", {
  long <- data.frame(vertex = c("a", "b", "c", "z"),
                     edge = c("e1", "e1", "e1", "solo"), w = 1)
  hg <- group_hypergraph(long, actor = "vertex", group = "edge",
                                    weight = "w")
  got <- hg_project(hg, method = "association")
  expect_true(all(is.finite(got$weight)))
  expect_false("z" %in% c(got$from, got$to))
  expect_equal(got$weight, c(0.5, 0.5, 0.5))
})

test_that("the matrix hand-off is symmetric, zero-diagonal and named", {
  hg <- fixture()
  w <- hg_project(hg, method = "association", what = "matrix")
  expect_equal(unname(as.matrix(w)), unname(t(as.matrix(w))))
  expect_equal(as.numeric(diag(as.matrix(w))), rep(0, 4))
  expect_identical(rownames(w), c("a", "b", "c", "d"))
  line <- hg_line_graph(hg, what = "matrix")
  expect_identical(rownames(line), c("e1", "e2", "e3", "e4"))
  expect_equal(as.numeric(diag(as.matrix(line))), rep(0, 4))
})

test_that("bad input raises classed conditions", {
  expect_error(hg_project(42), class = "honets_bad_input")
  expect_error(hg_line_graph(list()), class = "honets_bad_input")
  expect_error(
    hg_project(fixture(), method = "association", weighted = FALSE),
    class = "honets_bad_input"
  )
  expect_error(
    hg_project(fixture(), method = "association", weighted = TRUE),
    class = "honets_bad_input"
  )
  expect_error(hg_project(fixture(), weighted = NA), "TRUE or FALSE")
  expect_error(hg_line_graph(fixture(), s = 0), "at least 1")
  expect_error(hg_line_graph(fixture(), s = 1.5), "whole number")
  expect_error(hg_line_graph(fixture(), s = c(1, 2)), "single")
})
