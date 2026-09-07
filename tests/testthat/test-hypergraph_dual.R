# Dual hypergraph: transpose identity, involution, and equality with the
# opposite-orientation construction.

dual_corpus <- c(a = "salt and soup", b = "soup and stars",
                 c = "stars and salt soup")

test_that("the dual incidence is exactly the transpose", {
  hg <- text_hypergraph(dual_corpus, weight = "tfidf")
  dual <- dual_hypergraph(hg)
  expect_identical(dual$incidence, t(hg$incidence))
  expect_identical(dual$n_nodes, hg$n_hyperedges)
  expect_identical(dual$n_hyperedges, hg$n_nodes)
})

test_that("the dual of the dual restores the original incidence", {
  hg <- text_hypergraph(dual_corpus, weight = "tfidf")
  dual <- dual_hypergraph(hg)
  back <- dual_hypergraph(dual)
  expect_identical(back$incidence, hg$incidence)
})

test_that("the dual of a bag hypergraph equals the opposite orientation", {
  doc_hg <- text_hypergraph(dual_corpus, nodes = "doc", weight = "tfidf")
  word_hg <- text_hypergraph(dual_corpus, nodes = "word", weight = "tfidf")
  dual <- dual_hypergraph(doc_hg)
  expect_identical(dual$incidence, word_hg$incidence)
  expect_identical(dual$text$nodes, "word")
  expect_s3_class(dual, "text_hypergraph")
  # and every verb sees the same hypergraph
  dual_measures <- hg_measures(dual)
  word_measures <- hg_measures(word_hg)
  expect_identical(dual_measures, word_measures)
})

test_that("duals of non-bag constructions drop the text layer", {
  win <- text_hypergraph(c(d = "a b c a b"), construction = "window",
                         window = 2)
  dual <- dual_hypergraph(win)
  expect_false(inherits(dual, "text_hypergraph"))
  expect_identical(dual$incidence, t(win$incidence))
})

test_that("dual rejects non-hypergraphs with a classed error", {
  expect_error(dual_hypergraph("x"), class = "hypernets_bad_input")
})
