skip_on_cran()

.cc_docs <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  cooking_3 = "salt the onions and simmer the soup",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night",
  space_3 = "stars and a galaxy through the telescope"
)
.cc_stops <- c("the", "with", "and", "a", "this", "at", "on", "all",
               "through")

test_that("hg_cocluster returns one tidy row per document and word", {
  hg <- text_hypergraph(.cc_docs, stop_words = .cc_stops)
  out <- hg_cocluster(hg, k = 2, seed = 1)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("node", "role", "cluster"))
  expect_identical(out$node, c(hg$nodes, colnames(hg$incidence)))
  expect_identical(out$role, rep(c("node", "hyperedge"),
                                 c(hg$n_nodes, hg$n_hyperedges)))
  expect_identical(out$cluster[[1L]], "Cluster 1")
  docs <- stats::setNames(out$cluster, out$node)
  expect_identical(unname(docs[c("cooking_1", "cooking_2", "cooking_3")]),
                   rep("Cluster 1", 3))
  expect_identical(unname(docs[c("space_1", "space_2", "space_3")]),
                   rep("Cluster 2", 3))
  # each topic's own words land with its documents
  expect_identical(unname(docs[c("soup", "onions", "simmer")]),
                   rep("Cluster 1", 3))
  expect_identical(unname(docs[c("telescope", "stars", "galaxy")]),
                   rep("Cluster 2", 3))
  words <- hg_cocluster(hg, k = 2, seed = 1, role = "hyperedge")
  expect_identical(unique(words$role), "hyperedge")
  expect_identical(nrow(words), hg$n_hyperedges)
  embedding <- hg_cocluster(hg, k = 3, seed = 1, what = "embedding")
  expect_named(embedding, c("node", "role", "cluster", "dim1", "dim2"))
})

test_that("the embedding is Dhillon's scaled singular vectors (hand check, k = 2)", {
  hg <- text_hypergraph(.cc_docs, stop_words = .cc_stops)
  a <- as.matrix(hg$incidence)
  d1 <- rowSums(a)
  d2 <- colSums(a)
  s <- svd(a / sqrt(d1) / rep(sqrt(d2), each = nrow(a)))
  expected <- unname(c(s$u[, 2] / sqrt(d1), s$v[, 2] / sqrt(d2)))
  got <- hg_cocluster(hg, k = 2, seed = 1, what = "embedding")$dim1
  expect_equal(abs(got), abs(expected))
  expect_equal(abs(sum(got * expected)), sum(expected^2))
  # the leading singular value of the scaled matrix is 1 (the trivial pair)
  expect_equal(s$d[[1L]], 1)
})

test_that("co-clustering invariants: weight scale, sparse storage, seed", {
  hg <- text_hypergraph(.cc_docs, stop_words = .cc_stops)
  base <- hg_cocluster(hg, k = 2, seed = 1, what = "embedding")
  scaled <- hg
  scaled$incidence <- scaled$incidence * 7
  # A_n is scale-free; Z = D^-1/2 U shrinks by 1/sqrt(7), the partition stays
  rescaled <- hg_cocluster(scaled, k = 2, seed = 1, what = "embedding")
  expect_identical(rescaled$cluster, base$cluster)
  expect_equal(rescaled$dim1, base$dim1 / sqrt(7))
  sparse <- text_hypergraph(.cc_docs, stop_words = .cc_stops, sparse = TRUE)
  from_sparse <- hg_cocluster(sparse, k = 2, seed = 1, what = "embedding")
  expect_identical(from_sparse$cluster, base$cluster)
  expect_equal(from_sparse$dim1, base$dim1)
  set.seed(9)
  before <- stats::runif(1)
  set.seed(9)
  again <- hg_cocluster(hg, k = 2, seed = 1)
  expect_identical(stats::runif(1), before)
  expect_identical(again$cluster, base$cluster)
})

test_that("hg_cocluster raises classed errors", {
  hg <- text_hypergraph(.cc_docs, stop_words = .cc_stops)
  expect_error(hg_cocluster(hg, k = 1), class = "hypernets_bad_input")
  expect_error(hg_cocluster(hg, k = 1000), class = "hypernets_bad_input")
  empty <- hg
  empty$incidence[, 1L] <- 0
  expect_error(hg_cocluster(empty, k = 2), class = "hypernets_bad_input")
  expect_error(hg_cocluster(list(), k = 2), class = "hypernets_bad_input")
})
