# kNN embedding hypergraph: hand-computed neighbors, invariants, the
# text_hypergraph() knn construction, and classed error paths.

cosine <- function(u, v) sum(u * v) / (sqrt(sum(u^2)) * sqrt(sum(v^2)))

two_pairs <- matrix(
  c(1, 0.2,  0.9, 0.3,  0.2, 1,  0.3, 0.9),
  nrow = 4, byrow = TRUE,
  dimnames = list(c("a", "b", "c", "d"), NULL)
)

test_that("neighbors and weights are hand-computed correctly", {
  hg <- knn_hypergraph(two_pairs, k = 1)
  expect_identical(hg$n_hyperedges, 4L)
  expect_identical(hg$n_nodes, 4L)
  # nearest neighbor of a is b (and vice versa); of c is d (and vice versa)
  expect_equal(
    hg$incidence["b", "a"],
    cosine(two_pairs["a", ], two_pairs["b", ])
  )
  expect_equal(
    hg$incidence["d", "c"],
    cosine(two_pairs["c", ], two_pairs["d", ])
  )
  expect_identical(hg$incidence["c", "a"], 0)
  # the center sits in its own edge with weight exactly 1
  expect_identical(hg$incidence["a", "a"], 1)
  expect_identical(hg$incidence["d", "d"], 1)
})

test_that("every hyperedge has exactly k + 1 members", {
  hg <- knn_hypergraph(two_pairs, k = 2)
  edges <- hg_measures(hg, what = "edges")
  expect_true(all(edges$size == 3L))
  expect_identical(nrow(edges), 4L)
})

test_that("binary weighting gives an unweighted membership incidence", {
  hg <- knn_hypergraph(two_pairs, k = 1, weight = "binary")
  expect_true(all(hg$incidence %in% c(0, 1)))
  expect_identical(hg$knn$weight, "binary")
})

test_that("construction is deterministic, with alphabetical tie-breaks", {
  first <- knn_hypergraph(two_pairs, k = 2)
  second <- knn_hypergraph(two_pairs, k = 2)
  expect_identical(first, second)
  ties <- matrix(c(1, 0, 1, 0, 1, 0), nrow = 3, byrow = TRUE,
                 dimnames = list(c("z", "m", "a"), NULL))
  hg <- knn_hypergraph(ties, k = 1)
  # all three rows identical: the neighbor is always the alphabetically
  # first other item -- a for z and m, m for a
  expect_equal(hg$incidence["a", "z"], 1)
  expect_equal(hg$incidence["a", "m"], 1)
  expect_equal(hg$incidence["m", "a"], 1)
  expect_identical(hg$incidence["z", "m"], 0)
})

test_that("contract violations raise classed errors", {
  no_names <- matrix(seq_len(6) + 0.5, nrow = 3)
  expect_error(knn_hypergraph(no_names, k = 1), class = "honets_bad_input")
  dup <- matrix(seq_len(4) + 0.5, nrow = 2,
                dimnames = list(c("a", "a"), NULL))
  expect_error(knn_hypergraph(dup, k = 1), class = "honets_bad_input")
  expect_error(knn_hypergraph(two_pairs, k = 0))
  expect_error(knn_hypergraph(two_pairs, k = 4))
  zero_row <- matrix(c(1, 0, 0, 0), nrow = 2,
                     dimnames = list(c("a", "b"), NULL))
  expect_error(knn_hypergraph(zero_row, k = 1), class = "honets_bad_input")
})

test_that("orthogonal neighbors are refused under cosine weighting", {
  ortho <- matrix(c(1, 0, 0, 1), nrow = 2,
                  dimnames = list(c("a", "b"), NULL))
  expect_error(knn_hypergraph(ortho, k = 1),
               class = "honets_nonpositive_similarity")
  hg <- knn_hypergraph(ortho, k = 1, weight = "binary")
  expect_identical(hg$n_hyperedges, 2L)
})

test_that("text_hypergraph knn construction matches knn_hypergraph", {
  corpus <- c(a = "one", b = "two", c = "three", d = "four")
  direct <- knn_hypergraph(two_pairs, k = 2)
  thg <- text_hypergraph(corpus, construction = "knn", k = 2,
                         embeddings = two_pairs)
  expect_identical(thg$incidence, direct$incidence)
  expect_s3_class(thg, "text_hypergraph")
  expect_identical(thg$text$construction, "knn")
  vocab <- as.data.frame(thg, what = "vocabulary")
  expect_identical(nrow(vocab), 0L)
  tab <- as.data.frame(thg)
  expect_identical(names(tab), c("doc", "edge", "weight"))
  expect_identical(nrow(tab), 12L)
})

test_that("embeddings rownames are matched to document IDs in any order", {
  corpus <- c(a = "one", b = "two", c = "three", d = "four")
  scrambled <- two_pairs[c("d", "b", "a", "c"), , drop = FALSE]
  thg <- text_hypergraph(corpus, construction = "knn", k = 1,
                         embeddings = scrambled)
  reference <- text_hypergraph(corpus, construction = "knn", k = 1,
                               embeddings = two_pairs)
  expect_identical(thg$incidence, reference$incidence)
})

test_that("knn construction carries metadata and prints distinctly", {
  articles <- data.frame(
    key = c("a", "b", "c", "d"),
    txt = c("one", "two", "three", "four"),
    year = c(2020L, 2021L, 2022L, 2023L)
  )
  thg <- text_hypergraph(articles, column = "txt", id = "key",
                         construction = "knn", k = 1,
                         embeddings = two_pairs)
  docs <- as.data.frame(thg, what = "documents")
  expect_identical(docs$doc, c("a", "b", "c", "d"))
  expect_identical(docs$year, c(2020L, 2021L, 2022L, 2023L))
  expect_output(print(thg), "kNN embedding hyperedges: k = 1")
})

test_that("knn construction rejects token-layer arguments and bad shapes", {
  corpus <- c(a = "one", b = "two", c = "three", d = "four")
  expect_error(
    text_hypergraph(corpus, construction = "knn", embeddings = two_pairs,
                    stop_words = "the"),
    class = "honets_bad_input"
  )
  expect_error(
    text_hypergraph(corpus, construction = "knn", embeddings = two_pairs,
                    min_count = 2L),
    class = "honets_bad_input"
  )
  wrong_rows <- two_pairs[seq_len(3), , drop = FALSE]
  rownames(wrong_rows) <- NULL
  expect_error(
    text_hypergraph(corpus, construction = "knn", k = 1,
                    embeddings = wrong_rows),
    class = "honets_bad_input"
  )
  misnamed <- two_pairs
  rownames(misnamed) <- c("a", "b", "c", "x")
  expect_error(
    text_hypergraph(corpus, construction = "knn", k = 1,
                    embeddings = misnamed),
    class = "honets_bad_input"
  )
})

test_that("knn hypergraph clusters the two embedding groups", {
  corpus <- c(a = "one", b = "two", c = "three", d = "four")
  thg <- text_hypergraph(corpus, construction = "knn", k = 2,
                         embeddings = two_pairs)
  clusters <- hg_cluster(thg, k = 2, type = "random_walk", seed = 1)
  merged <- merge(clusters, data.frame(node = c("a", "b", "c", "d"),
                                       truth = c("p1", "p1", "p2", "p2")))
  agreement <- aggregate(cluster ~ truth, data = merged,
                         FUN = \(x) length(unique(x)))
  expect_true(all(agreement$cluster == 1L))
  expect_identical(length(unique(merged$cluster)), 2L)
})

test_that("covid_embeddings dataset is intact and aligned", {
  expect_identical(dim(covid_embeddings), c(165L, 384L))
  expect_identical(rownames(covid_embeddings), covid_abstracts$doc)
  expect_false(anyNA(covid_embeddings))
  expect_lt(max(abs(sqrt(rowSums(covid_embeddings^2)) - 1)), 1e-6)
})
