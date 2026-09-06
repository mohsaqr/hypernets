# Sparse core: parity with the dense engines (the in-package oracle chain:
# dense == HyperNetX), plus scale behavior and guards.

sparse_corpus <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
)
sparse_stop <- c("the", "with", "and", "a", "this", "at", "on", "all")

both <- function(...) {
  list(
    dense = text_hypergraph(sparse_corpus, stop_words = sparse_stop, ...),
    sparse = text_hypergraph(sparse_corpus, stop_words = sparse_stop,
                             sparse = TRUE, ...)
  )
}

test_that("sparse construction reproduces the dense incidence exactly", {
  for (w in c("n", "tfidf")) {
    for (nd in c("doc", "word")) {
      hg <- both(weight = w, nodes = nd)
      expect_s4_class(hg$sparse$incidence, "sparseMatrix")
      expect_identical(as.matrix(hg$sparse$incidence), hg$dense$incidence)
      expect_identical(hg$sparse$nodes, hg$dense$nodes)
      expect_identical(hg$sparse$size_distribution,
                       hg$dense$size_distribution)
      sparse_table <- as.data.frame(hg$sparse)
      dense_table <- as.data.frame(hg$dense)
      expect_identical(sparse_table, dense_table)
    }
  }
})

test_that("sparse measures match the dense engine", {
  hg <- both(weight = "tfidf")
  for (what in c("nodes", "edges", "summary", "overlap")) {
    sparse_measures <- hg_measures(hg$sparse, what = what)
    dense_measures <- hg_measures(hg$dense, what = what)
    expect_equal(sparse_measures, dense_measures, tolerance = 1e-12)
  }
})

test_that("sparse pagerank matches dense to machine precision", {
  hg <- both(weight = "tfidf")
  for (d in c(0.85, 1)) {
    sparse_pr <- hg_pagerank(hg$sparse, damping = d, tol = 1e-14,
                             max_iter = 20000L)
    dense_pr <- hg_pagerank(hg$dense, damping = d, tol = 1e-14,
                            max_iter = 20000L)
    expect_equal(sparse_pr, dense_pr, tolerance = 1e-10)
  }
  sparse_personalized <- hg_pagerank(hg$sparse, personalized = c(cooking_1 = 1))
  dense_personalized <- hg_pagerank(hg$dense, personalized = c(cooking_1 = 1))
  expect_equal(sparse_personalized, dense_personalized, tolerance = 1e-10)
})

test_that("sparse transduction (CG) matches the dense closed form", {
  hg <- both(weight = "tfidf")
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  for (type in c("zhou", "random_walk")) {
    sparse_fit <- hg_classify(hg$sparse, labels = labels, type = type)
    dense_fit <- hg_classify(hg$dense, labels = labels, type = type)
    expect_identical(sparse_fit$predicted, dense_fit$predicted)
    expect_equal(sparse_fit$score, dense_fit$score, tolerance = 1e-8)
    expect_equal(sparse_fit$margin, dense_fit$margin, tolerance = 1e-8)
  }
})

test_that("sparse class_mass transduction matches the dense engine", {
  hg <- both(weight = "tfidf")
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  sparse_fit <- hg_classify(hg$sparse, labels = labels,
                            normalization = "class_mass")
  dense_fit <- hg_classify(hg$dense, labels = labels,
                           normalization = "class_mass")
  expect_identical(sparse_fit$predicted, dense_fit$predicted)
  expect_equal(sparse_fit$score, dense_fit$score, tolerance = 1e-8)
  expect_equal(sparse_fit$margin, dense_fit$margin, tolerance = 1e-8)
})

test_that("sparse clustering recovers the same planted partition", {
  events <- data.frame(
    person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
               "d", "e", "f", "c", "d"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
                "m4", "m4", "m4", "m5", "m5"),
    w = 1
  )
  dense_hg <- group_hypergraph(events, actor = "person",
                                          group = "meeting", weight = "w")
  sparse_hg <- honets:::.thg_sparse_bipartite(
    events, actor = "person", group = "meeting", weight = "w"
  )
  dense_cl <- hg_cluster(dense_hg, k = 2, seed = 1)
  sparse_cl <- hg_cluster(sparse_hg, k = 2, seed = 1)
  expect_identical(sparse_cl, dense_cl)

  fit <- honets:::.thg_sparse_cluster(
    sparse_hg, k = 2, type = "zhou", edge_weights = NULL, nstart = 25L,
    seed = 1
  )
  dense_fit <- hypergraph_cluster(dense_hg, k = 2, seed = 1)
  expect_equal(fit$eigenvalues[seq_len(3)],
               dense_fit$eigenvalues[seq_len(3)], tolerance = 1e-8)
  expect_equal(fit$pi, dense_fit$pi, tolerance = 1e-10)
})

test_that("sparse duals round-trip and stay sparse", {
  hg <- both(weight = "tfidf")
  dual <- dual_hypergraph(hg$sparse)
  expect_s4_class(dual$incidence, "sparseMatrix")
  expect_identical(as.matrix(dual$incidence), t(hg$dense$incidence))
  back <- dual_hypergraph(dual)
  expect_identical(as.matrix(back$incidence), hg$dense$incidence)
})

test_that("unsupported sparse paths refuse with classed errors", {
  hg <- both()
  expect_error(hg_centrality(hg$sparse), class = "honets_sparse_unsupported")
  expect_error(hg_null_test(hg$sparse), class = "honets_sparse_unsupported")
  expect_error(
    text_hypergraph(sparse_corpus, construction = "window", sparse = TRUE),
    class = "honets_bad_input"
  )
  big_long <- data.frame(v = rep(c("x", "y"), each = 2100),
                         e = paste0("e", c(seq_len(2100), seq_len(2100))),
                         w = 1)
  big <- honets:::.thg_sparse_bipartite(big_long, actor = "v",
                                                group = "e", weight = "w")
  expect_error(hg_measures(big, what = "overlap"),
               class = "honets_sparse_too_large")
})

test_that("sparse scale: thousands of documents classify in seconds", {
  skip_on_cran()
  set.seed(42)
  # letters-only tokens: the tokenizer treats digits as separators
  vocab <- head(apply(expand.grid(letters, letters, letters), 1, paste0,
                      collapse = ""), 1405)
  # two disjoint vocabulary blocks; only the first 100 documents carry two
  # tokens from a tiny shared bridge pool, so the corpus is connected but
  # the bridge hyperedges stay small (the zhou Laplacian is binary -- a
  # bridge word present in every document would wash all labels into ties)
  make_doc <- \(i) {
    pool <- if (i %% 2L == 0L) seq_len(700) else 701:1400
    body <- sample(vocab[pool], 30, replace = TRUE)
    if (i <= 100L) {
      body <- c(body, sample(vocab[1401:1405], 2))
    }
    paste(body, collapse = " ")
  }
  docs <- vapply(seq_len(3000), make_doc, character(1))
  names(docs) <- sprintf("d%04d", seq_len(3000))
  elapsed <- system.time({
    hg <- text_hypergraph(docs, weight = "tfidf", sparse = TRUE)
    labels <- c(d0002 = "even", d0004 = "even", d0006 = "even",
                d0001 = "odd", d0003 = "odd", d0005 = "odd")
    fit <- hg_classify(hg, labels = labels, type = "zhou")
    pr <- hg_pagerank(hg)
  })[["elapsed"]]
  expect_identical(nrow(fit), 3000L)
  expect_identical(nrow(pr), 3000L)
  expect_lt(elapsed, 60)
  # the two vocabulary populations are recovered from one label each
  truth <- rep(c("odd", "even"), length.out = 3000)
  accuracy <- mean(fit$predicted == truth)
  expect_gt(accuracy, 0.95)
})
