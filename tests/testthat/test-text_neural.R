# HGNN neural tier: the propagation operator is verified against the
# oracle-verified Zhou spectral core (machine precision), the factorization
# against the hand formula, and training end-to-end for shape, determinism
# and fit. Forward-pass parity with the official implementation (DHG) lives
# in local_testing_and_equivalence/test-equiv-hgnn-dhg.R.

neural_corpus <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
)
neural_stop <- c("the", "with", "and", "a", "this", "at", "on", "all")

test_that("the HGNN propagation equals the Zhou similarity operator", {
  hg <- text_hypergraph(neural_corpus, stop_words = neural_stop)
  m <- hg$n_hyperedges
  w <- rep(1, m)
  parts <- .thg_hgnn_factor(hg, edge_weights = w)
  g <- as.matrix(parts$a %*% Matrix::t(parts$a))
  zhou <- .hl_build(hg, type = "zhou", edge_weights = w)
  s <- diag(hg$n_nodes) - (zhou$L + t(zhou$L)) / 2
  expect_equal(unname(g), unname(s), tolerance = 1e-12)
})

test_that("the factorization matches the hand formula, weighted", {
  events <- data.frame(
    person = c("a", "b", "c", "a", "b", "c", "d"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m3", "m3"),
    w = 1
  )
  hg <- group_hypergraph(events, actor = "person",
                                    group = "meeting", weight = "w")
  w <- c(2, 1, 3)
  parts <- .thg_hgnn_factor(hg, edge_weights = w)
  h <- (hg$incidence != 0) * 1
  de <- colSums(h)
  dv <- as.numeric(h %*% w)
  g_hand <- diag(1 / sqrt(dv)) %*% h %*% diag(w) %*% diag(1 / de) %*%
    t(h) %*% diag(1 / sqrt(dv))
  expect_equal(unname(as.matrix(parts$a %*% Matrix::t(parts$a))),
               unname(g_hand), tolerance = 1e-12)
})

test_that("hg_neural trains, predicts, and is deterministic", {
  skip_if_not_installed("torch")
  hg <- text_hypergraph(neural_corpus, stop_words = neural_stop)
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  fit <- hg_neural(hg, labels = labels, hidden = 8, epochs = 80,
                   validation = 0, seed = 1)
  expect_s3_class(fit, "data.frame")
  expect_identical(nrow(fit), 4L)
  expect_named(fit, c("node", "label", "predicted", "score", "margin"))
  expect_identical(subset(fit, node == "cooking_2")$predicted, "cooking")
  expect_identical(subset(fit, node == "space_2")$predicted, "space")
  history <- attr(fit, "history")
  expect_identical(nrow(history), 80L)
  expect_true(all(is.finite(history$loss)))
  expect_lt(history$loss[80L], history$loss[1L])
  refit <- hg_neural(hg, labels = labels, hidden = 8, epochs = 80,
                     validation = 0, seed = 1)
  expect_identical(fit$predicted, refit$predicted)
  expect_equal(fit$score, refit$score, tolerance = 1e-12)
})

test_that("sparse hypergraphs and explicit feature matrices work", {
  skip_if_not_installed("torch")
  hg_sparse <- text_hypergraph(neural_corpus, stop_words = neural_stop,
                               sparse = TRUE)
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  fit_sparse <- hg_neural(hg_sparse, labels = labels, hidden = 8,
                          epochs = 80, validation = 0, seed = 1)
  expect_identical(nrow(fit_sparse), 4L)
  expect_identical(subset(fit_sparse, node == "cooking_2")$predicted,
                   "cooking")

  hg <- text_hypergraph(neural_corpus, stop_words = neural_stop)
  feats <- diag(4)
  rownames(feats) <- hg$nodes
  fit_feats <- hg_neural(hg, labels = labels, features = feats,
                         hidden = 8, epochs = 80, validation = 0, seed = 1)
  expect_identical(nrow(fit_feats), 4L)
})

test_that("hg_neural argument contracts are enforced", {
  skip_if_not_installed("torch")
  hg <- text_hypergraph(neural_corpus, stop_words = neural_stop)
  expect_error(hg_neural(hg, labels = c(cooking_1 = "x")),
               class = "honets_bad_input")
  expect_error(hg_neural(hg, labels = c(zz = "x", cooking_1 = "y")),
               "Unknown node names")
  bad_feats <- diag(4)
  expect_error(
    hg_neural(hg, labels = c(cooking_1 = "x", space_1 = "y"),
              features = bad_feats),
    "rownames"
  )
  expect_error(
    hg_neural(hg, labels = c(cooking_1 = "x", space_1 = "y"),
              edge_weights = c(1, 2)),
    "one positive number per hyperedge"
  )
})
