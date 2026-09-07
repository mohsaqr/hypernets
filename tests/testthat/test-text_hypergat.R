# HyperGAT: the layer is checked against an independent plain-R
# re-implementation of the official math (double precision), the verb
# end-to-end for shape, determinism, sentence-set invariance and fit.
# Forward parity with the official PyTorch layer lives in
# local_testing_and_equivalence/test-equiv-hypergat-official.R.

hypergat_docs <- c(
  cooking_1 = "Simmer the soup slowly. Add onions and carrots to the pot.",
  cooking_2 = "This soup recipe needs salt. Serve the dish on a cold night.",
  cooking_3 = "Chop onions for the stew. Salt the broth and simmer.",
  space_1 = "The telescope revealed a distant galaxy. Stars everywhere.",
  space_2 = "Astronomers aimed the telescope. The stars were sharp.",
  space_3 = "A comet passed the galaxy. Astronomers watched the stars."
)
hypergat_labels <- c(cooking_1 = "cooking", cooking_2 = "cooking",
                     space_1 = "space", space_2 = "space")

test_that("native online LDA is deterministic and uses training documents only", {
  corpus <- .thg_hypergat_corpus(
    hypergat_docs, names(hypergat_docs), stop_words_en(), 1L, TRUE
  )
  args <- list(
    sentences = corpus$sentences, vocab = corpus$vocab,
    train_idx = match(names(hypergat_labels), corpus$doc_id),
    n_topics = 2L, top_n = 4L, max_iter = 2L, batch_size = 2L,
    seed = 0L
  )
  a <- do.call(.thg_hypergat_lda, args)
  b <- do.call(.thg_hypergat_lda, args)
  expect_equal(a$components, b$components, tolerance = 0)
  expect_identical(a$keywords, b$keywords)
  expect_identical(a$n_topics, 2L)
  expect_identical(a$top_n, 4L)
  expect_identical(a$training_documents, 4L)
  expect_true(all(lengths(a$keywords) == 4L))
  expect_true(all(unlist(a$keywords) %in% corpus$vocab))
})

test_that("semantic topics append official get_slice-style document edges", {
  vocab <- c("alpha", "beta", "gamma", "delta")
  sentences <- list(
    list(c(2L, 3L), c(3L, 4L)),
    list(c(3L, 5L))
  )
  keywords <- list(topic_1 = c("alpha", "gamma"),
                   topic_2 = c("beta", "delta"))
  docs <- .thg_hypergat_semantic_docs(sentences, vocab, keywords)
  expect_identical(docs[[1]][1:2], sentences[[1]])
  expect_identical(docs[[1]][[3]], c(2L, 4L))
  expect_identical(docs[[1]][[4]], 3L)
  expect_identical(docs[[2]][[2]], integer(0))
  expect_identical(docs[[2]][[3]], c(3L, 5L))
  if (requireNamespace("torch", quietly = TRUE)) {
    packed <- .thg_hypergat_batch(docs)
    adj <- as.array(packed$adj)
    expect_equal(adj[1, 1:4, 1:3], rbind(
      c(1, 1, 0), c(0, 1, 1), c(1, 0, 1), c(0, 1, 0)
    ))
    # The empty first semantic topic in document 2 remains an empty edge.
    expect_equal(adj[2, 2, ], c(0, 0, 0))
  }
})

# Plain-R double-precision reference for one dual-attention layer
# (official math, B = 1): an oracle independent of torch.
.ref_hypergat_layer <- function(x, adj, w2, w3, a, a2, ctx, w = NULL,
                                alpha = 0.2, concat = TRUE) {
  lrelu <- \(z) ifelse(z > 0, z, alpha * z)
  elu <- \(z) ifelse(z > 0, z, exp(z) - 1)
  softmax_rows <- \(m) t(apply(m, 1, \(r) exp(r - max(r)) /
                                 sum(exp(r - max(r)))))
  out_f <- ncol(w2)
  x4 <- x %*% w2
  val <- if (is.null(w)) x else x %*% w
  a_top <- a[seq_len(out_f), , drop = FALSE]
  a_bot <- a[out_f + seq_len(out_f), , drop = FALSE]
  s_node <- lrelu(as.numeric(x4 %*% a_bot) + as.numeric(ctx %*% a_top))
  e_scores <- matrix(s_node, nrow = nrow(adj), ncol = nrow(x),
                     byrow = TRUE)
  att_edge <- softmax_rows(ifelse(adj > 0, e_scores, -9e15))
  edge <- att_edge %*% val
  e4 <- edge %*% w3
  a2_top <- a2[seq_len(out_f), , drop = FALSE]
  a2_bot <- a2[out_f + seq_len(out_f), , drop = FALSE]
  s_pair <- lrelu(outer(as.numeric(e4 %*% a2_bot),
                        rep(1, nrow(x))) +
                    outer(rep(1, nrow(adj)), as.numeric(x4 %*% a2_top)))
  att_node <- softmax_rows(t(ifelse(adj > 0, s_pair, -9e15)))
  node <- att_node %*% edge
  if (isTRUE(concat)) elu(node) else node
}

test_that("the torch layer matches the plain-R reference math", {
  skip_if_not_installed("torch")
  set.seed(3)
  n <- 5L; e <- 3L; in_f <- 4L; out_f <- 4L
  x <- matrix(rnorm(n * in_f), n)
  adj <- matrix(0, e, n)
  adj[1, c(1, 2, 3)] <- 1
  adj[2, c(3, 4)] <- 1
  adj[3, c(1, 5)] <- 1
  w2 <- matrix(rnorm(in_f * out_f), in_f)
  w3 <- matrix(rnorm(out_f * out_f), out_f)
  a <- matrix(rnorm(2 * out_f), ncol = 1)
  a2 <- matrix(rnorm(2 * out_f), ncol = 1)
  ctx <- matrix(rnorm(out_f), nrow = 1)
  ref <- .ref_hypergat_layer(x, adj, w2, w3, a, a2, ctx, concat = TRUE)

  layer <- .thg_hypergat_layer(in_f, out_f, dropout = 0, alpha = 0.2,
                               transfer = FALSE, concat = TRUE)
  tt <- \(m) torch::torch_tensor(m, dtype = torch::torch_float())
  torch::with_no_grad({
    layer$weight2$copy_(tt(w2))
    layer$weight3$copy_(tt(w3))
    layer$a$copy_(tt(a))
    layer$a2$copy_(tt(a2))
    layer$word_context$copy_(tt(ctx))
  })
  layer$eval()
  got <- torch::with_no_grad(
    as.array(layer(tt(array(x, c(1, n, in_f))),
                   tt(array(adj, c(1, e, n)))))
  )
  expect_equal(got[1, , ], ref, tolerance = 1e-5)
})

test_that("hg_hypergat trains, predicts every document, deterministic", {
  skip_if_not_installed("torch")
  fit <- hg_hypergat(hypergat_docs, labels = hypergat_labels,
                     embed_dim = 16, hidden = 8, epochs = 30, lr = 0.05,
                     validation = 0, seed = 1)
  expect_s3_class(fit, "data.frame")
  expect_identical(nrow(fit), 6L)
  expect_named(fit, c("node", "label", "predicted", "score", "margin"))
  # enough capacity and epochs to fit the four labeled documents
  labeled <- subset(fit, !is.na(label))
  expect_identical(labeled$predicted, unname(hypergat_labels[labeled$node]))
  history <- attr(fit, "history")
  expect_identical(nrow(history), 30L)
  expect_lt(history$loss[30L], history$loss[1L])
  refit <- hg_hypergat(hypergat_docs, labels = hypergat_labels,
                       embed_dim = 16, hidden = 8, epochs = 30, lr = 0.05,
                       validation = 0, seed = 1)
  expect_identical(fit$predicted, refit$predicted)
  expect_equal(fit$score, refit$score, tolerance = 1e-12)
})

test_that("hg_hypergat exposes the full LDA semantic-hyperedge path", {
  skip_if_not_installed("torch")
  fit <- hg_hypergat(
    hypergat_docs, labels = hypergat_labels, semantic = "lda",
    lda_topics = 2, lda_top_n = 3, lda_max_iter = 2,
    embed_dim = 8, hidden = 4, epochs = 2, validation = 0, seed = 2
  )
  info <- attr(fit, "semantic")
  expect_identical(info$method, "lda")
  expect_identical(info$n_topics, 2L)
  expect_identical(info$top_n, 3L)
  expect_identical(info$training_documents, 4L)
  expect_length(info$keywords, 2L)

  supplied <- list(food = c("soup", "onions", "salt"),
                   sky = c("stars", "galaxy", "telescope"))
  replay <- hg_hypergat(
    hypergat_docs, labels = hypergat_labels, semantic = "lda",
    lda_keywords = supplied, embed_dim = 8, hidden = 4, epochs = 1,
    validation = 0, seed = 2
  )
  replay_info <- attr(replay, "semantic")
  expect_identical(replay_info$method, "precomputed")
  expect_identical(replay_info$keywords, supplied)
})

test_that("predictions are invariant to word order within sentences", {
  skip_if_not_installed("torch")
  shuffled <- c(
    cooking_1 = "the Simmer slowly soup. carrots and onions Add pot to the.",
    cooking_2 = "salt needs recipe soup This. cold a on dish the Serve night.",
    cooking_3 = "the for stew onions Chop. simmer and broth the Salt.",
    space_1 = "galaxy distant a revealed telescope The. everywhere Stars.",
    space_2 = "telescope the aimed Astronomers. sharp were stars The.",
    space_3 = "galaxy the passed comet A. stars the watched Astronomers."
  )
  fit <- hg_hypergat(hypergat_docs, labels = hypergat_labels,
                     embed_dim = 16, hidden = 8, epochs = 10,
                     validation = 0, seed = 1)
  refit <- hg_hypergat(shuffled, labels = hypergat_labels,
                       embed_dim = 16, hidden = 8, epochs = 10,
                       validation = 0, seed = 1)
  expect_identical(fit$predicted, refit$predicted)
  expect_equal(fit$score, refit$score, tolerance = 1e-5)
})

test_that("hg_hypergat argument contracts are enforced", {
  skip_if_not_installed("torch")
  expect_error(
    hg_hypergat(hypergat_docs, labels = c(cooking_1 = "x")),
    class = "hypernets_bad_input"
  )
  expect_error(
    hg_hypergat(hypergat_docs, labels = c(zz = "x", cooking_1 = "y")),
    "Unknown or dropped"
  )
  expect_error(
    hg_hypergat(data.frame(txt = hypergat_docs),
                labels = hypergat_labels),
    "`column` must name"
  )
  expect_warning(
    hg_hypergat(c(hypergat_docs, empty_1 = "the and of"),
                labels = hypergat_labels, embed_dim = 8, hidden = 4,
                epochs = 2, validation = 0),
    class = "hypernets_dropped_documents"
  )
})

test_that("hg_hypergat(what = 'attention') returns per-edge-normalised word attention", {
  skip_if_not_installed("torch")
  att <- hg_hypergat(hypergat_docs, labels = hypergat_labels,
                     embed_dim = 16, hidden = 8, epochs = 5, lr = 0.05,
                     validation = 0, seed = 1, what = "attention")
  expect_named(att, c("node", "word", "attention", "attention_2", "n_edges"))
  expect_true(all(att$attention > 0 & att$attention <= att$n_edges))
  expect_true(all(att$attention_2 > 0 & att$attention_2 <= att$n_edges))
  expect_true(all(att$n_edges >= 1L))
  expect_false(anyDuplicated(att[, c("node", "word")]) > 0)
  expect_identical(att, att[order(att$node, att$word), ])
  # each hyperedge's node-level weights sum to one, so a document's total
  # attention equals its hyperedge (sentence) count
  sentences <- lengths(strsplit(hypergat_docs, "[.!?;]+"))
  totals <- tapply(att$attention, att$node, sum)
  expect_equal(as.numeric(totals[names(sentences)]),
               as.numeric(sentences), tolerance = 1e-6)
  totals_2 <- tapply(att$attention_2, att$node, sum)
  expect_equal(as.numeric(totals_2[names(sentences)]),
               as.numeric(sentences), tolerance = 1e-6)
  # structural: in a one-sentence document every word's layer-1 output is
  # the same edge vector, so layer-2 attention is exactly uniform there,
  # while layer-1 attention (over the word embeddings) is not
  one <- hg_hypergat(c(a = "simmer the soup with onions and carrots",
                       b = "the telescope revealed a distant galaxy"),
                     labels = c(a = "x", b = "y"), embed_dim = 16, hidden = 8,
                     epochs = 5, lr = 0.05, validation = 0, seed = 1,
                     what = "attention")
  uniform <- 1 / as.numeric(table(one$node)[one$node])
  expect_equal(one$attention_2, uniform, tolerance = 1e-6)
  expect_false(isTRUE(all.equal(one$attention, uniform, tolerance = 1e-3)))
  # deterministic under a seed
  again <- hg_hypergat(hypergat_docs, labels = hypergat_labels,
                       embed_dim = 16, hidden = 8, epochs = 5, lr = 0.05,
                       validation = 0, seed = 1, what = "attention")
  expect_equal(att, again)
  # feeds hg_keywords(type = "attention")
  hg <- text_hypergraph(hypergat_docs, stop_words = stop_words_en())
  kw <- hg_keywords(hg, hypergat_labels, n = 3, scores = att)
  expect_named(kw, c("type", "cluster", "size", "rank", "word", "score",
                     "share", "n_docs"))
  expect_identical(unique(kw$type), "attention")
  expect_setequal(unique(kw$cluster), unique(hypergat_labels))
})
