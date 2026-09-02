# HyperGAT (Ding et al. 2020) natively in R via {torch}: inductive text
# classification on DOCUMENT-LEVEL hypergraphs -- each document is its own
# hypergraph with its unique words as vertices and its sentences as
# hyperedges -- through two dual-attention layers, masked mean pooling and
# a linear classifier.
#
# Semantics mirror the official implementation (kaize0409/
# HyperGAT_TextClassification, verified against layers.py/model.py/
# utils.py, which differ from the paper's equations in two places the
# code settles): node-level attention scores depend on the node and a
# GLOBAL trainable context vector (not on the edge), and the default
# construction uses sentence hyperedges only (LDA topic hyperedges are
# behind a flag upstream and are not implemented here; the corresponding
# ablation is "w/o semantic" in the paper's Table 4). Forward parity with
# the official layer is asserted in local_testing_and_equivalence/
# test-equiv-hypergat-official.R. The attention layer computes
# [q || y] %*% a as q %*% a_top + y %*% a_bottom, which avoids the
# official code's 4-D pair tensor without changing the math.

# One dual-attention layer on a dense batch: x [B, N, F_in], adj
# [B, E, N] binary. Returns [B, N, F_out].
.thg_hypergat_layer <- torch::nn_module(
  initialize = function(in_features, out_features, dropout, alpha,
                        transfer, concat) {
    self$out_features <- out_features
    self$p_drop <- dropout
    self$alpha <- alpha
    self$transfer <- transfer
    self$concat <- concat
    stdv <- 1 / sqrt(out_features)
    runif_par <- function(...) {
      torch::nn_parameter(torch::torch_empty(...)$uniform_(-stdv, stdv))
    }
    if (isTRUE(transfer)) {
      self$weight <- runif_par(in_features, out_features)
    }
    self$weight2 <- runif_par(in_features, out_features)
    self$weight3 <- runif_par(out_features, out_features)
    self$a <- runif_par(2L * out_features, 1L)
    self$a2 <- runif_par(2L * out_features, 1L)
    self$word_context <- runif_par(1L, out_features)
  },
  forward = function(x, adj) {
    out_f <- self$out_features
    x_att <- torch::torch_matmul(x, self$weight2)          # [B, N, out]
    if (isTRUE(self$transfer)) {
      x <- torch::torch_matmul(x, self$weight)
    }
    # node-level attention: score([context || x_att_n] a) per node,
    # broadcast over edges, masked softmax over each edge's members
    a_top <- self$a$narrow(1L, 1L, out_f)
    a_bot <- self$a$narrow(1L, out_f + 1L, out_f)
    ctx <- torch::torch_matmul(self$word_context, a_top)$reshape(1L)
    s_node <- torch::nnf_leaky_relu(
      torch::torch_matmul(x_att, a_bot)$squeeze(3L) + ctx, self$alpha
    )                                                       # [B, N]
    e <- s_node$unsqueeze(2L)$expand(c(-1L, adj$shape[2L], -1L))
    e <- torch::nnf_dropout(e, self$p_drop, training = self$training)
    neg <- torch::torch_full_like(e, -9e15)
    att_edge <- torch::nnf_softmax(torch::torch_where(adj > 0, e, neg),
                                   dim = 3L)
    edge <- torch::torch_matmul(att_edge, x)                # [B, E, out]
    edge <- torch::nnf_dropout(edge, self$p_drop, training = self$training)
    # edge-level attention: score([x_att_n || edge_att_e] a2), masked
    # softmax over each node's edges
    edge_att <- torch::torch_matmul(edge, self$weight3)     # [B, E, out]
    a2_top <- self$a2$narrow(1L, 1L, out_f)
    a2_bot <- self$a2$narrow(1L, out_f + 1L, out_f)
    s_pair <- torch::nnf_leaky_relu(
      torch::torch_matmul(x_att, a2_top)$squeeze(3L)$unsqueeze(2L) +
        torch::torch_matmul(edge_att, a2_bot)$squeeze(3L)$unsqueeze(3L),
      self$alpha
    )                                                       # [B, E, N]
    s_pair <- torch::nnf_dropout(s_pair, self$p_drop,
                                 training = self$training)
    neg2 <- torch::torch_full_like(s_pair, -9e15)
    att_node <- torch::nnf_softmax(
      torch::torch_where(adj > 0, s_pair, neg2)$transpose(2L, 3L), dim = 3L
    )                                                       # [B, N, E]
    node <- torch::torch_matmul(att_node, edge)             # [B, N, out]
    if (isTRUE(self$concat)) {
      node <- torch::nnf_elu(node)
    }
    node
  }
)

# The document classifier: embedding -> two dual-attention layers ->
# masked mean pool -> LayerNorm -> linear. Mirrors DocumentGraph +
# HGNN_ATT (embedding, norm and output layer initialized uniform with
# stdv = 1/sqrt(hidden), the official reset_parameters quirk).
.thg_hypergat_net <- torch::nn_module(
  initialize = function(vocab_size, embed_dim, hidden, n_class, dropout,
                        pretrained = NULL) {
    self$p_drop <- dropout
    self$embedding <- torch::nn_embedding(vocab_size + 1L, embed_dim,
                                          padding_idx = 1L)
    self$gat1 <- .thg_hypergat_layer(embed_dim, embed_dim, dropout,
                                     alpha = 0.2, transfer = FALSE,
                                     concat = TRUE)
    self$gat2 <- .thg_hypergat_layer(embed_dim, hidden, dropout,
                                     alpha = 0.2, transfer = TRUE,
                                     concat = FALSE)
    self$norm <- torch::nn_layer_norm(hidden, eps = 1e-6)
    self$out <- torch::nn_linear(hidden, n_class)
    stdv <- 1 / sqrt(hidden)
    torch::with_no_grad({
      self$embedding$weight$uniform_(-stdv, stdv)
      self$out$weight$uniform_(-stdv, stdv)
      self$out$bias$uniform_(-stdv, stdv)
      if (!is.null(pretrained)) {
        self$embedding$weight$copy_(
          torch::torch_tensor(pretrained, dtype = torch::torch_float())
        )
      }
    })
  },
  forward = function(items, adj, mask) {
    h <- self$embedding(items)                              # [B, N, emb]
    h <- self$gat1(h, adj)
    h <- torch::nnf_dropout(h, self$p_drop, training = self$training)
    h <- self$gat2(h, adj)                                  # [B, N, hid]
    m <- mask$unsqueeze(3L)
    pooled <- torch::torch_sum(h * m, dim = 2L) /
      torch::torch_sum(mask, dim = 2L, keepdim = TRUE)
    self$out(self$norm(pooled))
  }
)

# Tokenize documents into sentences of word ids. Returns a list per doc:
# integer-id sentences (pad id is 1; word ids start at 2), plus the vocab.
.thg_hypergat_corpus <- function(text, doc_id, stop_words, min_count,
                                 lowercase) {
  sentences <- lapply(strsplit(text, "[.!?;]+"), \(sents) {
    toks <- .thg_tokenize(sents, lowercase)
    toks <- lapply(toks, \(s) setdiff(s, stop_words))
    toks[lengths(toks) > 0L]
  })
  counts <- table(unlist(sentences))
  vocab <- sort(names(counts)[counts >= min_count])
  sentences <- lapply(sentences, \(doc) {
    kept <- lapply(doc, \(s) match(intersect(s, vocab), vocab) + 1L)
    kept[lengths(kept) > 0L]
  })
  keep <- lengths(sentences) > 0L
  if (!all(keep)) {
    warning(warningCondition(
      sprintf("%d document(s) had no usable tokens and were dropped",
              sum(!keep)),
      class = "honets_dropped_documents"
    ))
  }
  list(sentences = sentences[keep], doc_id = doc_id[keep], vocab = vocab)
}

# Pack a set of documents (list of integer-id sentence lists) into the
# dense batch tensors the layers consume.
.thg_hypergat_batch <- function(docs) {
  node_sets <- lapply(docs, \(doc) sort(unique(unlist(doc))))
  n_nodes <- max(lengths(node_sets))
  n_edges <- max(lengths(docs))
  b <- length(docs)
  items <- matrix(1L, nrow = b, ncol = n_nodes)   # pad id 1
  mask <- matrix(0, nrow = b, ncol = n_nodes)
  adj <- array(0, dim = c(b, n_edges, n_nodes))
  # sequential packing: each document writes its own slab of the batch
  # tensors; vectorizing across documents would need ragged indexing
  for (i in seq_along(docs)) {
    nodes_i <- node_sets[[i]]
    items[i, seq_along(nodes_i)] <- nodes_i
    mask[i, seq_along(nodes_i)] <- 1
    hits <- do.call(rbind, lapply(seq_along(docs[[i]]), \(s)
      cbind(s, match(docs[[i]][[s]], nodes_i))))
    adj[cbind(i, hits[, 1L], hits[, 2L])] <- 1
  }
  list(
    items = torch::torch_tensor(items, dtype = torch::torch_int64()),
    adj = torch::torch_tensor(adj, dtype = torch::torch_float()),
    mask = torch::torch_tensor(mask, dtype = torch::torch_float())
  )
}

#' HyperGAT document classifier
#'
#' Trains the dual-attention hypergraph network of Ding et al. (2020) --
#' every document becomes its own hypergraph (its unique words as
#' vertices, its sentences as hyperedges), two attention layers aggregate
#' words into sentences and sentences back into words, and a masked mean
#' pool feeds a linear classifier. Inductive: only labeled documents are
#' trained on, and every document (labeled or not) is scored. Needs the
#' suggested \pkg{torch} package. Semantics follow the official
#' implementation with sentence hyperedges (the paper's "w/o semantic"
#' construction; LDA topic hyperedges are not implemented).
#'
#' @param x A character vector of documents (names become ids) or a
#'   data.frame with a text column.
#' @param labels The known labels: a named character vector (names are
#'   document ids, values class labels) or a tidy data.frame with a
#'   `node` column and a `label`, `cluster` or `predicted` column. At
#'   least two classes.
#' @param column,id When `x` is a data.frame: the text column and the
#'   optional id column, as in [text_hypergraph()].
#' @param stop_words Words removed before building sentences (default
#'   [stop_words_en()]).
#' @param min_count Minimum corpus frequency for a word to become a
#'   vertex.
#' @param lowercase Lowercase the text first.
#' @param embed_dim,hidden Embedding and hidden width (official defaults
#'   300 and 100).
#' @param epochs,lr,dropout,batch_size,weight_decay Training
#'   hyperparameters (official defaults: 10, 0.001, 0.3, 8, 1e-6).
#' @param lr_decay,lr_step StepLR schedule: multiply the learning rate by
#'   `lr_decay` every `lr_step` epochs (official 0.1 every 3).
#' @param validation Fraction of the labeled documents held out
#'   (stratified) to pick the best epoch; `0` keeps the final weights.
#' @param class_weights `"balanced"` (default, inverse-frequency loss
#'   weights as in the official run script) or `"none"`.
#' @param embeddings Optional pretrained word-vector matrix (rownames are
#'   words, `embed_dim` columns); words not covered keep their random
#'   initialization.
#' @param seed Integer seed (R and torch); results are deterministic
#'   given a seed.
#' @param verbose Message the loss each epoch.
#' @return A base `data.frame`, one row per (kept) document: `node`,
#'   `label` (the given label or `NA`), `predicted`, `score` (softmax
#'   probability of the winning class), `margin` (winner minus runner-up).
#'   The training history is attached as attribute `"history"`
#'   (`epoch`, `loss`, `val_accuracy`).
#' @references
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with
#' less: Hypergraph attention networks for inductive text classification.
#' \emph{EMNLP 2020}.
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE)) {
#'   docs <- c(
#'     cooking_1 = "Simmer the soup. Add onions and carrots.",
#'     cooking_2 = "This soup recipe needs salt. Serve on a cold night.",
#'     space_1 = "The telescope revealed a galaxy. Stars everywhere.",
#'     space_2 = "Astronomers aimed the telescope. The stars were sharp."
#'   )
#'   hg_hypergat(docs, labels = c(cooking_1 = "cooking", space_1 = "space"),
#'               embed_dim = 16, hidden = 8, epochs = 5, validation = 0)
#' }
#' }
#' @export
hg_hypergat <- function(x, labels, column = NULL, id = NULL,
                        stop_words = stop_words_en(), min_count = 1L,
                        lowercase = TRUE, embed_dim = 300L, hidden = 100L,
                        epochs = 10L, lr = 0.001, dropout = 0.3,
                        batch_size = 8L, weight_decay = 1e-6,
                        lr_decay = 0.1, lr_step = 3L, validation = 0.1,
                        class_weights = c("balanced", "none"),
                        embeddings = NULL, seed = 1L, verbose = FALSE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(errorCondition(
      "hg_hypergat() needs the torch package: install.packages(\"torch\")",
      class = "honets_missing_torch", call = NULL
    ))
  }
  class_weights <- match.arg(class_weights)
  labels <- .thg_labels_input(labels)
  stopifnot(
    "`x` must be a character vector or a data.frame" =
      is.character(x) || is.data.frame(x),
    "`labels` must be a named character vector" =
      is.character(labels) && !is.null(names(labels)),
    "`embed_dim` must be a single positive integer" =
      length(embed_dim) == 1L && is.finite(embed_dim) && embed_dim >= 1,
    "`hidden` must be a single positive integer" =
      length(hidden) == 1L && is.finite(hidden) && hidden >= 1,
    "`epochs` must be a single positive integer" =
      length(epochs) == 1L && is.finite(epochs) && epochs >= 1,
    "`batch_size` must be a single positive integer" =
      length(batch_size) == 1L && is.finite(batch_size) && batch_size >= 1,
    "`dropout` must be a single number in [0, 1)" =
      is.numeric(dropout) && length(dropout) == 1L &&
      dropout >= 0 && dropout < 1,
    "`validation` must be a single number in [0, 1)" =
      is.numeric(validation) && length(validation) == 1L &&
      validation >= 0 && validation < 1,
    "`seed` must be a single integer" =
      length(seed) == 1L && is.finite(seed)
  )

  if (is.data.frame(x)) {
    stopifnot(
      "when `x` is a data.frame, `column` must name one of its columns" =
        is.character(column) && length(column) == 1L && column %in% names(x)
    )
    text <- as.character(x[[column]])
    doc_id <- if (is.null(id)) {
      sprintf("doc_%d", seq_len(nrow(x)))
    } else {
      stopifnot("`id` must name a column of `x`" =
                  is.character(id) && length(id) == 1L && id %in% names(x))
      as.character(x[[id]])
    }
  } else {
    text <- x
    doc_id <- names(x) %||% sprintf("doc_%d", seq_along(x))
  }
  stopifnot(
    "document ids must be unique and non-missing" =
      !anyNA(doc_id) && anyDuplicated(doc_id) == 0L
  )

  corpus <- .thg_hypergat_corpus(text, doc_id, stop_words, min_count,
                                 lowercase)
  unknown <- setdiff(names(labels), corpus$doc_id)
  if (length(unknown) > 0L) {
    stop("Unknown or dropped document ids in `labels`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  classes <- sort(unique(as.character(labels)))
  if (length(classes) < 2L) {
    stop(errorCondition(
      "`labels` must contain at least two distinct classes.",
      class = "honets_bad_input", call = NULL
    ))
  }

  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
  on.exit(if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE, after = FALSE)
  set.seed(seed)
  torch::torch_manual_seed(seed)

  pretrained <- NULL
  if (!is.null(embeddings)) {
    stopifnot(
      "`embeddings` needs rownames (words) and `embed_dim` columns" =
        !is.null(rownames(embeddings)) && ncol(embeddings) == embed_dim
    )
    pretrained <- matrix(stats::runif((length(corpus$vocab) + 1L) *
                                        embed_dim,
                                      -1 / sqrt(hidden), 1 / sqrt(hidden)),
                         nrow = length(corpus$vocab) + 1L)
    covered <- intersect(corpus$vocab, rownames(embeddings))
    pretrained[match(covered, corpus$vocab) + 1L, ] <-
      as.matrix(embeddings)[match(covered, rownames(embeddings)), ]
    pretrained[1L, ] <- 0
  }

  model <- .thg_hypergat_net(
    vocab_size = length(corpus$vocab), embed_dim = embed_dim,
    hidden = hidden, n_class = length(classes), dropout = dropout,
    pretrained = pretrained
  )

  target_all <- match(as.character(labels)[match(corpus$doc_id,
                                                 names(labels))], classes)
  labeled <- which(!is.na(target_all))
  val_idx <- integer(0)
  if (validation > 0) {
    # stratified holdout; classes too small to split stay fully in train
    val_pick <- unlist(lapply(split(labeled, target_all[labeled]), \(i) {
      k <- floor(validation * length(i))
      if (k >= 1L) sample(i, k) else integer(0)
    }))
    val_idx <- sort(val_pick)
  }
  train_idx <- setdiff(labeled, val_idx)

  loss_w <- if (identical(class_weights, "balanced")) {
    counts <- tabulate(target_all[train_idx], nbins = length(classes))
    torch::torch_tensor(length(train_idx) / (length(classes) *
                                               pmax(counts, 1L)),
                        dtype = torch::torch_float())
  } else {
    NULL
  }

  score_docs <- function(idx) {
    # sequential mini-batching over a fixed document list (eval mode)
    starts <- seq(1L, length(idx), by = 16L)
    probs <- lapply(starts, \(s) {
      chunk <- idx[s:min(s + 15L, length(idx))]
      batch <- .thg_hypergat_batch(corpus$sentences[chunk])
      torch::with_no_grad(
        as.matrix(torch::nnf_softmax(
          model(batch$items, batch$adj, batch$mask), dim = 2L
        ))
      )
    })
    do.call(rbind, probs)
  }

  optimizer <- torch::optim_adam(model$parameters, lr = lr,
                                 weight_decay = weight_decay)
  scheduler <- torch::lr_step(optimizer, step_size = lr_step,
                              gamma = lr_decay)
  best_val <- -Inf
  best_state <- NULL
  history <- vapply(seq_len(epochs), \(epoch) {
    # sequential training: epochs and mini-batches are ordered by Adam
    scheduler$step()
    model$train()
    order_idx <- sample(train_idx)
    starts <- seq(1L, length(order_idx), by = batch_size)
    losses <- vapply(starts, \(s) {
      chunk <- order_idx[s:min(s + batch_size - 1L, length(order_idx))]
      batch <- .thg_hypergat_batch(corpus$sentences[chunk])
      optimizer$zero_grad()
      out <- model(batch$items, batch$adj, batch$mask)
      y <- torch::torch_tensor(target_all[chunk],
                               dtype = torch::torch_int64())
      loss <- torch::nnf_cross_entropy(out, y, weight = loss_w)
      loss$backward()
      optimizer$step()
      as.numeric(loss)
    }, numeric(1))
    val_acc <- NA_real_
    if (length(val_idx) > 0L) {
      model$eval()
      val_probs <- score_docs(val_idx)
      val_pred <- max.col(val_probs, ties.method = "first")
      val_acc <- mean(val_pred == target_all[val_idx])
      if (val_acc > best_val) {
        best_val <<- val_acc
        best_state <<- lapply(model$state_dict(), \(p) p$clone())
      }
    }
    if (isTRUE(verbose)) {
      message(sprintf("epoch %d loss %.4f val %.3f", epoch,
                      mean(losses), val_acc))
    }
    c(mean(losses), val_acc)
  }, numeric(2))
  if (!is.null(best_state)) {
    model$load_state_dict(best_state)
  }

  model$eval()
  probs <- score_docs(seq_along(corpus$doc_id))
  dimnames(probs) <- list(corpus$doc_id, classes)
  lab_full <- rep(NA_character_, length(corpus$doc_id))
  lab_full[labeled] <- classes[target_all[labeled]]
  out <- .hl_score_predictions(probs, lab_full, "none")
  rownames(out) <- NULL
  attr(out, "history") <- data.frame(
    epoch = seq_len(epochs), loss = history[1L, ],
    val_accuracy = history[2L, ]
  )
  out
}

# The model consumes raw text rather than a pre-built hypergraph, so its
# descriptive alias belongs to the text family.
#' @rdname hg_hypergat
#' @export
text_hypergat <- hg_hypergat
