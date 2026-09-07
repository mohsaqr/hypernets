# ---- HyperGCN (Yadati et al. 2019) -------------------------------------

# Shared validation and training utilities for node-classification models.
.thg_neural_problem <- function(hg, labels, features) {
  .thg_check_hg(hg)
  labels <- .thg_labels_input(labels)
  stopifnot("`labels` must be a named character vector" =
              is.character(labels) && !is.null(names(labels)))
  unknown <- setdiff(names(labels), hg$nodes)
  if (length(unknown)) {
    stop("Unknown node names in `labels`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  classes <- sort(unique(as.character(labels)))
  if (length(classes) < 2L) {
    stop(errorCondition("`labels` must contain at least two distinct classes.",
                        class = "hypernets_bad_input", call = NULL))
  }
  x <- if (identical(features, "incidence")) {
    hg$incidence
  } else {
    stopifnot(
      "`features` must be \"incidence\" or a numeric matrix" =
        (is.matrix(features) && is.numeric(features)) ||
          methods::is(features, "sparseMatrix"),
      "`features` needs rownames matching the hypergraph nodes" =
        !is.null(rownames(features)) && all(hg$nodes %in% rownames(features))
    )
    features[match(hg$nodes, rownames(features)), , drop = FALSE]
  }
  list(labels = labels, classes = classes, features = x)
}

.thg_neural_fit <- function(make_model, forward, nodes, labels, classes,
                            epochs, lr, weight_decay, validation, seed,
                            verbose = FALSE) {
  stopifnot(
    "`epochs` must be a positive whole number" =
      is.numeric(epochs) && length(epochs) == 1L && is.finite(epochs) &&
        epochs >= 1 && epochs == round(epochs),
    "`validation` must be in [0, 1)" =
      is.numeric(validation) && length(validation) == 1L &&
        validation >= 0 && validation < 1,
    "`seed` must be a finite whole number" =
      is.numeric(seed) && length(seed) == 1L && is.finite(seed) &&
        seed == round(seed)
  )
  had_seed <- exists(".Random.seed", envir = globalenv())
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit(if (had_seed) {
    assign(".Random.seed", old_seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv())) {
    rm(".Random.seed", envir = globalenv())
  }, add = TRUE, after = FALSE)
  set.seed(as.integer(seed))
  torch::torch_manual_seed(as.integer(seed))
  model <- make_model()
  seed_idx <- match(names(labels), nodes)
  target <- match(as.character(labels), classes)
  val_idx <- integer(0)
  if (validation > 0) {
    val_idx <- sort(unlist(lapply(split(seq_along(target), target), function(i) {
      k <- floor(validation * length(i))
      if (k >= 1L) sample(i, k) else integer(0)
    })))
  }
  train_pos <- setdiff(seq_along(target), val_idx)
  train_t <- torch::torch_tensor(seed_idx[train_pos],
                                 dtype = torch::torch_int64())
  y_train <- torch::torch_tensor(target[train_pos],
                                 dtype = torch::torch_int64())
  optimizer <- torch::optim_adam(model$parameters, lr = lr,
                                 weight_decay = weight_decay)
  best_val <- -Inf
  best_state <- NULL
  history <- vapply(seq_len(as.integer(epochs)), function(epoch) {
    model$train()
    optimizer$zero_grad()
    out <- forward(model)
    loss <- torch::nnf_cross_entropy(out[train_t, ], y_train)
    loss$backward()
    optimizer$step()
    val_acc <- NA_real_
    if (length(val_idx)) {
      model$eval()
      pred <- torch::with_no_grad(
        as.integer(torch::torch_argmax(forward(model), dim = 2L))
      )
      val_acc <- mean(pred[seed_idx[val_idx]] == target[val_idx])
      if (val_acc > best_val) {
        best_val <<- val_acc
        best_state <<- lapply(model$state_dict(), function(p) p$clone())
      }
    }
    if (isTRUE(verbose) && epoch %% 20L == 0L) {
      message(sprintf("epoch %d loss %.4f", epoch, as.numeric(loss)))
    }
    c(as.numeric(loss), val_acc)
  }, numeric(2L))
  if (!is.null(best_state)) model$load_state_dict(best_state)
  model$eval()
  probs <- torch::with_no_grad(
    as.matrix(torch::nnf_softmax(forward(model), dim = 2L))
  )
  dimnames(probs) <- list(nodes, classes)
  lab_full <- rep(NA_character_, length(nodes))
  lab_full[seed_idx] <- as.character(labels)
  out <- .hl_score_predictions(probs, lab_full, "none")
  rownames(out) <- NULL
  attr(out, "history") <- data.frame(
    epoch = seq_len(as.integer(epochs)), loss = history[1L, ],
    val_accuracy = history[2L, ]
  )
  out
}

# Mediator graph from Algorithm 1/2. The farthest pair is selected by the
# transformed signal; deterministic lexicographic order resolves exact ties.
.thg_hypergcn_adjacency <- function(hg, signal, edge_weights = NULL,
                                    mediators = TRUE) {
  n <- hg$n_nodes
  stopifnot("`signal` must have one finite row per node" =
              (is.matrix(signal) || is.data.frame(signal)) &&
              nrow(signal) == n && all(is.finite(as.matrix(signal))))
  w <- edge_weights %||% rep(1, hg$n_hyperedges)
  stopifnot("`edge_weights` must be positive, one per hyperedge" =
              is.numeric(w) && length(w) == hg$n_hyperedges &&
              all(is.finite(w)) && all(w > 0))
  H <- hg$incidence != 0
  A <- diag(1, n)
  for (edge in seq_len(hg$n_hyperedges)) {
    members <- which(H[, edge])
    if (length(members) < 2L) next
    pairs <- utils::combn(members, 2L)
    distances <- rowSums((as.matrix(signal)[pairs[1L, ], , drop = FALSE] -
                            as.matrix(signal)[pairs[2L, ], , drop = FALSE])^2)
    pair <- pairs[, which.max(distances)]
    others <- setdiff(members, pair)
    links <- if (isTRUE(mediators) && length(others)) {
      cbind(c(pair[1L], rep(pair[1L], length(others)),
              rep(pair[2L], length(others))),
            c(pair[2L], others, others))
    } else {
      matrix(pair, nrow = 1L)
    }
    contribution <- if (isTRUE(mediators)) {
      w[edge] / (2 * length(members) - 3)
    } else w[edge]
    for (link in seq_len(nrow(links))) {
      i <- links[link, 1L]
      j <- links[link, 2L]
      A[i, j] <- A[i, j] + contribution
      A[j, i] <- A[j, i] + contribution
    }
  }
  degree <- rowSums(A)
  A <- A / sqrt(outer(degree, degree))
  dimnames(A) <- list(hg$nodes, hg$nodes)
  A
}

#' HyperGCN node classifier
#'
#' Implements Yadati et al.'s (2019) mediator expansion and two-layer GCN.
#' Every hyperedge selects the pair of vertices farthest apart in the current
#' transformed representation. With `method = "hypergcn"`, the expansion is
#' recomputed at both layers and every epoch (Algorithm 1). `"fast"` fixes it
#' once from the input features (Algorithm 2), and `"one"` uses only the
#' farthest pair without mediators (Algorithm 3). Exact distance ties are
#' resolved lexicographically for reproducibility.
#'
#' @inheritParams hg_neural
#' @param method `"hypergcn"`, `"fast"`, or `"one"` for Algorithms 1--3.
#' @param hidden Hidden width (default 128).
#' @param epochs Training epochs (default 200).
#' @param lr,weight_decay Adam learning rate and L2 penalty.
#' @param dropout Hidden-layer dropout.
#' @param validation Fraction of labeled nodes held out to select an epoch.
#' @param edge_weights Optional positive hyperedge weights.
#' @param seed Reproducible R and torch seed.
#' @param verbose Message loss every 20 epochs.
#' @return A prediction data.frame in the same format as [hg_neural()], with
#'   training history attached.
#' @references
#' Yadati, N., Nimishakavi, M., Yadav, P., Nitin, V., Louis, A., & Talukdar,
#' P. (2019). HyperGCN: A new method of training graph convolutional
#' networks on hypergraphs. \emph{NeurIPS 32}.
#' @export
hypergraph_hypergcn <- function(
    hg, labels, features = "incidence",
    method = c("hypergcn", "fast", "one"), hidden = 128L,
    epochs = 200L, lr = 0.01, weight_decay = 5e-4, dropout = 0.5,
    validation = 0.1, edge_weights = NULL, seed = 1L, verbose = FALSE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(errorCondition(
      "hypergraph_hypergcn() needs the torch package",
      class = "hypernets_missing_torch", call = NULL
    ))
  }
  method <- match.arg(method)
  problem <- .thg_neural_problem(hg, labels, features)
  stopifnot(
    "`hidden` must be a positive whole number" =
      is.numeric(hidden) && length(hidden) == 1L && hidden >= 1 &&
        hidden == round(hidden),
    "`dropout` must be in [0, 1)" =
      is.numeric(dropout) && length(dropout) == 1L &&
        dropout >= 0 && dropout < 1
  )
  x_matrix <- as.matrix(problem$features)
  x_t <- torch::torch_tensor(x_matrix, dtype = torch::torch_float())
  fixed <- if (method == "fast") {
    .thg_torch_sparse(.thg_hypergcn_adjacency(
      hg, x_matrix, edge_weights, mediators = TRUE
    ))
  } else NULL
  make_model <- function() {
    torch::nn_module(
      initialize = function() {
        self$lin1 <- torch::nn_linear(ncol(x_matrix), as.integer(hidden),
                                      bias = FALSE)
        self$lin2 <- torch::nn_linear(as.integer(hidden),
                                      length(problem$classes), bias = FALSE)
        self$drop <- torch::nn_dropout(dropout)
      },
      forward = function(x) {
        z1 <- self$lin1(x)
        a1 <- if (!is.null(fixed)) fixed else {
          signal <- as.matrix(z1$detach()$cpu())
          .thg_torch_sparse(.thg_hypergcn_adjacency(
            hg, signal, edge_weights, mediators = method != "one"
          ))
        }
        h <- self$drop(torch::nnf_relu(torch::torch_mm(a1, z1)))
        z2 <- self$lin2(h)
        a2 <- if (!is.null(fixed)) fixed else {
          signal <- as.matrix(z2$detach()$cpu())
          .thg_torch_sparse(.thg_hypergcn_adjacency(
            hg, signal, edge_weights, mediators = method != "one"
          ))
        }
        torch::torch_mm(a2, z2)
      }
    )()
  }
  out <- .thg_neural_fit(
    make_model, function(model) model(x_t), hg$nodes, problem$labels,
    problem$classes, epochs, lr, weight_decay, validation, seed, verbose
  )
  attr(out, "method") <- method
  out
}

#' @rdname hypergraph_hypergcn
#' @export
hg_hypergcn <- hypergraph_hypergcn
