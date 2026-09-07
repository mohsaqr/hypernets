# Neural tier: HGNN (Feng et al. 2019) natively in R via {torch}.
#
# The HGNN propagation matrix G = Dv^-1/2 H W De^-1 H^T Dv^-1/2 is exactly
# the Zhou (2006) similarity operator S = I - L_zhou that
# hypergraph_transduction() spreads with -- HGNN is learned filters on the
# operator the closed-form classifier uses fixed. That identity is asserted
# in the tests (machine precision against .hl_build()), tying the neural
# code to the oracle-verified spectral core. Layer semantics (bias applied
# BEFORE propagation, ReLU/dropout skipped on the last layer) match the
# official implementation (iMoonLab; DHG HGNNConv), verified against its
# source and by forward-pass parity in local_testing_and_equivalence/.

# `self` is injected by torch::nn_module at runtime
utils::globalVariables("self")

# Factor G = A %*% t(A) with A = Dv^-1/2 H diag(sqrt(w/de)): each layer
# then needs two sparse-dense products, never the dense n x n operator.
# Returns the sparse factor as a base dgCMatrix plus the degree vectors.
.thg_hgnn_factor <- function(hg, edge_weights = NULL) {
  membership <- Matrix::Matrix((hg$incidence != 0) * 1, sparse = TRUE)
  m <- ncol(membership)
  w <- edge_weights %||% rep(1, m)
  stopifnot(
    "`edge_weights` must be one positive number per hyperedge" =
      is.numeric(w) && length(w) == m && all(is.finite(w)) && all(w > 0)
  )
  de <- as.numeric(Matrix::colSums(membership))
  dv <- as.numeric(membership %*% w)
  stopifnot(
    "every hyperedge must contain at least one vertex" = all(de > 0),
    "every vertex must sit in at least one hyperedge" = all(dv > 0)
  )
  a <- Matrix::Diagonal(x = 1 / sqrt(dv)) %*% membership %*%
    Matrix::Diagonal(x = sqrt(w / de))
  list(a = a, w = w, de = de, dv = dv)
}

# dgCMatrix -> torch sparse COO tensor.
.thg_torch_sparse <- function(x) {
  triplet <- methods::as(x, "TsparseMatrix")
  torch::torch_sparse_coo_tensor(
    indices = torch::torch_tensor(rbind(triplet@i + 1L, triplet@j + 1L),
                                  dtype = torch::torch_int64()),
    values = torch::torch_tensor(triplet@x, dtype = torch::torch_float()),
    size = dim(x)
  )$coalesce()
}

#' Hypergraph neural network classifier (HGNN)
#'
#' Trains the two-layer hypergraph convolutional network of Feng et al.
#' (2019) on a hypergraph: each layer propagates linearly transformed
#' vertex features through \eqn{G = D_v^{-1/2} H W D_e^{-1} H^T
#' D_v^{-1/2}} (the Zhou similarity operator), with ReLU and dropout
#' between layers and a cross-entropy loss on the labeled vertices.
#' Needs the suggested \pkg{torch} package.
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hypergraph`),
#'   dense or sparse.
#' @param labels The known labels: a named character vector (names are
#'   node identifiers, values class labels) or a tidy data.frame with a
#'   `node` column and a `label`, `cluster` or `predicted` column. At
#'   least two classes.
#' @param features Vertex feature matrix, one row per node in `hg` node
#'   order (rownames must match the node names), e.g. sbert embeddings.
#'   The default `"incidence"` uses the hypergraph's own weighted
#'   incidence rows (tf-idf bag-of-words features for a `nodes = "doc"`
#'   text hypergraph).
#' @param hidden Width of the hidden layer.
#' @param epochs Training epochs (Adam).
#' @param lr,weight_decay Adam learning rate and L2 penalty. The defaults
#'   (`lr = 0.01`, `epochs = 600`) were selected by validation accuracy on
#'   held-out seeds for high-dimensional sparse text features; the original
#'   paper's `lr = 0.001`, 200 epochs underfits that regime.
#' @param dropout Dropout rate after the hidden layer.
#' @param validation Fraction of the labeled seeds held out (stratified)
#'   to pick the best epoch; `0` trains on all seeds for `epochs` and
#'   keeps the final weights.
#' @param edge_weights Optional positive hyperedge weights `W` (one per
#'   hyperedge); default all 1, the paper's initialization.
#' @param seed Integer seed for weight initialization, dropout and the
#'   validation split (torch and R RNGs); results are deterministic
#'   given a seed.
#' @param verbose If `TRUE`, message the loss every 20 epochs.
#' @return A base `data.frame`, one row per node, with columns `node`,
#'   `label` (the given label or `NA`), `predicted`, `score` (softmax
#'   probability of the winning class) and `margin` (winner minus
#'   runner-up probability). The training history is attached as
#'   attribute `"history"` (data.frame: `epoch`, `loss`,
#'   `val_accuracy`).
#' @references
#' Feng, Y., You, H., Zhang, Z., Ji, R., & Gao, Y. (2019). Hypergraph
#' neural networks. \emph{AAAI 33}.
#' @examples
#' \donttest{
#' if (requireNamespace("torch", quietly = TRUE)) {
#'   hg <- text_hypergraph(c(
#'     cooking_1 = "simmer the soup with onions and carrots",
#'     cooking_2 = "this soup recipe needs salt on a cold night",
#'     space_1 = "the telescope revealed a distant galaxy and stars",
#'     space_2 = "astronomers aimed the telescope at the stars all night"
#'   ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#'   hg_neural(hg, labels = c(cooking_1 = "cooking", space_1 = "space"),
#'             hidden = 8, epochs = 50, validation = 0)
#' }
#' }
#' @export
hg_neural <- function(hg, labels, features = "incidence", hidden = 128L,
                      epochs = 600L, lr = 0.01, weight_decay = 5e-4,
                      dropout = 0.5, validation = 0.1,
                      edge_weights = NULL, seed = 1L, verbose = FALSE) {
  .thg_check_hg(hg)
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(errorCondition(
      "hg_neural() needs the torch package: install.packages(\"torch\")",
      class = "hypernets_missing_torch", call = NULL
    ))
  }
  labels <- .thg_labels_input(labels)
  stopifnot(
    "`labels` must be a named character vector" =
      is.character(labels) && !is.null(names(labels)),
    "`hidden` must be a single positive integer" =
      length(hidden) == 1L && is.finite(hidden) && hidden >= 1,
    "`epochs` must be a single positive integer" =
      length(epochs) == 1L && is.finite(epochs) && epochs >= 1,
    "`dropout` must be a single number in [0, 1)" =
      is.numeric(dropout) && length(dropout) == 1L &&
      dropout >= 0 && dropout < 1,
    "`validation` must be a single number in [0, 1)" =
      is.numeric(validation) && length(validation) == 1L &&
      validation >= 0 && validation < 1,
    "`seed` must be a single integer" =
      length(seed) == 1L && is.finite(seed)
  )
  nodes <- hg$nodes
  unknown <- setdiff(names(labels), nodes)
  if (length(unknown) > 0L) {
    stop("Unknown node names in `labels`: ",
         paste(unknown, collapse = ", "), call. = FALSE)
  }
  classes <- sort(unique(as.character(labels)))
  if (length(classes) < 2L) {
    stop(errorCondition("`labels` must contain at least two distinct classes.",
                        class = "hypernets_bad_input", call = NULL))
  }

  # features: the incidence rows, or a caller matrix aligned by rownames
  x_mat <- if (identical(features, "incidence")) {
    hg$incidence
  } else {
    stopifnot(
      "`features` must be \"incidence\" or a numeric matrix" =
        (is.matrix(features) && is.numeric(features)) ||
        methods::is(features, "sparseMatrix"),
      "`features` needs rownames matching the hypergraph nodes" =
        !is.null(rownames(features)) && all(nodes %in% rownames(features))
    )
    features[match(nodes, rownames(features)), , drop = FALSE]
  }

  factor_parts <- .thg_hgnn_factor(hg, edge_weights)

  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
  on.exit(if (!is.null(old_seed)) {
    assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE, after = FALSE)
  set.seed(seed)
  torch::torch_manual_seed(seed)

  a_sp <- .thg_torch_sparse(factor_parts$a)
  at_sp <- .thg_torch_sparse(Matrix::t(factor_parts$a))
  prop <- function(x) torch::torch_mm(a_sp, torch::torch_mm(at_sp, x))

  x_t <- if (methods::is(x_mat, "sparseMatrix")) {
    .thg_torch_sparse(x_mat)
  } else {
    torch::torch_tensor(x_mat, dtype = torch::torch_float())
  }
  in_dim <- ncol(x_mat)
  n_class <- length(classes)

  # DHG/iMoonLab layer semantics: Linear (with bias) FIRST, then propagate
  model <- torch::nn_module(
    initialize = function() {
      self$lin1 <- torch::nn_linear(in_dim, hidden)
      self$lin2 <- torch::nn_linear(hidden, n_class)
      self$drop <- torch::nn_dropout(dropout)
    },
    forward = function(x) {
      h <- prop(torch::torch_mm(x, torch::torch_t(self$lin1$weight)) +
                  self$lin1$bias)
      h <- self$drop(torch::nnf_relu(h))
      prop(torch::torch_mm(h, torch::torch_t(self$lin2$weight)) +
             self$lin2$bias)
    }
  )()

  seed_idx <- match(names(labels), nodes)
  target <- match(as.character(labels), classes)
  val_idx <- integer(0)
  if (validation > 0) {
    # stratified holdout; classes too small to split stay fully in train
    val_pick <- unlist(lapply(split(seq_along(target), target), \(i) {
      k <- floor(validation * length(i))
      if (k >= 1L) sample(i, k) else integer(0)
    }))
    val_idx <- sort(val_pick)
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
  history <- vapply(seq_len(epochs), \(epoch) {
    # sequential training loop: each Adam step depends on the previous one
    model$train()
    optimizer$zero_grad()
    out <- model(x_t)
    loss <- torch::nnf_cross_entropy(out[train_t, ], y_train)
    loss$backward()
    optimizer$step()
    val_acc <- NA_real_
    if (length(val_idx) > 0L) {
      model$eval()
      pred <- torch::with_no_grad(
        as.integer(torch::torch_argmax(model(x_t), dim = 2L))
      )
      val_acc <- mean(pred[seed_idx[val_idx]] == target[val_idx])
      if (val_acc > best_val) {
        best_val <<- val_acc
        best_state <<- lapply(model$state_dict(), \(p) p$clone())
      }
    }
    if (isTRUE(verbose) && epoch %% 20L == 0L) {
      message(sprintf("epoch %d loss %.4f", epoch, as.numeric(loss)))
    }
    c(as.numeric(loss), val_acc)
  }, numeric(2))
  if (!is.null(best_state)) {
    model$load_state_dict(best_state)
  }

  model$eval()
  probs <- torch::with_no_grad(
    as.matrix(torch::nnf_softmax(model(x_t), dim = 2L))
  )
  dimnames(probs) <- list(nodes, classes)
  lab_full <- rep(NA_character_, length(nodes))
  lab_full[seed_idx] <- as.character(labels)
  out <- .hl_score_predictions(probs, lab_full, "none")
  rownames(out) <- NULL
  attr(out, "history") <- data.frame(
    epoch = seq_len(epochs), loss = history[1L, ],
    val_accuracy = history[2L, ]
  )
  out
}

#' @rdname hg_neural
#' @export
hypergraph_neural <- hg_neural
