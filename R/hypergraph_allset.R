# ---- AllSet (Chien et al. 2022) ----------------------------------------

# Built on demand: torch is in Suggests, so no torch call may run when
# the namespace is loaded.
.thg_deepset_module <- function() torch::nn_module(
  initialize = function(in_dim, out_dim) {
    self$phi <- torch::nn_sequential(
      torch::nn_linear(in_dim, out_dim), torch::nn_relu(),
      torch::nn_linear(out_dim, out_dim)
    )
    self$rho <- torch::nn_sequential(
      torch::nn_linear(out_dim, out_dim), torch::nn_relu(),
      torch::nn_linear(out_dim, out_dim)
    )
  },
  forward = function(x, incidence) {
    self$rho(torch::torch_mm(incidence, self$phi(x)))
  }
)
.thg_deepset <- function(...) .thg_deepset_module()(...)

# Built on demand: torch is in Suggests, so no torch call may run when
# the namespace is loaded.
.thg_pma_module <- function() torch::nn_module(
  initialize = function(dim, heads) {
    if (dim %% heads != 0L) stop("`hidden` must be divisible by `heads`.")
    self$dim <- dim
    self$heads <- heads
    self$head_dim <- as.integer(dim / heads)
    self$key <- torch::nn_linear(dim, dim, bias = FALSE)
    self$value <- torch::nn_linear(dim, dim, bias = FALSE)
    self$theta <- torch::nn_parameter(
      torch::torch_empty(1L, heads, self$head_dim)$uniform_(
        -1 / sqrt(dim), 1 / sqrt(dim)
      )
    )
    self$norm1 <- torch::nn_layer_norm(dim)
    self$mlp <- torch::nn_sequential(
      torch::nn_linear(dim, dim), torch::nn_relu(),
      torch::nn_linear(dim, dim)
    )
    self$norm2 <- torch::nn_layer_norm(dim)
  },
  forward = function(x, sets) {
    outputs <- lapply(sets, function(index) {
      idx <- torch::torch_tensor(index, dtype = torch::torch_int64())
      values <- x[idx, , drop = FALSE]
      n <- length(index)
      key <- self$key(values)$reshape(c(n, self$heads, self$head_dim))
      value <- self$value(values)$reshape(c(n, self$heads, self$head_dim))
      score <- torch::torch_sum(key * self$theta, dim = 3L) /
        sqrt(self$head_dim)
      attention <- torch::nnf_softmax(score, dim = 1L)$unsqueeze(3L)
      pooled <- torch::torch_sum(attention * value, dim = 1L)$reshape(
        c(1L, self$dim)
      )
      query <- self$theta$reshape(c(1L, self$dim))
      y <- self$norm1(query + pooled)
      self$norm2(y + self$mlp(y))
    })
    torch::torch_cat(outputs, dim = 1L)
  }
)
.thg_pma <- function(...) .thg_pma_module()(...)

.thg_allset_memberships <- function(hg) {
  H <- hg$incidence != 0
  if (any(Matrix::rowSums(H) == 0) || any(Matrix::colSums(H) == 0)) {
    stop("AllSet requires every vertex and hyperedge to have a membership.",
         call. = FALSE)
  }
  list(
    vertices_by_edge = lapply(seq_len(ncol(H)), function(e) {
      unname(which(H[, e]))
    }),
    edges_by_vertex = lapply(seq_len(nrow(H)), function(v) {
      unname(which(H[v, ]))
    }),
    H = Matrix::Matrix(H * 1, sparse = TRUE)
  )
}

#' AllSet hypergraph node classifier
#'
#' Implements the AllSet framework of Chien et al. (2022) as one complete
#' node-to-edge-to-node propagation layer. `model = "deepsets"` uses Eq. 7,
#' an MLP before and after each multiset sum. `model = "transformer"` uses
#' Eq. 8: pooling by multihead attention from a learned query, followed by
#' residual MLPs and layer normalization. Both directions have independent
#' learned multiset functions and are permutation invariant.
#'
#' @inheritParams hg_neural
#' @param model `"deepsets"` for AllDeepSets or `"transformer"` for
#'   AllSetTransformer.
#' @param heads Number of attention heads for AllSetTransformer. `hidden`
#'   must be divisible by `heads`.
#' @param hidden Hidden width (default 128).
#' @param epochs Training epochs (default 200).
#' @param lr,weight_decay Adam learning rate and L2 penalty.
#' @param dropout Hidden-layer dropout.
#' @param validation Fraction of labeled nodes held out to select an epoch.
#' @param seed Reproducible R and torch seed.
#' @param verbose Message loss every 20 epochs.
#' @return A prediction data.frame in the same format as [hg_neural()], with
#'   training history and the AllSet model name attached.
#' @references
#' Chien, E., Pan, C., Peng, J., & Milenkovic, O. (2022). You are AllSet: A
#' multiset function framework for hypergraph neural networks. \emph{ICLR
#' 2022}.
#' @export
hypergraph_allset <- function(
    hg, labels, features = "incidence",
    model = c("deepsets", "transformer"), hidden = 128L, heads = 4L,
    epochs = 200L, lr = 0.001, weight_decay = 5e-4, dropout = 0.5,
    validation = 0.1, seed = 1L, verbose = FALSE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(errorCondition("hypergraph_allset() needs the torch package",
                        class = "hypernets_missing_torch", call = NULL))
  }
  model <- match.arg(model)
  problem <- .thg_neural_problem(hg, labels, features)
  stopifnot(
    "`hidden` and `heads` must be positive whole numbers" =
      all(vapply(list(hidden, heads), function(z) {
        is.numeric(z) && length(z) == 1L && z >= 1 && z == round(z)
      }, logical(1L))),
    "`dropout` must be in [0, 1)" =
      is.numeric(dropout) && length(dropout) == 1L &&
        dropout >= 0 && dropout < 1
  )
  if (model == "transformer" && hidden %% heads != 0L) {
    stop("`hidden` must be divisible by `heads`.", call. = FALSE)
  }
  memberships <- .thg_allset_memberships(hg)
  x_matrix <- problem$features
  x_t <- if (methods::is(x_matrix, "sparseMatrix")) {
    .thg_torch_sparse(x_matrix)
  } else torch::torch_tensor(x_matrix, dtype = torch::torch_float())
  H <- .thg_torch_sparse(memberships$H)
  Ht <- .thg_torch_sparse(Matrix::t(memberships$H))
  make_model <- if (model == "deepsets") {
    function() {
      torch::nn_module(
        initialize = function() {
          self$input <- torch::nn_linear(ncol(x_matrix), as.integer(hidden))
          self$v_to_e <- .thg_deepset(as.integer(hidden), as.integer(hidden))
          self$e_to_v <- .thg_deepset(as.integer(hidden), as.integer(hidden))
          self$drop <- torch::nn_dropout(dropout)
          self$out <- torch::nn_linear(as.integer(hidden),
                                       length(problem$classes))
        },
        forward = function(x) {
          node <- torch::nnf_relu(self$input(x))
          edge <- self$v_to_e(node, Ht)
          node <- self$e_to_v(self$drop(edge), H)
          self$out(self$drop(torch::nnf_relu(node)))
        }
      )()
    }
  } else {
    function() {
      torch::nn_module(
        initialize = function() {
          self$input <- torch::nn_linear(ncol(x_matrix), as.integer(hidden))
          self$v_to_e <- .thg_pma(as.integer(hidden), as.integer(heads))
          self$e_to_v <- .thg_pma(as.integer(hidden), as.integer(heads))
          self$drop <- torch::nn_dropout(dropout)
          self$out <- torch::nn_linear(as.integer(hidden),
                                       length(problem$classes))
        },
        forward = function(x) {
          node <- torch::nnf_relu(self$input(x))
          edge <- self$v_to_e(node, memberships$vertices_by_edge)
          node <- self$e_to_v(self$drop(edge), memberships$edges_by_vertex)
          self$out(self$drop(torch::nnf_relu(node)))
        }
      )()
    }
  }
  out <- .thg_neural_fit(
    make_model, function(fit) fit(x_t), hg$nodes, problem$labels,
    problem$classes, epochs, lr, weight_decay, validation, seed, verbose
  )
  attr(out, "model") <- if (model == "deepsets") {
    "AllDeepSets"
  } else "AllSetTransformer"
  out
}

#' @rdname hypergraph_allset
#' @export
hg_allset <- hypergraph_allset

#' @rdname hypergraph_allset
#' @param ... Arguments passed to [hypergraph_allset()].
#' @export
hypergraph_alldeepsets <- function(...) {
  hypergraph_allset(..., model = "deepsets")
}

#' @rdname hypergraph_allset
#' @export
hypergraph_allset_transformer <- function(...) {
  hypergraph_allset(..., model = "transformer")
}
