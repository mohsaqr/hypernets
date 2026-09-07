# ---- HNHN (Dong et al. 2020) -------------------------------------------

# Paper Section 2.3 normalization operators. v_to_e is m x n and e_to_v is
# n x m; every row sums to one. beta controls vertex-degree weighting during
# V->E propagation, alpha hyperedge-cardinality weighting during E->V.
.thg_hnhn_operators <- function(hg, alpha = 0, beta = 0) {
  stopifnot(
    "`alpha` and `beta` must be finite numbers" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
        is.numeric(beta) && length(beta) == 1L && is.finite(beta)
  )
  H <- Matrix::Matrix((hg$incidence != 0) * 1, sparse = TRUE)
  vertex_degree <- as.numeric(Matrix::rowSums(H))
  edge_degree <- as.numeric(Matrix::colSums(H))
  stopifnot("every vertex and hyperedge must have positive degree" =
              all(vertex_degree > 0) && all(edge_degree > 0))
  vertex_weight <- vertex_degree^beta
  edge_weight <- edge_degree^alpha
  v_to_e_den <- as.numeric(Matrix::t(H) %*% vertex_weight)
  e_to_v_den <- as.numeric(H %*% edge_weight)
  v_to_e <- Matrix::Diagonal(x = 1 / v_to_e_den) %*%
    Matrix::t(H) %*% Matrix::Diagonal(x = vertex_weight)
  e_to_v <- Matrix::Diagonal(x = 1 / e_to_v_den) %*%
    H %*% Matrix::Diagonal(x = edge_weight)
  list(v_to_e = v_to_e, e_to_v = e_to_v,
       vertex_degree = vertex_degree, edge_degree = edge_degree,
       alpha = alpha, beta = beta)
}

#' HNHN node classifier
#'
#' Trains Hypergraph Networks with Hyperedge Neurons (Dong et al. 2020).
#' A full layer sends node features to learned hyperedge representations and
#' then sends those representations back to nodes. The paper's normalization
#' is implemented directly: `alpha` changes the relative contribution of
#' large hyperedges in edge-to-node messages, while `beta` changes the
#' contribution of high-degree nodes in node-to-edge messages. Zero recovers
#' ordinary neighborhood means in each direction.
#'
#' @inheritParams hg_neural
#' @param alpha Hyperedge-cardinality normalization exponent.
#' @param beta Vertex-degree normalization exponent.
#' @param hidden Hidden width (default 128).
#' @param epochs Training epochs (default 200).
#' @param lr,weight_decay Adam learning rate and L2 penalty.
#' @param dropout Hidden-layer dropout.
#' @param validation Fraction of labeled nodes held out to select an epoch.
#' @param seed Reproducible R and torch seed.
#' @param verbose Message loss every 20 epochs.
#' @return A prediction data.frame in the same format as [hg_neural()], with
#'   training history and the normalization exponents attached.
#' @references
#' Dong, Y., Sawin, W., & Bengio, Y. (2020). HNHN: Hypergraph networks with
#' hyperedge neurons. \emph{ICML Graph Representation Learning and Beyond
#' Workshop}.
#' @export
hypergraph_hnhn <- function(
    hg, labels, features = "incidence", hidden = 128L,
    alpha = 0, beta = 0, epochs = 200L, lr = 0.01,
    weight_decay = 5e-4, dropout = 0.5, validation = 0.1,
    seed = 1L, verbose = FALSE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(errorCondition("hypergraph_hnhn() needs the torch package",
                        class = "hypernets_missing_torch", call = NULL))
  }
  problem <- .thg_neural_problem(hg, labels, features)
  stopifnot(
    "`hidden` must be a positive whole number" =
      is.numeric(hidden) && length(hidden) == 1L && hidden >= 1 &&
        hidden == round(hidden),
    "`dropout` must be in [0, 1)" =
      is.numeric(dropout) && length(dropout) == 1L &&
        dropout >= 0 && dropout < 1
  )
  operators <- .thg_hnhn_operators(hg, alpha, beta)
  p_ve <- .thg_torch_sparse(operators$v_to_e)
  p_ev <- .thg_torch_sparse(operators$e_to_v)
  x_matrix <- problem$features
  x_t <- if (methods::is(x_matrix, "sparseMatrix")) {
    .thg_torch_sparse(x_matrix)
  } else torch::torch_tensor(x_matrix, dtype = torch::torch_float())
  make_model <- function() {
    torch::nn_module(
      initialize = function() {
        self$input <- torch::nn_linear(ncol(x_matrix), as.integer(hidden))
        self$edge <- torch::nn_linear(as.integer(hidden), as.integer(hidden))
        self$node <- torch::nn_linear(as.integer(hidden),
                                      length(problem$classes))
        self$drop <- torch::nn_dropout(dropout)
      },
      forward = function(x) {
        node <- self$drop(torch::nnf_relu(self$input(x)))
        edge <- torch::nnf_relu(self$edge(torch::torch_mm(p_ve, node)))
        self$node(torch::torch_mm(p_ev, self$drop(edge)))
      }
    )()
  }
  out <- .thg_neural_fit(
    make_model, function(model) model(x_t), hg$nodes, problem$labels,
    problem$classes, epochs, lr, weight_decay, validation, seed, verbose
  )
  attr(out, "normalization") <- c(alpha = alpha, beta = beta)
  out
}

#' @rdname hypergraph_hnhn
#' @export
hg_hnhn <- hypergraph_hnhn
