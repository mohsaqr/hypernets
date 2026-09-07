# ---- Heterogeneous Graph Attention (Linmei et al. 2019) -----------------

# Built on demand: torch is in Suggests, so no torch call may run when
# the namespace is loaded.
.thg_hgat_layer_module <- function() torch::nn_module(
  initialize = function(hidden, type_names, adjacency, type_index) {
    self$hidden <- hidden
    self$type_names <- type_names
    self$adjacency <- adjacency
    self$type_index <- type_index
    self$transform <- torch::nn_module_list(lapply(type_names, function(x) {
      torch::nn_linear(hidden, hidden, bias = FALSE)
    }))
    self$type_attention <- torch::nn_module_list(lapply(type_names, function(x) {
      torch::nn_linear(2L * hidden, 1L, bias = FALSE)
    }))
    self$node_attention <- torch::nn_linear(2L * hidden, 1L, bias = FALSE)
  },
  forward = function(h) {
    rows <- lapply(seq_len(nrow(self$adjacency)), function(v) {
      hv <- h[v, , drop = FALSE]
      active <- which(vapply(self$type_index, function(idx) {
        any(self$adjacency[v, idx] > 0)
      }, logical(1L)))
      type_scores <- lapply(active, function(type) {
        idx <- self$type_index[[type]]
        idx <- idx[self$adjacency[v, idx] > 0]
        weights <- torch::torch_tensor(self$adjacency[v, idx],
                                       dtype = torch::torch_float())
        neigh <- h[torch::torch_tensor(idx, dtype = torch::torch_int64()), ,
                   drop = FALSE]
        h_type <- torch::torch_sum(neigh * weights$unsqueeze(2L),
                                   dim = 1L, keepdim = TRUE)
        torch::nnf_leaky_relu(
          self$type_attention[[type]](torch::torch_cat(list(hv, h_type),
                                                        dim = 2L)),
          negative_slope = 0.2
        )
      })
      alpha <- torch::nnf_softmax(torch::torch_cat(type_scores, dim = 1L),
                                  dim = 1L)
      neighbor_scores <- list()
      neighbor_rows <- list()
      neighbor_types <- integer(0)
      cursor <- 0L
      for (position in seq_along(active)) {
        type <- active[position]
        idx <- self$type_index[[type]]
        idx <- idx[self$adjacency[v, idx] > 0]
        transformed <- self$transform[[type]](
          h[torch::torch_tensor(idx, dtype = torch::torch_int64()), ,
            drop = FALSE]
        )
        for (j in seq_along(idx)) {
          cursor <- cursor + 1L
          hu <- h[idx[j], , drop = FALSE]
          pair <- alpha[position, , drop = FALSE] *
            torch::torch_cat(list(hv, hu), dim = 2L)
          neighbor_scores[[cursor]] <- torch::nnf_leaky_relu(
            self$node_attention(pair), negative_slope = 0.2
          )
          neighbor_rows[[cursor]] <- transformed[j, , drop = FALSE]
          neighbor_types[cursor] <- type
        }
      }
      beta <- torch::nnf_softmax(torch::torch_cat(neighbor_scores, dim = 1L),
                                 dim = 1L)
      messages <- torch::torch_cat(neighbor_rows, dim = 1L)
      torch::nnf_relu(torch::torch_sum(messages * beta, dim = 1L,
                                       keepdim = TRUE))
    })
    torch::torch_cat(rows, dim = 1L)
  }
)
.thg_hgat_layer <- function(...) .thg_hgat_layer_module()(...)

.thg_hgat_inputs <- function(adjacency, node_types, features) {
  stopifnot(
    "`adjacency` must be a finite non-negative square numeric matrix" =
      is.matrix(adjacency) && is.numeric(adjacency) &&
        nrow(adjacency) == ncol(adjacency) &&
        all(is.finite(adjacency)) && all(adjacency >= 0),
    "`adjacency` must be symmetric" =
      isTRUE(all.equal(adjacency, t(adjacency), tolerance = 1e-12)),
    "`adjacency` needs unique row and column names in the same order" =
      !is.null(rownames(adjacency)) && !is.null(colnames(adjacency)) &&
        identical(rownames(adjacency), colnames(adjacency)) &&
        anyDuplicated(rownames(adjacency)) == 0L,
    "`node_types` must be named for every adjacency node" =
      is.character(node_types) && !is.null(names(node_types)) &&
        setequal(names(node_types), rownames(adjacency)),
    "`features` must be a named list, one matrix per node type" =
      is.list(features) && !is.null(names(features)) &&
        setequal(names(features), unique(node_types))
  )
  nodes <- rownames(adjacency)
  node_types <- node_types[nodes]
  type_names <- sort(unique(node_types))
  type_index <- lapply(type_names, function(type) which(node_types == type))
  names(type_index) <- type_names
  aligned <- lapply(type_names, function(type) {
    value <- features[[type]]
    stopifnot("each feature block must be a finite numeric matrix" =
                is.matrix(value) && is.numeric(value) &&
                all(is.finite(value)) && !is.null(rownames(value)))
    expected <- nodes[type_index[[type]]]
    if (!setequal(rownames(value), expected)) {
      stop(sprintf("feature rows for type `%s` must match its node names", type),
           call. = FALSE)
    }
    value[expected, , drop = FALSE]
  })
  names(aligned) <- type_names
  A <- adjacency
  diag(A) <- diag(A) + 1
  degree <- rowSums(A)
  A <- A / sqrt(outer(degree, degree))
  list(nodes = nodes, node_types = node_types, type_names = type_names,
       type_index = type_index, features = aligned, adjacency = A)
}

#' Heterogeneous graph attention classifier (HGAT)
#'
#' Implements Linmei et al.'s (2019) heterogeneous graph convolution and
#' dual-level attention. Node types retain separate input feature spaces and
#' transformation matrices. For every target node, type-level attention first
#' weighs the available neighboring information types; node-level attention
#' then weighs individual neighbors before their type-specific transformed
#' representations are combined. This is HGAT for heterogeneous information
#' networks, distinct from Ding et al.'s document-hypergraph [text_hypergat()].
#'
#' @param adjacency Symmetric non-negative adjacency matrix with identical
#'   unique row and column names.
#' @param node_types Named character vector assigning every node a type.
#' @param features Named list with one numeric feature matrix per type. Rows
#'   are named nodes of that type; feature widths may differ across types.
#' @param labels Named character labels for a subset of nodes.
#' @param hidden Common hidden dimension.
#' @param layers Number of dual-attention layers (paper default 2).
#' @param epochs,lr,weight_decay,dropout,validation,seed,verbose Training
#'   controls as in [hg_neural()].
#' @return A prediction data.frame for every heterogeneous-network node, with
#'   training history and node types attached as attributes.
#' @references
#' Linmei, H., Yang, T., Shi, C., Ji, H., & Li, X. (2019). Heterogeneous graph
#' attention networks for semi-supervised short text classification.
#' \emph{EMNLP-IJCNLP 2019}, 4821--4830.
#' @export
heterogeneous_hgat <- function(
    adjacency, node_types, features, labels, hidden = 128L, layers = 2L,
    epochs = 200L, lr = 0.005, weight_decay = 5e-4, dropout = 0.5,
    validation = 0.1, seed = 1L, verbose = FALSE) {
  if (!requireNamespace("torch", quietly = TRUE)) {
    stop(errorCondition("heterogeneous_hgat() needs the torch package",
                        class = "hypernets_missing_torch", call = NULL))
  }
  input <- .thg_hgat_inputs(adjacency, node_types, features)
  labels <- .thg_labels_input(labels)
  stopifnot(
    "`labels` must be a named character vector" =
      is.character(labels) && !is.null(names(labels)),
    "`hidden` and `layers` must be positive whole numbers" =
      all(vapply(list(hidden, layers), function(z) {
        is.numeric(z) && length(z) == 1L && z >= 1 && z == round(z)
      }, logical(1L))),
    "`dropout` must be in [0, 1)" =
      is.numeric(dropout) && length(dropout) == 1L &&
        dropout >= 0 && dropout < 1
  )
  unknown <- setdiff(names(labels), input$nodes)
  if (length(unknown)) {
    stop("Unknown node names in `labels`: ", paste(unknown, collapse = ", "),
         call. = FALSE)
  }
  classes <- sort(unique(as.character(labels)))
  if (length(classes) < 2L) {
    stop(errorCondition("`labels` must contain at least two distinct classes.",
                        class = "hypernets_bad_input", call = NULL))
  }
  x <- lapply(input$features, function(value) {
    torch::torch_tensor(value, dtype = torch::torch_float())
  })
  make_model <- function() {
    torch::nn_module(
      initialize = function() {
        self$input <- torch::nn_module_list(lapply(input$type_names,
                                                    function(type) {
          torch::nn_linear(ncol(input$features[[type]]), as.integer(hidden))
        }))
        self$attention <- torch::nn_module_list(lapply(seq_len(layers),
                                                        function(i) {
          .thg_hgat_layer(as.integer(hidden), input$type_names,
                          input$adjacency, input$type_index)
        }))
        self$drop <- torch::nn_dropout(dropout)
        self$out <- torch::nn_linear(as.integer(hidden), length(classes))
      },
      forward = function(values) {
        projected <- lapply(seq_along(input$type_names), function(type) {
          torch::nnf_relu(self$input[[type]](values[[type]]))
        })
        position <- integer(length(input$nodes))
        for (type in seq_along(input$type_index)) {
          position[input$type_index[[type]]] <- seq_along(input$type_index[[type]])
        }
        rows <- lapply(seq_along(input$nodes), function(node) {
          type <- match(input$node_types[node], input$type_names)
          projected[[type]][position[node], , drop = FALSE]
        })
        h <- torch::torch_cat(rows, dim = 1L)
        for (layer in seq_len(layers)) {
          h <- self$drop(self$attention[[layer]](h))
        }
        self$out(h)
      }
    )()
  }
  out <- .thg_neural_fit(
    make_model, function(model) model(x), input$nodes, labels, classes,
    epochs, lr, weight_decay, validation, seed, verbose
  )
  attr(out, "node_types") <- input$node_types
  out
}
