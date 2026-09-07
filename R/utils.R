# ---- Shared internal helpers ----
#
# Shared by all three structure families (memory networks, simplicial
# complexes, hypergraphs). Copied verbatim from Nestimate (R/utils.R,
# R/estimate_network.R, R/mcml.R) as part of the hypernets delegation (see
# Nestimate's HYPERNETS-DELEGATION-PLAN.md). Nestimate keeps its own copies of
# .coerce_sequence_input / .as_netobject / .extract_edges_from_matrix /
# .validate_mcml_matrix / .wrap_netobject (used elsewhere there);
# .ho_cograph_fields moves here outright (its only callers are the
# memory-family builders).

#' Coerce tna or netobject to labeled sequence data.frame
#'
#' When \code{data} is a \code{tna} or \code{netobject}, extracts the
#' sequence data and converts numeric state IDs to label names. This
#' allows \code{build_hon()}, \code{build_hypa()}, and other pathway
#' functions to accept model objects directly.
#'
#' @param data Input: data.frame, list, tna, or netobject.
#' @return A data.frame or list suitable for \code{.hon_parse_input()}.
#' @noRd
.coerce_sequence_input <- function(data) {
  if (inherits(data, "tna")) {
    if (is.null(data$data)) { # nocov start
      stop("tna object has no sequence data ($data). ",
           "Build the tna from sequence data, not a raw matrix.",
           call. = FALSE)
    } # nocov end
    df <- as.data.frame(data$data, stringsAsFactors = FALSE) # nocov start
    lbl <- attr(data$data, "labels") %||% data$labels
    if (!is.null(lbl) && length(lbl) > 0L &&
        (is.integer(df[[1]]) || is.numeric(df[[1]]))) {
      df[] <- lapply(df, function(col) {
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(lbl),
               NA_character_, lbl[idx])
      })
    }
    return(df) # nocov end
  }
  if (inherits(data, "cograph_network") && !inherits(data, "netobject")) {
    data <- .as_netobject(data)
  }
  if (inherits(data, "netobject")) {
    if (is.null(data$data)) { # nocov start
      stop("netobject has no sequence data ($data). ",
           "Build the network from sequence data.",
           call. = FALSE) # nocov end
    }
    df <- as.data.frame(data$data, stringsAsFactors = FALSE)
    lbl <- rownames(data$weights)
    if (!is.null(lbl) && length(lbl) > 0L &&
        (is.integer(df[[1]]) || is.numeric(df[[1]]))) {
      df[] <- lapply(df, function(col) { # nocov start
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(lbl),
               NA_character_, lbl[idx])
      }) # nocov end
    }
    return(df)
  }
  ## Bare sequence matrix (character / logical) -> wide data.frame.
  if (is.matrix(data) && !is.numeric(data)) {
    return(as.data.frame(data, stringsAsFactors = FALSE))
  }
  data
}

#' Convert pure cograph_network to dual-class netobject/cograph_network
#'
#' Internal converter so that the builders can accept either
#' \code{netobject} or \code{cograph_network} inputs transparently.
#' Objects that already have the \code{"netobject"} class are returned
#' unchanged.
#'
#' @param x A \code{netobject} (returned unchanged) or \code{cograph_network}.
#' @return A dual-class \code{c("netobject", "cograph_network")} object.
#' @noRd
.as_netobject <- function(x) {
  if (inherits(x, "netobject")) return(x)
  if (!inherits(x, "cograph_network")) {
    stop("Expected a netobject or cograph_network.", call. = FALSE)
  }

  mat <- x$weights
  if (is.null(mat)) {
    stop("cograph_network has no $weights matrix.", call. = FALSE)
  }
  if (!is.matrix(mat)) mat <- as.matrix(mat)
  if (!is.numeric(mat)) storage.mode(mat) <- "double"
  nodes_df <- x$nodes
  states <- nodes_df$label
  raw_data <- x$data
  directed <- x$directed %||% TRUE

  # Infer method from tna metadata or matrix symmetry
  tna_meta <- x$meta$tna
  method <- if (!is.null(tna_meta$method)) {
    tna_meta$method
  } else if (is.matrix(mat) && isSymmetric(mat)) {
    "co_occurrence"
  } else {
    "relative"
  }

  is_sequence_method <- method %in% c(
    "relative", "frequency", "co_occurrence", "attention"
  )

  # Decode integer-encoded tna data -> character labels
  # Only for sequence methods; association methods keep numeric data as-is
  if (!is.null(raw_data)) {
    raw_data <- as.data.frame(raw_data, stringsAsFactors = FALSE)
    if (is_sequence_method &&
        (is.integer(raw_data[[1]]) || is.numeric(raw_data[[1]]))) {
      raw_data[] <- lapply(raw_data, function(col) {
        idx <- as.integer(col)
        ifelse(is.na(idx) | idx < 1L | idx > length(states),
               NA_character_, states[idx])
      })
    }
  }

  edges <- .extract_edges_from_matrix(mat, directed = directed)

  structure(list(
    data = raw_data, weights = mat, nodes = nodes_df,
    edges = edges, directed = directed, method = method,
    params = list(), scaling = NULL, threshold = 0,
    n_nodes = length(states), n_edges = nrow(edges),
    level = NULL,
    meta = x$meta %||% list(source = "cograph", layout = NULL,
                            tna = list(method = method)),
    node_groups = x$node_groups
  ), class = c("netobject", "cograph_network"))
}

#' Extract a tidy edge list from a weight matrix
#'
#' @param mat Numeric weight matrix.
#' @param directed Logical. Directed network?
#' @return data.frame with integer \code{from}/\code{to} and \code{weight}.
#' @noRd
.extract_edges_from_matrix <- function(mat, directed = FALSE) {
  if (directed) {
    # Keep self-loops too: every non-zero entry is a real edge.
    idx <- which(mat != 0, arr.ind = TRUE)
  } else {
    # row <= col keeps the upper triangle PLUS the diagonal (one row
    # per self-loop, no double-count for undirected networks).
    idx <- which(mat != 0 & row(mat) <= col(mat), arr.ind = TRUE)
  }

  if (nrow(idx) == 0) {
    return(data.frame(
      from = integer(0), to = integer(0),
      weight = numeric(0), stringsAsFactors = FALSE
    ))
  }

  data.frame(
    from   = as.integer(idx[, 1]),
    to     = as.integer(idx[, 2]),
    weight = mat[idx],
    stringsAsFactors = FALSE
  )
}

#' Add cograph_network fields to a higher-order network object
#'
#' @param mat Square weight matrix with named rows/columns.
#' @param node_names Character vector of node names.
#' @param method Character. Method label for metadata.
#' @return Named list with \code{weights}, \code{nodes} (data.frame),
#'   \code{edges}, \code{directed}, \code{meta} fields.
#' @noRd
.ho_cograph_fields <- function(mat, node_names, method = "hon") {
  nodes_df <- data.frame(
    id = seq_along(node_names),
    label = node_names,
    name = node_names,
    stringsAsFactors = FALSE
  )
  edges <- .extract_edges_from_matrix(mat, directed = TRUE)
  list(
    weights = mat,
    nodes = nodes_df,
    edges = edges,
    directed = TRUE,
    n_nodes = length(node_names),
    n_edges = nrow(edges),
    meta = list(
      source = "nestimate",
      layout = NULL,
      tna = list(method = method)
    ),
    node_groups = NULL
  )
}

# =========================================================================
# Matrix validation + minimal netobject wrapper (hypergraph family)
# =========================================================================

.validate_mcml_matrix <- function(mat) {
  if (any(is.na(mat) | !is.finite(mat))) {
    stop("x matrix must contain finite non-missing weights.", call. = FALSE)
  }
  row_names <- rownames(mat)
  col_names <- colnames(mat)
  has_row_names <- !is.null(row_names)
  has_col_names <- !is.null(col_names)
  if (has_row_names && (any(is.na(row_names)) || any(!nzchar(row_names)))) {
    stop("x matrix row names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (has_col_names && (any(is.na(col_names)) || any(!nzchar(col_names)))) {
    stop("x matrix column names must not contain missing or empty values.",
         call. = FALSE)
  }
  if (has_row_names && anyDuplicated(row_names)) {
    stop("x matrix row names must be unique.", call. = FALSE)
  }
  if (has_col_names && anyDuplicated(col_names)) {
    stop("x matrix column names must be unique.", call. = FALSE)
  }
  if (has_row_names && has_col_names && !identical(row_names, col_names)) {
    stop("x matrix row and column names must be identical and in the same order.",
         call. = FALSE)
  }
  invisible(TRUE)
}



#' Wrap a weight matrix + optional data into a minimal netobject
#' @noRd
.wrap_netobject <- function(mat, data = NULL, method = "relative",
                            directed = TRUE, inits = NULL) {
  if (!is.matrix(mat) || !is.numeric(mat)) {
    stop("'mat' must be a numeric matrix.", call. = FALSE)
  }
  if (nrow(mat) != ncol(mat)) {
    stop("'mat' must be a square matrix.", call. = FALSE)
  }
  .validate_mcml_matrix(mat)
  if (!is.character(method) || length(method) != 1L || is.na(method) ||
      !nzchar(method)) {
    stop("'method' must be a single non-missing character value.",
         call. = FALSE)
  }
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    stop("'directed' must be TRUE or FALSE.", call. = FALSE)
  }
  states <- rownames(mat)
  if (is.null(states)) {
    states <- as.character(seq_len(nrow(mat)))
    dimnames(mat) <- list(states, states)
  }
  edges <- .extract_edges_from_matrix(mat, directed = directed)
  nodes_df <- data.frame(
    id = seq_along(states), label = states, name = states,
    x = NA_real_, y = NA_real_, stringsAsFactors = FALSE
  )
  if (!is.null(inits)) {
    if (!is.numeric(inits) || length(inits) != length(states) ||
        any(is.na(inits) | !is.finite(inits))) {
      stop("'inits' must be a finite numeric vector with one value per state.",
           call. = FALSE)
    }
    if (!is.null(names(inits))) {
      missing_inits <- setdiff(states, names(inits))
      extra_inits <- setdiff(names(inits), states)
      if (length(missing_inits) > 0L || length(extra_inits) > 0L) {
        stop("'inits' names must match matrix state names.", call. = FALSE)
      }
      inits <- inits[states]
    } else {
      names(inits) <- states
    }
  }

  structure(
    list(
      data       = data,
      weights    = mat,
      inits      = inits,
      nodes      = nodes_df,
      edges      = edges,
      directed   = directed,
      method     = method,
      params     = list(),
      scaling    = NULL,
      threshold  = 0,
      n_nodes    = length(states),
      n_edges    = nrow(edges),
      level      = NULL,
      meta       = list(source = "nestimate", layout = NULL,
                        tna = list(method = method)),
      node_groups = NULL
    ),
    class = c("netobject", "cograph_network")
  )
}
