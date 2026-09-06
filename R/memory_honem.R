# ---- HONEM: Higher-Order Network Embedding ----
#
# Implements HONEM (Saebi et al. 2020) for learning low-dimensional embeddings
# of Higher-Order Networks. Uses neighborhood matrix construction with
# exponentially decaying powers of the transition matrix, then truncated SVD.
#
# NO random walks or skip-gram - direct matrix factorization.

# ---------------------------------------------------------------------------
# Internal: Build transition matrix from HON adjacency
# ---------------------------------------------------------------------------

#' Build row-normalized transition matrix from HON output
#'
#' @param mat Square weighted adjacency matrix from build_hon().
#' @return Row-stochastic transition matrix.
#' @noRd
.honem_transition_matrix <- function(mat) {
  row_sums <- rowSums(mat)
  nonzero <- row_sums > 0
  result <- mat
  result[nonzero, ] <- result[nonzero, ] / row_sums[nonzero]
  result
}

# ---------------------------------------------------------------------------
# Internal: Compute neighborhood matrix
# ---------------------------------------------------------------------------

#' Compute HONEM neighborhood matrix
#'
#' S = (1/Z) * sum_{k=0}^{L} exp(-k) * D^{k+1}
#' where D is the row-normalized transition matrix and Z is the normalization.
#'
#' @param D Row-stochastic transition matrix.
#' @param max_power Integer. Maximum power L (default: 10 or diameter).
#' @return Dense neighborhood matrix S.
#' @noRd
.honem_neighborhood_matrix <- function(D, max_power = 10L) {
  n <- nrow(D)

  # Normalization constant: Z = sum_{k=0}^{L} exp(-k)
  weights <- exp(-(0L:max_power))
  Z <- sum(weights)

  # Accumulate weighted matrix powers
  S <- matrix(0, nrow = n, ncol = n)
  D_power <- D  # D^1

  vapply(seq_along(weights), function(k_idx) {
    S <<- S + weights[k_idx] * D_power
    D_power <<- D_power %*% D  # D^{k+2}
    0
  }, numeric(1L))

  S / Z
}

# ---------------------------------------------------------------------------
# Internal: Truncated SVD for embeddings
# ---------------------------------------------------------------------------

#' Compute embeddings via truncated SVD
#'
#' @param S Neighborhood matrix.
#' @param dim Integer. Embedding dimension.
#' @return List with embeddings (n x dim matrix), singular_values, explained.
#' @noRd
.honem_svd <- function(S, dim) {
  n <- nrow(S)
  dim <- min(dim, n - 1L)

  sv <- svd(S, nu = dim, nv = dim)

  # Embedding: U * sqrt(Sigma)
  sigma_sqrt <- sqrt(sv$d[seq_len(dim)])
  embeddings <- sweep(sv$u[, seq_len(dim), drop = FALSE], 2, sigma_sqrt, `*`)
  rownames(embeddings) <- rownames(S)
  colnames(embeddings) <- paste0("dim_", seq_len(dim))

  total_var <- sum(sv$d^2)
  explained <- sum(sv$d[seq_len(dim)]^2) / total_var

  list(
    embeddings = embeddings,
    singular_values = sv$d[seq_len(dim)],
    explained_variance = explained
  )
}

# ---------------------------------------------------------------------------
# Main function: build_honem
# ---------------------------------------------------------------------------

#' Build HONEM Embeddings for Higher-Order Networks
#'
#' Constructs low-dimensional embeddings from a Higher-Order Network (HON)
#' that preserve higher-order dependencies. Uses exponentially-decaying matrix
#' powers of the HON transition matrix followed by truncated SVD.
#'
#' HONEM is parameter-free and scalable - no random walks, skip-gram, or
#' hyperparameter tuning required.
#'
#' @param hon A \code{net_hon} object from \code{\link{build_hon}}, or a
#'   square weighted adjacency matrix.
#' @param dim Integer. Embedding dimension (default 32).
#' @param max_power Integer. Maximum walk length for neighborhood computation
#'   (default 10). Higher values capture longer-range structure.
#' @return An object of class \code{net_honem} with components:
#'   \describe{
#'     \item{embeddings}{Numeric matrix (n_nodes x dim) of node embeddings.}
#'     \item{nodes}{Character vector of node names.}
#'     \item{singular_values}{Numeric vector of top singular values.}
#'     \item{explained_variance}{Proportion of variance explained.}
#'     \item{dim}{Embedding dimension used.}
#'     \item{max_power}{Maximum power used.}
#'     \item{n_nodes}{Number of nodes embedded.}
#'   }
#'
#' @references
#' Saebi, M., Ciampaglia, G. L., Kaplan, L. M., & Chawla, N. V. (2020).
#' HONEM: Learning Embedding for Higher Order Networks. \emph{Big Data},
#' 8(4), 255-269. \doi{10.1089/big.2019.0169}
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon_2 <- build_hon(seqs, max_order = 2)
#' hem <- build_honem(hon_2, dim = 2)
#'
#' \donttest{
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
#'               c("B","C","D","A"), c("C","D","A","B"))
#' hon <- build_hon(trajs, max_order = 2)
#' emb <- build_honem(hon, dim = 4)
#' print(emb)
#' plot(emb)
#' }
#'
#' @export
build_honem <- function(hon, dim = 32L, max_power = 10L) {
  dim <- as.integer(dim)
  max_power <- as.integer(max_power)
  stopifnot(
    "'dim' must be >= 1" = dim >= 1L,
    "'max_power' must be >= 1" = max_power >= 1L
  )

  # Extract adjacency matrix
  if (inherits(hon, "net_hon")) {
    mat <- hon$matrix
  } else if (is.matrix(hon) && nrow(hon) == ncol(hon)) {
    mat <- hon
  } else {
    stop("'hon' must be a net_hon object or a square matrix")
  }

  n <- nrow(mat)
  if (n < 2L) stop("Need at least 2 nodes for embedding")

  dim <- min(dim, n - 1L)

  # Build transition matrix
  D <- .honem_transition_matrix(mat)

  # Compute neighborhood matrix
  S <- .honem_neighborhood_matrix(D, max_power)
  rownames(S) <- rownames(mat)
  colnames(S) <- colnames(mat)

  # Truncated SVD
  svd_result <- .honem_svd(S, dim)

  result <- list(
    embeddings = svd_result$embeddings,
    nodes = rownames(mat),
    singular_values = svd_result$singular_values,
    explained_variance = svd_result$explained_variance,
    dim = dim,
    max_power = max_power,
    n_nodes = n
  )

  class(result) <- "net_honem"
  result
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_honem
#'
#' @param x A \code{net_honem} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon_2 <- build_hon(seqs, max_order = 2)
#' hem <- build_honem(hon_2, dim = 2)
#' print(hem)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' hon   <- build_hon(seqs, max_order = 2L)
#' honem <- build_honem(hon, dim = 2L)
#' print(honem)
#' }
#'
#' @export
print.net_honem <- function(x, ...) {
  cat("HONEM: Higher-Order Network Embedding\n")
  cat(sprintf("  Nodes:      %d\n", x$n_nodes))
  cat(sprintf("  Dimensions: %d\n", x$dim))
  cat(sprintf("  Max power:  %d\n", x$max_power))
  cat(sprintf("  Variance explained: %.1f%%\n", x$explained_variance * 100))
  invisible(x)
}

#' Summary Method for net_honem
#'
#' @param object A \code{net_honem} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A data.frame with one row per node: column \code{node} (node
#'   label) followed by \code{dim1}, \code{dim2}, ..., \code{dim}\emph{d}
#'   embedding coordinates, returned visibly; the summary text is printed as
#'   a side effect.
#'
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon_2 <- build_hon(seqs, max_order = 2)
#' hem <- build_honem(hon_2, dim = 2)
#' summary(hem)
#'
#' \donttest{
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon <- build_hon(seqs, max_order = 3)
#' he <- build_honem(hon, dim = 2)
#' summary(he)
#' }
#'
#' @export
summary.net_honem <- function(object, ...) {
  cat("HONEM Summary\n\n")
  cat(sprintf("  Nodes: %d | Dimensions: %d\n", object$n_nodes, object$dim))
  cat(sprintf("  Variance explained: %.1f%%\n",
              object$explained_variance * 100))
  cat(sprintf("  Top singular values: %s\n",
              paste(round(object$singular_values[seq_len(min(5, object$dim))],
                          3), collapse = ", ")))
  cat(sprintf("  Embedding range: [%.3f, %.3f]\n",
              min(object$embeddings), max(object$embeddings)))

  emb <- object$embeddings
  node_labels <- rownames(emb)
  if (is.null(node_labels)) node_labels <- paste0("n", seq_len(nrow(emb)))
  dim_cols <- as.data.frame(emb)
  colnames(dim_cols) <- paste0("dim", seq_len(ncol(emb)))
  invisible(data.frame(node = node_labels, dim_cols,
                       stringsAsFactors = FALSE, row.names = NULL))
}

#' Plot Method for net_honem
#'
#' @param x A \code{net_honem} object.
#' @param dims Integer vector of length 2. Dimensions to plot (default: \code{c(1, 2)}).
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon_2 <- build_hon(seqs, max_order = 2)
#' hem <- build_honem(hon_2, dim = 2)
#' plot(hem)
#'
#' \donttest{
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' hon <- build_hon(seqs, max_order = 3)
#' he <- build_honem(hon, dim = 2)
#' plot(he)
#' }
#'
#' @export
plot.net_honem <- function(x, dims = c(1L, 2L), ...) {
  if (x$dim < 2L) {
    message("Need at least 2 dimensions to plot")
    return(invisible(x))
  }

  emb <- x$embeddings
  d1 <- dims[1L]
  d2 <- dims[2L]

  plot(emb[, d1], emb[, d2], pch = 19, col = "steelblue",
       xlab = sprintf("Dimension %d", d1),
       ylab = sprintf("Dimension %d", d2),
       main = sprintf("HONEM Embedding (%d nodes, %.0f%% var)",
                       x$n_nodes, x$explained_variance * 100),
       ...)

  if (x$n_nodes <= 50) {
    graphics::text(emb[, d1], emb[, d2], labels = x$nodes,
                   pos = 3, cex = 0.7)
  }

  invisible(x)
}

#' Coerce a net_honem to a tidy table
#'
#' @param x A `net_honem` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"embeddings"` (default) for the node coordinates, or
#'   `"variance"` for the per-dimension singular values.
#' @return A data.frame. For `what = "embeddings"`, one row per higher-order
#'   node: `node` followed by one column per embedding dimension (`dim1`,
#'   `dim2`, ...). For `what = "variance"`, one row per dimension: `dim`,
#'   `singular_value`, `proportion` (that dimension's share of the total
#'   squared singular value).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_honem <- function(x, row.names = NULL, optional = FALSE,
                                    ..., what = c("embeddings", "variance"),
                                    top = NULL) {
  what <- match.arg(what)
  if (what == "variance") {
    return(.ho_top(data.frame(
      dim            = seq_along(x$singular_values),
      singular_value = as.numeric(x$singular_values),
      proportion     = as.numeric(x$singular_values^2 /
                                    sum(x$singular_values^2)),
      stringsAsFactors = FALSE
    ), top))
  }
  emb <- x$embeddings
  out <- data.frame(node = x$nodes, stringsAsFactors = FALSE)
  emb_df <- as.data.frame(emb, stringsAsFactors = FALSE)
  names(emb_df) <- paste0("dim", seq_len(ncol(emb)))
  out <- cbind(out, emb_df)
  rownames(out) <- NULL
  .ho_top(out, top)
}
