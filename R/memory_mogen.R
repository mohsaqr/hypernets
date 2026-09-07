# ---- Multi-Order Generative Model (MOGen) ----
#
# Implements multi-order De Bruijn graph models for sequential data.
# Based on Scholtes (2017) "When is a Network a Network?" and
# Gote & Scholtes (2023) "Predicting variable-length paths using MOGen".
#
# Key idea: Build higher-order De Bruijn graphs at orders k=0,1,...,K,
# compute transition matrices, and use AIC/BIC/LRT to select the optimal
# Markov order for the data.

# ---------------------------------------------------------------------------
# Internal: Count k-grams and transitions
# ---------------------------------------------------------------------------

#' Count k-grams and transitions between consecutive k-grams
#'
#' At order k, nodes are k-tuples of states from trajectories.
#' Edges connect consecutive k-tuples (overlapping by k-1 states).
#'
#' @param trajectories List of character vectors.
#' @param k Integer order (k >= 1).
#' @return List with nodes, node_counts, edges (data.frame from/to/weight).
#' @noRd
.mogen_count_kgrams <- function(trajectories, k) {
  stopifnot(k >= 1L)

  results <- lapply(trajectories, function(traj) {
    n <- length(traj)
    if (n < k) return(list(nodes = character(0L), from = character(0L),
                            to = character(0L)))

    # k-grams via sliding window
    kgrams <- vapply(seq_len(n - k + 1L), function(i) {
      paste(traj[i:(i + k - 1L)], collapse = .HON_SEP)
    }, character(1L))

    # Transitions between consecutive k-grams
    if (length(kgrams) > 1L) {
      list(nodes = kgrams, from = kgrams[-length(kgrams)],
           to = kgrams[-1L])
    } else {
      list(nodes = kgrams, from = character(0L), to = character(0L))
    }
  })

  all_nodes <- unlist(lapply(results, `[[`, "nodes"))
  all_from <- unlist(lapply(results, `[[`, "from"))
  all_to <- unlist(lapply(results, `[[`, "to"))

  node_tab <- table(all_nodes)
  nodes <- names(node_tab)
  node_counts <- as.integer(node_tab)
  names(node_counts) <- nodes

  if (length(all_from) == 0L) {
    edges <- data.frame(from = character(0L), to = character(0L),
                        weight = integer(0L), stringsAsFactors = FALSE)
  } else {
    edge_keys <- paste(all_from, all_to, sep = "\x02")
    edge_tab <- table(edge_keys)
    edge_split <- strsplit(names(edge_tab), "\x02", fixed = TRUE)
    edges <- data.frame(
      from = vapply(edge_split, `[`, character(1L), 1L),
      to = vapply(edge_split, `[`, character(1L), 2L),
      weight = as.integer(edge_tab),
      stringsAsFactors = FALSE
    )
  }

  list(nodes = nodes, node_counts = node_counts, edges = edges)
}

# ---------------------------------------------------------------------------
# Internal: Build transition matrix
# ---------------------------------------------------------------------------

#' Build row-stochastic transition matrix from k-gram edge counts
#'
#' @param nodes Character vector of node names.
#' @param edges Data frame with from, to, weight columns.
#' @return Named square matrix (row-stochastic).
#' @noRd
.mogen_transition_matrix <- function(nodes, edges) {
  n <- length(nodes)
  mat <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))

  if (nrow(edges) > 0L) {
    idx <- cbind(match(edges$from, nodes), match(edges$to, nodes))
    mat[idx] <- edges$weight
  }

  row_sums <- rowSums(mat)
  nonzero <- row_sums > 0
  mat[nonzero, ] <- mat[nonzero, ] / row_sums[nonzero]
  mat
}

# ---------------------------------------------------------------------------
# Internal: Order-0 marginal distribution
# ---------------------------------------------------------------------------

#' Compute order-0 marginal distribution over states
#'
#' @param trajectories List of character vectors.
#' @return Named numeric vector of state probabilities.
#' @noRd
.mogen_marginal <- function(trajectories) {
  all_states <- unlist(trajectories)
  tab <- table(all_states)
  probs <- as.numeric(tab) / sum(tab)
  names(probs) <- names(tab)
  probs
}

# ---------------------------------------------------------------------------
# Internal: Log-likelihood computation
# ---------------------------------------------------------------------------

#' Compute log-likelihood of trajectories under k-th order model
#'
#' Uses hierarchical decomposition:
#' - Step 1: order-0 probability (state marginal)
#' - Steps 2..k: use increasing orders (order j for step j)
#' - Steps k+1..l: use order-k transitions
#'
#' @param trajectories List of character vectors.
#' @param k Integer model order.
#' @param trans_mats List of transition matrices (index j+1 = order j).
#'   Index 1 is order-0 marginal (named numeric vector).
#' @return Numeric scalar (log-likelihood).
#' @noRd
.mogen_log_likelihood <- function(trajectories, k, trans_mats) {
  log_eps <- log(.Machine$double.eps)
  p0 <- trans_mats[[1L]]

  sum(vapply(trajectories, function(traj) {
    n <- length(traj)
    if (n == 0L) return(0)

    # Order 0: initial state probability
    ll <- if (traj[1L] %in% names(p0) && p0[traj[1L]] > 0)
      log(p0[traj[1L]]) else log_eps

    if (n < 2L) return(ll)

    # Steps 2..n: use order min(step-1, k)
    step_lls <- vapply(2L:n, function(step) {
      order_used <- min(step - 1L, k)

      if (order_used == 0L) {
        # Order 0: each state drawn from marginal (no transition context)
        if (traj[step] %in% names(p0) && p0[traj[step]] > 0) # nocov start
          log(p0[traj[step]]) else log_eps # nocov end
      } else {
        src_start <- step - order_used
        src_key <- paste(traj[src_start:(step - 1L)], collapse = .HON_SEP)
        tgt_key <- paste(traj[(src_start + 1L):step], collapse = .HON_SEP)

        tm <- trans_mats[[order_used + 1L]]
        if (src_key %in% rownames(tm) && tgt_key %in% colnames(tm)) {
          p <- tm[src_key, tgt_key]
          if (p > 0) log(p) else log_eps
        } else {
          log_eps # nocov
        }
      }
    }, numeric(1L))

    ll + sum(step_lls)
  }, numeric(1L)))
}

# ---------------------------------------------------------------------------
# Internal: Degrees of freedom
# ---------------------------------------------------------------------------

#' Compute degrees of freedom for a row-stochastic transition matrix
#'
#' For each row with m non-zero entries: m - 1 free parameters.
#'
#' @param trans_mat Square transition matrix.
#' @return Integer (total DOF).
#' @noRd
.mogen_layer_dof <- function(trans_mat) {
  row_nonzero <- rowSums(trans_mat > 0)
  as.integer(sum(pmax(row_nonzero - 1L, 0L)))
}

# ---------------------------------------------------------------------------
# Main function: build_mogen
# ---------------------------------------------------------------------------

#' Build Multi-Order Generative Model (MOGen)
#'
#' Constructs higher-order De Bruijn graphs from sequential trajectory data and
#' selects the optimal Markov order using AIC, BIC, or likelihood ratio tests.
#'
#' At order k, nodes are k-tuples of states and edges represent transitions
#' between overlapping k-tuples. The model tests increasingly complex Markov
#' orders and selects the one that best balances fit and parsimony.
#'
#' @param data A data.frame (rows = trajectories, columns = time points) or
#'   a list of character/numeric vectors (one per trajectory).
#' @param max_order Integer. Maximum Markov order to test (default 5).
#'   Must be a whole number; a non-integer value (e.g. \code{2.7}) is an
#'   error rather than being silently truncated.
#' @param criterion Character. Model selection criterion: \code{"aic"}
#'   (default), \code{"bic"}, or \code{"lrt"} (likelihood ratio test).
#' @param lrt_alpha Numeric. Significance threshold for LRT (default 0.01).
#' @return An object of class \code{net_mogen} with components:
#'   \describe{
#'     \item{optimal_order}{Selected optimal Markov order.}
#'     \item{criterion}{Which criterion was used for selection.}
#'     \item{orders}{Integer vector of tested orders (0 to max_order).}
#'     \item{aic}{Named numeric vector of AIC values per order.}
#'     \item{bic}{Named numeric vector of BIC values per order.}
#'     \item{log_likelihood}{Named numeric vector of log-likelihoods.}
#'     \item{dof}{Named integer vector of cumulative DOF per model.}
#'     \item{layer_dof}{Named integer vector of per-layer DOF.}
#'     \item{transition_matrices}{List of transition matrices (index 1 = order 0).}
#'     \item{states}{Unique first-order states.}
#'     \item{n_paths}{Number of trajectories.}
#'     \item{n_observations}{Total number of state observations.}
#'   }
#'
#' @references
#' Scholtes, I. (2017). When is a Network a Network? Multi-Order Graphical
#' Model Selection in Pathways and Temporal Networks. \emph{KDD 2017}.
#'
#' Gote, C. & Scholtes, I. (2023). Predicting variable-length paths in
#' networked systems using multi-order generative models. \emph{Applied
#' Network Science}, 8, 62.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#'
#' \donttest{
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
#'               c("B","C","D","A"), c("C","D","A","B"))
#' m <- build_mogen(trajs, max_order = 3)
#' print(m)
#' plot(m)
#' }
#'
#' @export
build_mogen <- function(data, max_order = 5L, criterion = c("aic", "bic", "lrt"),
                        lrt_alpha = 0.01) {
  data <- .coerce_sequence_input(data)
  criterion <- match.arg(criterion)
  stopifnot(
    "'data' must be a data.frame or list" =
      is.data.frame(data) || is.list(data),
    "'max_order' must be a whole number" =
      is.numeric(max_order) && length(max_order) == 1L &&
      max_order == round(max_order),
    "'max_order' must be >= 1" = max_order >= 1L,
    "'lrt_alpha' must be in (0, 1)" = lrt_alpha > 0 && lrt_alpha < 1
  )
  max_order <- as.integer(max_order)

  trajectories <- .hon_parse_input(data, collapse_repeats = FALSE)
  if (length(trajectories) == 0L) {
    stop("No valid trajectories (each must have at least 2 states)")
  }

  # Cap max_order at max path length - 1
  max_path_len <- max(vapply(trajectories, length, integer(1L)))
  if (max_order >= max_path_len) {
    max_order <- max_path_len - 1L
    message(sprintf("max_order capped at %d (longest path = %d)", max_order,
                    max_path_len))
  }

  n_obs <- sum(vapply(trajectories, length, integer(1L)))

  # --- Order 0: marginal distribution ---
  marginal <- .mogen_marginal(trajectories)
  trans_mats <- vector("list", max_order + 1L)
  count_mats <- vector("list", max_order + 1L)
  trans_mats[[1L]] <- marginal
  count_mats[[1L]] <- marginal * sum(vapply(trajectories, length, integer(1L)))

  layer_dofs <- integer(max_order + 1L)
  layer_dofs[1L] <- length(marginal) - 1L

  # --- Orders 1..max_order: De Bruijn graphs ---
  lapply(seq_len(max_order), function(k) {
    kg <- .mogen_count_kgrams(trajectories, k)
    tm <- .mogen_transition_matrix(kg$nodes, kg$edges)
    # Also store raw count matrix
    cm <- matrix(0L, nrow = length(kg$nodes), ncol = length(kg$nodes),
                 dimnames = list(kg$nodes, kg$nodes))
    if (nrow(kg$edges) > 0L) {
      idx <- cbind(match(kg$edges$from, kg$nodes),
                   match(kg$edges$to, kg$nodes))
      cm[idx] <- kg$edges$weight
    }
    trans_mats[[k + 1L]] <<- tm
    count_mats[[k + 1L]] <<- cm
    layer_dofs[k + 1L] <<- .mogen_layer_dof(tm)
    NULL
  })

  # --- Cumulative DOF and log-likelihoods ---
  cum_dofs <- cumsum(layer_dofs)

  logliks <- vapply(0L:max_order, function(k) {
    .mogen_log_likelihood(trajectories, k, trans_mats[seq_len(k + 1L)])
  }, numeric(1L))

  # --- Information criteria ---
  aics <- 2 * cum_dofs - 2 * logliks
  bics <- log(n_obs) * cum_dofs - 2 * logliks

  order_names <- paste0("order_", 0L:max_order)
  names(aics) <- order_names
  names(bics) <- order_names
  names(logliks) <- order_names
  names(cum_dofs) <- order_names
  names(layer_dofs) <- order_names

  # --- Select optimal order ---
  if (criterion == "aic") {
    optimal_order <- which.min(aics) - 1L
  } else if (criterion == "bic") {
    optimal_order <- which.min(bics) - 1L
  } else {
    # Likelihood ratio test: sequential testing
    optimal_order <- 0L
    for (k in seq_len(max_order)) {
      x <- -2 * (logliks[k] - logliks[k + 1L])
      df_diff <- layer_dofs[k + 1L]
      if (df_diff > 0L && x > 0) {
        p_val <- stats::pchisq(x, df = df_diff, lower.tail = FALSE)
        if (p_val < lrt_alpha) {
          optimal_order <- k
        } else {
          break
        }
      }
    }
  }

  # Extract optimal-order transition matrix for cograph compatibility
  opt_mat <- trans_mats[[optimal_order + 1L]]
  if (is.null(dim(opt_mat))) {
    # Order 0: marginal vector -> 1xn matrix
    opt_mat <- matrix(opt_mat, nrow = 1,
                      dimnames = list("marginal", names(opt_mat)))
  }
  opt_nodes <- rownames(opt_mat) %||% names(marginal)
  cg <- .ho_cograph_fields(opt_mat, opt_nodes, method = "mogen")

  result <- list(
    optimal_order = as.integer(optimal_order),
    criterion = criterion,
    orders = 0L:max_order,
    aic = aics,
    bic = bics,
    log_likelihood = logliks,
    dof = cum_dofs,
    layer_dof = layer_dofs,
    transition_matrices = trans_mats,
    count_matrices = count_mats,
    states = names(marginal),
    n_paths = length(trajectories),
    n_observations = n_obs,
    weights = cg$weights,
    nodes = cg$nodes,
    edges = cg$edges,
    directed = TRUE,
    n_nodes = cg$n_nodes,
    n_edges = cg$n_edges,
    meta = cg$meta,
    node_groups = NULL
  )

  class(result) <- c("net_mogen", "cograph_network")
  result
}

# ---------------------------------------------------------------------------
# Utility: Extract transitions at any order
# ---------------------------------------------------------------------------

#' Extract Transition Table from a MOGen Model
#'
#' Returns a data frame of all transitions at a given Markov order,
#' sorted by count (descending). Each row shows the full path as a readable
#' sequence of states, along with the observed count and transition probability.
#'
#' At order k, each edge in the De Bruijn graph represents a (k+1)-step path.
#' For example, at order 2, the edge from node "AI -> FAIL" to node
#' "FAIL -> SOLVE" represents the three-step path AI -> FAIL -> SOLVE.
#' The \code{path} column reconstructs this full sequence for readability.
#'
#' @param x A \code{net_mogen} object from \code{build_mogen()}.
#' @param order Integer. Which order's transitions to extract. Must be a
#'   whole number; a non-integer value is an error rather than being
#'   silently truncated. Defaults to the optimal order selected by the model.
#' @param min_count Integer. Minimum observed count to include (default 1).
#'   Use this to filter out rare transitions that have unreliable probabilities.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return A data frame with columns:
#'   \describe{
#'     \item{path}{The full state sequence (e.g., "AI -> FAIL -> SOLVE").}
#'     \item{count}{Number of times this transition was observed.}
#'     \item{probability}{Transition probability P(to | from).}
#'     \item{from}{The context / conditioning states (k-gram source node).}
#'     \item{to}{The predicted next state.}
#'   }
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#' mogen_transitions(mg, order = 1)
#'
#' \donttest{
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
#'               c("B","C","D","A"), c("C","D","A","B"))
#' m <- build_mogen(trajs, max_order = 3)
#' mogen_transitions(m, order = 1)
#' }
#'
#' @export
mogen_transitions <- function(x, order = NULL, min_count = 1L, top = NULL) {
  stopifnot(inherits(x, "net_mogen"))
  if (is.null(order)) order <- x$optimal_order
  stopifnot(
    "'order' must be a whole number" =
      is.numeric(order) && length(order) == 1L &&
      order == round(order)
  )
  order <- as.integer(order)
  min_count <- as.integer(min_count)
  stopifnot(
    "'order' must be >= 1" = order >= 1L,
    "'order' exceeds max tested order" = order <= max(x$orders)
  )

  tm <- x$transition_matrices[[order + 1L]]
  cm <- x$count_matrices[[order + 1L]]

  # Find edges with sufficient count
  idx <- which(cm >= min_count, arr.ind = TRUE)

  if (nrow(idx) == 0L) {
    return(data.frame(path = character(0L), count = integer(0L),
                      probability = numeric(0L), from = character(0L),
                      to = character(0L), stringsAsFactors = FALSE))
  }

  from_raw <- rownames(tm)[idx[, 1L]]
  to_raw <- colnames(tm)[idx[, 2L]]

  # Reconstruct full path and extract context/next_state
  parsed <- lapply(seq_len(nrow(idx)), function(i) {
    from_parts <- strsplit(from_raw[i], .HON_SEP, fixed = TRUE)[[1L]]
    to_parts <- strsplit(to_raw[i], .HON_SEP, fixed = TRUE)[[1L]]
    next_state <- to_parts[length(to_parts)]
    list(
      path = paste(c(from_parts, next_state), collapse = " -> "),
      from = paste(from_parts, collapse = " -> "),
      to = next_state
    )
  })

  result <- data.frame(
    path = vapply(parsed, `[[`, character(1L), "path"),
    count = as.integer(cm[idx]),
    probability = round(tm[idx], 4),
    from = vapply(parsed, `[[`, character(1L), "from"),
    to = vapply(parsed, `[[`, character(1L), "to"),
    stringsAsFactors = FALSE
  )
  result <- result[order(-result$count), ]
  rownames(result) <- NULL
  .ho_top(result, top)
}

# ---------------------------------------------------------------------------
# Utility: Count path frequencies
# ---------------------------------------------------------------------------

#' Count Path Frequencies in Trajectory Data
#'
#' Counts the frequency of k-step paths (k-grams) across all trajectories.
#' Useful for understanding which sequences dominate the data before applying
#' formal models.
#'
#' @param data A list of character vectors (trajectories) or a data.frame
#'   (rows = trajectories, columns = time points).
#' @param k Integer. Length of the path / n-gram (default 2). A k of 2 counts
#'   individual transitions; k of 3 counts two-step paths, etc. Must be a
#'   whole number; a non-integer value is an error rather than being
#'   silently truncated.
#' @param top Integer or `NULL`. Return only the first `top` rows of the
#'   count-sorted table (default `NULL` returns every path).
#' @return A data frame with columns: \code{path}, \code{count},
#'   \code{proportion}.
#'
#' @examples
#' trajs <- list(c("A","B","C","D"), c("A","B","D","C"))
#' path_counts(trajs, k = 2)
#'
#' \donttest{
#' path_counts(trajs, k = 3, top = 10)
#' }
#'
#' @export
path_counts <- function(data, k = 2L, top = NULL) {
  data <- .coerce_sequence_input(data)
  stopifnot(
    "'k' must be a whole number" =
      is.numeric(k) && length(k) == 1L && k == round(k)
  )
  k <- as.integer(k)
  stopifnot("'k' must be >= 2" = k >= 2L)

  # Normalize input
  if (is.data.frame(data)) {
    trajectories <- lapply(seq_len(nrow(data)), function(i) {
      as.character(unlist(data[i, ], use.names = FALSE))
    })
  } else {
    trajectories <- lapply(data, as.character)
  }

  all_grams <- unlist(lapply(trajectories, function(traj) {
    traj <- traj[!is.na(traj)]
    n <- length(traj)
    if (n < k) return(character(0L))
    vapply(seq_len(n - k + 1L), function(i) {
      paste(traj[i:(i + k - 1L)], collapse = " -> ")
    }, character(1L))
  }))

  tbl <- sort(table(all_grams), decreasing = TRUE)
  result <- data.frame(
    path = names(tbl),
    count = as.integer(tbl),
    proportion = round(as.numeric(tbl) / sum(tbl), 4),
    stringsAsFactors = FALSE
  )
  rownames(result) <- NULL

  .ho_top(result, top)
}

# NOTE: state_frequencies() is deliberately NOT here. It is defined in
# Nestimate's R/mogen.R, is on the htna downstream contract
# (Nestimate/tests/testthat/test-contract-htna.R), and was carved out at the
# hypernets T0 move on purpose -- see "Carve-out discovered during T0" in
# Nestimate/HYPERNETS-DELEGATION-PLAN.md. Exporting it here would mask
# Nestimate's when both packages are attached. At T2, Nestimate relocates it
# out of mogen.R rather than delegating it.

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_mogen
#'
#' @param x A \code{net_mogen} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#' print(mg)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' mog <- build_mogen(seqs, max_order = 2L)
#' print(mog)
#' }
#'
#' @export
print.net_mogen <- function(x, ...) {
  cat("Multi-Order Generative Model (MOGen)\n")
  cat(sprintf("  Optimal order:  %d (by %s)\n", x$optimal_order, x$criterion))
  cat(sprintf("  Orders tested:  0 to %d\n", max(x$orders)))
  cat(sprintf("  States:         %d\n", length(x$states)))
  cat(sprintf("  Paths:          %d (%d observations)\n",
              x$n_paths, x$n_observations))

  aic_min <- which.min(x$aic) - 1L
  bic_min <- which.min(x$bic) - 1L

  fmt_ic <- function(vals, opt_idx) {
    parts <- sprintf("%.1f", vals)
    parts[opt_idx] <- paste0(parts[opt_idx], "*")
    paste(parts, collapse = " | ")
  }

  cat(sprintf("  AIC:            %s\n", fmt_ic(x$aic, aic_min + 1L)))
  cat(sprintf("  BIC:            %s\n", fmt_ic(x$bic, bic_min + 1L)))
  if (aic_min == bic_min) {
    cat(sprintf("  (* = minimum;  AIC and BIC agree on order %d)\n", aic_min))
  } else {
    cat(sprintf("  (* = minimum;  AIC favors order %d, BIC favors order %d)\n",
                aic_min, bic_min))
  }
  invisible(x)
}

#' Summary Method for net_mogen
#'
#' @param object A \code{net_mogen} object.
#' @param ... Additional arguments (ignored).
#'
#' @return A per-order model-selection data.frame with columns \code{order},
#'   \code{layer_dof}, \code{cum_dof}, \code{loglik}, \code{aic}, \code{bic},
#'   \code{best} (\code{"AIC"}/\code{"BIC"}/\code{"AIC+BIC"} marker) and
#'   \code{selected} (\code{"<--"} on the chosen order), returned visibly;
#'   the summary text is printed as a side effect.
#'
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#' summary(mg)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' mog <- build_mogen(seqs, max_order = 2L)
#' summary(mog)
#' }
#'
#' @export
summary.net_mogen <- function(object, ...) {
  cat("Multi-Order Generative Model (MOGen) Summary\n\n")
  cat(sprintf("  States: %s\n", paste(object$states, collapse = ", ")))
  cat(sprintf("  Paths: %d | Observations: %d\n",
              object$n_paths, object$n_observations))

  aic_min <- which.min(object$aic) - 1L
  bic_min <- which.min(object$bic) - 1L
  cat(sprintf("  Best by AIC: order %d  |  Best by BIC: order %d\n",
              aic_min, bic_min))
  cat(sprintf("  Selected:    order %d (by %s)\n\n",
              object$optimal_order, object$criterion))

  res <- data.frame(
    order     = object$orders,
    layer_dof = as.integer(object$layer_dof),
    cum_dof   = as.integer(object$dof),
    loglik    = round(object$log_likelihood, 2),
    aic       = round(object$aic, 2),
    bic       = round(object$bic, 2),
    stringsAsFactors = FALSE
  )

  res$best <- ""
  res$best[aic_min + 1L] <- "AIC"
  res$best[bic_min + 1L] <- if (res$best[bic_min + 1L] == "") "BIC" else
    paste0(res$best[bic_min + 1L], "+BIC")
  sel_idx <- object$optimal_order + 1L
  res$selected <- ""
  res$selected[sel_idx] <- "<--"

  invisible(res)
}

#' Plot Method for net_mogen
#'
#' @param x A \code{net_mogen} object.
#' @param type Character. Plot type: \code{"ic"} (default) or \code{"likelihood"}.
#' @param ... Additional arguments passed to \code{\link[graphics]{plot}}.
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
#' mg <- build_mogen(seqs, max_order = 2)
#' plot(mg)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' mog <- build_mogen(seqs, max_order = 2L)
#' plot(mog, type = "ic")
#' }
#'
#' @export
plot.net_mogen <- function(x, type = c("ic", "likelihood"), ...) {
  type <- match.arg(type)

  if (type == "ic") {
    orders <- x$orders
    aic_vals <- x$aic
    bic_vals <- x$bic

    ylim <- range(c(aic_vals, bic_vals))
    plot(orders, aic_vals, type = "b", pch = 19, col = "steelblue",
         xlab = "Markov Order", ylab = "Information Criterion",
         main = "MOGen: Order Selection", ylim = ylim, xaxt = "n", ...)
    graphics::axis(1, at = orders)
    graphics::lines(orders, bic_vals, type = "b", pch = 17, col = "coral")
    graphics::legend("topright", legend = c("AIC", "BIC"),
                     col = c("steelblue", "coral"), pch = c(19, 17),
                     lty = 1, bty = "n")
    graphics::abline(v = x$optimal_order, lty = 2, col = "gray40")
    graphics::text(x$optimal_order, ylim[1], sprintf("k*=%d", x$optimal_order),
                   pos = 4, col = "gray40")
  } else {
    orders <- x$orders
    plot(orders, x$log_likelihood, type = "b", pch = 19, col = "steelblue",
         xlab = "Markov Order", ylab = "Log-Likelihood",
         main = "MOGen: Log-Likelihood by Order", xaxt = "n", ...)
    graphics::axis(1, at = orders)
    graphics::abline(v = x$optimal_order, lty = 2, col = "gray40")
  }

  invisible(x)
}

#' Coerce a net_mogen to a tidy table
#'
#' @param x A `net_mogen` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"orders"` (default) for the per-order model-selection table,
#'   or `"transitions"` for the fitted transition probabilities (the same
#'   table [mogen_transitions()] returns).
#' @param order Integer or `NULL`. Only for `what = "transitions"`: which
#'   layer's transitions to return. Default `NULL` uses the selected order.
#' @return A data.frame. For `what = "orders"`, one row per candidate order:
#'   `order`, `log_likelihood`, `aic`, `bic`, `dof`, `layer_dof`, and
#'   `optimal` (logical, `TRUE` for the selected order). For
#'   `what = "transitions"`, one row per transition.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_mogen <- function(x, row.names = NULL, optional = FALSE,
                                    ..., what = c("orders", "transitions"),
                                    order = NULL, top = NULL) {
  what <- match.arg(what)
  if (what == "transitions") {
    return(mogen_transitions(x, order = order %||% x$optimal_order,
                             top = top))
  }
  out <- data.frame(
    order          = as.integer(x$orders),
    log_likelihood = as.numeric(x$log_likelihood),
    aic            = as.numeric(x$aic),
    bic            = as.numeric(x$bic),
    dof            = as.integer(x$dof),
    layer_dof      = as.integer(x$layer_dof),
    stringsAsFactors = FALSE
  )
  out$optimal <- out$order == x$optimal_order
  rownames(out) <- NULL
  .ho_top(out, top)
}
