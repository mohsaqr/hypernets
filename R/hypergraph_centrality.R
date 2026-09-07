# ---- Hypergraph eigenvector centralities (HON-5) -------------------------
# Three variants from Benson (2019, arXiv:1807.09644):
#   * clique-motif ("CEC") - standard eigenvector centrality on the
#                            clique-expanded pairwise graph
#   * Z-eigenvector ("Z")  - linear tensor eigenvector
#   * H-eigenvector ("H")  - H-eigenvector (power-k-1 recurrence)
#
# All three solved by power iteration on the hyperedge list.

#' Hypergraph eigenvector centralities
#'
#' Computes one or more eigenvector-style centralities on a
#' [net_hypergraph][build_hypergraph]: *clique-motif* (CEC),
#' *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC). Each variant
#' captures influence differently - CEC flattens group structure via
#' clique expansion, while ZEC and HEC propagate through the
#' higher-order groups directly.
#'
#' @param hg A `net_hypergraph` (from [build_hypergraph()],
#'   [group_hypergraph()], or [window_hypergraph()]).
#' @param type Character vector, any subset of
#'   `c("clique", "Z", "H", "pagerank", "subhypergraph")`. The default computes the three
#'   eigenvector variants; request `"pagerank"` explicitly.
#' @param max_iter Maximum number of power-iteration steps. Default
#'   `1000`.
#' @param tol Convergence tolerance on the L1 change between successive
#'   iterates. Default `1e-8`.
#' @param normalize Logical. If `TRUE` (default), each returned
#'   centrality vector is L2-normalized to unit norm (compatible with
#'   `igraph::eigen_centrality()`'s scale for type `"clique"`). Does not
#'   apply to `"pagerank"`, which always sums to 1, or
#'   `"subhypergraph"`, which is returned on its natural log scale.
#' @param damping Single numeric in (0, 1). PageRank damping factor
#'   (probability of following the walk rather than teleporting).
#'   Default `0.85`. Only used by `type = "pagerank"`.
#' @param edge_weights NULL, a single positive number (recycled to every
#'   hyperedge, e.g. `1` for unit weights), or a positive numeric vector,
#'   one per hyperedge. Hyperedge weights of the EDVW random walk behind
#'   `type = "pagerank"`; `NULL` defaults to the window counts for
#'   hypergraphs built by [window_hypergraph()], else the Hayashi et al.
#'   dispersion heuristic (unit weights on a binary incidence). Only used
#'   by `type = "pagerank"`.
#'
#' @param sort_by `NULL` (node order, default) or one of the requested
#'   `type`s - return the table ranked by that centrality, largest first,
#'   ties broken by node name so the order is deterministic.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#'
#' @return A data.frame, one row per node: `node`, followed by one column
#'   per requested `type`, in the order requested. Compare variants by
#'   reading across the columns; rank by one with `sort_by =`.
#'
#' @details
#' **Clique-motif eigenvector centrality (CEC)**: forms the
#' clique-expanded pairwise graph \eqn{W} where
#' \eqn{W_{ij} = |\{e : i, j \in e\}|} and returns the leading
#' eigenvector of \eqn{W}. Equivalent to running
#' `igraph::eigen_centrality()` on [clique_expansion()] output.
#'
#' **Z-eigenvector centrality (ZEC)**: solves the linear
#' eigen-equation on the hyperedge tensor,
#' \deqn{\lambda\, x_i \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j,}
#' via power iteration. Works for hypergraphs with mixed edge sizes.
#'
#' **H-eigenvector centrality (HEC)**: solves the power-k-1
#' eigen-equation,
#' \deqn{\lambda\, x_i^{k-1} \;=\; \sum_{e \ni i}\; \prod_{j \in e,\; j \neq i} x_j.}
#' For uniform hypergraphs (all hyperedges of size \eqn{k}), this is
#' equivalent to normalizing the ZEC update by the geometric-mean
#' exponent \eqn{1/(k-1)}. For mixed sizes, the effective exponent is
#' taken from the largest hyperedge; expect slightly different rankings
#' from ZEC in the mixed case.
#'
#' **Hypergraph PageRank** (`"pagerank"`): the stationary distribution of
#' the damped EDVW random walk of Chitra & Raphael (2019): from node
#' \eqn{v}, pick a hyperedge \eqn{e \ni v} with probability proportional
#' to its weight \eqn{w(e)}, then a node \eqn{u \in e} with probability
#' proportional to its edge-dependent vertex weight \eqn{\gamma_e(u)}
#' (the incidence cell, i.e. occurrence totals for
#' [window_hypergraph()]); with probability \eqn{1 - damping} teleport
#' uniformly. Their collapse theorem: with edge-*independent* vertex
#' weights (a binary incidence) the walk is equivalent to PageRank on the
#' weighted clique expansion with edge weights
#' \eqn{\sum_{e \ni u,v} w(e)/\delta(e)} -- the hypergraph adds
#' information exactly when \eqn{\gamma} is edge-dependent. Nodes left in
#' no hyperedge (possible after `min_weight`/`min_size` filtering)
#' teleport from every step and receive only teleportation mass. The
#' undamped stationary distribution of the same walk is the `pi` column
#' reported by [hypergraph_cluster()].
#'
#' **Subhypergraph centrality** (`"subhypergraph"`): the logarithm of the
#' diagonal of the matrix exponential of the clique adjacency derived from
#' binary incidence (so entries count shared hyperedges),
#' \eqn{\log[\exp(W)]_{ii}}. It counts closed walks based at each node with a
#' factorial penalty for length and matches the implementation used by
#' HypergraphX 1.5 in the legal-hypergraphs analysis.
#'
#' @seealso [build_hypergraph()], [clique_expansion()],
#'   [hypergraph_measures()].
#'
#' @examples
#' df <- data.frame(
#'   member  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
#' )
#' hg <- group_hypergraph(df, "member", "session")
#'
#' # One row per node, one column per variant - compare across the columns
#' hypergraph_centrality(hg)
#'
#' # Ranked by one of them
#' hypergraph_centrality(hg, sort_by = "clique")
#'
#' # PageRank of the EDVW random walk
#' hypergraph_centrality(hg, type = "pagerank")
#'
#' @references
#' Benson, A. R. (2019). Three hypergraph eigenvector centralities.
#' \emph{SIAM Journal on Mathematics of Data Science} 1(2), 293-312.
#' arXiv:1807.09644.
#'
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. \emph{Proceedings of the 36th
#' International Conference on Machine Learning}, PMLR 97, 1172-1181.
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{Proceedings of CIKM
#' 2020}, 495-504. \doi{10.1145/3340531.3412034}
#'
#' Estrada, E., & Rodriguez-Velazquez, J. A. (2005). Complex networks as
#' hypergraphs. *arXiv preprint physics/0505137*.
#'
#' @note The `"clique"`, `"Z"`, and `"H"` variants have an independent
#'   XGI equivalence runner in `local_testing_and_equivalence/`; clique is
#'   also validated against `igraph::eigen_centrality` (cosine ~ 1).
#'   `"pagerank"` is checked against
#'   `igraph::page_rank` on the collapsed graph of the edge-independent
#'   case plus a dense linear-system solve. Clean-room list-based tensor
#'   equations remain in the unit suite for deterministic validation without
#'   Python.
#'
#' @export
hypergraph_centrality <- function(hg,
                                   type     = c("clique", "Z", "H"),
                                   max_iter = 1000L,
                                   tol      = 1e-8,
                                   normalize = TRUE,
                                   damping = 0.85,
                                   edge_weights = NULL,
                                   sort_by = NULL,
                                   top = NULL) {
  stopifnot(
    inherits(hg, "net_hypergraph"),
    is.numeric(max_iter), length(max_iter) == 1L, max_iter > 0,
    is.numeric(tol), length(tol) == 1L, tol > 0,
    is.logical(normalize), length(normalize) == 1L,
    "`damping` must be a single number strictly between 0 and 1" =
      is.numeric(damping) && length(damping) == 1L && is.finite(damping) &&
      damping > 0 && damping < 1
  )
  type <- match.arg(type, choices = c("clique", "Z", "H", "pagerank",
                                      "subhypergraph"),
                    several.ok = TRUE)

  n     <- hg$n_nodes
  nodes <- hg$nodes

  # Degenerate: no nodes -> a zero-row table with the requested columns
  if (n == 0L) {
    res <- data.frame(node = character(0L), stringsAsFactors = FALSE)
    res[type] <- rep(list(numeric(0L)), length(type))
    return(res)
  }

  # Degenerate: no hyperedges -> zero eigenvector centralities; PageRank
  # is pure teleportation, so it is uniform (it must still sum to 1)
  if (hg$n_hyperedges == 0L) {
    res <- data.frame(node = nodes, stringsAsFactors = FALSE)
    res[type] <- rep(list(rep(0, n)), length(type))
    if ("pagerank" %in% type) res$pagerank <- rep(1 / n, n)
    return(.ho_top(res, top))
  }

  # Shared initialization: uniform positive vector
  x0 <- rep(1 / sqrt(n), n)

  # Hyperedges as integer vector lists; pre-compute for Z/H
  hyperedges <- hg$hyperedges
  edge_sizes <- vapply(hyperedges, length, integer(1L))
  k_max      <- max(edge_sizes)

  out <- list()

  # ---- CEC: power iteration on clique-expansion W ----
  if ("clique" %in% type) {
    B_bin <- (hg$incidence > 0) * 1.0
    W <- tcrossprod(B_bin)
    diag(W) <- 0
    x <- x0
    for (iter in seq_len(max_iter)) {
      y <- as.numeric(W %*% x)
      nrm <- sqrt(sum(y^2))
      if (nrm == 0) break
      y <- y / nrm
      if (sum(abs(y - x)) < tol) break
      x <- y
    }
    # Sign convention: positive entries (W is non-negative so Perron vec is positive)
    if (any(x != 0) && sum(x) < 0) x <- -x
    if (!normalize && any(x != 0)) x <- x / max(abs(x))
    out$clique <- stats::setNames(x, nodes)
  }

  # ---- ZEC: lambda x = sum_{e contains i} prod_{j in e,j!=i} x_j ----
  if ("Z" %in% type) {
    out$Z <- stats::setNames(
      .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                            exponent = 1L, max_iter = max_iter, tol = tol,
                            x0 = x0, normalize = normalize),
      nodes
    )
  }

  # ---- HEC: lambda x^{k-1} = sum_{e contains i} prod_{j in e,j!=i} x_j ----
  if ("H" %in% type) {
    out$H <- stats::setNames(
      .hg_tensor_power_iter(hyperedges, edge_sizes, n,
                            exponent = k_max - 1L,
                            max_iter = max_iter, tol = tol,
                            x0 = x0, normalize = normalize),
      nodes
    )
  }

  # ---- PageRank: damped EDVW random walk (Chitra & Raphael 2019) ----
  if ("pagerank" %in% type) {
    rw <- .hl_rw_transition(hg, edge_weights)
    P <- rw$P
    # Nodes in no hyperedge (d_v = 0) have undefined walk rows: they
    # teleport uniformly from every step (standard dangling-node fix)
    dangling <- rw$d_v == 0
    if (any(dangling)) P[dangling, ] <- 1 / n
    x <- rep(1 / n, n)
    converged <- FALSE
    # power iteration: an inherently sequential fixed-point loop
    for (iter in seq_len(max_iter)) {
      y <- damping * as.numeric(crossprod(P, x)) + (1 - damping) / n
      if (sum(abs(y - x)) < tol) {
        x <- y
        converged <- TRUE
        break
      }
      x <- y
    }
    if (!converged) {
      warning(warningCondition(
        sprintf(
          "PageRank did not converge in %d iterations (L1 change > %g).",
          as.integer(max_iter), tol),
        class = "hypernets_no_converge"
      ))
    }
    out$pagerank <- stats::setNames(x, nodes)
  }

  # ---- Log subhypergraph centrality (Estrada & Rodriguez-Velazquez) ----
  if ("subhypergraph" %in% type) {
    b <- (hg$incidence != 0) * 1
    adjacency <- as.matrix(tcrossprod(b))
    diag(adjacency) <- 0
    eig <- eigen(adjacency, symmetric = TRUE)
    shift <- max(eig$values)
    terms <- sweep(eig$vectors^2, 2L, exp(eig$values - shift), `*`)
    value <- shift + log(rowSums(terms))
    out$subhypergraph <- stats::setNames(value, nodes)
  }

  # One row per node, one column per requested type, in the order the user
  # asked for them.
  res <- data.frame(node = nodes, stringsAsFactors = FALSE)
  res[type] <- lapply(out[type], unname)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, type)
    res <- res[order(-res[[sort_by]], res$node), , drop = FALSE]
    rownames(res) <- NULL
  }
  .ho_top(res, top)
}

# Shared tensor-power-iteration kernel.
# Uses Kolda-Mayo SSHOPM shift: x_{k+1} ~ f(x_k) + shift * x_k
# which guarantees monotone convergence for non-negative tensors
# (Chang, Pearson & Zhang 2009 / Kolda & Mayo 2011).
#
# exponent = 1    => Z-eigenvector (no post-root)
# exponent = k-1  => H-eigenvector (k-1-root)
#' @noRd
.hg_tensor_power_iter <- function(hyperedges, edge_sizes, n, exponent,
                                   max_iter, tol, x0, normalize,
                                   shift = 1) {
  x <- x0
  for (iter in seq_len(max_iter)) {
    y <- numeric(n)
    for (e_idx in seq_along(hyperedges)) {
      e  <- hyperedges[[e_idx]]
      ke <- edge_sizes[e_idx]
      if (ke < 2L) next
      x_e <- x[e]
      # Product over all members excluding each i, done in O(k) per edge
      # via total-product / x_j (handling zeros explicitly)
      total <- prod(x_e)
      if (any(x_e == 0)) {
        nz_idx <- which(x_e != 0)
        if (length(nz_idx) == ke - 1L) {
          zero_pos <- setdiff(seq_len(ke), nz_idx)
          y[e[zero_pos]] <- y[e[zero_pos]] + prod(x_e[nz_idx])
        }
        # With >=2 zeros in x_e, every leave-one-out product is zero.
      } else {
        for (idx in seq_len(ke)) {
          y[e[idx]] <- y[e[idx]] + total / x_e[idx]
        }
      }
    }

    if (exponent > 1L) {
      # k-1-root, preserving sign (for non-negative hypergraphs all >= 0)
      y <- sign(y) * abs(y)^(1 / exponent)
    }

    # SSHOPM shift: bias towards x to stabilize oscillating tensor iterates
    y <- y + shift * x

    nrm <- sqrt(sum(y^2))
    if (nrm == 0) {
      x <- y
      break
    }
    y <- y / nrm

    if (sum(abs(y - x)) < tol) {
      x <- y
      break
    }
    x <- y
  }

  if (any(x != 0) && sum(x) < 0) x <- -x
  if (!normalize && any(x != 0)) x <- x / max(abs(x))
  x
}
