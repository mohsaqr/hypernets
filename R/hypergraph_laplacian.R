# ---- Hypergraph Laplacians, Spectral Clustering, Transduction ----------
#
# Zhou, Huang & Scholkopf (2006): normalized hypergraph Laplacian
#   L = I - D_v^{-1/2} H W D_e^{-1} H^T D_v^{-1/2}
# on the BINARY incidence pattern H with hyperedge weights w(e).
#
# Hayashi, Aksoy, Park & Park (2020): random-walk Laplacian with
# edge-dependent vertex weights (EDVW). The WEIGHTED incidence cells
# gamma_e(v) drive the walk (pick e ~ w(e), then v ~ gamma_e(v)); the
# non-reversible walk is symmetrized through its stationary distribution
# with Chung's (2005) directed Laplacian. With a binary incidence and unit
# hyperedge weights the two Laplacians coincide (the walk is reversible
# with pi ~ d(v)) - pinned by an invariant test.
#
# Reference implementation used as the equivalence oracle: HyperNetX
# `hypernetx/algorithms/clustering/laplacians_clustering.py` (PNNL;
# author-adjacent). Its hyperedge-weight default - the population standard
# deviation of the edge's vertex weights plus one - is reproduced here as
# `edge_weights = NULL` for `type = "random_walk"`.

#' Normalized hypergraph Laplacian
#'
#' Computes the normalized Laplacian of a hypergraph, either the classic
#' Zhou-Huang-Scholkopf form on the binary incidence pattern
#' (`type = "zhou"`) or the random-walk form with edge-dependent vertex
#' weights (`type = "random_walk"`), in which the weighted incidence cells
#' (e.g. the summed weights produced by [group_hypergraph()]) determine
#' where a random walker lands inside a hyperedge, and the resulting
#' non-reversible walk is symmetrized through its stationary distribution
#' (Chung 2005). Both Laplacians are symmetric positive semi-definite with
#' eigenvalues in `[0, 2]`; for a binary incidence with unit hyperedge
#' weights they coincide.
#'
#' @param hg A `net_hypergraph` from [build_hypergraph()] or
#'   [group_hypergraph()]. Must be connected and have at least one
#'   hyperedge.
#' @param type Character. `"zhou"` (default) for the Zhou et al. (2006)
#'   normalized Laplacian on the binary incidence pattern, or
#'   `"random_walk"` for the Hayashi et al. (2020) EDVW random-walk
#'   Laplacian on the weighted incidence.
#' @param edge_weights A single positive number (recycled) or a numeric vector of positive hyperedge weights
#'   (length `hg$n_hyperedges`), or `NULL` for the default. Hypergraphs
#'   built by [window_hypergraph()] default to their window counts (for
#'   both types). Otherwise the default is type-specific: unit weights for
#'   `"zhou"`; for `"random_walk"` the Hayashi et al. heuristic - the
#'   population standard deviation of each hyperedge's (non-zero) vertex
#'   weights plus one - which reduces to unit weights on a binary
#'   incidence.
#'
#' @return A symmetric `n_nodes` x `n_nodes` numeric matrix (node names as
#'   dimnames) with attributes `type` (the Laplacian type),
#'   `pi` (named stationary distribution of the underlying random walk)
#'   and `edge_weights` (the hyperedge weights actually used).
#'
#' @references
#' Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
#' Clustering, classification, and embedding. \emph{NeurIPS 19}.
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{CIKM 2020}, 495-504.
#' \doi{10.1145/3340531.3412034}
#'
#' Chung, F. (2005). Laplacians and the Cheeger inequality for directed
#' graphs. \emph{Annals of Combinatorics}, 9(1), 1-19.
#'
#' @examples
#' events <- data.frame(
#'   person = c("a", "b", "c", "a", "b", "d", "c", "d", "e", "e", "a"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
#'               "m4", "m4"),
#'   hours = c(2, 1, 1, 3, 2, 1, 2, 2, 4, 1, 1)
#' )
#' hg <- group_hypergraph(events, actor = "person", group = "meeting",
#'                        weight = "hours")
#' L <- hypergraph_laplacian(hg, type = "random_walk")
#' range(eigen(L, symmetric = TRUE, only.values = TRUE)$values)
#'
#' @export
hypergraph_laplacian <- function(hg,
                                 type = c("zhou", "random_walk"),
                                 edge_weights = NULL) {
  type <- match.arg(type)
  .hl_validate_hg(hg)
  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  L <- parts$L
  attr(L, "type") <- type
  attr(L, "pi") <- parts$pi
  attr(L, "edge_weights") <- parts$w
  L
}

#' Spectral or symmetric-NMF clustering of hypergraph vertices
#'
#' Partitions the nodes of a hypergraph into `k` clusters with either of
#' Hayashi et al.'s (2020) representative-digraph algorithms. `algorithm =
#' "spectral"` (RDC-Spec) row-normalizes the `k` smallest Laplacian
#' eigenvectors and applies k-means. `algorithm = "symnmf"` (RDC-Sym)
#' computes a rank-`k` non-negative factorization `T ~= U U'` of the
#' normalized similarity `T = I - L`, then assigns each vertex to the
#' largest entry in its row of `U`, exactly as Algorithm 2 specifies.
#' With `type = "random_walk"` and a weighted
#' incidence (e.g. from [group_hypergraph()] with `weight =`), the
#' edge-dependent vertex weights genuinely change the partition - with
#' edge-independent weights the walk collapses to a graph random walk
#' (Chitra & Raphael 2019).
#'
#' Both solvers are stochastic: `nstart` initializations are used and a
#' `seed` fixes the result. Report stability across seeds for consequential
#' results.
#'
#' @param hg A connected `net_hypergraph`.
#' @param k Integer number of clusters, `2 <= k <= n_nodes - 1`.
#' @param type,edge_weights Passed to [hypergraph_laplacian()].
#' @param algorithm Character. `"spectral"` (default; RDC-Spec) or
#'   `"symnmf"` (RDC-Sym).
#' @param nstart Integer. Random restarts (default 25).
#' @param seed Optional integer seed for initialization.
#' @param max_iter Maximum multiplicative-update iterations for
#'   `algorithm = "symnmf"`.
#' @param tol Relative objective tolerance for `algorithm = "symnmf"`.
#'
#' @return An object of class `net_hypergraph_cluster`: a list with
#'   `$clusters` (data.frame, one row per node: `node`, `cluster` - labels
#'   `"Cluster 1"`, `"Cluster 2"`, ... ordered by first appearance),
#'   `$embedding` (node x k row-normalized spectral embedding used by
#'   k-means, dims `dim1..dimk`), `$k`, `$type`, `$eigenvalues` (full
#'   Laplacian spectrum, increasing), `$eigengap` (gap after the k-th
#'   eigenvalue), `$sizes` (data.frame `cluster`/`size`), `$pi`
#'   (stationary distribution) and `$params`. Has `print`, `summary`,
#'   `plot` and `as.data.frame` methods; `as.data.frame()` returns one row
#'   per node with `node`, `cluster`, the stationary probability `pi`, and
#'   the embedding coordinates.
#'
#' @references
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. \emph{CIKM 2020}, 495-504.
#' \doi{10.1145/3340531.3412034}
#'
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. \emph{ICML 2019}.
#'
#' @examples
#' events <- data.frame(
#'   person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
#'              "d", "e", "f", "c", "d"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
#'               "m4", "m4", "m4", "m5", "m5")
#' )
#' hg <- group_hypergraph(events, actor = "person", group = "meeting")
#' cl <- hypergraph_cluster(hg, k = 2, seed = 1)
#' cl
#' as.data.frame(cl)
#'
#' @export
hypergraph_cluster <- function(hg, k,
                               type = c("zhou", "random_walk"),
                               edge_weights = NULL,
                               nstart = 25L,
                               seed = NULL,
                               algorithm = c("spectral", "symnmf"),
                               max_iter = 500L,
                               tol = 1e-6) {
  type <- match.arg(type)
  algorithm <- match.arg(algorithm)
  .hl_validate_hg(hg)
  stopifnot(
    "`k` must be a single whole number" =
      is.numeric(k) && length(k) == 1L && is.finite(k) && k == round(k),
    "`nstart` must be a single positive whole number" =
      is.numeric(nstart) && length(nstart) == 1L && nstart >= 1 &&
        nstart == round(nstart),
    "`max_iter` must be a single positive whole number" =
      is.numeric(max_iter) && length(max_iter) == 1L && max_iter >= 1 &&
        max_iter == round(max_iter),
    "`tol` must be a single positive finite number" =
      is.numeric(tol) && length(tol) == 1L && is.finite(tol) && tol > 0
  )
  k <- as.integer(k)
  n <- hg$n_nodes
  if (k < 2L || k > n - 1L) {
    stop(sprintf("`k` must be between 2 and n_nodes - 1 (= %d), got %d.",
                 n - 1L, k), call. = FALSE)
  }
  if (!is.null(seed)) set.seed(as.integer(seed))

  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  eig <- eigen((parts$L + t(parts$L)) / 2, symmetric = TRUE)
  # eigen() returns decreasing order; take the k SMALLEST, increasing
  ord <- rev(seq_len(n))
  values <- eig$values[ord]
  if (algorithm == "spectral") {
    U <- eig$vectors[, ord[seq_len(k)], drop = FALSE]
    # Row-normalize the spectral embedding to unit length (zero rows kept)
    row_norm <- sqrt(rowSums(U^2))
    nz <- row_norm > 0
    U[nz, ] <- U[nz, , drop = FALSE] / row_norm[nz]
    km <- stats::kmeans(U, centers = k, nstart = as.integer(nstart),
                        iter.max = 100L)
    assignment <- km$cluster
    diagnostics <- list(tot_withinss = km$tot.withinss)
  } else {
    T <- diag(n) - parts$L
    # T is non-negative analytically; remove only floating-point undershoot.
    if (min(T) < -1e-10) {
      stop("`I - L` has negative entries and cannot be factorized by SymNMF.",
           call. = FALSE)
    }
    T <- pmax((T + t(T)) / 2, 0)
    fit <- .hl_symnmf(T, k = k, nstart = as.integer(nstart), seed = seed,
                      max_iter = as.integer(max_iter), tol = tol)
    U <- fit$factor
    assignment <- max.col(U, ties.method = "first")
    diagnostics <- fit[c("objective", "objective_history", "iterations",
                         "converged", "restart")]
  }

  # Deterministic labels: "Cluster 1" = first node's cluster, etc.
  relabel <- match(assignment, unique(assignment))
  cluster_lab <- paste("Cluster", relabel)
  clusters <- data.frame(node = hg$nodes, cluster = cluster_lab,
                         stringsAsFactors = FALSE)
  dimnames(U) <- list(hg$nodes, paste0("dim", seq_len(k)))
  sizes <- as.data.frame(table(cluster = cluster_lab),
                         stringsAsFactors = FALSE)
  names(sizes) <- c("cluster", "size")
  sizes <- sizes[order(sizes$cluster), , drop = FALSE]
  rownames(sizes) <- NULL

  structure(
    list(
      clusters    = clusters,
      embedding   = U,
      k           = k,
      type        = type,
      algorithm   = algorithm,
      eigenvalues = values,
      eigengap    = if (k < n) values[k + 1L] - values[k] else NA_real_,
      sizes       = sizes,
      pi          = parts$pi,
      n_nodes     = n,
      n_hyperedges = hg$n_hyperedges,
      params = c(list(edge_weights = parts$w, nstart = as.integer(nstart),
                      seed = seed, max_iter = as.integer(max_iter), tol = tol),
                 diagnostics)
    ),
    class = "net_hypergraph_cluster"
  )
}

#' Transductive label spreading on a hypergraph
#'
#' Semi-supervised classification of hypergraph nodes by the regularization
#' framework of Zhou et al. (2006): given labels for a subset of nodes, the
#' scores `F = (1 - xi) * (I - xi * S)^{-1} Y` spread the labels over the
#' hypergraph, where `S = I - L` is the normalized similarity operator of
#' the chosen Laplacian and `Y` is the label indicator matrix. Each node is
#' assigned the class with the highest score. This is the non-neural
#' ancestor of hypergraph-attention text classifiers: with documents as
#' hyperedges over words (or vice versa) it classifies unlabeled nodes from
#' a handful of labeled ones.
#'
#' @param hg A connected `net_hypergraph`.
#' @param labels Node labels. Either a named vector (names = node names,
#'   values = class labels) covering a subset of nodes, or a full-length
#'   vector aligned with `hg$nodes` with `NA` for unlabeled nodes. At least
#'   two distinct classes must be labeled.
#' @param xi Numeric in `(0, 1)`. Spreading coefficient (default `0.99`);
#'   larger values weight the hypergraph structure more relative to the
#'   initial labels.
#' @param type,edge_weights Passed to [hypergraph_laplacian()].
#' @param normalization Decision rule applied to the score matrix before
#'   the argmax. `"none"` (default) is the raw Zhou (2006) rule.
#'   `"class_mass"` divides each class column by its total spread mass
#'   (class-mass normalization, Zhu et al. 2003) before the argmax; use it
#'   when the labeled seeds are class-imbalanced, where the raw rule can
#'   collapse every prediction onto the majority class.
#'
#' @return An object of class `net_hypergraph_transduction`: a list with
#'   `$predictions` (data.frame, one row per node: `node`, `label` (given,
#'   `NA` if unlabeled), `predicted`, `score` (winning class score),
#'   `margin` (winning minus runner-up score)), `$classes`, `$scores`
#'   (node x class score matrix), `$xi`, `$type`, `$normalization`,
#'   `$n_labeled` and `$params`. Has `print`, `summary`, `plot` and `as.data.frame` methods;
#'   `as.data.frame(x, what = "scores")` returns the tidy long score table.
#'
#' @references
#' Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
#' Clustering, classification, and embedding. \emph{NeurIPS 19}.
#'
#' Zhu, X., Ghahramani, Z., & Lafferty, J. (2003). Semi-supervised learning
#' using Gaussian fields and harmonic functions. \emph{ICML 20}.
#'
#' @examples
#' events <- data.frame(
#'   person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
#'              "d", "e", "f", "c", "d"),
#'   meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
#'               "m4", "m4", "m4", "m5", "m5")
#' )
#' hg <- group_hypergraph(events, actor = "person", group = "meeting")
#' tr <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"))
#' tr
#' as.data.frame(tr)
#'
#' @export
hypergraph_transduction <- function(hg, labels, xi = 0.99,
                                    type = c("zhou", "random_walk"),
                                    edge_weights = NULL,
                                    normalization = c("none", "class_mass")) {
  type <- match.arg(type)
  normalization <- match.arg(normalization)
  .hl_validate_hg(hg)
  stopifnot(
    "`xi` must be a single number in (0, 1)" =
      is.numeric(xi) && length(xi) == 1L && xi > 0 && xi < 1
  )
  n <- hg$n_nodes
  nodes <- hg$nodes

  # Normalize `labels` to a full-length character vector with NAs
  lab <- if (!is.null(names(labels))) {
    unknown <- setdiff(names(labels), nodes)
    if (length(unknown) > 0L) {
      stop("Unknown node names in `labels`: ",
           paste(unknown, collapse = ", "), call. = FALSE)
    }
    full <- rep(NA_character_, n)
    full[match(names(labels), nodes)] <- as.character(labels)
    full
  } else {
    if (length(labels) != n) {
      stop(sprintf(paste0("Unnamed `labels` must have length n_nodes ",
                          "(= %d), got %d."), n, length(labels)),
           call. = FALSE)
    }
    as.character(labels)
  }
  classes <- sort(unique(lab[!is.na(lab)]))
  if (length(classes) < 2L) {
    stop("`labels` must contain at least two distinct classes.",
         call. = FALSE)
  }

  parts <- .hl_build(hg, type = type, edge_weights = edge_weights)
  S <- diag(n) - (parts$L + t(parts$L)) / 2
  Y <- vapply(classes, function(cl) as.numeric(!is.na(lab) & lab == cl),
              numeric(n))
  F_scores <- (1 - xi) * solve(diag(n) - xi * S, Y)
  dimnames(F_scores) <- list(nodes, classes)

  predictions <- .hl_score_predictions(F_scores, lab, normalization)

  structure(
    list(
      predictions = predictions,
      classes     = classes,
      scores      = F_scores,
      xi          = xi,
      type        = type,
      normalization = normalization,
      n_labeled   = sum(!is.na(lab)),
      n_nodes     = n,
      params      = list(edge_weights = parts$w)
    ),
    class = "net_hypergraph_transduction"
  )
}

# ---- Internal machinery -------------------------------------------------

# Decision rule on a node x class score matrix. "none" is the Zhou (2006)
# argmax on the raw spread scores; "class_mass" first divides each class
# column by its total spread mass (class-mass normalization, Zhu et al.
# 2003), which prevents an imbalanced seed set from collapsing every
# prediction onto the majority class. `score` and `margin` are reported on
# the matrix the argmax actually used. Shared by hypergraph_transduction(),
# its sparse counterpart, and the neural classifiers (hg_neural, hg_hypergat).
.hl_score_predictions <- function(F_scores, lab, normalization) {
  decision <- if (identical(normalization, "class_mass")) {
    mass <- colSums(F_scores)
    if (any(mass <= 0)) {
      stop(errorCondition(
        "a class received zero total spread mass; class_mass normalization is undefined",
        class = "honets_bad_input", call = NULL
      ))
    }
    sweep(F_scores, 2L, mass, "/")
  } else {
    F_scores
  }
  n <- nrow(F_scores)
  classes <- colnames(F_scores)
  win <- max.col(decision, ties.method = "first")
  score <- decision[cbind(seq_len(n), win)]
  runner <- vapply(seq_len(n), function(i) {
    max(decision[i, -win[i]])
  }, numeric(1L))
  data.frame(
    node      = rownames(F_scores),
    label     = lab,
    predicted = classes[win],
    score     = score,
    margin    = score - runner,
    stringsAsFactors = FALSE
  )
}

#' Validate a net_hypergraph for spectral work (connected, non-degenerate)
#' @noRd
.hl_validate_hg <- function(hg) {
  stopifnot(
    "`hg` must be a net_hypergraph (build_hypergraph/group_hypergraph)" =
      inherits(hg, "net_hypergraph"),
    "`hg` must have at least 2 nodes" = hg$n_nodes >= 2L,
    "`hg` must have at least 1 hyperedge" = hg$n_hyperedges >= 1L
  )
  if (!.hl_is_connected(hg$incidence > 0)) {
    stop(errorCondition(
      paste0("The hypergraph is not connected; the random walk has no ",
             "unique stationary distribution. Analyze components ",
             "separately."),
      class = "honets_hypergraph_disconnected", call = NULL
    ))
  }
  invisible(TRUE)
}

#' Connectivity of the hypergraph via BFS on the node co-membership graph
#' @noRd
.hl_is_connected <- function(pattern) {
  n <- nrow(pattern)
  if (n == 1L) return(TRUE)
  adj <- tcrossprod(pattern * 1.0) > 0
  reached <- c(TRUE, rep(FALSE, n - 1L))
  # BFS frontier expansion; bounded by n iterations (diameter <= n - 1)
  repeat {
    nxt <- reached | (as.vector(adj %*% reached) > 0)
    if (identical(nxt, reached)) break
    reached <- nxt
  }
  all(reached)
}

#' Population standard deviation (ddof = 0, numpy convention)
#' @noRd
.hl_pop_sd <- function(x) sqrt(mean((x - mean(x))^2))

#' Build Laplacian + stationary distribution + edge weights for a type
#' @return list(L, pi, w)
#' @noRd
.hl_build <- function(hg, type, edge_weights) {
  gamma <- hg$incidence * 1.0
  pattern <- (gamma > 0) * 1.0
  n <- nrow(gamma)
  m <- ncol(gamma)
  nodes <- rownames(gamma) %||% hg$nodes

  edge_weights <- .hl_check_edge_weights(edge_weights, m)

  if (type == "zhou") {
    # Binary incidence pattern; default hyperedge weights are the window
    # counts when present (window_hypergraph), else unit
    w <- as.numeric(edge_weights %||% hg$window_counts %||% rep(1, m))
    delta_e <- colSums(pattern)
    d_v <- as.vector(pattern %*% w)
    Hs <- pattern / sqrt(d_v)
    Theta <- Hs %*% ((w / delta_e) * t(Hs))
    L <- diag(n) - (Theta + t(Theta)) / 2
    pi_v <- d_v / sum(d_v)
  } else {
    # Hayashi EDVW random walk on the weighted incidence; the transition
    # matrix is shared with hypergraph_centrality(type = "pagerank")
    rw <- .hl_rw_transition(hg, edge_weights)
    w <- rw$w
    P <- rw$P

    # Stationary distribution: dominant left eigenvector of P
    eg <- eigen(t(P))
    lead <- which.max(Mod(eg$values))
    pi_v <- Re(eg$vectors[, lead])
    pi_v <- pi_v / sum(pi_v)
    if (min(pi_v) < -1e-8) {
      stop("Stationary distribution has negative mass; the walk appears ",
           "reducible despite the connectivity check.", call. = FALSE)
    }
    pi_v <- pmax(pi_v, 0)
    pi_v <- pi_v / sum(pi_v)

    # Chung (2005) directed Laplacian, symmetrized via pi
    G <- (P * sqrt(pi_v))                       # scale rows by sqrt(pi)
    G <- sweep(G, 2L, sqrt(pi_v), "/")          # scale cols by 1/sqrt(pi)
    L <- diag(n) - (G + t(G)) / 2
  }

  dimnames(L) <- list(nodes, nodes)
  names(pi_v) <- nodes
  list(L = L, pi = pi_v, w = w)
}

#' EDVW random-walk transition matrix of a net_hypergraph
#'
#' `P[v, u] = sum_e w(e) 1[v in e] / d(v) * gamma_e(u) / delta(e)`
#' (Chitra & Raphael 2019; Hayashi et al. 2020). Shared by .hl_build's
#' random_walk branch and hypergraph_centrality(type = "pagerank").
#' Default hyperedge weights: window counts when present
#' (window_hypergraph), else the Hayashi dispersion heuristic. Rows of
#' nodes that sit in no hyperedge (d_v = 0) come back NaN; callers that
#' tolerate such nodes (PageRank teleportation) must replace them.
#'
#' @param hg A `net_hypergraph`.
#' @param edge_weights NULL or positive numeric vector, one per hyperedge.
#' @return list(P, w, d_v, nodes).
#' @noRd
.hl_rw_transition <- function(hg, edge_weights = NULL) {
  gamma <- hg$incidence * 1.0
  pattern <- (gamma > 0) * 1.0
  m <- ncol(gamma)
  edge_weights <- .hl_check_edge_weights(edge_weights, m)
  w <- as.numeric(edge_weights %||% hg$window_counts %||%
                    vapply(seq_len(m), function(j) {
                      .hl_pop_sd(gamma[gamma[, j] > 0, j]) + 1
                    }, numeric(1L)))
  delta_e <- colSums(gamma)
  d_v <- as.vector(pattern %*% w)
  A <- sweep(pattern, 2L, w, "*") / d_v      # A[v,e] = w(e) 1[v in e]/d(v)
  B <- t(gamma) / delta_e                     # B[e,u] = gamma_e(u)/delta(e)
  list(P = A %*% B, w = w, d_v = d_v,
       nodes = rownames(gamma) %||% hg$nodes)
}

#' Validate edge_weights; a positive scalar recycles to all hyperedges
#' @noRd
.hl_check_edge_weights <- function(edge_weights, m) {
  if (is.null(edge_weights)) return(NULL)
  stopifnot(
    "`edge_weights` must be positive numeric, length 1 or one per hyperedge" =
      is.numeric(edge_weights) &&
      (length(edge_weights) == 1L || length(edge_weights) == m) &&
      all(is.finite(edge_weights)) && all(edge_weights > 0)
  )
  if (length(edge_weights) == 1L) {
    edge_weights <- rep(as.numeric(edge_weights), m)
  }
  edge_weights
}

# ---- S3: net_hypergraph_cluster ----------------------------------------

#' Print method for net_hypergraph_cluster
#'
#' @param x A `net_hypergraph_cluster` object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.net_hypergraph_cluster <- function(x, ...) {
  algorithm <- x$algorithm %||% "spectral"
  label <- switch(algorithm, spectral = "RDC-Spec spectral",
                  symnmf = "RDC-Sym symmetric-NMF", algorithm)
  cat("Hypergraph ", label, " clustering (", x$type, " Laplacian)\n",
      sep = "")
  cat(sprintf("  Nodes: %d | Hyperedges: %d | k: %d\n",
              x$n_nodes, x$n_hyperedges, x$k))
  cat(sprintf("  Cluster sizes: %s\n",
              paste(sprintf("%s = %d", x$sizes$cluster, x$sizes$size),
                    collapse = ", ")))
  if (!identical(algorithm, "spectral")) {
    cat(sprintf("  Objective: %.6g | Iterations: %d | Converged: %s\n",
                x$params$objective, x$params$iterations,
                if (isTRUE(x$params$converged)) "yes" else "no"))
  } else {
    cat(sprintf("  Eigengap after k: %.4f\n", x$eigengap))
  }
  invisible(x)
}

#' Summary method for net_hypergraph_cluster
#'
#' @param object A `net_hypergraph_cluster` object.
#' @param ... Additional arguments (ignored).
#' @return A data.frame, one row per cluster: `cluster`, `size`, `share`.
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @export
summary.net_hypergraph_cluster <- function(object, ...) {
  out <- object$sizes
  out$share <- out$size / sum(out$size)
  invisible(out)
}

#' Coerce a net_hypergraph_cluster to a data.frame
#'
#' @param x A `net_hypergraph_cluster` object.
#' @param ... Additional arguments (ignored).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return The tidy assignment table: one row per node, columns `node`,
#'   `cluster`, `pi` (stationary probability of the node under the
#'   Laplacian's random walk) and the spectral-embedding or NMF-factor coordinates
#'   `dim1..dimk`.
#' @export
as.data.frame.net_hypergraph_cluster <- function(x, ..., top = NULL) {
  out <- x$clusters
  out$pi <- as.numeric(x$pi)
  .ho_top(cbind(out, as.data.frame(x$embedding), row.names = NULL), top)
}

#' Plot method for net_hypergraph_cluster
#'
#' Two diagnostic panels. `"spectrum"`: scree plot of the Laplacian
#' spectrum with the k used for clustering marked - the eigengap after k
#' supports (or questions) the choice of k. `"embedding"`: the nodes in
#' the first two spectral-embedding dimensions, labelled, coloured and
#' shaped by cluster, sized by stationary probability - the geometry
#' k-means actually clustered. `"both"` (default) arranges the two side
#' by side (via gridExtra when available, base grid viewports otherwise).
#'
#' @param x A `net_hypergraph_cluster` object.
#' @param what Character. `"both"` (default), `"spectrum"`, or
#'   `"embedding"`.
#' @param n_values Integer. How many smallest eigenvalues to show in the
#'   spectrum panel (default: `min(3 * k, n_nodes)`).
#' @param ... Additional arguments (ignored).
#' @return For `"spectrum"`/`"embedding"`, the ggplot object. For
#'   `"both"`, the arranged gtable when gridExtra is installed (drawn on
#'   the current device), otherwise the two panels are drawn via grid
#'   viewports and the list of the two ggplots is returned invisibly.
#' @export
plot.net_hypergraph_cluster <- function(x,
                                        what = c("both", "spectrum",
                                                 "embedding"),
                                        n_values = NULL, ...) {
  what <- match.arg(what)
  has_spectrum <- any(is.finite(x$eigenvalues))
  if (!has_spectrum && what == "spectrum") {
    stop("This clustering objective has no Laplacian spectrum; use `what = \"embedding\"`.",
         call. = FALSE)
  }
  if (!has_spectrum && what == "both") what <- "embedding"
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#999999", "#000000")

  n_show <- as.integer(n_values %||% min(3L * x$k, x$n_nodes))
  df_s <- data.frame(
    index = seq_len(n_show),
    value = x$eigenvalues[seq_len(n_show)],
    used  = ifelse(seq_len(n_show) <= x$k, "used", "unused")
  )
  p_spec <- ggplot2::ggplot(df_s, ggplot2::aes(x = .data$index,
                                               y = .data$value,
                                               color = .data$used,
                                               shape = .data$used)) +
    ggplot2::geom_line(color = "#999999", linewidth = 0.4) +
    ggplot2::geom_point(size = 2.5) +
    ggplot2::geom_vline(xintercept = x$k + 0.5, linetype = "dashed",
                        color = "#999999") +
    ggplot2::scale_color_manual(
      values = c(used = "#0072B2", unused = "#E69F00")
    ) +
    ggplot2::scale_shape_manual(values = c(used = 16, unused = 17)) +
    ggplot2::labs(
      title = sprintf("Laplacian spectrum (k = %d)", x$k),
      subtitle = sprintf("Eigengap after k: %.4f", x$eigengap),
      x = "Eigenvalue index (increasing)", y = "Eigenvalue",
      color = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (what == "spectrum") return(p_spec)

  df_e <- as.data.frame(x)
  p_emb <- ggplot2::ggplot(df_e, ggplot2::aes(x = .data$dim1,
                                              y = .data$dim2,
                                              color = .data$cluster,
                                              shape = .data$cluster)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$pi), alpha = 0.85) +
    ggplot2::geom_text(ggplot2::aes(label = .data$node),
                       vjust = -1.1, size = 3.5, show.legend = FALSE) +
    ggplot2::scale_color_manual(values = okabe) +
    ggplot2::scale_shape_manual(
      values = rep(c(16, 17, 15, 18, 8, 7, 3, 4, 6), length.out = x$k)
    ) +
    ggplot2::scale_size_continuous(range = c(2, 6),
                                   guide = "none") +
    ggplot2::expand_limits(
      x = range(df_e$dim1) + c(-0.15, 0.15) * diff(range(df_e$dim1)),
      y = range(df_e$dim2) + c(-0.1, 0.25) * diff(range(df_e$dim2))
    ) +
    ggplot2::labs(
      title = switch(x$algorithm %||% "spectral",
                     spectral = "Spectral embedding (first two dimensions)",
                     symnmf = "SymNMF factors (first two dimensions)",
                     joint = "Joint-NMF vertex factors (first two dimensions)",
                     joint_symmetric = paste0(
                       "Joint-Symmetric-NMF vertex factors ",
                       "(first two dimensions)")),
      subtitle = "Point size = stationary probability of the node",
      x = "dim1", y = "dim2", color = NULL, shape = NULL
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
  if (what == "embedding") return(p_emb)

  if (requireNamespace("gridExtra", quietly = TRUE)) {
    return(gridExtra::grid.arrange(p_spec, p_emb, ncol = 2))
  }
  grid::grid.newpage()
  grid::pushViewport(grid::viewport(layout = grid::grid.layout(1L, 2L)))
  print(p_spec, vp = grid::viewport(layout.pos.row = 1L,
                                    layout.pos.col = 1L))
  print(p_emb, vp = grid::viewport(layout.pos.row = 1L,
                                   layout.pos.col = 2L))
  grid::popViewport()
  invisible(list(spectrum = p_spec, embedding = p_emb))
}

# ---- S3: net_hypergraph_transduction -----------------------------------

#' Print method for net_hypergraph_transduction
#'
#' @param x A `net_hypergraph_transduction` object.
#' @param ... Additional arguments (ignored).
#' @return The input object, invisibly.
#' @export
print.net_hypergraph_transduction <- function(x, ...) {
  cat("Hypergraph transductive label spreading (", x$type,
      " Laplacian, xi = ", format(x$xi), ")\n", sep = "")
  tab <- table(x$predictions$predicted)
  cat(sprintf("  Nodes: %d (%d labeled) | Classes: %s\n",
              x$n_nodes, x$n_labeled, paste(x$classes, collapse = ", ")))
  cat(sprintf("  Predicted: %s\n",
              paste(sprintf("%s = %d", names(tab), as.integer(tab)),
                    collapse = ", ")))
  invisible(x)
}

#' Summary method for net_hypergraph_transduction
#'
#' @param object A `net_hypergraph_transduction` object.
#' @param ... Additional arguments (ignored).
#' @return A data.frame, one row per class: `class`, `n_labeled`,
#'   `n_predicted`, `mean_margin` (mean winning margin among the nodes
#'   predicted into the class).
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @export
summary.net_hypergraph_transduction <- function(object, ...) {
  p <- object$predictions
  out <- do.call(rbind, lapply(object$classes, function(cl) {
    sel <- p$predicted == cl
    data.frame(
      class       = cl,
      n_labeled   = sum(!is.na(p$label) & p$label == cl),
      n_predicted = sum(sel),
      mean_margin = if (any(sel)) mean(p$margin[sel]) else NA_real_,
      stringsAsFactors = FALSE
    )
  }))
  rownames(out) <- NULL
  invisible(out)
}

#' Coerce a net_hypergraph_transduction to a data.frame
#'
#' @param x A `net_hypergraph_transduction` object.
#' @param row.names Ignored (present for S3 consistency with the generic).
#' @param optional Ignored (present for S3 consistency with the generic).
#' @param ... Additional arguments (ignored).
#' @param what Character. `"predictions"` (default) for the one-row-per-node
#'   table, `"scores"` for the tidy long score table (one row per node x
#'   class: `node`, `class`, `score`).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @return A data.frame as selected by `what`.
#' @export
as.data.frame.net_hypergraph_transduction <- function(
    x, row.names = NULL, optional = FALSE, ...,
    what = c("predictions", "scores"), top = NULL) {
  what <- match.arg(what)
  if (what == "predictions") return(.ho_top(x$predictions, top))
  .ho_top(data.frame(
    node  = rep(rownames(x$scores), times = ncol(x$scores)),
    class = rep(colnames(x$scores), each = nrow(x$scores)),
    score = as.vector(x$scores),
    stringsAsFactors = FALSE
  ), top)
}

#' Plot method for net_hypergraph_transduction
#'
#' Heatmap of the full node-by-class score matrix: rows are nodes (grouped
#' by predicted class), columns are classes, tile shading and printed
#' values are the spreading scores. Seed nodes (given labels) carry a
#' black tile border, and each node's winning class is marked with a dot,
#' so agreement between seeds, scores, and decisions is visible in one
#' panel. Rows whose winning and runner-up scores are close (small
#' `margin`) are the assignments to distrust.
#'
#' @param x A `net_hypergraph_transduction` object.
#' @param ... Additional arguments (ignored).
#' @return A ggplot object, invisibly printable.
#' @export
plot.net_hypergraph_transduction <- function(x, ...) {
  sc <- as.data.frame(x, what = "scores")
  p <- x$predictions
  sc$seeded <- !is.na(p$label[match(sc$node, p$node)]) &
    p$label[match(sc$node, p$node)] == sc$class
  sc$winner <- p$predicted[match(sc$node, p$node)] == sc$class
  # order rows by predicted class, then margin (most confident on top)
  ord <- p$node[order(p$predicted, -p$margin)]
  sc$node <- factor(sc$node, levels = rev(ord))

  ggplot2::ggplot(sc, ggplot2::aes(x = .data$class, y = .data$node)) +
    ggplot2::geom_tile(ggplot2::aes(fill = .data$score),
                       color = "white", linewidth = 0.6) +
    ggplot2::geom_tile(
      data = sc[sc$seeded, , drop = FALSE],
      fill = NA, color = "#000000", linewidth = 1.1
    ) +
    ggplot2::geom_point(
      data = sc[sc$winner, , drop = FALSE],
      ggplot2::aes(x = .data$class, y = .data$node),
      shape = 21, size = 2.6, stroke = 0.9,
      fill = "white", color = "#D55E00",
      inherit.aes = FALSE
    ) +
    ggplot2::geom_text(ggplot2::aes(label = sprintf("%.3f", .data$score)),
                       size = 3, vjust = -0.9) +
    ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#4A6FE3") +
    ggplot2::labs(
      title = "Transductive label spreading: score matrix",
      subtitle = sprintf(paste0("%d labeled of %d nodes, xi = %s. Black ",
                                "border = seed; dot = winning class."),
                         x$n_labeled, x$n_nodes, format(x$xi)),
      x = "Class", y = NULL, fill = "Score"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom",
                   panel.grid = ggplot2::element_blank())
}
