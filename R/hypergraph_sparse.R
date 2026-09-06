# Sparse core (v0.4 scale mandate). Large corpora as Matrix::dgCMatrix
# incidences, with every heavy computation expressed as sparse mat-vec
# operators: PageRank by power iteration, transduction by conjugate
# gradient on the similarity operator, spectral clustering by RSpectra
# partial eigenpairs. The shipped dense engines
# (HyperNetX-parity-verified) are the equivalence oracle: sparse and dense
# must agree on every small corpus (tested to 1e-8 or better).

.thg_is_sparse <- function(hg) methods::is(hg$incidence, "sparseMatrix")

# Sparse bipartite constructor mirroring group_hypergraph():
# sorted vertex/edge names, summed weights, same top-level fields.
.thg_sparse_bipartite <- function(long, member, group, weight) {
  vertices <- sort(unique(long[[member]]))
  edges <- sort(unique(long[[group]]))
  incidence <- Matrix::sparseMatrix(
    i = match(long[[member]], vertices),
    j = match(long[[group]], edges),
    x = as.numeric(long[[weight]]),
    dims = c(length(vertices), length(edges)),
    dimnames = list(vertices, edges)
  )
  sizes <- Matrix::colSums(incidence > 0)
  size_tab <- table(sizes)
  size_distribution <- as.integer(size_tab)
  names(size_distribution) <- paste0("size_", names(size_tab))
  structure(
    list(
      incidence = incidence,
      nodes = vertices,
      n_nodes = length(vertices),
      n_hyperedges = length(edges),
      size_distribution = size_distribution,
      params = list(sparse = TRUE)
    ),
    class = "net_hypergraph"
  )
}

# Hayashi SD+1 default edge weights, computed on the dgCMatrix slots
# without densifying: population SD of each column's non-zeros, plus one.
.thg_sparse_edge_weights <- function(incidence) {
  csp <- methods::as(incidence, "CsparseMatrix")
  nnz <- diff(csp@p)
  s1 <- Matrix::colSums(csp)
  s2 <- Matrix::colSums(csp^2)
  mu <- s1 / nnz
  unname(sqrt(pmax(s2 / nnz - mu^2, 0)) + 1)
}

# Connectivity by alternating vertex -> edge -> vertex BFS, never forming
# the vertex-vertex co-membership matrix (which densifies at scale).
.thg_sparse_connected <- function(membership) {
  n <- nrow(membership)
  if (n == 1L) return(TRUE)
  reached <- c(TRUE, rep(FALSE, n - 1L))
  # frontier expansion, bounded by n iterations; sequential by nature
  repeat {
    edge_hit <- as.numeric(Matrix::crossprod(membership, reached)) > 0
    nxt <- reached | (as.numeric(membership %*% edge_hit) > 0)
    if (identical(nxt, reached)) break
    reached <- nxt
  }
  all(reached)
}

# EDVW walk operators over a (sparse or dense) incidence: left action
# v |-> vP, right actions P v and t(P) v, all as mat-vec chains.
.thg_walk_operators <- function(hg, edge_weights = NULL) {
  incidence <- hg$incidence
  membership <- (incidence > 0) * 1
  w <- edge_weights %||% .thg_sparse_edge_weights(incidence)
  stopifnot(
    "`edge_weights` must be positive and one per hyperedge" =
      is.numeric(w) && length(w) == ncol(incidence) &&
      all(is.finite(w)) && all(w > 0)
  )
  delta <- as.numeric(Matrix::colSums(incidence))
  d_v <- as.numeric(membership %*% w)
  if (any(d_v <= 0)) {
    stop(errorCondition(
      "every vertex must belong to at least one hyperedge",
      class = "honets_bad_input", call = NULL
    ))
  }
  scale_e <- w / delta
  list(
    left = function(v) {
      as.numeric(Matrix::crossprod(membership, v / d_v) * scale_e) |>
        (\(z) as.numeric(incidence %*% z))()
    },
    right = function(v) {
      as.numeric(membership %*% (scale_e *
        as.numeric(Matrix::crossprod(incidence, v)))) / d_v
    },
    d_v = d_v, delta = delta, w = w,
    membership = membership, incidence = incidence
  )
}

# Stationary distribution of the EDVW walk by power iteration on the left
# action (sequential fixed point; bounded, convergence surfaced).
.thg_sparse_stationary <- function(ops, tol = 1e-13, max_iter = 20000L) {
  n <- length(ops$d_v)
  pi_vec <- rep(1 / n, n)
  iter <- 0L
  # power iteration: sequential dependence between iterates
  while (iter < max_iter) {
    iter <- iter + 1L
    nxt <- ops$left(pi_vec)
    nxt <- nxt / sum(nxt)
    if (sum(abs(nxt - pi_vec)) < tol) {
      return(nxt)
    }
    pi_vec <- nxt
  }
  warning(warningCondition(
    sprintf("stationary distribution not converged in %d iterations",
            max_iter),
    class = "honets_no_converge"
  ))
  pi_vec
}

# Similarity operator S for the chosen Laplacian type (S = I - L), as a
# symmetric mat-vec, plus the stationary distribution it implies.
.thg_similarity_operator <- function(hg, type, edge_weights = NULL) {
  incidence <- hg$incidence
  membership <- (incidence > 0) * 1
  if (!.thg_sparse_connected(membership)) {
    stop(errorCondition(
      paste0("The hypergraph is not connected; the random walk has no ",
             "unique stationary distribution. Analyze components ",
             "separately."),
      class = "honets_hypergraph_disconnected", call = NULL
    ))
  }
  if (identical(type, "zhou")) {
    w <- edge_weights %||% rep(1, ncol(incidence))
    delta <- as.numeric(Matrix::colSums(membership))
    d_v <- as.numeric(membership %*% w)
    root_d <- sqrt(d_v)
    scale_e <- w / delta
    pi_v <- d_v / sum(d_v)
    smult <- function(v) {
      z <- as.numeric(Matrix::crossprod(membership, v / root_d)) * scale_e
      as.numeric(membership %*% z) / root_d
    }
  } else {
    ops <- .thg_walk_operators(hg, edge_weights)
    w <- ops$w
    pi_v <- .thg_sparse_stationary(ops)
    root_pi <- sqrt(pi_v)
    smult <- function(v) {
      g <- root_pi * ops$right(v / root_pi)
      gt <- ops$left(root_pi * v) / root_pi
      (g + gt) / 2
    }
  }
  names(pi_v) <- rownames(incidence)
  list(smult = smult, pi = pi_v, w = w, n = nrow(incidence))
}

# Conjugate gradient for (I - xi * S) x = b; SPD because the spectrum of S
# lies in [-1, 1] and xi < 1.
.thg_cg <- function(smult, b, xi, tol = 1e-12, max_iter = 10000L) {
  x <- numeric(length(b))
  r <- b
  p <- r
  rs <- sum(r^2)
  iter <- 0L
  # CG: each step depends on the previous residual; sequential by nature
  while (iter < max_iter && sqrt(rs) > tol) {
    iter <- iter + 1L
    ap <- p - xi * smult(p)
    alpha <- rs / sum(p * ap)
    x <- x + alpha * p
    r <- r - alpha * ap
    rs_new <- sum(r^2)
    p <- r + (rs_new / rs) * p
    rs <- rs_new
  }
  if (sqrt(rs) > tol) {
    warning(warningCondition(
      sprintf("conjugate gradient not converged in %d iterations", max_iter),
      class = "honets_no_converge"
    ))
  }
  x
}

# Sparse spectral clustering: top-k eigenpairs of S (= bottom-k of L)
# through RSpectra's function interface, then the same row-normalized
# k-means as the dense engine.
.thg_sparse_cluster <- function(hg, k, type, edge_weights, nstart, seed) {
  sim <- .thg_similarity_operator(hg, type, edge_weights)
  n <- sim$n
  k_ask <- min(k + 1L, n - 1L)
  eig <- RSpectra::eigs_sym(
    \(x, args) sim$smult(x), k = k_ask, which = "LA", n = n
  )
  if (length(eig$values) < k) {
    stop(errorCondition(
      "partial eigendecomposition did not converge; try the dense engine",
      class = "honets_no_converge", call = NULL
    ))
  }
  lap_values <- 1 - eig$values                 # increasing: smallest of L
  U <- eig$vectors[, seq_len(k), drop = FALSE]
  row_norm <- sqrt(rowSums(U^2))
  nz <- row_norm > 0
  U[nz, ] <- U[nz, , drop = FALSE] / row_norm[nz]

  if (!is.null(seed)) set.seed(as.integer(seed))
  km <- stats::kmeans(U, centers = k, nstart = as.integer(nstart),
                      iter.max = 100L)
  relabel <- match(km$cluster, unique(km$cluster))
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
      clusters = clusters,
      embedding = U,
      k = as.integer(k),
      type = type,
      algorithm = "spectral",
      eigenvalues = lap_values,
      eigengap = if (length(lap_values) > k) {
        lap_values[k + 1L] - lap_values[k]
      } else {
        NA_real_
      },
      sizes = sizes,
      pi = sim$pi,
      n_nodes = n,
      n_hyperedges = hg$n_hyperedges,
      params = list(edge_weights = sim$w, nstart = as.integer(nstart),
                    seed = seed, tot_withinss = km$tot.withinss,
                    sparse = TRUE)
    ),
    class = "net_hypergraph_cluster"
  )
}

# Sparse transduction: the dense engine's closed form, solved per class by
# conjugate gradient on the similarity operator.
.thg_sparse_transduction <- function(hg, labels, xi, type, edge_weights,
                                     normalization = "none") {
  n <- hg$n_nodes
  nodes <- hg$nodes
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
      stop(sprintf("Unnamed `labels` must have length n_nodes (= %d).", n),
           call. = FALSE)
    }
    as.character(labels)
  }
  classes <- sort(unique(lab[!is.na(lab)]))
  if (length(classes) < 2L) {
    stop("`labels` must contain at least two distinct classes.",
         call. = FALSE)
  }
  stopifnot(
    "`xi` must be a single number in (0, 1)" =
      is.numeric(xi) && length(xi) == 1L && xi > 0 && xi < 1
  )
  sim <- .thg_similarity_operator(hg, type, edge_weights)
  F_scores <- vapply(classes, \(cl) {
    y <- as.numeric(!is.na(lab) & lab == cl)
    (1 - xi) * .thg_cg(sim$smult, y, xi = xi)
  }, numeric(n))
  dimnames(F_scores) <- list(nodes, classes)

  predictions <- .hl_score_predictions(F_scores, lab, normalization)

  structure(
    list(predictions = predictions, classes = classes, scores = F_scores,
         xi = xi, type = type, normalization = normalization,
         n_labeled = sum(!is.na(lab)), n_nodes = n,
         params = list(edge_weights = sim$w, sparse = TRUE)),
    class = "net_hypergraph_transduction"
  )
}

# Sparse structural measures (nodes/edges/summary); the pairwise tables are
# guarded, since they are quadratic in the hyperedge count.
.thg_sparse_measures <- function(hg, what) {
  incidence <- hg$incidence
  membership <- (incidence > 0) * 1
  sizes <- as.numeric(Matrix::colSums(membership))
  if (identical(what, "nodes")) {
    triplet <- methods::as(membership, "TsparseMatrix")
    max_size <- tapply(sizes[triplet@j + 1L], triplet@i + 1L, max)
    co <- Matrix::tcrossprod(membership)
    Matrix::diag(co) <- 0
    # engine definition: strength(v) = sum of the SIZES of v's edges
    return(data.frame(
      node = hg$nodes,
      hyperdegree = as.integer(Matrix::rowSums(membership)),
      strength = as.numeric(membership %*% sizes),
      max_edge_size = as.integer(max_size[as.character(seq_len(hg$n_nodes))]),
      n_neighbors = as.integer(Matrix::rowSums(co > 0)),
      row.names = NULL
    ))
  }
  if (identical(what, "edges")) {
    return(data.frame(edge = colnames(incidence), size = as.integer(sizes),
                      row.names = NULL))
  }
  if (identical(what, "overlap")) {
    if (hg$n_hyperedges > 2000L) {
      stop(errorCondition(
        sprintf(
          "the overlap table is quadratic in hyperedges (%d here); compute it on a subset",
          hg$n_hyperedges
        ),
        class = "honets_sparse_too_large", call = NULL
      ))
    }
    co <- as.matrix(Matrix::crossprod(membership))
    pair <- which(upper.tri(co), arr.ind = TRUE)
    s1 <- sizes[pair[, "row"]]
    s2 <- sizes[pair[, "col"]]
    ov <- co[pair]
    return(data.frame(
      edge_1 = colnames(incidence)[pair[, "row"]],
      edge_2 = colnames(incidence)[pair[, "col"]],
      overlap = as.numeric(ov),
      overlap_coefficient = as.numeric(ov / pmin(s1, s2)),
      jaccard = as.numeric(ov / (s1 + s2 - ov)),
      row.names = NULL
    ))
  }
  # summary; pairwise participation only while the vertex count keeps the
  # co-membership crossprod affordable
  participation <- if (hg$n_nodes <= 5000L) {
    co_nodes <- Matrix::tcrossprod(membership) > 0
    (Matrix::nnzero(co_nodes) - hg$n_nodes) /
      (hg$n_nodes * (hg$n_nodes - 1))
  } else {
    NA_real_
  }
  data.frame(
    measure = c("n_nodes", "n_hyperedges", "density", "avg_edge_size",
                "pairwise_participation"),
    value = c(hg$n_nodes, hg$n_hyperedges,
              Matrix::nnzero(membership) / (hg$n_nodes * hg$n_hyperedges),
              mean(sizes), participation),
    row.names = NULL
  )
}
