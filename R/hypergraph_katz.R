# Katz centrality of a hypergraph (ROADMAP Phase 2, item 5): the Katz (1953)
# status index on the hypergraph adjacency of Estrada & Rodriguez-Velazquez
# (2006), A_ij = number of hyperedges containing both i and j (i != j) --
# the weighted clique expansion, which is also XGI's clique-motif matrix.
# Battiston et al. (2020, Sec. III.B.3) state the same series on a
# higher-order adjacency.
#
#   katz = sum_{k >= 1} alpha^k A^k 1 = (I - alpha A)^{-1} alpha A 1,
#   0 < alpha < 1 / lambda_max(A).
#
# The engine sits behind hg_centrality(type = "katz", alpha = ...); the patch
# to R/text_verbs.R is recorded in the Phase 2 report. Oracles
# (local_testing_and_equivalence/test-oracle-katz.R): igraph::alpha_centrality
# on the weighted clique expansion (Bonacich & Lloyd's alpha centrality is
# katz + 1) and XGI 0.10.2 katz_centrality (fixed alpha = 2^-n, L1 scale).

# Weighted clique-expansion adjacency with zero diagonal, dense or sparse.
.hg_katz_adjacency <- function(hg) {
  membership <- .hg_membership(hg)
  adjacency <- Matrix::tcrossprod(membership)
  Matrix::diag(adjacency) <- 0
  adjacency
}

# Largest eigenvalue of a symmetric non-negative adjacency (its spectral
# radius, by Perron-Frobenius).
.hg_katz_lambda <- function(adjacency) {
  if (methods::is(adjacency, "sparseMatrix") && nrow(adjacency) > 10L) {
    general <- methods::as(adjacency, "generalMatrix")
    return(RSpectra::eigs_sym(general, k = 1L, which = "LA")$values[1L])
  }
  eigen(as.matrix(adjacency), symmetric = TRUE,
        only.values = TRUE)$values[1L]
}

# Katz centrality: one value per node, on its natural scale (the
# alpha-discounted count of walks of length >= 1 leaving the node).
.hg_katz_fit <- function(hg, alpha) {
  .thg_check_hg(hg)
  if (!(is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
        alpha > 0)) {
    stop(errorCondition(
      "`alpha` must be a single positive number below 1 / lambda_max",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  adjacency <- .hg_katz_adjacency(hg)
  lambda <- .hg_katz_lambda(adjacency)
  # the eigenvalue carries rounding error: treat the boundary as outside
  if (alpha * lambda >= 1 - 1e-10) {
    stop(errorCondition(
      sprintf(paste0("`alpha` = %g does not make the Katz series converge: ",
                     "it must be below 1 / lambda_max = %g"),
              alpha, 1 / lambda),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  n <- nrow(adjacency)
  rhs <- alpha * as.numeric(adjacency %*% rep(1, n))
  # I - alpha A is symmetric positive definite for alpha < 1 / lambda_max,
  # so a Cholesky solve is exact and stable.
  if (methods::is(adjacency, "sparseMatrix")) {
    system <- Matrix::Diagonal(n) - alpha * adjacency
    factor <- Matrix::Cholesky(methods::as(system, "symmetricMatrix"))
    katz <- as.numeric(Matrix::solve(factor, rhs))
  } else {
    system <- diag(n) - alpha * as.matrix(adjacency)
    factor <- chol(system)
    katz <- backsolve(factor, forwardsolve(t(factor), rhs))
  }
  data.frame(node = rownames(hg$incidence), katz = as.numeric(katz),
             row.names = NULL, stringsAsFactors = FALSE)
}
