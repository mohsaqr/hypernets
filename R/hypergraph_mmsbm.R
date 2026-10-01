# Hy-MMSBM: mixed-membership stochastic block model for hypergraphs
# (Ruggeri, Contisciani, Battiston & De Bacco 2023). Ported from the authors'
# reference implementation, github.com/nickruggeri/Hy-MMSBM (MIT),
# src/model/model.py: the EM updates `_w_update()` / `_u_update()`, the
# Poisson parameters `poisson_params()`, the log-likelihood and the
# initialisation `_init_w()` / `_init_u()`. Local oracle:
# local_testing_and_equivalence/test-oracle-mmsbm.R.
#
# Notation (paper, Eqs. 1-9): B is the binary incidence (nodes x hyperedges),
# A_e the hyperedge weight (multiplicity), u the N x K memberships, w the
# symmetric K x K affinity, s_e = sum_{i in e} u_i and
#   lambda_e = sum_{i<j in e} u_i' w u_j = (s_e' w s_e - sum_{i in e} u_i' w u_i) / 2
# (Appendix C, Eq. C2), which is linear in |e|: every sum below is exact for
# hyperedges of any size, so no size truncation is needed.

# Binary incidence pattern as a sparse dgCMatrix (nodes x hyperedges).
.mmsbm_pattern <- function(incidence) {
  pattern <- if (methods::is(incidence, "sparseMatrix")) {
    Matrix::drop0(incidence)
  } else {
    Matrix::Matrix(incidence != 0, sparse = TRUE)
  }
  pattern <- methods::as(methods::as(methods::as(pattern, "CsparseMatrix"),
                                     "generalMatrix"), "dMatrix")
  pattern@x <- rep(1, length(pattern@x))
  pattern
}

# lambda_e for every hyperedge (Eq. C2), plus the hyperedge sums s_e (E x K).
.mmsbm_lambda <- function(B, u, w) {
  s <- as.matrix(Matrix::crossprod(B, u))
  self_term <- as.numeric(Matrix::crossprod(B, rowSums((u %*% w) * u)))
  list(lambda = 0.5 * (rowSums((s %*% w) * s) - self_term), s = s)
}

# sum_{i<j in V} u_i' w u_j, over every node pair.
.mmsbm_pair_sum <- function(u, w) {
  u_sum <- colSums(u)
  0.5 * (sum((u_sum %*% w) * u_sum) - sum((u %*% w) * u))
}

# Log-likelihood, paper Eq. (5), up to additive constants free of (u, w):
#   -C sum_{i<j} u_i' w u_j + sum_e A_e log(lambda_e).
# With C = 1 this is the authors' HyMMSBM.log_likelihood().
.mmsbm_loglik <- function(B, A, u, w, C = 1) {
  lambda <- .mmsbm_lambda(B, u, w)$lambda
  -C * .mmsbm_pair_sum(u, w) + sum(A * log(lambda))
}

# EM / MAP update of w (Eq. 9 / A2, C absorbed into w as in model.py).
.mmsbm_w_update <- function(B, A, u, w, w_prior) {
  lp <- .mmsbm_lambda(B, u, w)
  multiplier <- A / lp$lambda
  first <- crossprod(lp$s, lp$s * multiplier)
  weighting <- as.numeric(B %*% multiplier)
  second <- crossprod(u, u * weighting)
  numerator <- 0.5 * w * (first - second)
  u_sum <- colSums(u)
  denominator <- 0.5 * (outer(u_sum, u_sum) - crossprod(u))
  numerator / (denominator + w_prior)
}

# EM / MAP update of u (Eq. 8 / A1, C absorbed into w as in model.py).
.mmsbm_u_update <- function(B, A, u, w, u_prior) {
  lp <- .mmsbm_lambda(B, u, w)
  multiplier <- A / lp$lambda
  first <- as.matrix(B %*% (lp$s * multiplier))
  weighting_sum <- as.numeric(B %*% multiplier)
  numerator <- u * ((first - weighting_sum * u) %*% w)
  denominator <- matrix(as.numeric(w %*% colSums(u)), nrow(u), ncol(u),
                        byrow = TRUE) - u %*% w
  numerator / (denominator + u_prior)
}

# Random initialisation, as model.py `_init_w()` then `_init_u()`: w from
# Exp(w_prior) (uniform when the prior is 0), symmetrised from its upper
# triangle and diagonal when assortative; u from Exp(u_prior) (uniform when 0).
.mmsbm_init <- function(n, k, assortative, u_prior, w_prior) {
  draw_w <- if (w_prior == 0) stats::runif(k * k) else stats::rexp(k * k, w_prior)
  w <- matrix(draw_w, k, k)
  w[lower.tri(w)] <- t(w)[lower.tri(w)]
  if (assortative) w <- diag(diag(w), k, k)
  draw_u <- if (u_prior == 0) stats::runif(n * k) else stats::rexp(n * k, u_prior)
  list(u = matrix(draw_u, n, k), w = w)
}

# Has the fit stopped moving between two consecutive steps?
# "parameters" is HyMMSBM.fit()'s rule: ||w - w_old||_F / K < tol and
# ||u - u_old||_F / N < tol. It is not invariant to the rescaling
# u -> c u, w -> w / c^2 that leaves the likelihood unchanged, and under the
# default w prior (with no u prior) the MAP objective keeps gaining by
# shrinking w and growing u, so the rule never fires there.
# "membership" compares the row-normalised memberships instead (largest
# absolute change), which that rescaling leaves unchanged.
.mmsbm_settled <- function(u, w, old_u, old_w, tol, criterion) {
  if (criterion == "parameters") {
    return(sqrt(sum((w - old_w)^2)) / ncol(u) < tol &&
             sqrt(sum((u - old_u)^2)) / nrow(u) < tol)
  }
  change <- abs(.mmsbm_share(u) - .mmsbm_share(old_u))
  max(change[is.finite(change)], 0) < tol
}

# One EM run from (u, w), as HyMMSBM.fit(): w is updated before u in every
# step; with a tolerance, the change between consecutive steps is checked
# every `check_every` steps (0-based step index divisible by it, not the
# first) with .mmsbm_settled(). At the end C is removed from w (or from u,
# as sqrt(C), when w is fixed).
.mmsbm_em <- function(B, A, u, w, max_iter, tol, check_every, u_prior,
                      w_prior, C, criterion = "membership", fixed_u = FALSE,
                      fixed_w = FALSE) {
  old_u <- u
  old_w <- w
  converged <- FALSE
  finite <- TRUE
  iterations <- 0L
  # EM is a fixed-point iteration: step t starts from the parameters of step
  # t - 1 and may stop early, so the steps cannot be vectorised.
  for (it in seq_len(max_iter) - 1L) {
    if (!fixed_w) w <- .mmsbm_w_update(B, A, u, w, w_prior)
    if (!fixed_u) u <- .mmsbm_u_update(B, A, u, w, u_prior)
    iterations <- it + 1L
    if (!all(is.finite(u)) || !all(is.finite(w))) {
      finite <- FALSE
      break
    }
    if (!is.null(tol) && it %% check_every == 0L && it > 0L &&
        .mmsbm_settled(u, w, old_u, old_w, tol, criterion)) {
      converged <- TRUE
      break
    }
    old_u <- u
    old_w <- w
  }
  if (!fixed_w) {
    w <- w / C
  } else if (!fixed_u) {
    u <- u / sqrt(C)
  }
  list(u = u, w = w, iterations = iterations, converged = converged,
       finite = finite,
       loglik = if (finite) .mmsbm_loglik(B, A, u, w, C) else NA_real_)
}

# Validate a single number against a predicate, raising hypernets_bad_input.
.mmsbm_check_number <- function(x, arg, ok, what) {
  if (!(is.numeric(x) && length(x) == 1L && !is.na(x) && is.finite(x) &&
        isTRUE(ok(x)))) {
    .thg_bad_input(sprintf("`%s` must be %s", arg, what))
  }
  invisible(x)
}

.mmsbm_whole <- function(x) isTRUE(all.equal(x, round(x)))

#' Mixed-membership communities of a hypergraph (Hy-MMSBM)
#'
#' Fits the hypergraph mixed-membership stochastic block model of Ruggeri et
#' al. (2023), which gives every node a membership vector over `k`
#' communities from a likelihood model of the whole hypergraph. The weight
#' \eqn{A_e} of a hyperedge \eqn{e} is Poisson with mean
#' \eqn{\lambda_e / \kappa_{|e|}}, where
#' \deqn{\lambda_e = \sum_{i < j \in e} u_i^\top w\, u_j,}
#' \eqn{u} is the non-negative \eqn{N \times k} membership matrix, \eqn{w} the
#' symmetric non-negative \eqn{k \times k} affinity between communities and
#' \eqn{\kappa_n = \frac{n(n-1)}{2}\binom{N-2}{n-2}} normalises for the size
#' \eqn{n} of the hyperedge. Summed over every possible hyperedge up to the
#' largest size \eqn{D}, the normalisers give the constant
#' \eqn{C = \sum_{n=2}^{D} 2 / (n(n-1)) = 2(1 - 1/D)} and the log-likelihood
#' (Eq. 5, up to additive constants that do not depend on \eqn{u}, \eqn{w})
#' \deqn{\mathcal{L} = -C \sum_{i<j} u_i^\top w\, u_j + \sum_e A_e \log \lambda_e.}
#'
#' The fit is the authors' expectation-maximisation (Eqs. 8-9, with the
#' exponential priors of their Appendix A), ported from their reference code:
#' in each step \eqn{w} is updated, then \eqn{u}; \eqn{C} is absorbed into
#' \eqn{w} during the iterations and divided out at the end. Every quantity
#' goes through \eqn{\lambda_e = (s_e^\top w s_e - \sum_{i\in e} u_i^\top w
#' u_i)/2} with \eqn{s_e = \sum_{i \in e} u_i} and sparse products with the
#' incidence matrix, so one step costs \eqn{O(\mathrm{nnz}\,k + E k^2)} and
#' hyperedges of any size (a word used in every document) are handled
#' exactly. EM finds a local optimum only: the fit is run from `nstart`
#' random starts and the one with the highest \eqn{\mathcal{L}} is kept
#' (Algorithm 1). The restarts table reports the log-likelihood, the
#' iterations, the convergence flag of every start, and the adjusted Rand
#' index of its hard partition with the kept one, so the stability of the
#' solution is visible.
#'
#' Hyperedges of size one carry no pair and are outside the model (their
#' \eqn{\lambda_e} is zero); they are left out and counted. `max_size`
#' optionally leaves out larger hyperedges too, as the authors'
#' `--max_hye_size`; the default `NULL` keeps them all. A node left with no
#' hyperedge has \eqn{u_i = 0} after the first step and no membership (`NA`),
#' with a `hypernets_isolated_nodes` warning. A node that does lie in fitted
#' hyperedges can still have its memberships driven to zero (or to subnormal
#' remnants) by the multiplicative EM updates; any row whose total is below
#' the precision of the largest row is treated as having no membership
#' (`NA`) and reported with a `hypernets_collapsed_membership` warning,
#' rather than normalising remnants into spurious exact mixtures.
#'
#' The membership of node \eqn{i} in community \eqn{k} is
#' \eqn{u_{ik} / \sum_q u_{iq}}; the hard `community` is its largest entry.
#' The size of \eqn{u_i} reflects how active the node is, its direction the
#' communities it mixes; at the maximum of the likelihood many memberships
#' sit at or near 0 and 1, and mixing is read from the nodes that do not.
#' Communities are numbered by decreasing total membership. On a
#' [text_hypergraph()] with documents as nodes, the communities are topics
#' and the memberships the documents' topic shares.
#'
#' @param hg A hypernets `net_hg` ([text_hypergraph()], [group_hypergraph()],
#'   [window_hypergraph()], ...).
#' @param k Number of communities, a whole number between 1 and the number of
#'   nodes.
#' @param assortative If `TRUE`, `w` is diagonal (communities interact only
#'   within themselves); default `FALSE`, the full affinity.
#' @param nstart Number of random starts (default `10L`, the authors'
#'   `--training_rounds`).
#' @param max_iter Maximum EM steps per start (default `5000L`).
#' @param tol Convergence tolerance (default `1e-5`), compared with the
#'   change between two consecutive steps every `check_every` steps. A start
#'   that hits `max_iter` first is reported as not converged; if the kept
#'   start did not converge, a `hypernets_no_converge` warning is raised.
#' @param check_every Steps between convergence checks (default `10L`, as
#'   the authors').
#' @param criterion What must stop changing: `"membership"` (default), the
#'   largest absolute change of any membership; or `"parameters"`, the
#'   authors' rule, the Frobenius change of `w` divided by `k` and of `u`
#'   divided by the number of nodes. The likelihood is unchanged when `u` is
#'   multiplied by \eqn{c} and `w` divided by \eqn{c^2}; with the default
#'   `w_prior` and no `u_prior` the fit keeps drifting along that direction
#'   (`w` shrinks, `u` grows) while the memberships and the log-likelihood
#'   have settled, so `"parameters"` rarely fires there.
#' @param w_prior,u_prior Rates of the exponential priors on `w` and `u`
#'   (MAP estimation); `0` means none (maximum likelihood). The defaults
#'   `w_prior = 1`, `u_prior = 0` are the authors' "half-Bayesian" setting
#'   that fixes the scale of the parameters. As in the authors' code, the
#'   `w` prior acts on \eqn{C w} (the parametrisation used during the
#'   iterations).
#' @param max_size `NULL` (default: every hyperedge) or the largest
#'   hyperedge size kept; it is also the \eqn{D} in \eqn{C}. With `NULL`,
#'   \eqn{D} is the largest size in the data.
#' @param edge_weights `NULL` or positive hyperedge weights \eqn{A_e}, one
#'   per hyperedge or one recycled. `NULL` uses the window counts of a
#'   [window_hypergraph()] and 1 otherwise. The incidence is read as binary:
#'   the stored cell weights (word counts, tf-idf) do not enter the model.
#' @param seed `NULL` (default: the current random stream) or a whole number
#'   for the random starts; the caller's stream is restored on exit.
#' @return An object of class `net_hg_mmsbm`, a list read through
#'   `hg_get(x, what = )`:
#'   * `"membership"` (default): one row per node and community, columns
#'     `node`, `community`, `membership` (sums to one over a node's rows) and
#'     `membership_weight` (the fitted, unnormalised membership parameter
#'     `u`);
#'   * `"nodes"`: one row per node, `node`, `community` (the largest
#'     membership) and `membership` (its value);
#'   * `"affinity"`: one row per community pair with `from <= to`, columns
#'     `from`, `to`, `affinity` (the fitted `w`);
#'   * `"restarts"`: one row per start, `run`, `log_likelihood`,
#'     `iterations`, `converged`, `best`, `ari_to_best` (adjusted Rand index
#'     of the start's
#'     hard partition with the kept one).
#'
#'   `print()` gives the fit, `summary()` one row per community (`community`,
#'   `size` hard-assigned nodes, `mass` total membership, `mean_membership`
#'   of its assigned nodes, `w_within`), and `plot()` the memberships, the
#'   affinity or the restarts. Raises `hypernets_bad_input` for a non-hypergraph,
#'   an invalid argument or no hyperedge of size two or more.
#' @references
#' Ruggeri, N., Contisciani, M., Battiston, F., & De Bacco, C. (2023).
#' Community detection in large hypergraphs. *Science Advances*, 9(28),
#' eadg9159. \doi{10.1126/sciadv.adg9159}
#'
#' Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
#' Classification*, 2, 193--218. \doi{10.1007/BF01908075}
#' @seealso [hg_membership()] for soft membership read from a hard spectral
#'   partition.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots and salt",
#'   cooking_2 = "this soup recipe needs onions, carrots and salt",
#'   cooking_3 = "carrots and onions make a sweet soup",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars",
#'   space_3 = "a distant galaxy full of stars",
#'   both = "astronomers eat soup with carrots under the stars"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "of", "under"))
#' fit <- hg_mmsbm(hg, k = 2, seed = 1)
#' fit
#' hg_get(fit, what = "nodes")
#' hg_get(fit, what = "restarts")
#' @export
hg_mmsbm <- function(hg, k, assortative = FALSE, nstart = 10L,
                     max_iter = 5000L, tol = 1e-5, check_every = 10L,
                     criterion = c("membership", "parameters"),
                     w_prior = 1, u_prior = 0, max_size = NULL,
                     edge_weights = NULL, seed = NULL) {
  .thg_check_hg(hg)
  criterion <- match.arg(criterion)
  n <- hg$n_nodes
  whole <- function(lo, hi = Inf) \(x) .mmsbm_whole(x) && x >= lo && x <= hi
  .mmsbm_check_number(k, "k", whole(1, n),
                      sprintf("a whole number between 1 and %d (the nodes)", n))
  .mmsbm_check_number(nstart, "nstart", whole(1), "a whole number >= 1")
  .mmsbm_check_number(max_iter, "max_iter", whole(1), "a whole number >= 1")
  .mmsbm_check_number(check_every, "check_every", whole(1),
                      "a whole number >= 1")
  .mmsbm_check_number(tol, "tol", \(x) x > 0, "a positive number")
  .mmsbm_check_number(w_prior, "w_prior", \(x) x >= 0, "a number >= 0")
  .mmsbm_check_number(u_prior, "u_prior", \(x) x >= 0, "a number >= 0")
  if (!is.null(max_size)) {
    .mmsbm_check_number(max_size, "max_size", whole(2),
                        "NULL or a whole number >= 2")
  }
  if (!is.null(seed)) {
    .mmsbm_check_number(seed, "seed", .mmsbm_whole, "NULL or a whole number")
  }
  if (!(is.logical(assortative) && length(assortative) == 1L &&
        !is.na(assortative))) {
    .thg_bad_input("`assortative` must be TRUE or FALSE")
  }
  k <- as.integer(k)

  B_all <- .mmsbm_pattern(hg$incidence)
  m <- ncol(B_all)
  A_all <- edge_weights %||% hg$window_counts %||% rep(1, m)
  if (!(is.numeric(A_all) && length(A_all) %in% c(1L, m) &&
        all(is.finite(A_all)) && all(A_all > 0))) {
    .thg_bad_input(sprintf(
      "`edge_weights` must be positive numbers, one or one per hyperedge (%d)", m))
  }
  A_all <- rep_len(as.numeric(A_all), m)

  sizes <- Matrix::colSums(B_all)
  singleton <- sizes < 2
  oversize <- !singleton & !is.null(max_size) & sizes > (max_size %||% Inf)
  keep <- !singleton & !oversize
  if (!any(keep)) {
    .thg_bad_input("the hypergraph has no hyperedge of size 2 or more to fit")
  }
  B <- B_all[, keep, drop = FALSE]
  A <- A_all[keep]
  D <- as.integer(max_size %||% max(sizes[keep]))
  C <- 2 * (1 - 1 / D)
  nodes <- rownames(hg$incidence) %||% hg$nodes
  isolated <- Matrix::rowSums(B) == 0

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    saved_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
    on.exit(.thg_rng_restore(had_seed, saved_seed), add = TRUE, after = FALSE)
    set.seed(as.integer(seed))
  }
  # Starts drawn serially, one after the other, from one stream.
  fits <- lapply(seq_len(nstart), \(start) {
    init <- .mmsbm_init(n, k, assortative, u_prior, w_prior)
    .mmsbm_em(B, A, init$u, init$w, max_iter = as.integer(max_iter), tol = tol,
              check_every = as.integer(check_every), u_prior = u_prior,
              w_prior = w_prior, C = C, criterion = criterion)
  })
  finite <- vapply(fits, \(f) f$finite, logical(1L))
  if (!any(finite)) {
    stop(errorCondition(
      "every Hy-MMSBM start produced non-finite parameters",
      class = "hypernets_no_converge", call = NULL))
  }
  if (!all(finite)) {
    warning(warningCondition(sprintf(
      "%d of %d Hy-MMSBM starts produced non-finite parameters and were discarded",
      sum(!finite), nstart), class = "hypernets_no_converge", call = NULL))
  }
  loglik <- vapply(fits, \(f) f$loglik, numeric(1L))
  best <- which.max(replace(loglik, !finite, -Inf))
  fit <- fits[[best]]
  converged <- vapply(fits, \(f) f$converged, logical(1L))
  if (!converged[best]) {
    warning(warningCondition(sprintf(
      "the kept Hy-MMSBM start did not converge in %d steps (tol = %g)",
      as.integer(max_iter), tol), class = "hypernets_no_converge", call = NULL))
  }

  # Communities numbered by decreasing total membership (ties: fitted order).
  share <- .mmsbm_share(fit$u)
  mass <- colSums(share, na.rm = TRUE)
  ord <- order(-mass, seq_len(k))
  labels <- paste("Community", seq_len(k))
  u <- fit$u[, ord, drop = FALSE]
  w <- fit$w[ord, ord, drop = FALSE]
  share <- share[, ord, drop = FALSE]
  dimnames(u) <- dimnames(share) <- list(nodes, labels)
  dimnames(w) <- list(labels, labels)

  hard_of <- \(f) .mmsbm_hard(.mmsbm_share(f$u))
  best_hard <- hard_of(fit)
  ari <- vapply(seq_len(nstart), \(j) {
    if (!finite[j]) return(NA_real_)
    .thg_ari(best_hard[!isolated], hard_of(fits[[j]])[!isolated])
  }, numeric(1L))

  if (any(isolated)) {
    warning(warningCondition(sprintf(
      "%d node(s) belong to no hyperedge of size 2 or more and have no membership (NA): %s",
      sum(isolated), paste(utils::head(nodes[isolated], 5L), collapse = ", ")),
      class = "hypernets_isolated_nodes", call = NULL))
  }
  collapsed <- !isolated & rowSums(is.na(share)) > 0L
  if (any(collapsed)) {
    warning(warningCondition(sprintf(paste0(
      "%d node(s) lie in fitted hyperedges but their memberships collapsed ",
      "to zero under EM (no membership, NA): %s. Try more starts, another ",
      "k or `assortative`."),
      sum(collapsed), paste(utils::head(nodes[collapsed], 5L), collapse = ", ")),
      class = "hypernets_collapsed_membership", call = NULL))
  }
  # subnormal affinities and memberships are zero to working precision
  w[abs(w) < .Machine$double.xmin] <- 0
  u[abs(u) < .Machine$double.xmin] <- 0

  structure(list(
    u = u, w = w, membership = share, loglik = fit$loglik, k = k,
    restarts = data.frame(run = seq_len(nstart), log_likelihood = loglik,
                          iterations = vapply(fits, \(f) f$iterations,
                                              integer(1L)),
                          converged = converged,
                          best = seq_len(nstart) == best, ari_to_best = ari),
    n_nodes = n, n_hyperedges = sum(keep),
    dropped = c(singleton = sum(singleton), oversize = sum(oversize)),
    isolated = nodes[isolated], collapsed = nodes[collapsed],
    max_size = D, C = C,
    params = list(assortative = assortative, nstart = as.integer(nstart),
                  max_iter = as.integer(max_iter), tol = tol,
                  check_every = as.integer(check_every),
                  criterion = criterion, w_prior = w_prior,
                  u_prior = u_prior, max_size = max_size, seed = seed)
  ), class = "net_hg_mmsbm")
}

# Row-normalised memberships; NA for a node whose u is all zero.
.mmsbm_share <- function(u) {
  total <- rowSums(u)
  out <- u / total
  # A row whose total is below the precision of the largest row has
  # collapsed (EM's multiplicative updates drove it to zero or to subnormal
  # remnants); dividing the remnants would report arbitrary exact fractions
  # (1/3, 1/2, ...) as if they were mixtures. It has no membership.
  collapsed <- total <= .Machine$double.eps * max(total, 0)
  out[collapsed, ] <- NA_real_
  out
}

# Column of the largest membership of every node; NA for a node with none.
.mmsbm_hard <- function(share) {
  vapply(seq_len(nrow(share)), \(i) {
    if (anyNA(share[i, ])) NA_integer_ else which.max(share[i, ])
  }, integer(1L))
}

# The hard community table, one row per node.
.mmsbm_nodes <- function(x) {
  share <- x$membership
  top <- .mmsbm_hard(share)
  data.frame(
    node = rownames(share),
    community = colnames(share)[top],
    membership = share[cbind(seq_len(nrow(share)), top)],
    stringsAsFactors = FALSE
  )
}

#' @rdname hg_mmsbm
#' @param x A `net_hg_mmsbm` object.
#' @param what Which table (see Value).
#' @param top `NULL` (default, every row) or the number of first rows to
#'   return.
#' @param ... Unused; for S3 consistency.
#' @export
hg_get.net_hg_mmsbm <- function(x, what = c("membership", "nodes",
                                            "affinity", "restarts"), ...,
                                top = NULL) {
  what <- match.arg(what)
  out <- switch(what,
    membership = data.frame(
      node = rep(rownames(x$u), times = x$k),
      community = rep(colnames(x$u), each = nrow(x$u)),
      membership = as.numeric(x$membership),
      membership_weight = as.numeric(x$u),
      stringsAsFactors = FALSE
    )[order(rep(seq_len(nrow(x$u)), times = x$k),
            rep(seq_len(x$k), each = nrow(x$u))), , drop = FALSE],
    nodes = .mmsbm_nodes(x),
    affinity = {
      pairs <- which(upper.tri(x$w, diag = TRUE), arr.ind = TRUE)
      pairs <- pairs[order(pairs[, "row"], pairs[, "col"]), , drop = FALSE]
      data.frame(from = rownames(x$w)[pairs[, "row"]],
                 to = colnames(x$w)[pairs[, "col"]],
                 affinity = x$w[pairs], stringsAsFactors = FALSE)
    },
    restarts = x$restarts
  )
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' @rdname hg_mmsbm
#' @param n Number of rows of the default table to print. Default `10`.
#' @export
print.net_hg_mmsbm <- function(x, n = 10L, ...) {
  cat(sprintf("Hy-MMSBM mixed-membership communities (%s affinity), k = %d\n",
              if (isTRUE(x$params$assortative)) "assortative" else "full",
              x$k))
  cat(sprintf(paste0("  Nodes: %d | Hyperedges fitted: %d (left out: %d of ",
                     "size 1, %d above max_size) | D = %d\n"),
              x$n_nodes, x$n_hyperedges, x$dropped[["singleton"]],
              x$dropped[["oversize"]], x$max_size))
  r <- x$restarts
  kept <- r[r$best, , drop = FALSE]
  cat(sprintf(paste0("  Log-likelihood: %.6g | Starts: %d (%d converged) | ",
                     "Kept start: %d, %d steps, %s\n"),
              x$loglik, nrow(r), sum(r$converged), kept$run,
              kept$iterations,
              if (isTRUE(kept$converged)) "converged" else "not converged"))
  n_none <- length(x$isolated) + length(x$collapsed)
  if (n_none > 0L) {
    cat(sprintf("  No membership (NA): %d node(s) -- %d isolated, %d collapsed\n",
                n_none, length(x$isolated), length(x$collapsed)))
  }
  .ho_print_table(x, n)
  invisible(x)
}

#' @rdname result-summary
#' @export
summary.net_hg_mmsbm <- function(object, ...) {
  .ho_summary(object, list(communities = .mmsbm_communities(object)))
}

# One row per community: hard size, total membership mass, mean membership
# of its hard members, and the within-community affinity.
.mmsbm_communities <- function(object) {
  hard <- .mmsbm_nodes(object)
  labels <- colnames(object$u)
  data.frame(
    community = labels,
    size = vapply(labels, \(l) sum(hard$community == l, na.rm = TRUE),
                  integer(1L), USE.NAMES = FALSE),
    mass = colSums(object$membership, na.rm = TRUE),
    mean_membership = vapply(labels, \(l) {
      own <- hard$membership[hard$community %in% l]
      if (length(own)) mean(own) else NA_real_
    }, numeric(1L), USE.NAMES = FALSE),
    w_within = diag(object$w),
    row.names = NULL, stringsAsFactors = FALSE
  )
}

#' @rdname hg_mmsbm
#' @param type For `plot()`: `"membership"` (default) stacks each node's
#'   memberships in a bar, nodes grouped in panels by hard community;
#'   `"affinity"` draws `w` as a labelled tile map; `"restarts"` the
#'   log-likelihood of every start, shaped by convergence.
#' @param labels For `plot(type = "membership")`: show node names on the axis
#'   (default: when there are at most 40 nodes).
#' @export
plot.net_hg_mmsbm <- function(x, type = c("membership", "affinity",
                                          "restarts"),
                              labels = NULL, ...) {
  type <- match.arg(type)
  palette <- rep_len(c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                       "#D55E00", "#CC79A7", "#999999", "#000000"), x$k)
  if (type == "affinity") {
    d <- expand.grid(from = colnames(x$w), to = colnames(x$w),
                     stringsAsFactors = FALSE)
    d$w <- x$w[cbind(match(d$from, colnames(x$w)), match(d$to, colnames(x$w)))]
    d$from <- factor(d$from, levels = rev(colnames(x$w)))
    d$to <- factor(d$to, levels = colnames(x$w))
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data$to, y = .data$from,
                                      fill = .data$w)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = signif(.data$w, 3)), size = 3) +
        ggplot2::scale_fill_gradient(low = "#FFFFFF", high = "#4A6FE3",
                                     name = "affinity w") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::theme_minimal(base_size = 12)
    )
  }
  if (type == "restarts") {
    d <- x$restarts
    d$status <- ifelse(d$converged, "converged", "not converged")
    return(
      ggplot2::ggplot(d, ggplot2::aes(x = .data$run, y = .data$log_likelihood,
                                      shape = .data$status)) +
        ggplot2::geom_point(size = 2.5, colour = "#0072B2") +
        ggplot2::geom_point(data = d[d$best, , drop = FALSE], size = 5,
                            shape = 1, colour = "#D55E00") +
        ggplot2::scale_shape_manual(values = c(converged = 16,
                                               `not converged` = 4),
                                    name = NULL) +
        ggplot2::labs(x = "start", y = "log-likelihood",
                      caption = "open circle: the kept start") +
        ggplot2::theme_minimal(base_size = 12)
    )
  }
  memb <- hg_get(x, what = "membership")
  hard <- .mmsbm_nodes(x)
  hard <- hard[!is.na(hard$community), , drop = FALSE]
  node_order <- hard$node[order(match(hard$community, colnames(x$u)),
                                -hard$membership, hard$node)]
  memb <- memb[memb$node %in% node_order, , drop = FALSE]
  memb$node <- factor(memb$node, levels = node_order)
  memb$community <- factor(memb$community, levels = colnames(x$u))
  memb$panel <- factor(hard$community[match(memb$node, hard$node)],
                       levels = colnames(x$u))
  show_labels <- labels %||% (length(node_order) <= 40L)
  p <- ggplot2::ggplot(memb, ggplot2::aes(x = .data$node, y = .data$membership,
                                          fill = .data$community)) +
    ggplot2::geom_col(width = 1, colour = if (show_labels) "white" else NA,
                      linewidth = 0.2) +
    ggplot2::facet_grid(cols = ggplot2::vars(.data$panel), scales = "free_x",
                        space = "free_x") +
    ggplot2::scale_fill_manual(values = stats::setNames(palette,
                                                        colnames(x$u)),
                               name = "community") +
    ggplot2::scale_y_continuous(limits = c(0, 1 + 1e-9), expand = c(0, 0)) +
    ggplot2::labs(x = "node (panel: largest membership)", y = "membership") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.x = ggplot2::element_blank())
  if (isTRUE(show_labels)) {
    p + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90,
                                                           hjust = 1,
                                                           vjust = 0.5))
  } else {
    p + ggplot2::theme(axis.text.x = ggplot2::element_blank())
  }
}
