testthat::skip_on_cran()

# Hy-MMSBM (Ruggeri et al. 2023). The authors' code is the oracle in
# local_testing_and_equivalence/test-oracle-mmsbm.R; these tests check the
# formulas against brute-force sums, EM invariants, and the verb's contract.

# A hypergraph from a list of hyperedges (character node names).
.mm_hg <- function(edges) {
  group_hypergraph(
    data.frame(member = unlist(edges),
               edge = sprintf("e%03d", rep(seq_along(edges), lengths(edges)))),
    actor = "member", group = "edge")
}

# Two planted groups a1..a6 and b1..b6, dense inside, two bridges.
.mm_planted <- function() {
  set.seed(42)
  a <- paste0("a", 1:6)
  b <- paste0("b", 1:6)
  within <- c(lapply(1:25, \(i) sample(a, sample(2:4, 1L))),
              lapply(1:25, \(i) sample(b, sample(2:4, 1L))))
  .mm_hg(c(within, list(c("a1", "b1"), c("a2", "b2", "b3"))))
}

# Brute-force EM step from the paper's Eqs. (7)-(9) with C absorbed into w
# (C = 1) and exponential priors (Appendix A): explicit sums over every
# hyperedge, node pair and community pair.
.mm_brute_step <- function(edges, A, u, w, w_prior = 0, u_prior = 0) {
  k <- ncol(u)
  n <- nrow(u)
  pairs_of <- \(e) if (length(e) < 2L) matrix(integer(), 0L, 2L) else
    t(utils::combn(e, 2L))
  lambda_of <- \(u, w) vapply(edges, \(e) {
    p <- pairs_of(e)
    sum(vapply(seq_len(nrow(p)), \(r) sum(u[p[r, 1L], ] %*% w %*% u[p[r, 2L], ]),
               numeric(1L)))
  }, numeric(1L))
  lambda <- lambda_of(u, w)
  # w: numerator sum_e A_e sum_{i<j in e} w_kq (u_ik u_jq + u_iq u_jk) / 2 / lambda_e
  num_w <- Reduce(`+`, lapply(seq_along(edges), \(e) {
    p <- pairs_of(edges[[e]])
    Reduce(`+`, lapply(seq_len(nrow(p)), \(r) {
      ui <- u[p[r, 1L], ]
      uj <- u[p[r, 2L], ]
      (outer(ui, uj) + outer(uj, ui)) / 2
    })) * w * A[e] / lambda[e]
  }))
  all_pairs <- t(utils::combn(n, 2L))
  den_w <- Reduce(`+`, lapply(seq_len(nrow(all_pairs)), \(r) {
    ui <- u[all_pairs[r, 1L], ]
    uj <- u[all_pairs[r, 2L], ]
    (outer(ui, uj) + outer(uj, ui)) / 2
  }))
  w_new <- num_w / (den_w + w_prior)
  lambda <- lambda_of(u, w_new)
  # u: numerator sum_{e ni i} A_e sum_{j in e, j != i} sum_q u_ik w_kq u_jq / lambda_e
  u_new <- t(vapply(seq_len(n), \(i) {
    num <- Reduce(`+`, lapply(seq_along(edges), \(e) {
      if (!(i %in% edges[[e]])) return(numeric(k))
      others <- setdiff(edges[[e]], i)
      u[i, ] * as.numeric(w_new %*% colSums(u[others, , drop = FALSE])) *
        A[e] / lambda[e]
    }))
    den <- as.numeric(w_new %*% colSums(u[-i, , drop = FALSE]))
    num / (den + u_prior)
  }, numeric(k)))
  list(u = u_new, w = w_new)
}

.mm_sparse <- function(edges, n) {
  Matrix::sparseMatrix(i = unlist(edges),
                       j = rep(seq_along(edges), lengths(edges)),
                       x = 1, dims = c(n, length(edges)))
}

test_that("lambda_e via Eq. (C2) equals the explicit pair sum (Eq. 2)", {
  set.seed(1)
  n <- 7L
  edges <- list(1:3, c(2L, 5L), c(1L, 4L, 6L, 7L), 1:7)
  u <- matrix(stats::runif(n * 3L), n, 3L)
  w <- matrix(stats::runif(9L), 3L, 3L)
  w <- (w + t(w)) / 2
  fast <- hypernets:::.mmsbm_lambda(.mm_sparse(edges, n), u, w)$lambda
  slow <- vapply(edges, \(e) {
    p <- t(utils::combn(e, 2L))
    sum(vapply(seq_len(nrow(p)), \(r) sum(u[p[r, 1L], ] %*% w %*% u[p[r, 2L], ]),
               numeric(1L)))
  }, numeric(1L))
  expect_equal(fast, slow, tolerance = 1e-13)
})

test_that("one EM step equals the brute-force Eqs. (7)-(9), with and without priors", {
  set.seed(2)
  n <- 6L
  edges <- list(1:3, c(2L, 4L), c(3L, 5L, 6L), c(1L, 6L), 2:5)
  A <- c(1, 2, 1.5, 1, 3)
  B <- .mm_sparse(edges, n)
  u <- matrix(stats::runif(n * 2L), n, 2L)
  w <- matrix(c(1.2, 0.3, 0.3, 0.8), 2L, 2L)
  lapply(list(c(0, 0), c(1, 0), c(2.5, 0.5)), \(pr) {
    brute <- .mm_brute_step(edges, A, u, w, w_prior = pr[1L], u_prior = pr[2L])
    w_new <- hypernets:::.mmsbm_w_update(B, A, u, w, pr[1L])
    u_new <- hypernets:::.mmsbm_u_update(B, A, u, w_new, pr[2L])
    expect_equal(w_new, brute$w, tolerance = 1e-12)
    expect_equal(u_new, brute$u, tolerance = 1e-12)
  })
})

test_that("maximum-likelihood EM never decreases the log-likelihood", {
  set.seed(3)
  n <- 10L
  edges <- lapply(1:25, \(i) sort(sample(n, sample(2:5, 1L))))
  B <- .mm_sparse(edges, n)
  A <- rep(1, length(edges))
  init <- hypernets:::.mmsbm_init(n, 3L, FALSE, 0, 0)
  u <- init$u
  w <- init$w
  ll <- vapply(seq_len(60L), \(t) {
    w <<- hypernets:::.mmsbm_w_update(B, A, u, w, 0)
    u <<- hypernets:::.mmsbm_u_update(B, A, u, w, 0)
    hypernets:::.mmsbm_loglik(B, A, u, w, C = 1)
  }, numeric(1L))
  expect_true(all(diff(ll) > -1e-9 * abs(ll[-1L])))
})

test_that("the log-likelihood is invariant to u -> c u, w -> w / c^2", {
  set.seed(4)
  n <- 8L
  edges <- lapply(1:12, \(i) sort(sample(n, sample(2:4, 1L))))
  B <- .mm_sparse(edges, n)
  A <- stats::runif(length(edges), 1, 2)
  u <- matrix(stats::runif(n * 2L), n, 2L)
  w <- matrix(c(1, 0.2, 0.2, 0.5), 2L, 2L)
  base <- hypernets:::.mmsbm_loglik(B, A, u, w, C = 1.5)
  expect_equal(hypernets:::.mmsbm_loglik(B, A, 3 * u, w / 9, C = 1.5), base,
               tolerance = 1e-12)
})

test_that("a repeated hyperedge is the same as weight 2", {
  set.seed(5)
  n <- 7L
  edges <- lapply(1:10, \(i) sort(sample(n, sample(2:4, 1L))))
  init <- hypernets:::.mmsbm_init(n, 2L, FALSE, 0, 1)
  twice <- hypernets:::.mmsbm_em(
    .mm_sparse(c(edges, edges[1L]), n), rep(1, 11L), init$u, init$w,
    max_iter = 50L, tol = NULL, check_every = 10L, u_prior = 0, w_prior = 1,
    C = 1.5)
  weighted <- hypernets:::.mmsbm_em(
    .mm_sparse(edges, n), c(2, rep(1, 9L)), init$u, init$w,
    max_iter = 50L, tol = NULL, check_every = 10L, u_prior = 0, w_prior = 1,
    C = 1.5)
  expect_equal(twice$u, weighted$u, tolerance = 1e-12)
  expect_equal(twice$w, weighted$w, tolerance = 1e-12)
  expect_equal(twice$loglik, weighted$loglik, tolerance = 1e-12)
})

test_that("the reported log-likelihood is Eq. (5) at the rescaled w", {
  hg <- .mm_planted()
  fit <- hg_mmsbm(hg, k = 2, nstart = 2L, seed = 1)
  B <- hypernets:::.mmsbm_pattern(hg$incidence)
  A <- rep(1, ncol(B))
  expect_identical(fit$C, 2 * (1 - 1 / 4))
  eq5 <- -fit$C * hypernets:::.mmsbm_pair_sum(fit$u, fit$w) +
    sum(A * log(hypernets:::.mmsbm_lambda(B, fit$u, fit$w)$lambda))
  expect_equal(fit$loglik, eq5, tolerance = 1e-12)
  # = the C = 1 form at C w, minus sum(A) log C
  expect_equal(fit$loglik,
               hypernets:::.mmsbm_loglik(B, A, fit$u, fit$C * fit$w, C = 1) -
                 sum(A) * log(fit$C), tolerance = 1e-12)
})

test_that("hg_mmsbm recovers two planted groups", {
  hg <- .mm_planted()
  truth <- substr(hg$nodes, 1L, 1L)
  lapply(c(FALSE, TRUE), \(assortative) {
    fit <- hg_mmsbm(hg, k = 2, assortative = assortative, seed = 1)
    nodes <- hg_get(fit, what = "nodes")
    expect_identical(nodes$node, hg$nodes)
    expect_equal(hypernets:::.thg_ari(truth, nodes$community), 1)
    expect_true(all(fit$restarts$converged))
    if (assortative) {
      expect_identical(fit$w[1L, 2L], 0)
    }
  })
})

test_that("hg_get() returns the documented tidy tables", {
  hg <- .mm_planted()
  fit <- hg_mmsbm(hg, k = 3, nstart = 4L, seed = 2)
  expect_s3_class(fit, "net_hg_mmsbm")

  memb <- hg_get(fit)
  expect_identical(class(memb), "data.frame")
  expect_named(memb, c("node", "community", "membership",
                       "membership_weight"))
  expect_identical(nrow(memb), 12L * 3L)
  expect_identical(attr(memb, "row.names"), seq_len(nrow(memb)))
  sums <- vapply(split(memb$membership, memb$node), sum, numeric(1L))
  expect_equal(unname(sums), rep(1, 12L), tolerance = 1e-12)
  expect_true(all(memb$membership >= 0 & memb$membership_weight >= 0))

  nodes <- hg_get(fit, what = "nodes")
  expect_named(nodes, c("node", "community", "membership"))
  expect_identical(nrow(nodes), 12L)
  expect_true(all(nodes$membership >= 1 / 3 - 1e-12))

  aff <- hg_get(fit, what = "affinity")
  expect_named(aff, c("from", "to", "affinity"))
  expect_identical(nrow(aff), 6L)
  expect_true(all(aff$affinity >= 0))

  rs <- hg_get(fit, what = "restarts")
  expect_named(rs, c("run", "log_likelihood", "iterations", "converged",
                     "best", "ari_to_best"))
  expect_identical(nrow(rs), 4L)
  expect_identical(sum(rs$best), 1L)
  expect_identical(rs$log_likelihood[rs$best], max(rs$log_likelihood))
  expect_identical(rs$ari_to_best[rs$best], 1)
  expect_identical(fit$loglik, max(rs$log_likelihood))

  expect_identical(nrow(hg_get(fit, top = 5)), 5L)

  s <- summary(fit)$communities
  expect_named(s, c("community", "size", "mass", "mean_membership",
                    "w_within"))
  expect_identical(sum(s$size), 12L)
  expect_equal(sum(s$mass), 12, tolerance = 1e-12)
  # communities are numbered by decreasing total membership
  expect_false(is.unsorted(rev(s$mass)))
})

test_that("the seed makes fits reproducible and restores the caller's stream", {
  hg <- .mm_planted()
  set.seed(99)
  before <- stats::runif(1L)
  set.seed(99)
  a <- hg_mmsbm(hg, k = 2, nstart = 3L, seed = 7)
  after <- stats::runif(1L)
  b <- hg_mmsbm(hg, k = 2, nstart = 3L, seed = 7)
  expect_identical(before, after)
  expect_identical(a$u, b$u)
  expect_identical(a$restarts, b$restarts)
})

test_that("size-1 hyperedges are left out and max_size truncates", {
  hg <- .mm_hg(list(c("a", "b", "c"), c("a", "b"), "c", c("c", "d"),
                    c("a", "b", "c", "d", "e"), c("d", "e")))
  fit <- hg_mmsbm(hg, k = 2, nstart = 2L, seed = 1)
  expect_identical(fit$dropped[["singleton"]], 1L)
  expect_identical(fit$n_hyperedges, 5L)
  expect_identical(fit$max_size, 5L)
  small <- hg_mmsbm(hg, k = 2, nstart = 2L, seed = 1, max_size = 3)
  expect_identical(small$dropped[["oversize"]], 1L)
  expect_identical(small$n_hyperedges, 4L)
  expect_identical(small$C, 2 * (1 - 1 / 3))
})

test_that("a node with no hyperedge of size 2+ gets NA and a classed warning", {
  hg <- .mm_hg(list(c("a", "b", "c"), c("a", "b"), c("b", "c"), "z"))
  expect_warning(fit <- hg_mmsbm(hg, k = 2, nstart = 2L, seed = 1),
                 class = "hypernets_isolated_nodes")
  nodes <- hg_get(fit, what = "nodes")
  expect_true(is.na(nodes$community[nodes$node == "z"]))
  expect_identical(fit$isolated, "z")
})

test_that("non-convergence is surfaced as hypernets_no_converge", {
  hg <- .mm_planted()
  expect_warning(fit <- hg_mmsbm(hg, k = 2, nstart = 2L, max_iter = 3L,
                                 seed = 1),
                 class = "hypernets_no_converge")
  expect_false(any(fit$restarts$converged))
})

test_that("invalid input raises hypernets_bad_input", {
  hg <- .mm_planted()
  bad <- \(...) expect_error(hg_mmsbm(...), class = "hypernets_bad_input")
  bad(list(), k = 2)
  bad(hg, k = 0)
  bad(hg, k = 1.5)
  bad(hg, k = 13)
  bad(hg, k = 2, nstart = 0)
  bad(hg, k = 2, tol = -1)
  bad(hg, k = 2, w_prior = -1)
  bad(hg, k = 2, max_size = 1)
  bad(hg, k = 2, edge_weights = -1)
  bad(hg, k = 2, edge_weights = c(1, 2))
  bad(hg, k = 2, assortative = NA)
  bad(hg, k = 2, seed = "a")
  bad(.mm_hg(list("a", "b")), k = 1)
})

test_that("hg_mmsbm fits a text hypergraph with a word in every document", {
  hg <- text_hypergraph(c(
    cooking_1 = "today simmer the soup with onions and carrots and salt",
    cooking_2 = "today this soup recipe needs onions and salt on a cold night",
    cooking_3 = "today carrots and onions make a sweet soup",
    space_1 = "today the telescope revealed a distant galaxy and stars",
    space_2 = "today astronomers aimed the telescope at the stars all night",
    space_3 = "today a distant galaxy full of stars"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all",
                    "of"))
  # "today" is a hyperedge on every document
  expect_identical(max(Matrix::colSums(hg$incidence > 0)), 6)
  fit <- hg_mmsbm(hg, k = 2, seed = 1)
  expect_identical(fit$max_size, 6L)
  expect_true(all(fit$restarts$converged))
  expect_identical(nrow(hg_get(fit, what = "nodes")), 6L)
})

test_that("print, summary and plot work", {
  fit <- hg_mmsbm(.mm_planted(), k = 2, nstart = 3L, seed = 1)
  expect_output(print(fit), "Hy-MMSBM")
  expect_s3_class(plot(fit), "ggplot")
  expect_s3_class(plot(fit, type = "affinity"), "ggplot")
  expect_s3_class(plot(fit, type = "restarts"), "ggplot")
  expect_s3_class(plot(fit, labels = FALSE), "ggplot")
})

test_that("collapsed memberships are NA with a classed warning, not fake mixtures", {
  u <- rbind(c(0.2, 0.8, 0), c(4.9e-324, 4.9e-324, 4.9e-324), c(0, 0, 0))
  share <- .mmsbm_share(u)
  expect_equal(share[1L, ], c(0.2, 0.8, 0))
  expect_true(all(is.na(share[2L, ])))   # underflow remnants, not 1/3 each
  expect_true(all(is.na(share[3L, ])))
})

test_that("community-labelled tables feed hg_agreement(); NA labels are dropped", {
  a <- data.frame(node = c("p", "q", "r", "s"), community = c("A", "A", "B", NA))
  b <- data.frame(node = c("p", "q", "r", "s"), cluster = c("x", "x", "y", "y"))
  expect_warning(out <- hg_agreement(a, b), class = "hypernets_missing_labels")
  expect_identical(out$n, 3L)
  expect_equal(out$ari, 1)
})
