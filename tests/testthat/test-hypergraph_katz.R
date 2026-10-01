# Katz centrality engine (.hg_katz_fit) behind hg_centrality(type = "katz"):
# closed form on a complete graph, the truncated walk series (a different
# algorithm from the Cholesky solve), igraph's alpha centrality, dense/sparse
# parity, invariants and error classes. XGI oracle in
# local_testing_and_equivalence/test-oracle-katz.R.

.kz_hg <- function(edges, sparse = FALSE) {
  long <- data.frame(
    member = unlist(edges),
    group = rep(names(edges), lengths(edges)),
    stringsAsFactors = FALSE
  )
  group_hypergraph(long, actor = "member", group = "group", sparse = sparse)
}

set.seed(23)
kz_edges <- lapply(seq_len(30), \(i) paste0("n", sample(14, sample(2:4, 1))))
names(kz_edges) <- sprintf("e%02d", seq_along(kz_edges))

test_that("one hyperedge of size n is K_n: closed form alpha (n-1) / (1 - alpha (n-1))", {
  hg <- .kz_hg(list(e1 = c("a", "b", "c", "d", "e")))
  out <- .hg_katz_fit(hg, alpha = 0.1)
  expect_identical(out$node, c("a", "b", "c", "d", "e"))
  expect_equal(out$katz, rep(0.1 * 4 / (1 - 0.1 * 4), 5), tolerance = 1e-14)
})

test_that("Katz equals the truncated walk series sum_k alpha^k A^k 1", {
  hg <- .kz_hg(kz_edges)
  membership <- (hg$incidence > 0) * 1
  adjacency <- tcrossprod(membership)
  diag(adjacency) <- 0
  lambda <- max(eigen(adjacency, symmetric = TRUE, only.values = TRUE)$values)
  alpha <- 0.5 / lambda
  # 200 terms: the remainder is below (1/2)^200 relative
  terms <- Reduce(\(v, k) alpha * as.numeric(adjacency %*% v),
                  seq_len(200), accumulate = TRUE,
                  init = rep(1, nrow(adjacency)))
  series <- Reduce(`+`, terms[-1L])
  out <- .hg_katz_fit(hg, alpha = alpha)
  expect_equal(out$katz, series, tolerance = 1e-12)
})

test_that("Katz + 1 is igraph's alpha centrality on the weighted clique expansion", {
  skip_if_not_installed("igraph")
  hg <- .kz_hg(kz_edges)
  # build the clique expansion from the edge list, not from our adjacency
  pairs <- do.call(rbind, lapply(kz_edges, \(e) t(utils::combn(e, 2))))
  g <- igraph::graph_from_edgelist(pairs, directed = FALSE)
  igraph::E(g)$weight <- 1
  g <- igraph::simplify(g, edge.attr.comb = list(weight = "sum"))
  alpha <- 0.02
  ref <- igraph::alpha_centrality(g, alpha = alpha, exo = 1,
                                  weights = igraph::E(g)$weight)
  out <- .hg_katz_fit(hg, alpha = alpha)
  expect_equal(out$katz + 1, unname(ref[out$node]), tolerance = 1e-10)
})

test_that("dense and sparse incidences agree", {
  dense <- .hg_katz_fit(.kz_hg(kz_edges), alpha = 0.03)
  sparse <- .hg_katz_fit(.kz_hg(kz_edges, sparse = TRUE), alpha = 0.03)
  expect_equal(sparse, dense, tolerance = 1e-12)
})

test_that("invariants: positivity, relabelling, monotone in alpha", {
  base <- .hg_katz_fit(.kz_hg(kz_edges), alpha = 0.03)
  expect_true(all(base$katz > 0))
  relabel <- setNames(paste0("z", sample(14)), paste0("n", seq_len(14)))
  moved <- lapply(rev(kz_edges), \(e) unname(relabel[e]))
  names(moved) <- paste0("f", seq_along(moved))
  shuffled <- .hg_katz_fit(.kz_hg(moved), alpha = 0.03)
  expect_equal(shuffled$katz[match(relabel[base$node], shuffled$node)],
               base$katz, tolerance = 1e-12)
  larger <- .hg_katz_fit(.kz_hg(kz_edges), alpha = 0.04)
  expect_true(all(larger$katz > base$katz))
})

test_that("alpha outside (0, 1 / lambda_max) raises hypernets_bad_input", {
  hg <- .kz_hg(list(e1 = c("a", "b", "c")))    # K_3: lambda_max = 2
  expect_error(.hg_katz_fit(hg, alpha = 0.5), class = "hypernets_bad_input")
  expect_error(.hg_katz_fit(hg, alpha = 0), class = "hypernets_bad_input")
  expect_error(.hg_katz_fit(hg, alpha = c(0.1, 0.2)),
               class = "hypernets_bad_input")
  expect_error(.hg_katz_fit(list(), alpha = 0.1),
               class = "hypernets_bad_input")
  expect_silent(.hg_katz_fit(hg, alpha = 0.49))
})

test_that("hg_centrality(type = 'katz') delegates to the engine", {
  skip_if_not("alpha" %in% names(formals(hg_centrality.net_hg)),
              "hg_centrality() katz patch not applied yet")
  hg <- .kz_hg(kz_edges)
  out <- hg_centrality(hg, type = c("clique", "katz"), alpha = 0.03)
  expect_identical(names(out), c("node", "clique", "katz"))
  expect_identical(out$katz, .hg_katz_fit(hg, alpha = 0.03)$katz)
  expect_identical(out$clique, hg_centrality(hg, type = "clique")$clique)
  sparse <- hg_centrality(.kz_hg(kz_edges, sparse = TRUE), type = "katz",
                          alpha = 0.03)
  expect_equal(sparse$katz, out$katz, tolerance = 1e-12)
  expect_error(hg_centrality(hg, type = "katz"),
               class = "hypernets_bad_input")
  expect_error(hg_centrality(.kz_hg(kz_edges, sparse = TRUE),
                             type = c("katz", "clique"), alpha = 0.03),
               class = "hypernets_sparse_unsupported")
  top <- hg_centrality(hg, type = "katz", alpha = 0.03, sort_by = "katz",
                       n = 3)
  expect_identical(top$katz, sort(out$katz, decreasing = TRUE)[1:3])
})
