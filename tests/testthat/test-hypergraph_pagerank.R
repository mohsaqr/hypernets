# EDVW hypergraph PageRank: independent linear-solve reference, the
# Chitra & Raphael collapse theorem, hypergraph_centrality parity, and
# invariants.

pr_corpus <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
)
pr_stop <- c("the", "with", "and", "a", "this", "at", "on", "all")

# Independent reference: build the EDVW transition entry by entry from the
# published formula, then solve the PageRank linear system directly
# (pi = (1 - d) u (I - d P)^{-1}) -- a different algorithm and code path
# from the package's power iteration.
.reference_transition <- function(hg, edge_weights) {
  R <- hg$incidence
  ids <- rownames(R)
  delta <- colSums(R)
  entry <- function(v, w) {
    e_of_v <- which(R[v, ] > 0)
    d_v <- sum(edge_weights[e_of_v])
    sum(vapply(e_of_v, \(e) edge_weights[e] / d_v * R[w, e] / delta[e],
               numeric(1)))
  }
  outer(ids, ids, Vectorize(entry))
}

test_that("pagerank matches a direct linear solve of the published formula", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop, weight = "tfidf")
  omega <- rep(1, hg$n_hyperedges)
  p_ref <- .reference_transition(hg, omega)
  d <- 0.85
  u <- rep(1 / hg$n_nodes, hg$n_nodes)
  pi_ref <- as.numeric((1 - d) * u %*% solve(diag(hg$n_nodes) - d * p_ref))
  pi_ref <- pi_ref / sum(pi_ref)

  out <- hg_pagerank(hg, damping = d, edge_weights = omega, tol = 1e-14)
  expect_identical(out$node, hg$nodes)
  expect_equal(out$pagerank, pi_ref, tolerance = 1e-10)
})

test_that("damping = 1 equals the spectral engine's stationary distribution", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop, weight = "tfidf")
  out <- hg_pagerank(hg, damping = 1, tol = 1e-15, max_iter = 10000L)
  engine <- hypergraph_cluster(hg, k = 2, type = "random_walk",
                                          seed = 1)
  expect_equal(
    out$pagerank,
    unname(engine$pi[out$node]),
    tolerance = 1e-12
  )
})

test_that("edge-independent weights collapse to the graph walk (Chitra & Raphael)", {
  # gamma_e(v) = gamma(v) for every edge: each vertex carries one weight
  # everywhere it appears. Theorem: the hypergraph walk equals the random
  # walk on the clique-expansion graph with w(u, v) = sum_e omega_e *
  # gamma(u) gamma(v) / delta_e over shared edges (self-loops included), so
  # its stationary distribution is the normalized node strength.
  gamma <- c(a = 1, b = 2, c = 3, d = 1.5)
  long <- data.frame(
    vertex = c("a", "b", "c", "b", "c", "d", "a", "d"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
    w = c(1, 2, 3, 2, 3, 1.5, 1, 1.5)
  )
  hg <- group_hypergraph(long, actor = "vertex", group = "edge",
                                    weight = "w")
  omega <- c(1, 2, 0.5)

  R <- hg$incidence
  delta <- colSums(R)
  # w(u, v) = sum over shared edges of omega_e * gamma(u) * gamma(v) /
  # delta_e; with gamma_e(v) = gamma(v), R[v, e] IS gamma(v) on membership.
  strength_of <- function(u) {
    sum(vapply(seq_len(nrow(R)), \(v) {
      shared <- which(R[u, ] > 0 & R[v, ] > 0)
      sum(omega[shared] * R[u, shared] * R[v, shared] / delta[shared])
    }, numeric(1)))
  }
  strengths <- vapply(seq_len(nrow(R)), strength_of, numeric(1))
  pi_graph <- strengths / sum(strengths)

  out <- hg_pagerank(hg, damping = 1, edge_weights = omega, tol = 1e-15,
                     max_iter = 10000L)
  expect_equal(out$pagerank, pi_graph, tolerance = 1e-12)
})

test_that("pagerank is a probability vector and permutation-invariant", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop, weight = "tfidf")
  out <- hg_pagerank(hg)
  expect_equal(sum(out$pagerank), 1)
  expect_true(all(out$pagerank > 0))

  shuffled <- text_hypergraph(rev(pr_corpus), stop_words = pr_stop,
                              weight = "tfidf")
  out_shuffled <- hg_pagerank(shuffled)
  merged <- merge(out, out_shuffled, by = "node")
  expect_equal(merged$pagerank.x, merged$pagerank.y, tolerance = 1e-12)
})

test_that("personalization concentrates teleport mass", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop)
  uniform <- hg_pagerank(hg)
  focused <- hg_pagerank(hg, personalized = c(cooking_1 = 1))
  merged <- merge(uniform, focused, by = "node", suffixes = c("_u", "_p"))
  target <- subset(merged, node == "cooking_1")
  expect_gt(target$pagerank_p, target$pagerank_u)
  # a character vector of node names is uniform teleport over those nodes
  by_name <- hg_pagerank(hg, personalized = "cooking_1")
  expect_equal(by_name$pagerank, focused$pagerank, tolerance = 1e-12)
  expect_error(hg_pagerank(hg, personalized = "no_such_node"),
               class = "hypernets_bad_input")
})

test_that("sort_by and n select without user-side subsetting", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop)
  top <- hg_pagerank(hg, sort_by = "pagerank", n = 2)
  expect_identical(nrow(top), 2L)
  expect_true(all(diff(top$pagerank) <= 0))
})

test_that("contract violations and non-convergence are classed", {
  hg <- text_hypergraph(pr_corpus, stop_words = pr_stop)
  expect_error(hg_pagerank(hg, damping = 0))
  expect_error(hg_pagerank(hg, damping = 1.2))
  expect_error(hg_pagerank(hg, personalized = c(nope = 1)),
               class = "hypernets_bad_input")
  expect_error(hg_pagerank(hg, personalized = c(0.5)),
               class = "hypernets_bad_input")
  expect_error(hg_pagerank(hg, edge_weights = c(1, 2)))
  expect_error(hg_pagerank(42), class = "hypernets_bad_input")
  expect_warning(hg_pagerank(hg, max_iter = 1L), class = "hypernets_no_converge")
})

test_that("hg_pagerank equals hypergraph_centrality(type = \"pagerank\")", {
  # same walk, same damping, uniform teleport: the verb and the engine must
  # agree to solver tolerance on binary and on weighted incidences
  long <- data.frame(
    vertex = c("a", "b", "c", "b", "c", "d", "a", "d", "e", "c", "e"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e3", "e4", "e4"),
    w = c(2, 1, 3, 1, 4, 2, 1, 1, 5, 1, 2)
  )
  binary <- group_hypergraph(long, actor = "vertex", group = "edge")
  weighted <- group_hypergraph(long, actor = "vertex", group = "edge",
                               weight = "w")
  for (hg in list(binary, weighted)) {
    verb <- hg_pagerank(hg, damping = 0.85, tol = 1e-14)
    engine <- hypergraph_centrality(hg, type = "pagerank", damping = 0.85,
                                    tol = 1e-14)
    expect_identical(verb$node, engine$node)
    expect_equal(verb$pagerank, engine$pagerank, tolerance = 1e-10)
  }
})
