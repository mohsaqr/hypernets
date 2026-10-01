# hg_modularity() (Kaminski et al. 2019) and IRMM (Kumar et al. 2020).
# The HyperNetX oracle lives in
# local_testing_and_equivalence/test-equiv-hg-modularity-{score,irmm}.R.

.hgm_from_list <- function(edges) {
  long <- data.frame(member = unlist(edges),
                     edge = rep(names(edges), lengths(edges)))
  group_hypergraph(long, "member", "edge")
}

.hgm_planted <- function(seed, n_blocks = 3L, block_size = 8L, n_within = 10L,
                         n_cross = 3L) {
  set.seed(seed)
  nodes <- sprintf("v%02d", seq_len(n_blocks * block_size))
  block <- rep(seq_len(n_blocks), each = block_size)
  within <- unlist(lapply(seq_len(n_blocks), function(b) {
    lapply(seq_len(n_within), function(i) sample(nodes[block == b], sample(2:5, 1L)))
  }), recursive = FALSE)
  cross <- lapply(seq_len(n_cross), function(i) sample(nodes, sample(2:4, 1L)))
  members <- c(within, cross)
  names(members) <- sprintf("e%03d", seq_along(members))
  list(hg = .hgm_from_list(members),
       truth = data.frame(node = nodes, community = block))
}

test_that("hand-computed modularity: two disjoint pairs and one split triple", {
  skip_on_cran()
  h <- .hgm_from_list(list(e1 = c("a", "b"), e2 = c("c", "d")))
  p <- data.frame(node = c("a", "b", "c", "d"), community = c(1, 1, 2, 2))
  # EC = 2/2; each community has vol share 1/2, DT = 2 * 2 * (1/2)^2 / 2
  vals <- vapply(c("linear", "majority", "strict"), function(t) {
    hg_modularity(h, p, type = t)$modularity
  }, numeric(1L))
  expect_equal(unname(vals), rep(0.5, 3), tolerance = 1e-14)

  h3 <- .hgm_from_list(list(e1 = c("a", "b", "c")))
  p3 <- data.frame(node = c("a", "b", "c"), community = c("x", "x", "y"))
  # linear: EC = 2/3; DT = 2/3 (4/9 + 2/9) + (8/27 + 1/27) = 7/9
  expect_equal(hg_modularity(h3, p3)$modularity, 2 / 3 - 7 / 9,
               tolerance = 1e-14)
  # strict: EC = 0; DT = (2/3)^3 + (1/3)^3 = 1/3
  expect_equal(hg_modularity(h3, p3, type = "strict")$modularity, -1 / 3,
               tolerance = 1e-14)
  # majority: EC = 1; DT = 4/9 + 2/9 + 8/27 + 1/27 = 1
  expect_equal(hg_modularity(h3, p3, type = "majority")$modularity, 0,
               tolerance = 1e-14)
})

test_that("strict modularity equals the Kaminski et al. (2019) closed form", {
  skip_on_cran()
  fx <- .hgm_planted(3L)
  h <- fx$hg
  set.seed(4L)
  labels <- sample(3L, nrow(fx$truth), replace = TRUE)
  p <- data.frame(node = fx$truth$node, community = labels)
  inc <- (h$incidence > 0) * 1
  lab <- labels[match(rownames(inc), fx$truth$node)]
  sizes <- colSums(inc)
  inside <- vapply(seq_len(ncol(inc)), function(j) {
    length(unique(lab[inc[, j] > 0])) == 1L
  }, logical(1L))
  vol <- tapply(rowSums(inc), lab, sum) / sum(inc)
  tax <- sum(vapply(sizes, function(d) sum(vol^d), numeric(1L)))
  closed <- (sum(inside) - tax) / ncol(inc)
  expect_equal(hg_modularity(h, p, type = "strict")$modularity, closed,
               tolerance = 1e-12)
})

test_that("invariants: relabelling, node order, weight scale, decomposition", {
  skip_on_cran()
  fx <- .hgm_planted(8L)
  h <- fx$hg
  set.seed(9L)
  p <- data.frame(node = fx$truth$node,
                  community = sample(4L, nrow(fx$truth), replace = TRUE))
  w <- round(stats::runif(h$n_hyperedges, 0.2, 4), 3)
  types <- c("linear", "majority", "strict")
  base <- vapply(types, function(t) {
    hg_modularity(h, p, type = t, edge_weights = w)$modularity
  }, numeric(1L))
  relabelled <- transform(p, community = c("q", "r", "s", "t")[community])
  shuffled <- relabelled[sample(nrow(relabelled)), ]
  moved <- vapply(types, function(t) {
    hg_modularity(h, shuffled, type = t, edge_weights = w)$modularity
  }, numeric(1L))
  scaled <- vapply(types, function(t) {
    hg_modularity(h, p, type = t, edge_weights = 7.5 * w)$modularity
  }, numeric(1L))
  expect_equal(moved, base, tolerance = 1e-12)
  expect_equal(scaled, base, tolerance = 1e-12)

  per <- hg_modularity(h, p, edge_weights = w, what = "communities")
  expect_named(per, c("community", "n_nodes", "volume", "edge_contribution",
                      "degree_tax", "modularity"))
  expect_equal(sum(per$modularity), base[["linear"]], tolerance = 1e-12)
  expect_equal(sum(per$n_nodes), nrow(fx$truth))
  expect_false(is.unsorted(rev(per$modularity)))

  # one community owns every hyperedge and pays the full tax: q = 0
  one <- data.frame(node = fx$truth$node, community = 1L)
  zero <- vapply(types, function(t) hg_modularity(h, one, type = t)$modularity,
                 numeric(1L))
  expect_equal(unname(zero), rep(0, 3), tolerance = 1e-12)

  # the planted partition beats a random one
  expect_gt(hg_modularity(h, fx$truth)$modularity, base[["linear"]])
})

test_that("score table columns; sparse storage and input forms agree", {
  skip_on_cran()
  fx <- .hgm_planted(12L)
  h <- fx$hg
  out <- hg_modularity(h, fx$truth)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("type", "modularity", "edge_contribution", "degree_tax",
                      "n_communities", "n_edges", "total_weight"))
  expect_equal(out$modularity, out$edge_contribution - out$degree_tax,
               tolerance = 1e-14)
  expect_equal(out$n_communities, 3L)

  hs <- h
  hs$incidence <- Matrix::Matrix(h$incidence, sparse = TRUE)
  expect_equal(hg_modularity(hs, fx$truth), out, tolerance = 1e-14)

  named <- stats::setNames(fx$truth$community, fx$truth$node)
  expect_equal(hg_modularity(h, named), out)
  in_order <- fx$truth$community[match(h$nodes, fx$truth$node)]
  expect_equal(hg_modularity(h, in_order), out)
})

test_that("hg_modularity raises classed errors on bad input", {
  skip_on_cran()
  h <- .hgm_from_list(list(e1 = c("a", "b"), e2 = c("b", "c")))
  expect_error(hg_modularity(list(), c(a = 1)), class = "hypernets_bad_input")
  expect_error(hg_modularity(h, c(a = 1, b = 1)), class = "hypernets_bad_input")
  expect_error(hg_modularity(h, c(a = 1, b = NA, c = 2)),
               class = "hypernets_bad_input")
  expect_error(hg_modularity(h, c(a = 1, b = 1, c = 2), type = "nope"))
  expect_error(hg_modularity(h, c(a = 1, b = 1, c = 2), edge_weights = c(1, -1)))
})

test_that("IRMM reweighting follows Kumar et al. (2020) Eq. 6", {
  skip_on_cran()
  h <- .hgm_from_list(list(e1 = c("a", "b", "c", "d")))
  parts <- .hg_mod_parts(h)
  # 3:1 split, c = 2, m = 1: (1/4 + 1/2) * (4 + 2) = 4.5
  expect_equal(.hg_irmm_reweight(parts, c(1, 1, 1, 2)), 4.5)
  # 2:2 split: (1/3 + 1/3) * 6 = 4, the balanced minimum
  expect_equal(.hg_irmm_reweight(parts, c(1, 1, 2, 2)), 4)
  # uncut, c = 1: (1/5) * 5 = 1
  expect_equal(.hg_irmm_reweight(parts, c(1, 1, 1, 1)), 1)
  # the random-walk two-section preserves weighted degree
  h2 <- .hgm_planted(2L)$hg
  p2 <- .hg_mod_parts(h2)
  a <- .hg_irmm_two_section(p2, p2$weights)
  multi <- p2$sizes > 1
  expect_equal(unname(rowSums(a)),
               as.numeric((h2$incidence > 0)[, multi] %*% p2$weights[multi]))
  expect_true(isSymmetric(a))
})

test_that("IRMM recovers planted communities across seeds and restores the RNG", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  fx <- .hgm_planted(21L)
  set.seed(99L)
  before <- .Random.seed
  fit <- .hg_irmm_communities(fx$hg, n_runs = 5L, seeds = c(1L, 7L, 13L, 42L, 2026L))
  expect_identical(.Random.seed, before)
  expect_s3_class(fit, "hg_communities")
  runs <- hg_get(fit, what = "runs")
  expect_named(runs, c("run", "seed", "n_communities", "iterations",
                       "converged", "max_weight_change", "modularity"))
  expect_true(all(runs$converged))
  truth <- fx$truth$community[match(fit$medoid$node, fx$truth$node)]
  aris <- vapply(split(fit$partitions, fit$partitions$run), function(r) {
    .thg_ari(r$community[match(fit$medoid$node, r$node)], truth)
  }, numeric(1L))
  expect_equal(unname(aris), rep(1, 5), tolerance = 1e-12)
  expect_equal(unname(diag(fit$similarity$ami)), rep(1, 5))
  expect_equal(nrow(fit$weights), fx$hg$n_hyperedges)
  expect_equal(fit$projection, .hg_irmm_two_section(.hg_mod_parts(fx$hg),
                                                     fit$weights$weight))
  # the medoid's modularity is the hg_modularity() of its partition
  expect_equal(runs$modularity[fit$medoid_run],
               hg_modularity(fx$hg, fit$medoid)$modularity, tolerance = 1e-12)

  again <- .hg_irmm_communities(fx$hg, n_runs = 2L, seeds = c(1L, 7L))
  expect_identical(again$partitions,
                   subset(fit$partitions, fit$partitions$run <= 2L))

  # sparse incidence runs the same fit
  hs <- fx$hg
  hs$incidence <- Matrix::Matrix(fx$hg$incidence, sparse = TRUE)
  sparse_fit <- .hg_irmm_communities(hs, n_runs = 2L, seeds = c(1L, 7L))
  expect_equal(sparse_fit$partitions, again$partitions)
  expect_equal(sparse_fit$weights, again$weights, tolerance = 1e-12)
  expect_equal(sparse_fit$projection, again$projection, tolerance = 1e-12)
})

test_that("IRMM surfaces non-convergence and rejects bad controls", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  fx <- .hgm_planted(5L)
  expect_warning(
    fit <- .hg_irmm_communities(fx$hg, n_runs = 2L, delta = 1e-12, max_iter = 1L),
    class = "hypernets_no_converge"
  )
  expect_false(any(hg_get(fit, what = "runs")$converged))
  expect_error(.hg_irmm_communities(fx$hg, n_runs = 1L, delta = 0),
               class = "hypernets_bad_input")
  expect_error(.hg_irmm_communities(fx$hg, n_runs = 2L, seeds = 1L),
               class = "hypernets_bad_input")
  expect_error(.hg_irmm_communities(list(), n_runs = 1L),
               class = "hypernets_bad_input")
})

test_that("hg_communities(type = \"irmm\") dispatches to the IRMM engine", {
  skip_on_cran()
  skip_if_not_installed("igraph")
  skip_if_not("type" %in% names(formals(hg_communities.net_hg)),
              "hg_communities() has no `type =` yet (patch not applied)")
  fx <- .hgm_planted(21L)
  direct <- .hg_irmm_communities(fx$hg, n_runs = 3L, seeds = 1:3)
  via <- hg_communities(fx$hg, type = "irmm", n_runs = 3L, seeds = 1:3)
  expect_identical(via, direct)
})
