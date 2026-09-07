# ---- hon_centrality() tests ----------------------------------------------

.hc_cycle <- function(n = 4L) {
  # a -> b -> c -> a, first order only
  build_hon(replicate(n, rep(c("a", "b", "c"), 5), simplify = FALSE),
            max_order = 1L)
}

.hc_planted <- function(n = 6L) {
  build_hon(c(replicate(n, rep(c("a", "b", "c"), 4), simplify = FALSE),
              replicate(n, rep(c("x", "b", "d"), 4), simplify = FALSE)),
            max_order = 2L)
}

test_that("directed 3-cycle matches the hand-computed centralities", {
  # Nodes a, b, c with edges a->b, b->c, c->a.
  #   PageRank: symmetric cycle -> uniform 1/3.
  #   Betweenness: each ordered pair has one shortest path; the only
  #     two-hop paths are a->b->c, b->c->a, c->a->b, so every node is an
  #     intermediate exactly once.
  #   Harmonic closeness (incoming): distances into a are 1 (from c) and
  #     2 (from b) -> 1 + 1/2 = 1.5.
  cen <- hon_centrality(.hc_cycle())
  expect_identical(cen$state, c("a", "b", "c"))
  expect_equal(cen$pagerank, rep(1 / 3, 3), tolerance = 1e-10)
  expect_equal(cen$betweenness, rep(1, 3), tolerance = 1e-12)
  expect_equal(cen$closeness, rep(1.5, 3), tolerance = 1e-12)
})

test_that("a first-order network projects to itself", {
  hon <- .hc_cycle()
  projected <- hon_centrality(hon, project = TRUE)
  raw <- hon_centrality(hon, project = FALSE)
  expect_identical(projected$state, raw$node)
  expect_identical(raw$order, rep(1L, 3L))
  expect_equal(projected$pagerank, raw$pagerank, tolerance = 1e-12)
  expect_equal(projected$betweenness, raw$betweenness, tolerance = 1e-12)
  expect_equal(projected$closeness, raw$closeness, tolerance = 1e-12)
})

test_that("projected PageRank conserves mass except under 'all'", {
  hon <- .hc_planted()
  results <- lapply(c("scaled", "last", "first"), function(pj) {
    cen <- hon_centrality(hon, type = "pagerank", projection = pj)
    expect_equal(sum(cen$pagerank), 1, tolerance = 1e-9)
    cen
  })
  expect_length(results, 3L)
  # 'all' hands the full value to every state on the path, so a network
  # with higher-order nodes must exceed 1
  all_pr <- hon_centrality(hon, type = "pagerank", projection = "all")
  expect_gt(sum(all_pr$pagerank), 1)
  # the three mass-conserving projections are not all identical here
  expect_false(isTRUE(all.equal(results[[1L]]$pagerank,
                                results[[2L]]$pagerank)))
})

test_that("projections place higher-order mass on the intended states", {
  hon <- .hc_planted()
  raw <- hon_centrality(hon, type = "pagerank", project = FALSE)
  ho <- subset(raw, order == 2L)
  expect_identical(sort(ho$node), c("a -> b", "x -> b"))
  last_pr <- hon_centrality(hon, type = "pagerank", projection = "last")
  first_pr <- hon_centrality(hon, type = "pagerank", projection = "first")
  # "a -> b" and "x -> b" both end in b, so 'last' concentrates their
  # mass on b while 'first' splits it between a and x
  expect_gt(subset(last_pr, state == "b")$pagerank,
            subset(first_pr, state == "b")$pagerank)
  expect_gt(subset(first_pr, state == "a")$pagerank,
            subset(last_pr, state == "a")$pagerank)
})

test_that("memory changes the ranking relative to the first-order network", {
  seqs <- c(replicate(8, rep(c("a", "b", "c"), 4), simplify = FALSE),
            replicate(8, rep(c("x", "b", "d"), 4), simplify = FALSE))
  hon_ho <- build_hon(seqs, max_order = 2L)
  ho <- hon_centrality(hon_ho)
  hon_fo <- build_hon(seqs, max_order = 1L)
  fo <- hon_centrality(hon_fo)
  expect_identical(ho$state, fo$state)
  # b is the shared hub: with memory its flow splits across two contexts
  expect_false(isTRUE(all.equal(ho$pagerank, fo$pagerank)))
  expect_lt(subset(ho, state == "b")$pagerank,
            subset(fo, state == "b")$pagerank)
})

test_that("centralities are invariant under state relabeling", {
  relabel <- c(a = "p", b = "q", c = "r")
  seqs <- replicate(6, rep(c("a", "b", "c"), 5), simplify = FALSE)
  hon1 <- build_hon(seqs, max_order = 2L)
  cen1 <- hon_centrality(hon1)
  hon2 <- build_hon(lapply(seqs, function(s) unname(relabel[s])),
                    max_order = 2L)
  cen2 <- hon_centrality(hon2)
  expect_identical(unname(relabel[cen1$state]), cen2$state)
  expect_equal(cen1$pagerank, cen2$pagerank, tolerance = 1e-12)
  expect_equal(cen1$betweenness, cen2$betweenness, tolerance = 1e-12)
  expect_equal(cen1$closeness, cen2$closeness, tolerance = 1e-12)
})

test_that("results are non-negative, finite and deterministically ordered", {
  set.seed(31)
  results <- lapply(1:8, function(i) {
    seqs <- replicate(10, sample(letters[1:5], sample(10:25, 1),
                                 replace = TRUE), simplify = FALSE)
    cen <- hon_centrality(build_hon(seqs, max_order = 3L))
    expect_identical(cen$state, sort(cen$state))
    expect_true(all(is.finite(cen$pagerank)))
    expect_true(all(cen$betweenness >= 0))
    expect_true(all(cen$closeness >= 0))
    expect_equal(sum(cen$pagerank), 1, tolerance = 1e-9)
    cen
  })
  expect_length(results, 8L)
})

test_that("type selection returns exactly the requested columns", {
  hon <- .hc_planted()
  pr_only <- hon_centrality(hon, type = "pagerank")
  expect_identical(names(pr_only), c("state", "pagerank"))
  cl_bt <- hon_centrality(hon, type = c("closeness", "betweenness"))
  expect_identical(names(cl_bt), c("state", "betweenness", "closeness"))
  pr_raw <- hon_centrality(hon, project = FALSE, type = "pagerank")
  expect_identical(names(pr_raw), c("node", "order", "pagerank"))
})

test_that("damping and weighting change PageRank as expected", {
  hon <- .hc_planted()
  d_low <- hon_centrality(hon, type = "pagerank", damping = 0.2)
  d_high <- hon_centrality(hon, type = "pagerank", damping = 0.95)
  expect_false(isTRUE(all.equal(d_low$pagerank, d_high$pagerank)))
  # as damping -> 0 the walk is pure teleportation, so the higher-order
  # PageRank tends to the uniform distribution over higher-order nodes
  # (the projection then splits it across each node's states)
  raw_tiny <- hon_centrality(hon, type = "pagerank", project = FALSE,
                             damping = 1e-6)
  expect_equal(raw_tiny$pagerank, rep(1 / nrow(raw_tiny), nrow(raw_tiny)),
               tolerance = 1e-5)
  wt <- hon_centrality(hon, type = "pagerank", weighted = TRUE)
  expect_equal(sum(wt$pagerank), 1, tolerance = 1e-9)
})

test_that("error paths and the path-enumeration guard", {
  hon <- .hc_planted()
  expect_error(hon_centrality(list(a = 1)), "net_hon")
  expect_error(hon_centrality(hon, damping = 0), "damping")
  expect_error(hon_centrality(hon, damping = 1), "damping")
  expect_error(hon_centrality(hon, tol = 0), "tol")
  expect_error(hon_centrality(hon, project = NA), "project")
  expect_error(hon_centrality(hon, type = "nope"), "arg")
  expect_error(hon_centrality(hon, type = "betweenness", max_paths = 1),
               class = "hypernets_too_many_paths")
  expect_warning(
    hon_centrality(hon, type = "pagerank", max_iter = 1L),
    class = "hypernets_no_converge")
})

test_that("build_hon accepts long format identically to a list", {
  seqs <- list(s1 = rep(c("a", "b", "c"), 4), s2 = rep(c("a", "b", "d"), 4),
               s3 = rep(c("x", "b", "c"), 4))
  long <- data.frame(
    code = unlist(seqs, use.names = FALSE),
    id   = rep(names(seqs), times = lengths(seqs)),
    t    = unlist(lapply(lengths(seqs), seq_len), use.names = FALSE),
    stringsAsFactors = FALSE
  )
  h_list <- build_hon(seqs, max_order = 2L)
  h_long <- build_hon(long, action = "code", actor = "id", time = "t",
                      max_order = 2L)
  expect_identical(h_list$matrix, h_long$matrix)
  expect_identical(h_list$ho_edges, h_long$ho_edges)
  expect_identical(h_list$first_order_states, h_long$first_order_states)
  expect_error(build_hon(seqs, actor = "id"), "actor")
})

test_that("first-order distances are the minimum over all realizations", {
  # Regression guard: several higher-order node pairs can realize the
  # same (state, state) cell, so the projection must take the minimum
  # rather than whichever candidate is written last.
  set.seed(77)
  results <- lapply(1:6, function(i) {
    seqs <- replicate(10, sample(letters[1:4], sample(12:25, 1),
                                 replace = TRUE), simplify = FALSE)
    hon <- build_hon(seqs, max_order = 3L)
    nodes <- rownames(hon$matrix)
    node_paths <- hypernets:::.hoc_node_paths(nodes)
    succ <- hypernets:::.hoc_succ(hon$matrix)
    n <- length(nodes)
    dist_ho <- do.call(rbind, lapply(seq_len(n), function(s) {
      hypernets:::.hoc_bfs_dist(succ, s, n)
    }))
    states <- hon$first_order_states
    got <- hypernets:::.hoc_dist_first(dist_ho, node_paths, states)
    # independent slow reference: minimum over every (v, w) realization
    ref <- matrix(Inf, length(states), length(states),
                  dimnames = list(states, states))
    for (v in seq_len(n)) {
      for (w in seq_len(n)) {
        if (!is.finite(dist_ho[v, w])) next
        s1 <- node_paths[[v]][1L]
        d1 <- node_paths[[w]][length(node_paths[[w]])]
        cand <- dist_ho[v, w] + length(node_paths[[v]]) - 1L
        if (cand < ref[s1, d1]) ref[s1, d1] <- cand
      }
    }
    diag(ref) <- 0
    expect_identical(got, ref)
    got
  })
  expect_length(results, 6L)
})

test_that("sort_by orders by a computed centrality, largest first", {
  hon <- .hc_planted()
  by_pr <- hon_centrality(hon, type = "pagerank", project = FALSE,
                          sort_by = "pagerank")
  expect_true(all(diff(by_pr$pagerank) <= 1e-12))
  by_cl <- hon_centrality(hon, type = c("pagerank", "closeness"),
                          sort_by = "closeness")
  expect_true(all(diff(by_cl$closeness) <= 1e-12))
  # only a requested type may be used as the sort key
  expect_error(hon_centrality(hon, type = "pagerank", sort_by = "closeness"),
               "arg")
})
