testthat::skip_on_cran()

# ---- One verb per idea: dispatch on the input class ----
#
# hg_centrality(), hg_communities() and hypa() pick their estimator from
# the class of `x`. Each method keeps its own arguments and defaults; an
# argument only the other method takes raises hypernets_bad_input.

.dp_seqs <- function() {
  list(c("a", "h", "b", "a", "h", "b", "a"),
       c("c", "h", "d", "c", "h", "d", "c"),
       c("a", "h", "b", "a", "h", "b"),
       c("c", "h", "d", "c", "h", "d"))
}

.dp_hg <- function() {
  group_hypergraph(data.frame(
    member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z",
               "a", "x"),
    edge = c(rep(paste0("e", 1:4), each = 3), "e5", "e5")
  ), "member", "edge")
}

test_that("hg_centrality dispatches on memory networks and hypergraphs", {
  h <- hon(.dp_seqs(), max_order = 2L)
  mem <- hg_centrality(h)
  expect_named(mem, c("state", "pagerank", "betweenness", "closeness"))
  expect_identical(hg_centrality(h, type = "pagerank", project = FALSE),
                   hg_centrality.net_hon(h, type = "pagerank",
                                         project = FALSE))
  hg <- .dp_hg()
  expect_named(hg_centrality(hg, type = "clique"), c("node", "clique"))
})

test_that("hg_centrality refuses an argument of the other method", {
  h <- hon(.dp_seqs(), max_order = 2L)
  expect_error(hg_centrality(h, alpha = 0.1), class = "hypernets_bad_input")
  expect_error(hg_centrality(h, n = 3), class = "hypernets_bad_input")
  expect_error(hg_centrality(.dp_hg(), project = FALSE),
               class = "hypernets_bad_input")
  expect_error(hg_centrality(.dp_hg(), top = 2),
               class = "hypernets_bad_input")
  expect_error(hg_centrality(data.frame(a = 1)), class = "hypernets_bad_input")
})

test_that("hg_communities dispatches on memory networks and hypergraphs", {
  h <- hon(.dp_seqs(), max_order = 2L)
  comm <- hg_communities(h, trials = 2L)
  expect_s3_class(comm, "net_hon_communities")
  expect_identical(comm, hg_communities.net_hon(h, trials = 2L))
  skip_if_not_installed("igraph")
  fit <- hg_communities(.dp_hg(), n_runs = 2L, trials = 2L, seeds = 1:2)
  expect_s3_class(fit, "hg_communities")
})

test_that("hg_communities refuses an argument of the other method", {
  h <- hon(.dp_seqs(), max_order = 2L)
  expect_error(hg_communities(h, n_runs = 2L), class = "hypernets_bad_input")
  expect_error(hg_communities(h, type = "irmm"), class = "hypernets_bad_input")
  expect_error(hg_communities(.dp_hg(), teleportation = 0.1),
               class = "hypernets_bad_input")
  expect_error(hg_communities(.dp_hg(), partition = NULL),
               class = "hypernets_bad_input")
  expect_error(hg_communities(matrix(1)), class = "hypernets_bad_input")
})

test_that("hypa dispatches on hypergraphs and on sequences", {
  pairs <- hypa(.dp_hg(), min_count = 1L)
  expect_s3_class(pairs, "data.frame")
  expect_true(all(c("from", "to", "observed", "expected") %in% names(pairs)))
  paths <- hypa(.dp_seqs(), min_count = 1L)
  expect_s3_class(paths, "net_hypa")
  expect_identical(.unview(paths), Nestimate::build_hypa(.dp_seqs(), min_count = 1L))
})

test_that("hypa refuses an argument of the other method", {
  expect_error(hypa(.dp_hg(), order = 2L), class = "hypernets_bad_input")
  expect_error(hypa(.dp_hg(), action = "a"), class = "hypernets_bad_input")
  expect_error(hypa(.dp_seqs(), top = 2L), class = "hypernets_bad_input")
})

test_that("the dispatching verbs keep each method's defaults", {
  expect_identical(formals(hg_centrality.net_hg)$tol, 1e-8)
  expect_identical(formals(hg_centrality.net_hon)$tol, 1e-12)
  expect_identical(formals(hg_communities.net_hg)$trials, 100L)
  expect_identical(formals(hg_communities.net_hon)$trials, 10L)
  expect_identical(formals(hypa.net_hg)$min_count, 5L)
  expect_identical(formals(hypa.default)$min_count, 5L)
})
