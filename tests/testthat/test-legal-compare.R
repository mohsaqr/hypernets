testthat::skip_on_cran()

.compare_fixture <- function() {
  hg <- group_hypergraph(
    data.frame(
      member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
      edge = rep(paste0("e", 1:4), each = 3)
    ),
    "member", "edge"
  )
  hg$edge_data <- data.frame(edge = paste0("e", 1:4),
                             source = c("s1", "s1", "s2", "s2"),
                             stringsAsFactors = FALSE)
  hg
}

test_that("Infomap runs on the citation projection, directed or not", {
  skip_if_not_installed("igraph")
  hg <- .compare_fixture()
  undirected <- hg_communities(hg, n_runs = 2, trials = 2, seeds = 1:2,
                               method = "citation")
  expect_identical(undirected$params$method, "citation")
  expect_false(undirected$params$directed)
  expect_setequal(undirected$medoid$node, c("a", "b", "c", "x", "y", "z", "s1", "s2"))
  directed <- hg_communities(hg, n_runs = 2, trials = 2, seeds = 1:2,
                             method = "citation", directed = TRUE)
  expect_true(directed$params$directed)
  expect_false(isSymmetric(unname(directed$projection)))
  expect_error(hg_communities(hg, n_runs = 2, trials = 2, directed = TRUE),
               class = "hypernets_bad_input")
})

test_that("hg_compare_communities tabulates sizes, similarity and summaries", {
  skip_if_not_installed("igraph")
  hg <- .compare_fixture()
  mh <- hg_communities(hg, n_runs = 2, trials = 2, seeds = 1:2)
  bh <- hg_communities(hg, n_runs = 2, trials = 2, seeds = 1:2,
                       duplicate_edges = "collapse")
  cmp <- hg_compare_communities(mh = mh, bh = bh)
  expect_s3_class(cmp, "hypernets_community_comparison")
  expect_identical(cmp$models, c("mh", "bh"))
  s <- as.data.frame(cmp)
  expect_identical(names(s), c("model", "medoid_seed", "n_runs", "n_communities",
                               "n_singletons", "n_nontrivial", "largest",
                               "second", "balance"))
  expect_identical(s$n_communities, c(2L, 2L))
  expect_equal(s$balance, c(1, 1))
  sim <- as.data.frame(cmp, what = "similarity")
  expect_identical(names(sim), c("model_a", "model_b", "n_nodes", "ami", "ari", "nmi"))
  expect_equal(sim$ari, 1)
  expect_equal(sim$ami, 1)
  sizes <- as.data.frame(cmp, what = "sizes")
  expect_identical(sizes$n_nodes, c(3L, 3L, 3L, 3L))
  cmp_summary <- summary(cmp)
  expect_identical(cmp_summary, s)
  expect_output(print(cmp), "2 representations")
  sizes_plot <- plot(cmp)
  expect_s3_class(sizes_plot, "ggplot")
  similarity_plot <- plot(cmp, what = "similarity")
  expect_s3_class(similarity_plot, "ggplot")
  m <- as.data.frame(cmp, what = "matrix")
  expect_identical(dimnames(m), list(c("mh", "bh"), c("mh", "bh")))
  expect_equal(m["bh", "mh"], sim$ami)
  expect_equal(m["mh", "bh"], sim$ari)
  expect_true(all(is.na(diag(m))))
  # a list works too; one fit or unnamed fits do not
  from_list <- hg_compare_communities(list(mh = mh, bh = bh))
  expect_identical(from_list$summary, cmp$summary)
  expect_error(hg_compare_communities(mh = mh), class = "hypernets_bad_input")
  expect_error(hg_compare_communities(mh, bh), class = "hypernets_bad_input")
  expect_error(as.data.frame(cmp, what = "quality"), class = "hypernets_bad_input")
})

test_that("quality scores each medoid on its own projection", {
  skip_if_not_installed("igraph")
  hg <- .compare_fixture()
  mh <- hg_communities(hg, n_runs = 2, trials = 2, seeds = 1:2)
  mgu <- hg_communities(hg, n_runs = 2, trials = 2, seeds = 1:2,
                        method = "citation")
  cmp <- hg_compare_communities(mh = mh, mgu = mgu, hg = hg)
  q <- as.data.frame(cmp, what = "quality")
  expect_identical(q$model, c("mh", "mgu"))
  expect_identical(names(q)[1:3], c("model", "coverage", "weighted_coverage"))
  # two disjoint triangles: the association medoid is the perfect cut
  expect_equal(q$coverage[[1L]], 1)
  expect_equal(q$performance[[1L]], 1)
  direct <- hg_community_quality(hg, mgu, method = "citation")
  expect_equal(q$modularity[[2L]], direct$modularity)
})
