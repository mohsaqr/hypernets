test_that("G(n,p) generator follows Bernoulli incidence model", {
  h <- hg_sample_gnp(12, m = 7, p = 0.3, seed = 17)
  expect_s3_class(h, "net_hypergraph")
  expect_equal(dim(h$incidence), c(12, 7))
  expect_true(all(h$incidence %in% 0:1))
  expect_identical(h, hg_sample_gnp(12, m = 7, p = 0.3, seed = 17))
  expect_identical(hypergraph_sample_gnp, hg_sample_gnp)
})

test_that("G(n,p) keeps valid limiting and Poisson cases", {
  expect_true(all(hg_sample_gnp(5, m = 3, p = 0, seed = 1)$incidence == 0))
  expect_true(all(hg_sample_gnp(5, m = 3, p = 1, seed = 1)$incidence == 1))
  h <- hg_sample_gnp(5, p = 0.2, lambda = 0, seed = 1)
  expect_equal(h$n_hyperedges, 0L)
})

test_that("uniform generator fixes every hyperedge size", {
  h <- hg_sample_uniform(20, 30, k = 4, seed = 4)
  expect_equal(unname(colSums(h$incidence)), rep(4, 30))
  expect_equal(h$size_distribution, c(size_4 = 30L))
  expect_identical(hypergraph_sample_uniform, hg_sample_uniform)
})

test_that("regular generator fixes every node degree", {
  h <- hg_sample_regular(20, 9, k = 3, seed = 4)
  expect_equal(unname(rowSums(h$incidence)), rep(3, 20))
  expect_identical(hypergraph_sample_regular, hg_sample_regular)
})

test_that("SBM generator records planted blocks and fixed sizes", {
  P <- matrix(c(0.8, 0.05, 0.05, 0.8), 2, 2)
  h <- hg_sample_sbm(P = P, block_sizes = c(8, 8), d = 3, seed = 8)
  expect_s3_class(h, "net_hypergraph")
  expect_equal(as.integer(table(h$blocks)), c(8L, 8L))
  expect_true(all(colSums(h$incidence) == 3))
  expect_identical(h, hg_sample_sbm(P = P, block_sizes = c(8, 8),
                                    d = 3, seed = 8))
  expect_identical(hypergraph_sample_sbm, hg_sample_sbm)
})

test_that("SBM variable sizes are at least dyadic", {
  P <- matrix(1, 2, 2)
  h <- hg_sample_sbm(P = P, block_sizes = c(10, 10), d = 1,
                     variable_size = TRUE, seed = 2)
  expect_true(all(colSums(h$incidence) >= 2))
})

test_that("random generators restore caller RNG state", {
  set.seed(101)
  before <- .Random.seed
  hg_sample_gnp(10, 4, p = .2, seed = 9)
  expect_identical(.Random.seed, before)
  hg_sample_uniform(10, 4, 2, seed = 9)
  expect_identical(.Random.seed, before)
  hg_sample_regular(10, 4, 2, seed = 9)
  expect_identical(.Random.seed, before)
})

test_that("random generators reject impossible contracts", {
  expect_error(hg_sample_gnp(4, 2, p = 2), "probability")
  expect_error(hg_sample_uniform(4, 2, k = 5), "cannot exceed")
  expect_error(hg_sample_regular(4, 2, k = 3), "cannot exceed")
  expect_error(hg_sample_sbm(P = diag(2), block_sizes = c(2, 2), d = 3),
               "not enough nodes")
})
