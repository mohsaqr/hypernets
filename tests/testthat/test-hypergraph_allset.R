.allset_fixture <- function() {
  events <- data.frame(
    node = c("a", "b", "c", "b", "c", "d"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e2")
  )
  group_hypergraph(events, actor = "node", group = "edge")
}

test_that("AllSet memberships retain the two paper multisets", {
  hg <- .allset_fixture()
  sets <- .thg_allset_memberships(hg)
  expect_identical(sets$vertices_by_edge, list(1:3, 2:4))
  expect_identical(sets$edges_by_vertex,
                   list(1L, c(1L, 2L), c(1L, 2L), 2L))
})

test_that("AllDeepSets Eq. 7 is permutation invariant", {
  skip_if_not_installed("torch")
  torch::torch_manual_seed(1)
  f <- .thg_deepset(3L, 4L)
  x <- torch::torch_tensor(matrix(1:9, 3, 3),
                           dtype = torch::torch_float())
  one_set <- .thg_torch_sparse(Matrix::Matrix(matrix(1, 1, 3), sparse = TRUE))
  reverse_set <- .thg_torch_sparse(Matrix::Matrix(matrix(1, 1, 3),
                                                   sparse = TRUE))
  reverse_idx <- torch::torch_tensor(c(3L, 2L, 1L),
                                     dtype = torch::torch_int64())
  expect_equal(as.array(f(x, one_set)),
               as.array(f(x[reverse_idx, ], reverse_set)), tolerance = 1e-6)
})

test_that("AllSetTransformer Eq. 8 is permutation invariant", {
  skip_if_not_installed("torch")
  torch::torch_manual_seed(2)
  f <- .thg_pma(4L, 2L)
  x <- torch::torch_tensor(matrix(rnorm(12), 3, 4),
                           dtype = torch::torch_float())
  a <- f(x, list(1:3))
  b <- f(x, list(3:1))
  expect_equal(as.array(a), as.array(b), tolerance = 1e-6)
})

test_that("AllDeepSets and AllSetTransformer train end to end", {
  skip_if_not_installed("torch")
  corpus <- c(a1 = "apple pear fruit", a2 = "pear fruit sweet",
              b1 = "orbit star galaxy", b2 = "star galaxy space")
  hg <- text_hypergraph(corpus)
  labels <- c(a1 = "food", b1 = "space")
  for (model in c("deepsets", "transformer")) {
    fit <- hypergraph_allset(
      hg, labels, model = model, hidden = 8, heads = 2,
      epochs = 5, validation = 0, seed = 3
    )
    expect_s3_class(fit, "data.frame")
    expect_equal(nrow(fit), 4L)
    expect_true(all(is.finite(attr(fit, "history")$loss)))
  }
})
