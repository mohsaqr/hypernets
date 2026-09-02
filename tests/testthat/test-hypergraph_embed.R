test_that("embedding verb returns only coordinates and stationary mass", {
  h <- group_hypergraph(
    data.frame(node = c("a", "b", "c", "b", "c", "d", "a", "d"),
               edge = c(rep("e1", 3), rep("e2", 3), rep("e3", 2))),
    "node", "edge"
  )
  out <- hg_embed(h, dimensions = 2, seed = 4)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("node", "pi", "dim1", "dim2"))
  expect_equal(nrow(out), h$n_nodes)
  expect_equal(sum(out$pi), 1)
  expect_identical(attr(out, "method"), "spectral")
})

test_that("embedding verb exposes the shared engine coordinates", {
  h <- group_hypergraph(
    data.frame(node = c("a", "b", "c", "b", "c", "d", "a", "d"),
               edge = c(rep("e1", 3), rep("e2", 3), rep("e3", 2))),
    "node", "edge"
  )
  direct <- hg_embed(h, dimensions = 2, seed = 9)
  clustered <- hg_cluster(h, k = 2, seed = 9, what = "embedding")
  expect_equal(unclass(direct[, c("node", "pi", "dim1", "dim2")]),
               unclass(clustered[, c("node", "pi", "dim1", "dim2")]))
  expect_identical(hypergraph_embed, hg_embed)
})

test_that("embedding verb validates dimensions", {
  h <- group_hypergraph(
    data.frame(node = c("a", "b", "b", "c"),
               edge = rep(c("e1", "e2"), each = 2)),
    "node", "edge"
  )
  expect_error(hg_embed(h, dimensions = 1.5), "whole number")
  expect_error(hg_embed(h, dimensions = 1), "between 2")
})
