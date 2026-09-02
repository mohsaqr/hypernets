test_that("HyperGCN mediator expansion implements Algorithms 1 and 2", {
  events <- data.frame(node = letters[1:4], edge = "e")
  hg <- group_hypergraph(events, member = "node", group = "edge")
  signal <- matrix(c(0, 3, 1, 2), ncol = 1)
  A <- .thg_hypergcn_adjacency(hg, signal, mediators = TRUE)
  raw <- diag(1, 4)
  links <- rbind(c(1, 2), c(1, 3), c(2, 3), c(1, 4), c(2, 4))
  for (i in seq_len(nrow(links))) {
    raw[links[i, 1], links[i, 2]] <- raw[links[i, 1], links[i, 2]] + 1 / 5
    raw[links[i, 2], links[i, 1]] <- raw[links[i, 2], links[i, 1]] + 1 / 5
  }
  expected <- raw / sqrt(outer(rowSums(raw), rowSums(raw)))
  expect_equal(unname(A), expected, tolerance = 1e-12)
  expect_true(isSymmetric(unname(A)))
})

test_that("1-HyperGCN retains only the farthest pair", {
  events <- data.frame(node = letters[1:4], edge = "e")
  hg <- group_hypergraph(events, member = "node", group = "edge")
  A <- .thg_hypergcn_adjacency(hg, matrix(c(0, 3, 1, 2), ncol = 1),
                                mediators = FALSE)
  expect_gt(A[1, 2], 0)
  expect_equal(A[1, 3], 0)
  expect_equal(A[2, 4], 0)
})

test_that("HyperGCN variants train end to end", {
  skip_if_not_installed("torch")
  corpus <- c(
    cooking_1 = "simmer soup onions carrots",
    cooking_2 = "soup recipe salt carrots",
    space_1 = "telescope galaxy stars orbit",
    space_2 = "astronomers telescope stars orbit"
  )
  hg <- text_hypergraph(corpus)
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  for (method in c("hypergcn", "fast", "one")) {
    fit <- hypergraph_hypergcn(
      hg, labels, method = method, hidden = 8, epochs = 10,
      validation = 0, seed = 2
    )
    expect_s3_class(fit, "data.frame")
    expect_equal(nrow(fit), hg$n_nodes)
    expect_identical(attr(fit, "method"), method)
    expect_true(all(is.finite(attr(fit, "history")$loss)))
  }
})
