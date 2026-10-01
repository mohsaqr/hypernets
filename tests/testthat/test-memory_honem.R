# ===========================================================================
# Section 4: honem end-to-end
# ===========================================================================

test_that("honem returns net_honem class from matrix", {
  mat <- matrix(c(0, 2, 0, 0,
                  0, 0, 3, 0,
                  1, 0, 0, 2,
                  0, 1, 0, 0), 4, 4, byrow = TRUE,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  emb <- honem(mat, dim = 2L, max_power = 3L)

  expect_s3_class(emb, "net_honem")
  expect_equal(emb$n_nodes, 4L)
  expect_equal(emb$dim, 2L)
  expect_equal(nrow(emb$embeddings), 4L)
  expect_equal(ncol(emb$embeddings), 2L)
})

test_that("honem works with net_hon object", {
  trajs <- list(c("A", "B", "C", "D"), c("A", "B", "D", "C"),
                c("B", "C", "D", "A"), c("C", "D", "A", "B"))
  hon <- hon(trajs, max_order = 2L, method = "hon")
  emb <- honem(hon, dim = 3L)

  expect_s3_class(emb, "net_honem")
  expect_equal(emb$n_nodes, hon$n_nodes)
})

test_that("honem preserves node names", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2, byrow = TRUE,
                dimnames = list(c("X", "Y"), c("X", "Y")))
  emb <- honem(mat, dim = 1L)

  expect_equal(emb$nodes, c("X", "Y"))
  expect_equal(rownames(emb$embeddings), c("X", "Y"))
})

test_that("honem rejects invalid input", {
  expect_error(honem(42), "net_hon object or a square matrix")
  expect_error(honem(matrix(1, 1, 1)), "at least 2 nodes")
})

# ===========================================================================
# Section 5: Embedding quality
# ===========================================================================

test_that("honem embeddings reflect network structure", {
  # Create a network with two clusters
  mat <- matrix(0, 6, 6, dimnames = list(LETTERS[1:6], LETTERS[1:6]))
  # Cluster 1: A, B, C (strong connections)
  mat["A", "B"] <- 5; mat["B", "C"] <- 5; mat["C", "A"] <- 5
  # Cluster 2: D, E, F (strong connections)
  mat["D", "E"] <- 5; mat["E", "F"] <- 5; mat["F", "D"] <- 5
  # Weak cross-cluster link
  mat["C", "D"] <- 1

  emb <- honem(mat, dim = 2L, max_power = 5L)

  # Nodes within same cluster should be closer than across clusters
  d_within1 <- sqrt(sum((emb$embeddings["A", ] - emb$embeddings["B", ])^2))
  d_across <- sqrt(sum((emb$embeddings["A", ] - emb$embeddings["E", ])^2))

  expect_true(d_within1 < d_across)
})

# ===========================================================================
# Section 6: S3 methods
# ===========================================================================

test_that("print.net_honem works", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- honem(mat, dim = 1L)
  out <- capture.output(print(emb))
  expect_true(any(grepl("Higher-order network embedding", out)))
})

test_that("summary.net_honem works", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- honem(mat, dim = 1L)
  out <- capture.output(print(summary(emb)))
  expect_true(any(grepl("variance", out)))
})

test_that("plot.net_honem works", {
  mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3,
                dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
  emb <- honem(mat, dim = 2L)
  expect_no_error(plot(emb))
})

# ===========================================================================
# Section 7: Coverage for plot.net_honem with dim < 2
# ===========================================================================

test_that("plot.net_honem issues message and returns invisibly when dim < 2", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- honem(mat, dim = 1L)
  expect_equal(emb$dim, 1L)
  expect_message(result <- plot(emb), "at least 2 dimensions")
  expect_identical(result, emb)
})


# ===========================================================================
# Section 8: Input validation — additional edge cases
# ===========================================================================

test_that("honem rejects non-square matrix", {
  mat <- matrix(1:6, 2, 3)
  expect_error(honem(mat), "net_hon object or a square matrix")
})

test_that("honem rejects dim < 1", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(honem(mat, dim = 0L), "'dim' must be >= 1")
})

test_that("honem rejects max_power < 1", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  expect_error(honem(mat, max_power = 0L), "'max_power' must be >= 1")
})

test_that("honem handles all-zero matrix (degenerate embeddings)", {
  mat <- matrix(0, 3, 3, dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L)
  # Should not error — produces zero embeddings

  expect_s3_class(emb, "net_honem")
  expect_equal(nrow(emb$embeddings), 3L)
  expect_equal(ncol(emb$embeddings), 2L)
  # All embeddings should be zero (no information in a zero matrix)
  expect_true(all(emb$embeddings == 0))
})

test_that("honem rejects non-matrix non-hon input types", {
  expect_error(honem(data.frame(a = 1, b = 2)), "net_hon object or a square matrix")
  expect_error(honem("not a matrix"), "net_hon object or a square matrix")
  expect_error(honem(list(a = 1)), "net_hon object or a square matrix")
})


# ===========================================================================
# Section 11: End-to-end — field verification
# ===========================================================================

test_that("honem result has all expected fields", {
  mat <- matrix(c(0, 2, 1, 3, 0, 2, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L, max_power = 5L)

  # All required fields present
  expect_true("embeddings" %in% names(emb))
  expect_true("nodes" %in% names(emb))
  expect_true("singular_values" %in% names(emb))
  expect_true("explained_variance" %in% names(emb))
  expect_true("dim" %in% names(emb))
  expect_true("max_power" %in% names(emb))
  expect_true("n_nodes" %in% names(emb))
})

test_that("honem embeddings dimensions match n_nodes x dim", {
  mat <- matrix(c(0, 3, 0, 1,
                  2, 0, 1, 0,
                  0, 2, 0, 3,
                  1, 0, 2, 0), 4, 4, byrow = TRUE,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  emb <- honem(mat, dim = 3L)

  expect_equal(nrow(emb$embeddings), emb$n_nodes)
  expect_equal(ncol(emb$embeddings), emb$dim)
  expect_equal(emb$n_nodes, 4L)
  expect_equal(emb$dim, 3L)
})

test_that("honem singular_values are non-negative and sorted descending", {
  set.seed(12)
  mat <- matrix(sample(0:5, 25, replace = TRUE), 5, 5,
                dimnames = list(LETTERS[1:5], LETTERS[1:5]))
  diag(mat) <- 0
  emb <- honem(mat, dim = 4L)

  expect_true(all(emb$singular_values >= 0))
  # Check sorted descending (or equal)
  diffs <- diff(emb$singular_values)
  expect_true(all(diffs <= 1e-10))
})

test_that("honem explained_variance is in [0, 1]", {
  mat <- matrix(c(0, 5, 1, 3, 0, 2, 1, 4, 0), 3, 3, byrow = TRUE,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L)

  expect_gte(emb$explained_variance, 0)
  expect_lte(emb$explained_variance, 1)
})

test_that("honem max_power stored correctly", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- honem(mat, dim = 1L, max_power = 7L)
  expect_equal(emb$max_power, 7L)
})


# ===========================================================================
# Section 12: S3 methods — detailed output checks
# ===========================================================================

test_that("print.net_honem output contains key info", {
  mat <- matrix(c(0, 2, 1, 3, 0, 2, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L, max_power = 5L)
  out <- capture.output(print(emb))

  # a header with the size, then the embedding table
  expect_identical(out[1L], "Higher-order network embedding: 3 nodes in 2 dimensions")
  expect_match(out[2L], "node +dim1 +dim2")
  expect_identical(length(out), 2L + 3L)
})

test_that("print.net_honem returns invisibly", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- honem(mat, dim = 1L)
  out <- capture.output(result <- print(emb))
  expect_identical(result, emb)
})

test_that("summary.net_honem output contains variance and singular values info", {
  mat <- matrix(c(0, 2, 1, 3, 0, 2, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L, max_power = 5L)
  out <- capture.output(print(summary(emb)))
  expect_identical(out[1L], "Higher-order network embedding: 3 nodes in 2 dimensions")
  expect_true(any(grepl("Tables: embeddings, variance", out)))
  expect_true(any(grepl("singular_value", out)))
})

test_that("summary.net_honem returns tidy embedding data.frame", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  emb <- honem(mat, dim = 1L)
  result <- summary(emb)$embeddings
  expect_s3_class(result, "data.frame")
  expect_equal(result$node, c("A", "B"))
  expect_true("dim1" %in% names(result))
})

test_that("plot.net_honem works with dim >= 2 and labels <= 50 nodes", {
  mat <- matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), 3, 3,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L)
  expect_no_error(result <- plot(emb))
  expect_identical(result, emb)
})

test_that("plot.net_honem suppresses labels for > 50 nodes", {
  # Create a larger matrix — plot should skip text labels
  n <- 55
  mat <- matrix(sample(0:3, n * n, replace = TRUE), n, n)
  diag(mat) <- 0
  nms <- paste0("N", seq_len(n))
  dimnames(mat) <- list(nms, nms)
  emb <- honem(mat, dim = 2L)
  # Just verify it doesn't error with >50 nodes
  expect_no_error(plot(emb))
})

test_that("plot.net_honem respects custom dims argument", {
  mat <- matrix(c(0, 3, 1, 2,
                  1, 0, 2, 0,
                  0, 1, 0, 3,
                  2, 0, 1, 0), 4, 4, byrow = TRUE,
                dimnames = list(LETTERS[1:4], LETTERS[1:4]))
  emb <- honem(mat, dim = 3L)
  # Plot dims 1 and 3
  expect_no_error(plot(emb, dims = c(1L, 3L)))
})


# ===========================================================================
# Section 14: Dim capping in honem
# ===========================================================================

test_that("honem caps dim at n-1 for small matrices", {
  mat <- matrix(c(0, 1, 1, 0), 2, 2,
                dimnames = list(c("A", "B"), c("A", "B")))
  # Request dim=100 but only 2 nodes => capped at 1
  emb <- honem(mat, dim = 100L)
  expect_equal(emb$dim, 1L)
  expect_equal(ncol(emb$embeddings), 1L)
})

test_that("honem with dim equal to n-1 works", {
  mat <- matrix(c(0, 2, 1, 3, 0, 2, 1, 3, 0), 3, 3, byrow = TRUE,
                dimnames = list(LETTERS[1:3], LETTERS[1:3]))
  emb <- honem(mat, dim = 2L)  # n-1 = 2
  expect_equal(emb$dim, 2L)
  expect_equal(ncol(emb$embeddings), 2L)
})
