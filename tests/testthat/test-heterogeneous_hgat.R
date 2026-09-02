.hgat_fixture <- function() {
  nodes <- c("d1", "d2", "t1", "e1")
  A <- matrix(0, 4, 4, dimnames = list(nodes, nodes))
  A["d1", c("t1", "e1")] <- 1
  A["d2", "t1"] <- 1
  A <- A + t(A)
  types <- c(d1 = "document", d2 = "document", t1 = "topic", e1 = "entity")
  features <- list(
    document = matrix(c(1, 0, 0, 1), 2, 2, byrow = TRUE,
                      dimnames = list(c("d1", "d2"), NULL)),
    topic = matrix(c(0.8, 0.2, 0), 1, 3,
                   dimnames = list("t1", NULL)),
    entity = matrix(c(0.1, 0.9, 0.2, 0.8), 1, 4,
                    dimnames = list("e1", NULL))
  )
  list(A = A, types = types, features = features)
}

test_that("HGAT aligns heterogeneous feature spaces and normalizes A+I", {
  x <- .hgat_fixture()
  input <- .thg_hgat_inputs(x$A, x$types, x$features)
  expect_identical(input$type_names, c("document", "entity", "topic"))
  expect_equal(input$features$document,
               x$features$document[c("d1", "d2"), , drop = FALSE])
  expect_true(isSymmetric(input$adjacency))
  expect_true(all(diag(input$adjacency) > 0))
})

test_that("heterogeneous HGAT trains across unequal feature widths", {
  skip_if_not_installed("torch")
  x <- .hgat_fixture()
  fit <- heterogeneous_hgat(
    x$A, x$types, x$features, labels = c(d1 = "a", d2 = "b"),
    hidden = 6, layers = 2, epochs = 5, validation = 0, seed = 1
  )
  expect_s3_class(fit, "data.frame")
  expect_equal(nrow(fit), 4L)
  expect_identical(attr(fit, "node_types"), x$types)
  expect_true(all(is.finite(attr(fit, "history")$loss)))
})

test_that("HGAT validates type-specific feature rows", {
  x <- .hgat_fixture()
  rownames(x$features$topic) <- "wrong"
  expect_error(.thg_hgat_inputs(x$A, x$types, x$features),
               "feature rows for type")
})
