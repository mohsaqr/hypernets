.hnhn_fixture <- function() {
  events <- data.frame(
    node = c("a", "b", "c", "a", "b", "b", "c", "d"),
    edge = c("e1", "e1", "e1", "e2", "e2", "e3", "e3", "e3")
  )
  group_hypergraph(events, member = "node", group = "edge")
}

test_that("HNHN alpha=beta=0 gives directional neighborhood means", {
  hg <- .hnhn_fixture()
  op <- .thg_hnhn_operators(hg, alpha = 0, beta = 0)
  H <- (hg$incidence != 0) * 1
  expect_equal(unname(as.matrix(op$v_to_e)),
               unname(diag(1 / colSums(H)) %*% t(H)), tolerance = 1e-12)
  expect_equal(unname(as.matrix(op$e_to_v)),
               unname(diag(1 / rowSums(H)) %*% H), tolerance = 1e-12)
  expect_equal(as.numeric(Matrix::rowSums(op$v_to_e)),
               rep(1, hg$n_hyperedges))
  expect_equal(as.numeric(Matrix::rowSums(op$e_to_v)), rep(1, hg$n_nodes))
})

test_that("HNHN exponents implement the Section 2.3 equations", {
  hg <- .hnhn_fixture()
  alpha <- -0.5
  beta <- 1.5
  op <- .thg_hnhn_operators(hg, alpha, beta)
  H <- (hg$incidence != 0) * 1
  vd <- rowSums(H)
  ed <- colSums(H)
  ve <- diag(1 / as.numeric(t(H) %*% vd^beta)) %*%
    t(H) %*% diag(vd^beta)
  ev <- diag(1 / as.numeric(H %*% ed^alpha)) %*%
    H %*% diag(ed^alpha)
  expect_equal(as.matrix(op$v_to_e), ve, tolerance = 1e-12)
  expect_equal(as.matrix(op$e_to_v), ev, tolerance = 1e-12)
})

test_that("HNHN trains and reports normalization", {
  skip_if_not_installed("torch")
  corpus <- c(a1 = "apple pear fruit", a2 = "pear fruit sweet",
              b1 = "orbit star galaxy", b2 = "star galaxy space")
  hg <- text_hypergraph(corpus)
  fit <- hypergraph_hnhn(
    hg, c(a1 = "food", b1 = "space"), hidden = 8, alpha = -0.5,
    beta = 0.5, epochs = 10, validation = 0, seed = 1
  )
  expect_s3_class(fit, "data.frame")
  expect_equal(nrow(fit), 4L)
  expect_equal(attr(fit, "normalization"), c(alpha = -0.5, beta = 0.5))
  expect_true(all(is.finite(attr(fit, "history")$loss)))
})
