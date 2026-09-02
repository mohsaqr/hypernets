test_that("AMI and NMI match the paper's scikit-learn convention", {
  x <- data.frame(node = letters[1:6], label = c(1, 1, 1, 2, 2, 2))
  y <- data.frame(node = letters[1:6], label = c(1, 1, 2, 2, 3, 3))
  out <- hg_agreement(x, y, method = c("ari", "ami", "nmi"))
  expect_equal(out$ari, 0.242424242424242, tolerance = 1e-12)
  expect_equal(out$ami, 0.298792458170890, tolerance = 1e-12)
  expect_equal(out$nmi, 0.515803742979389, tolerance = 1e-12)
  expect_equal(hg_agreement(x, x, method = c("ami", "nmi"))[, c("ami", "nmi")],
               data.frame(ami = 1, nmi = 1))
})

test_that("Infomap ensemble returns AMI medoid and all representations", {
  skip_if_not_installed("igraph")
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "c", "x", "y", "z"),
               event = rep(c("left", "right"), each = 3)),
    "member", "event"
  )
  fit <- hg_communities(h, n_runs = 3, trials = 3, seeds = 1:3)
  expect_s3_class(fit, "hg_communities")
  expect_equal(nrow(fit$medoid), 6)
  expect_equal(nrow(fit$partitions), 18)
  expect_true(fit$medoid_run %in% 1:3)
  expect_equal(unname(diag(fit$similarity$ami)), rep(1, 3))
  expect_equal(as.data.frame(fit), fit$medoid)
  expect_equal(nrow(as.data.frame(fit, what = "ami")), 9)
})

test_that("paper partition quality is exact on disconnected cliques", {
  skip_if_not_installed("igraph")
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "c", "x", "y", "z"),
               event = rep(c("left", "right"), each = 3)),
    "member", "event"
  )
  labels <- c(a = 1, b = 1, c = 1, x = 2, y = 2, z = 2)
  q <- hg_community_quality(h, labels)
  expect_equal(q$coverage, 1)
  expect_equal(q$weighted_coverage, 1)
  expect_equal(q$performance, 1)
  expect_equal(q$modularity, 0.5)
  expect_equal(q$conductance, 0)
})

test_that("community conductance is exact on a path cut", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c", "c", "d"),
               event = rep(c("ab", "bc", "cd"), each = 2)),
    "member", "event"
  )
  labels <- c(a = "left", b = "left", c = "right", d = "right")
  q <- hg_community_quality(h, labels)
  # One crossing edge; each side has weighted volume 3.
  expect_equal(q$conductance, 1 / 3)
})

test_that("one-community conductance is undefined", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c"),
               event = rep(c("ab", "bc"), each = 2)),
    "member", "event"
  )
  expect_true(is.na(hg_community_quality(h, c(a = 1, b = 1, c = 1))$conductance))
})
