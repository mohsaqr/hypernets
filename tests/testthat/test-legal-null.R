testthat::skip_on_cran()

# {1, 2, 3} twice and {1, 2, 4}: one repeated set; 9 member pairs with
# multiplicity against 5 distinct pairs, so 4 repeated pairs.
.null_fixture <- function() {
  group_hypergraph(
    data.frame(member = c(1, 2, 3, 1, 2, 3, 1, 2, 4),
               event = rep(c("e1", "e2", "e3"), each = 3)),
    "member", "event"
  )
}

test_that("repeated-collaboration statistics are counted exactly", {
  m <- (.null_fixture()$incidence > 0) * 1L
  got <- honets:::.thg_null_statistics(m, c("repeated_edges", "repeated_pairs"))
  expect_equal(unname(got), c(1, 4))
  out <- hg_null_test(.null_fixture(), statistic = c("repeated_edges", "repeated_pairs"),
                      method = "assignment", n = 19, seed = 1)
  expect_identical(out$statistic, c("repeated_edges", "repeated_pairs"))
  expect_equal(out$observed, c(1, 4))
  expect_identical(out$method, rep("assignment", 2))
})

test_that("the degree-ordered assignment preserves both margins exactly", {
  set.seed(4)
  m <- matrix(0L, 8, 12)
  m[cbind(sample(8, 36, replace = TRUE), rep(1:12, each = 3))] <- 1L
  m <- (m > 0) * 1L
  for (i in seq_len(20)) {
    draw <- honets:::.thg_assignment_draw(m)
    expect_identical(rowSums(draw), rowSums(m))
    expect_identical(colSums(draw), colSums(m))
    expect_true(all(draw %in% c(0L, 1L)))
  }
})

test_that("hg_null_test with the assignment null is reproducible", {
  a <- hg_null_test(.null_fixture(), statistic = "repeated_pairs",
                    method = "assignment", n = 19, seed = 3)
  b <- hg_null_test(.null_fixture(), statistic = "repeated_pairs",
                    method = "assignment", n = 19, seed = 3)
  expect_equal(a, b)
})
