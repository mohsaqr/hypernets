skip_on_cran()

planted_hg <- function(seed = 1) {
  set.seed(seed)
  bg <- data.frame(actor = sample(letters[1:12], 400, TRUE),
                   group = paste0("g", sample(60, 400, TRUE)))
  planted <- data.frame(actor = rep(c("x", "y"), 25),
                        group = rep(paste0("p", 1:25), each = 2))
  group_hypergraph(rbind(bg, planted), actor = "actor", group = "group")
}

test_that("hg_hypa recovers a planted over-represented pair", {
  out <- hg_hypa(planted_hg(), min_count = 3L)
  hit <- out[out$from == "x" & out$to == "y", ]
  expect_identical(nrow(hit), 1L)
  expect_identical(hit$anomaly, "over")
  expect_gt(hit$ratio, 2)
  # and it does not cry wolf on the random background
  background <- out[!(out$from == "x" & out$to == "y"), ]
  expect_true(all(background$anomaly == "normal"))
})

test_that("the null is anchored to all pairs, not to the survivors", {
  hg <- planted_hg()
  low <- hg_hypa(hg, min_count = 2L)
  high <- hg_hypa(hg, min_count = 8L)
  a <- low[low$from == "x" & low$to == "y", "expected"]
  b <- high[high$from == "x" & high$to == "y", "expected"]
  # raising the floor changes which pairs are tested, never the null itself
  expect_equal(a, b)
  expect_gt(nrow(low), nrow(high))
})

test_that("hg_hypa returns a tidy frame with ratio = observed / expected", {
  out <- hg_hypa(planted_hg(), min_count = 3L)
  expect_s3_class(out, "data.frame")
  expect_identical(names(out),
    c("from", "to", "observed", "expected", "ratio", "p_under", "p_over",
      "p_adjusted_under", "p_adjusted_over", "anomaly"))
  expect_equal(out$ratio, out$observed / out$expected)
  # ordered by significance, not by ratio (ratio favours pairs at the floor)
  expect_false(is.unsorted(out$p_adjusted_over))
  expect_true(all(out$observed >= 3L))
})

test_that("hg_hypa honours top and p_adjust = 'none'", {
  hg <- planted_hg()
  expect_identical(nrow(hg_hypa(hg, min_count = 3L, top = 4L)), 4L)
  raw <- hg_hypa(hg, min_count = 3L, p_adjust = "none")
  expect_identical(raw$p_adjusted_over, raw$p_over)
})

test_that("hg_hypa rejects bad input by class", {
  hg <- planted_hg()
  expect_error(hg_hypa(hg, min_count = 0), "min_count")
  expect_error(hg_hypa(hg, alpha = 2), "alpha")
  expect_error(hg_hypa(hg, p_adjust = "nope"), class = "hypernets_bad_input")
  expect_error(hg_hypa(hg, min_count = 10000L),
               class = "hypernets_empty_result")
})

test_that("hg_hypa is invariant to relabelling the nodes", {
  hg <- planted_hg()
  out <- hg_hypa(hg, min_count = 3L)
  relabelled <- hg
  map <- stats::setNames(paste0("N", seq_along(hg$nodes)), hg$nodes)
  relabelled$nodes <- unname(map[hg$nodes])
  rownames(relabelled$incidence) <- relabelled$nodes
  out2 <- hg_hypa(relabelled, min_count = 3L)
  expect_equal(sort(out$ratio), sort(out2$ratio))
})
