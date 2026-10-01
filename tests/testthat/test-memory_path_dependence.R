testthat::skip_on_cran()

test_that("hg_get.net_path_dependence filters, sorts and truncates", {
  set.seed(3)
  # memory() (Nestimate's) reads a wide frame, one sequence per row
  wide <- as.data.frame(t(replicate(30, sample(c("a", "b", "c"), 15,
                                               replace = TRUE))),
                        stringsAsFactors = FALSE)
  pd <- memory(wide, order = 2L, min_count = 1L)
  expect_s3_class(pd, "net_path_dependence")
  full <- hg_get(pd)
  expect_identical(class(full), "data.frame")
  # the estimator's context table in the package vocabulary
  expect_identical(full, .ho_rename(pd$contexts, c(
    n = "count", H_order1 = "entropy_first_order",
    H_orderk = "entropy_order_k", H_drop = "entropy_drop", KL = "kl",
    top_o1 = "top_first_order", top_ok = "top_order_k")),
    ignore_attr = "row.names")

  kept <- hg_get(pd, min_count = 10L)
  expect_true(all(kept$count >= 10L))
  expect_identical(nrow(kept), sum(full$count >= 10L))

  by_kl <- hg_get(pd, sort_by = "kl")
  expect_false(is.unsorted(rev(by_kl$kl)))
  expect_identical(hg_get(pd, sort_by = "kl", top = 3L),
                   `rownames<-`(utils::head(by_kl, 3L), NULL))
  expect_error(hg_get(pd, min_count = 0L), "`min_count` must be")
})


test_that("hg_get(flips =) filters contexts that change the prediction", {
  seqs <- lapply(1:30, function(i) c("a", "b", if (i %% 2) "c" else "d",
                                     "a", "b", if (i %% 3) "c" else "d"))
  pd <- memory(seqs, order = 2, min_count = 1)
  everything <- hg_get(pd)
  flipped <- hg_get(pd, flips = TRUE)
  kept <- hg_get(pd, flips = FALSE)
  expect_true(all(flipped$flips))
  expect_false(any(kept$flips))
  expect_identical(nrow(flipped) + nrow(kept), nrow(everything))
  expect_error(hg_get(pd, flips = NA))
})
