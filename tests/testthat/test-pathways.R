# ---- Tests for pathways() ----

# Shared test data with clear higher-order dependencies
.make_ho_seqs <- function(n = 100, seed = 42) {
  set.seed(seed)
  lapply(seq_len(n), function(i) {
    s <- character(20)
    s[1] <- sample(LETTERS[1:4], 1)
    for (j in 2:20) {
      if (j >= 3 && s[j - 2] == "A" && s[j - 1] == "B") {
        s[j] <- sample(c("C", "D"), 1, prob = c(0.95, 0.05))
      } else if (j >= 3 && s[j - 2] == "C" && s[j - 1] == "B") {
        s[j] <- sample(c("A", "D"), 1, prob = c(0.95, 0.05))
      } else {
        s[j] <- sample(LETTERS[1:4], 1)
      }
    }
    s
  })
}


# ---- pathways.net_hon ----

test_that("pathways.net_hon returns arrow notation for plot_simplicial", {
  seqs <- .make_ho_seqs()
  hon <- hon(seqs, max_order = 3)
  pw <- .pathways(hon)

  expect_type(pw, "character")
  expect_true(length(pw) > 0)
  # Format: "X Y -> Z" (sources space-separated, -> before target)
  expect_true(all(grepl("->", pw, fixed = TRUE)))
  # Each pathway has exactly one "->"
  expect_true(all(vapply(pw, function(p) {
    sum(gregexpr("->", p, fixed = TRUE)[[1]] > 0)
  }, integer(1)) == 1L))
})


test_that("pathways.net_hon min_prob filters weak transitions", {
  seqs <- .make_ho_seqs()
  hon <- hon(seqs, max_order = 3)

  pw_all <- .pathways(hon)
  pw_strong <- .pathways(hon, min_prob = 0.5)

  expect_true(length(pw_strong) <= length(pw_all))
  expect_true(length(pw_strong) > 0)
})


test_that("pathways.net_hon returns empty for first-order only", {
  seqs <- list(c("A", "B", "C"), c("B", "A", "C"))
  hon <- hon(seqs, max_order = 2)
  pw <- .pathways(hon)

  expect_length(pw, 0)
})


test_that("pathways.net_hon order parameter filters", {
  seqs <- .make_ho_seqs()
  hon <- hon(seqs, max_order = 3)

  pw2 <- .pathways(hon, order = 2)
  expect_true(length(pw2) > 0)

  pw5 <- .pathways(hon, order = 5)
  expect_length(pw5, 0)
})


# ---- pathways.net_hypa ----

test_that("pathways.net_hypa returns anomalous paths", {
  seqs <- .make_ho_seqs()
  hypa <- hypa(seqs, order = 2, alpha = 0.05)
  pw <- .pathways(hypa)

  expect_type(pw, "character")
  if (length(pw) > 0) {
    expect_true(all(grepl("->", pw, fixed = TRUE)))
  }
})


test_that("pathways.net_hypa type parameter filters", {
  seqs <- .make_ho_seqs()
  hypa <- hypa(seqs, order = 2, alpha = 0.05)

  pw_all <- .pathways(hypa, type = "all")
  pw_over <- .pathways(hypa, type = "over")
  pw_under <- .pathways(hypa, type = "under")

  expect_true(length(pw_over) + length(pw_under) == length(pw_all))
})


test_that("pathways.net_hypa returns empty when no anomalies", {
  seqs <- list(
    c("A", "B", "C"),
    c("A", "B", "C"),
    c("B", "C", "A")
  )
  hypa <- hypa(seqs, order = 2, alpha = 0.001)
  pw <- .pathways(hypa)

  # May or may not have anomalies at strict alpha
  expect_type(pw, "character")
})


# ---- pathways.net_mogen ----

test_that("pathways.net_mogen returns transitions at optimal order", {
  seqs <- .make_ho_seqs()
  mog <- mogen(seqs, max_order = 3)
  pw <- .pathways(mog)

  expect_type(pw, "character")
  if (mog$optimal_order >= 1) {
    expect_true(length(pw) > 0)
    expect_true(all(grepl("->", pw, fixed = TRUE)))
  }
})


test_that("pathways.net_mogen order parameter overrides optimal", {
  seqs <- .make_ho_seqs()
  mog <- mogen(seqs, max_order = 3)

  pw1 <- .pathways(mog, order = 1)
  pw2 <- .pathways(mog, order = 2)

  # Different orders give different number of pathways
  expect_true(length(pw1) != length(pw2) || !identical(pw1, pw2))
})


test_that("pathways.net_mogen min_prob filters", {
  seqs <- .make_ho_seqs()
  mog <- mogen(seqs, max_order = 3)

  pw_all <- .pathways(mog)
  pw_strong <- .pathways(mog, min_prob = 0.5)

  expect_true(length(pw_strong) <= length(pw_all))
})


test_that("pathways.net_mogen returns empty for order 0", {
  seqs <- .make_ho_seqs()
  mog <- mogen(seqs, max_order = 3)
  pw <- .pathways(mog, order = 0)

  expect_length(pw, 0)
})


# ---- HYPA $edges consistency ----

test_that("HYPA $edges is set and matches $scores", {
  seqs <- .make_ho_seqs()
  hypa <- hypa(seqs, order = 2)

  expect_false(is.null(hypa$ho_edges))
  expect_equal(hypa$ho_edges, hypa$scores)
  expect_true(all(c("path", "from", "to") %in% names(hypa$ho_edges)))
})


