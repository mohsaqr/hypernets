testthat::skip_on_cran()

.motif_hg <- function(edges) {
  group_hypergraph(
    do.call(rbind, lapply(seq_along(edges), function(i) {
      data.frame(member = edges[[i]], event = paste0("e", i))
    })), "member", "event"
  )
}

# The pairwise loop the vectorised census replaced, kept as its oracle.
.motif_loop <- function(edges) {
  if (length(edges) < 2L) return(c(Y = 0L, T = 0L, O = 0L))
  four_sets <- character()
  for (i in seq_len(length(edges) - 1L)) {
    for (j in (i + 1L):length(edges)) {
      if (length(intersect(edges[[i]], edges[[j]])) == 2L) {
        union_nodes <- sort(union(edges[[i]], edges[[j]]))
        if (length(union_nodes) == 4L) {
          four_sets <- c(four_sets, paste(union_nodes, collapse = ","))
        }
      }
    }
  }
  if (!length(four_sets)) return(c(Y = 0L, T = 0L, O = 0L))
  pair_count <- as.integer(table(four_sets))
  c(Y = sum(pair_count == 1L), T = sum(pair_count == 3L),
    O = sum(pair_count == 6L))
}

test_that("INVARIANT: vectorised motif census equals the pairwise loop", {
  set.seed(11)
  for (r in seq_len(25)) {
    n <- sample(6:14, 1)
    m <- sample(4:24, 1)
    edges <- unique(lapply(seq_len(m), function(i) sort(sample(n, 3))))
    expect_identical(honets:::.thg_yto_counts(edges), .motif_loop(edges))
  }
})

test_that("motif tests carry their draws and plot", {
  h <- .motif_hg(list(c(1, 2, 3), c(1, 2, 4), c(1, 3, 5), c(2, 4, 5)))
  out <- hg_motifs(h, n = 9, seed = 7)
  expect_s3_class(out, "honets_motifs")
  expect_s3_class(out, "data.frame")
  test <- as.data.frame(out)
  expect_identical(class(test), "data.frame")
  expect_identical(test$motif, c("Y", "T", "O"))
  draws <- as.data.frame(out, what = "draws")
  expect_identical(names(draws), c("run", "motif", "count"))
  expect_identical(nrow(draws), 27L)
  expect_equal(mean(draws$count[draws$motif == "Y"]), test$null_mean[[1L]])
  expect_s3_class(plot(out), "ggplot")
  expect_s3_class(plot(out, motif = "T"), "ggplot")
  expect_identical(class(hg_motifs(h, what = "counts")), "data.frame")
})
