.legal_motif_hg <- function(edges) {
  group_hypergraph(
    do.call(rbind, lapply(seq_along(edges), function(i) {
      data.frame(member = edges[[i]], event = paste0("e", i))
    })), "member", "event"
  )
}

test_that("induced Y, T and O motifs are counted exactly", {
  triples <- list(c(1, 2, 3), c(1, 2, 4), c(1, 3, 4), c(2, 3, 4))
  expect_equal(hg_motifs(.legal_motif_hg(triples[1:2]), what = "counts")$count,
               c(1L, 0L, 0L))
  expect_equal(hg_motifs(.legal_motif_hg(triples[1:3]), what = "counts")$count,
               c(0L, 1L, 0L))
  expect_equal(hg_motifs(.legal_motif_hg(triples), what = "counts")$count,
               c(0L, 0L, 1L))
})

test_that("configuration MCMC preserves degree and cardinality before collapse", {
  edges <- list(c(1, 2, 3), c(1, 4, 5), c(2, 4, 6), c(3, 5, 6))
  set.seed(10)
  draw <- honets:::.thg_configuration_mcmc(edges, steps = 100,
                                            collapse = FALSE)
  expect_equal(sort(lengths(draw)), sort(lengths(edges)))
  expect_equal(tabulate(unlist(draw), 6), tabulate(unlist(edges), 6))
})

test_that("motif null profiles are reproducible and preserve caller RNG", {
  h <- .legal_motif_hg(list(c(1, 2, 3), c(1, 2, 4), c(1, 3, 5), c(2, 4, 5)))
  set.seed(99)
  before <- .Random.seed
  a <- hg_motifs(h, n = 9, seed = 7)
  after <- .Random.seed
  b <- hg_motifs(h, n = 9, seed = 7)
  expect_identical(before, after)
  expect_equal(a, b)
  expect_equal(sum(a$normalized_delta^2), 1, tolerance = 1e-12)
})

test_that("motif census applies to temporal snapshots", {
  dat <- data.frame(
    member = c(1, 2, 3, 1, 2, 4),
    event = rep(c("e1", "e2"), each = 3), time = rep(1:2, each = 3)
  )
  thg <- temporal_hypergraph(dat, "member", "event", time = "time")
  out <- hg_motifs(thg, what = "counts")
  expect_equal(out$count[out$time == "1"], c(0L, 0L, 0L))
  expect_equal(out$count[out$time == "2"], c(1L, 0L, 0L))
})

test_that("motifs require a 3-uniform hypergraph", {
  h <- .legal_motif_hg(list(c(1, 2), c(1, 2, 3)))
  expect_error(hg_motifs(h, what = "counts"), class = "honets_bad_input")
})
