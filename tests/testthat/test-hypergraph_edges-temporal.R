testthat::skip_on_cran()

.edges_thg <- function() {
  seats <- data.frame(
    case = rep(c("A", "B", "C"), each = 3),
    arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
    constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3)
  )
  temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                      start = "constituted", end = "concluded")
}

test_that("temporal hg_edges adds a time column per snapshot", {
  out <- hg_edges(.edges_thg())
  expect_identical(names(out), c("time", "edge", "size", "weight",
                                 "n_incident_edges", "n_neighbors"))
  expect_identical(out$time, c(1, 2, 2, 3, 3, 4, 4, 6))
  expect_identical(out$edge, c("A", "A", "B", "A", "B", "A", "C", "C"))
  # at time 2, A and B share a1: each meets the other and sees two outsiders
  expect_identical(out$n_incident_edges[out$time == 2], c(1L, 1L))
  expect_identical(out$n_neighbors[out$time == 2], c(2L, 2L))
})

test_that("hg_edges summary is a series over time", {
  s <- hg_edges(.edges_thg(), what = "summary", measure = "n_neighbors")
  expect_s3_class(s, "hypernets_series")
  expect_identical(names(s), c("time", "n_edges", "mean", "sd", "min", "q25",
                               "median", "q75", "max"))
  expect_identical(s$n_edges, c(1L, 2L, 2L, 2L, 1L))
  # A and B share a1 at times 2 and 3; A and C share p1 at time 4
  expect_equal(s$mean, c(0, 2, 2, 2, 0))
  summary_plot <- plot(s, columns = c("mean", "max"))
  expect_s3_class(summary_plot, "ggplot")
  snap_two <- hypergraph_snapshot(.edges_thg(), at = 2)
  one <- hg_edges(snap_two, what = "summary")
  expect_identical(class(one), "data.frame")
  expect_equal(one$mean, 3)
})

test_that("several thresholds give one block each, and distributions carry time", {
  snap <- hypergraph_snapshot(.edges_thg(), at = 2)
  multi <- hg_edges(snap, s = 1:2)
  expect_identical(names(multi)[1L], "s")
  expect_identical(multi$s, c(1L, 1L, 2L, 2L))
  expect_identical(multi$n_incident_edges, c(1L, 1L, 0L, 0L))
  single <- hg_edges(snap)
  expect_false("s" %in% names(single))
  d <- hg_edges(.edges_thg(), what = "distribution", at = c(2, 6))
  expect_s3_class(d, "hypernets_distribution")
  expect_identical(d$time, c(2, 6))
  expect_equal(d$ccdf, c(1, 1))
  distribution_plot <- plot(d)
  expect_s3_class(distribution_plot, "ggplot")
  expect_error(hg_edges(snap, s = 0), class = "hypernets_bad_input")
})
