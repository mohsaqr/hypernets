testthat::skip_on_cran()

.plot_fixture <- function() {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e", "f"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e4")
  )
  hg <- group_hypergraph(dat, "member", "event")
  hg$edge_data <- data.frame(
    edge = c("e1", "e2", "e3", "e4"),
    sector = c("x", "y", "x", "z"),
    stringsAsFactors = FALSE
  )
  hg
}

test_that("fill selectors map to Okabe-Ito colours or the ramp", {
  discrete <- honets:::.thg_fill_colours(c("x", "y", "x"))
  expect_identical(discrete$colours, unname(discrete$palette[c("x", "y", "x")]))
  expect_identical(unname(discrete$palette), honets:::.thg_okabe_ito[1:2])
  numeric <- honets:::.thg_fill_colours(c(2, 12, 7))
  expect_null(numeric$palette)
  expect_identical(numeric$colours[1L], toupper("#F0E442"))
  expect_identical(numeric$colours[2L], toupper("#0072B2"))
  constant <- honets:::.thg_fill_colours(c(3, 3))
  expect_identical(constant$colours[1L], constant$colours[2L])
})

test_that("plot.net_hypergraph draws through cograph and does not print twice", {
  hg <- .plot_fixture()
  # cograph::plot_simplicial() prints as it draws; the method must not leave
  # a page on the open device before the caller prints the result.
  grDevices::pdf(NULL)
  before <- grDevices::dev.cur()
  p <- plot(hg)
  expect_identical(grDevices::dev.cur(), before)
  grDevices::dev.off()
  expect_s3_class(p, "ggplot")
  # one hyperedge blob per hyperedge with two or more members (e4 is a
  # singleton), drawn largest first
  polygons <- Filter(function(l) inherits(l$geom, "GeomPolygon"), p$layers)
  expect_length(polygons, 3L)
})

test_that("plot.net_hypergraph returns a ggplot for every selector", {
  hg <- .plot_fixture()
  default_plot <- plot(hg)
  expect_s3_class(default_plot, "ggplot")
  padded_plot <- plot(hg, padding = 0.03)
  expect_s3_class(padded_plot, "ggplot")
  size_plot <- plot(hg, color_by = "size")
  expect_s3_class(size_plot, "ggplot")
  sector_plot <- plot(hg, color_by = "sector", linetype_by = "sector")
  expect_s3_class(sector_plot, "ggplot")
  named_plot <- plot(hg, color_by = c(e1 = "u", e2 = "v", e3 = "u", e4 = "v"))
  expect_s3_class(named_plot, "ggplot")
  numeric_plot <- plot(hg, color_by = 1:4, labels = FALSE)
  expect_s3_class(numeric_plot, "ggplot")
  circle_plot <- plot(hg, labels = c(a = "Alpha"), layout = "circle")
  expect_s3_class(circle_plot, "ggplot")
})

test_that("a layout table is reused and validated", {
  hg <- .plot_fixture()
  pos <- honets:::.thg_layout(hg, "spring", seed = 3L)
  expect_identical(names(pos), c("node", "x", "y"))
  expect_identical(pos$node, hg$nodes)
  expect_true(all(pos$x >= 0 & pos$x <= 1 & pos$y >= 0 & pos$y <= 1))
  again <- honets:::.thg_layout(hg, pos, seed = 1L)
  expect_equal(again$x, pos$x)
  layout_plot <- plot(hg, layout = pos)
  expect_s3_class(layout_plot, "ggplot")
  expect_error(plot(hg, layout = pos[-1, ]), class = "honets_bad_input")
})

test_that("plot.net_hypergraph rejects bad selectors", {
  hg <- .plot_fixture()
  expect_error(plot(hg, color_by = "nope"), class = "honets_bad_input")
  expect_error(plot(hg, color_by = 1:3), class = "honets_bad_input")
  expect_error(plot(hg, color_by = c(e1 = 1)), class = "honets_bad_input")
  expect_error(plot(hg, labels = c("A", "B")), class = "honets_bad_input")
  expect_error(plot(hg, alpha = 2))
  singletons <- group_hypergraph(
    data.frame(member = c("a", "b"), event = c("e1", "e2")), "member", "event"
  )
  expect_error(plot(singletons), class = "honets_bad_input")
})
