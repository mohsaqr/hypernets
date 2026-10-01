testthat::skip_on_cran()

# The blob look of plot.net_hg(): title boxes, notes, the point-size
# variant, units, the bottom legend, and pieces laid out in a row. Parity with
# the pipeline helper plot_blobs() on the real event data is checked in
# local_testing_and_equivalence/test-equiv-plot-blobs-eventdata.R.

.blob_members <- function() {
  data.frame(
    group = c(rep("Hint", 3), rep("Question", 3), rep("Exit", 3)),
    state = c("Wrong", "Hint", "Retry", "Wrong", "Question", "Retry",
              "Wrong", "Assurance", "Quit"),
    trials = c(rep(1200, 3), rep(800, 3), rep(95.5, 3)),
    stringsAsFactors = FALSE
  )
}

.blob_hg <- function(members = .blob_members()) {
  group_hypergraph(members, actor = "state", group = "group")
}

.layers_of <- function(p, geom) {
  Filter(function(l) inherits(l$geom, geom), p$layers)
}

test_that("titles write name, count and unit beside each pebble", {
  hg <- .blob_hg()
  p <- plot(hg, color_by = "trials", titles = "trials", unit = "trials",
            title_prefix = "Mostly ")
  boxes <- .layers_of(p, "GeomLabel")
  expect_length(boxes, 1L)
  box <- boxes[[1L]]$data
  # in order of decreasing count, whole numbers with a separator, others to
  # three significant digits
  expect_identical(box$label, c("Mostly Hint\n1,200 trials",
                                "Mostly Question\n800 trials",
                                "Mostly Exit\n95.5 trials"))
  # each box lies beyond its own pebble, seen from the figure's centre
  hulls <- p$layers[[1L]]$data
  centre <- c(mean(range(hulls$x)), mean(range(hulls$y)))
  beyond <- vapply(seq_len(nrow(box)), function(i) {
    ring <- hulls[hulls$hyperedge == c("Hint", "Question", "Exit")[i], ]
    heading <- c(mean(ring$x), mean(ring$y)) - centre
    heading <- heading / sqrt(sum(heading^2))
    reach <- max((ring$x - centre[1L]) * heading[1L] +
                   (ring$y - centre[2L]) * heading[2L])
    along <- (box$x[i] - centre[1L]) * heading[1L] +
      (box$y[i] - centre[2L]) * heading[2L]
    abs(along - (reach + 0.06)) < 1e-9
  }, logical(1L))
  expect_true(all(beyond))
  expect_true(all(box$hjust %in% c(0, 0.5, 1)))
  # TRUE writes names only; a named vector works like a column
  named <- plot(hg, titles = TRUE)
  expect_setequal(.layers_of(named, "GeomLabel")[[1L]]$data$label,
                  c("Hint", "Question", "Exit"))
  vector_titles <- plot(hg, titles = c(Hint = 1200, Question = 800, Exit = 95.5),
                        unit = "steps")
  expect_identical(.layers_of(vector_titles, "GeomLabel")[[1L]]$data$label,
                   c("Hint\n1,200 steps", "Question\n800 steps", "Exit\n95.5 steps"))
})

test_that("notes add a line; unknown names are skipped", {
  hg <- .blob_hg()
  as_table <- plot(hg, titles = "trials", unit = "trials",
                   notes = data.frame(group = c("Exit", "nope"),
                                      note = c("solved 6%", "x")))
  expect_identical(.layers_of(as_table, "GeomLabel")[[1L]]$data$label[3L],
                   "Exit\n95.5 trials\nsolved 6%")
  as_vector <- plot(hg, titles = "trials", unit = "trials",
                    notes = c(Exit = "solved 6%"))
  expect_identical(.layers_of(as_vector, "GeomLabel")[[1L]]$data,
                   .layers_of(as_table, "GeomLabel")[[1L]]$data)
})

test_that("a box that would sit on an earlier one moves below it", {
  # two hyperedges with the same members have the same outline and so the
  # same anchor: the second box drops by 8% of the figure's height
  hg <- group_hypergraph(
    data.frame(state = c("a", "b", "c", "a", "b", "c", "c", "d"),
               group = c("g1", "g1", "g1", "g2", "g2", "g2", "g3", "g3"),
               trials = c(5, 5, 5, 4, 4, 4, 1, 1)),
    actor = "state", group = "group")
  p <- plot(hg, titles = "trials")
  box <- .layers_of(p, "GeomLabel")[[1L]]$data
  spread <- diff(range(p$layers[[1L]]$data$y))
  expect_equal(box$x[1L], box$x[2L])
  expect_equal(box$y[1L] - box$y[2L], 0.08 * spread)
  # title_gap moves every box further out along its heading
  far <- .layers_of(plot(hg, titles = "trials", title_gap = 0.2), "GeomLabel")[[1L]]$data
  centre <- c(mean(range(p$layers[[1L]]$data$x)), mean(range(p$layers[[1L]]$data$y)))
  expect_true(all(sqrt((far$x - centre[1L])^2 + (far$y - centre[2L])^2) >
                    sqrt((box$x - centre[1L])^2 + (box$y - centre[2L])^2)))
  expect_error(plot(hg, title_gap = -1), class = "hypernets_bad_input")
})

test_that("a one-member hyperedge gets its box beside its node", {
  hg <- group_hypergraph(
    data.frame(state = c("a", "b", "c", "b", "c", "d", "d"),
               group = c("g1", "g1", "g1", "g2", "g2", "g2", "solo"),
               trials = c(3, 3, 3, 2, 2, 2, 1)),
    actor = "state", group = "group")
  p <- plot(hg, titles = "trials", layout = "circle")
  box <- .layers_of(p, "GeomLabel")[[1L]]$data
  expect_identical(nrow(box), 3L)
  expect_true(all(is.finite(box$x) & is.finite(box$y)))
})

test_that("node sizes without direction are points with an area legend", {
  hg <- .blob_hg()
  sizes <- data.frame(state = c("Wrong", "Hint", "Retry", "Question",
                                "Assurance", "Quit"),
                      trials = c(2095.5, 1200, 2000, 800, 95.5, 95.5))
  p <- plot(hg, node_sizes = sizes)
  points <- .layers_of(p, "GeomPoint")
  expect_length(points, 1L)
  expect_true("size" %in% names(points[[1L]]$mapping))
  expect_null(.layers_of(p, "GeomPolygon")[-1L][1L][[1L]])
  size_scale <- p$scales$get_scales("size")
  expect_identical(size_scale$name, "Trials with the event")
  built <- ggplot2::ggplot_build(p)
  expect_identical(size_scale$get_labels(),
                   .thg_comma(size_scale$get_breaks()))
  expect_null(p$labels$caption)
  # each label sits above its own dot, higher for larger dots
  texts <- .layers_of(p, "GeomText")
  expect_length(texts, 9L)
  value <- sizes$trials[match(hg$nodes, sizes$state)]
  expect_equal(texts[[9L]]$data$vjust, -0.4 - 1.8 * sqrt(value / max(value)))
  # a generic legend title without a unit or size_title
  unitless <- group_hypergraph(.blob_members()[c("group", "state")],
                               actor = "state", group = "group")
  expect_identical(plot(unitless, node_sizes = sizes)$scales$get_scales("size")$name,
                   "Event value")
  expect_identical(.thg_comma(c(1000, 2.5, 1e6)), c("1,000", "2.5", "1,000,000"))
  expect_identical(.thg_comma(c(3667, 95.5, 0.12345), count = TRUE),
                   c("3,667", "95.5", "0.123"))
})

test_that("the count attribute colours, titles and names the unit by itself", {
  hg <- .blob_hg()
  short <- plot(hg)
  long <- plot(hg, color_by = "trials", titles = "trials", unit = "trials")
  expect_identical(ggplot2::ggplot_build(short)$data,
                   ggplot2::ggplot_build(long)$data)
  expect_identical(short$scales$get_scales("fill")$name, "Trials")
  expect_identical(.layers_of(short, "GeomLabel")[[1L]]$data$label,
                   c("Hint\n1,200 trials", "Question\n800 trials",
                     "Exit\n95.5 trials"))
  # overrides win
  expect_identical(plot(hg, legend_title = "n")$scales$get_scales("fill")$name, "n")
  expect_identical(plot(hg, unit = "steps")$scales$get_scales("fill")$name, "Steps")
  expect_length(.layers_of(plot(hg, titles = FALSE), "GeomLabel"), 0L)
  # a numeric attribute that is the same on every hyperedge gives way to the
  # one that varies; two that vary leave the choice to the call
  with_share <- .blob_members()
  with_share$share <- 1
  expect_identical(
    ggplot2::ggplot_build(plot(.blob_hg(with_share)))$data,
    ggplot2::ggplot_build(short)$data)
  with_share$share <- rep(c(0.5, 0.7, 0.9), each = 3L)
  undecided <- plot(.blob_hg(with_share))
  expect_length(.layers_of(undecided, "GeomLabel"), 0L)
  expect_length(unique(undecided$layers[[1L]]$data$colour), 1L)
  # no titles on panels that are titled already
  expect_s3_class(plot(hg, dismantled = TRUE), "ggplot")
})

test_that("the default look: legends below, margins, halo, alpha 0.5", {
  hg <- .blob_hg()
  p <- plot(hg)
  expect_identical(p$theme$legend.position, "bottom")
  expect_equal(as.numeric(p$theme$plot.margin), c(40, 130, 30, 130))
  expect_equal(as.numeric(p$theme$legend.key.width), 2.5)
  expect_null(p$theme$legend.box)
  # the wide keys are for the colour bar; a discrete legend keeps default
  # key widths
  by_group <- plot(hg, color_by = c(Hint = "a", Question = "b", Exit = "c"))
  guide <- by_group$scales$get_scales("fill")$guide
  expect_equal(as.numeric(guide$params$theme$legend.key.width), 1.2)
  # haloed bold labels at size 4.2, pebbles at alpha 0.5
  texts <- .layers_of(p, "GeomText")
  expect_length(texts, 9L)
  expect_identical(texts[[9L]]$aes_params$fontface, "bold")
  expect_identical(texts[[9L]]$aes_params$size, 4.2)
  expect_identical(p$layers[[1L]]$aes_params$alpha, 0.5)
  # explicit arguments
  chosen <- plot(hg, alpha = 0.2, label_size = 3)
  expect_identical(chosen$layers[[1L]]$aes_params$alpha, 0.2)
  expect_identical(.layers_of(chosen, "GeomText")[[9L]]$aes_params$size, 3)
  # three legends stack
  moves <- data.frame(from = c("Wrong", "Hint"), to = c("Hint", "Retry"),
                      trials = c(10, 5))
  stacked <- plot(hg, transitions = moves)
  expect_identical(stacked$theme$legend.box, "vertical")
})

test_that("pieces = \"row\" sets disconnected pieces side by side", {
  members <- data.frame(
    state = c("a", "b", "c", "b", "c", "d", "x", "y", "z"),
    group = c("g1", "g1", "g1", "g2", "g2", "g2", "small", "small", "lone"),
    trials = c(9, 9, 9, 8, 8, 8, 20, 20, 1)
  )
  hg <- group_hypergraph(members, actor = "state", group = "group")
  precedence <- order(-c(g1 = 9, g2 = 8, lone = 1, small = 20)[colnames(hg$incidence)],
                      colnames(hg$incidence))
  row <- .thg_row_layout(hg, seed = 1L, center = NULL, padding = 0.045,
                         precedence = precedence)
  # most prominent piece first (x, y), then the g1/g2 piece, then the lone node
  expect_identical(sort(row$node), sort(hg$nodes))
  first <- row[row$node %in% c("x", "y"), ]
  second <- row[row$node %in% c("a", "b", "c", "d"), ]
  lone <- row[row$node == "z", ]
  expect_true(all(first$x <= 1 + 1e-9))
  expect_true(all(second$x >= 1.4 - 1e-9 & second$x <= 2.4 + 1e-9))
  expect_equal(c(lone$x, lone$y), c(0.5 + 2.8, 0.5))
  # each piece is laid out exactly as it is when drawn alone
  alone <- group_hypergraph(members[members$group %in% c("g1", "g2"), ],
                            actor = "state", group = "group")
  own <- .thg_positions(alone, "bipartite", 1L, NULL, 0.045)$nodes
  expect_equal(second$x[match(own$node, second$node)] - 1.4, own$x)
  expect_equal(second$y[match(own$node, second$node)], own$y)
  # the row is the default; "packed" packs; a connected hypergraph is the same
  p_row <- plot(hg, pieces = "row")
  p_default <- plot(hg)
  expect_equal(p_row$layers[[1L]]$data, p_default$layers[[1L]]$data)
  p_packed <- plot(hg, pieces = "packed")
  expect_false(isTRUE(all.equal(p_row$layers[[2L]]$data,
                                p_packed$layers[[2L]]$data)))
  connected <- .blob_hg()
  expect_identical(plot(connected, pieces = "row")$layers[[2L]]$data,
                   plot(connected)$layers[[2L]]$data)
})

test_that("the blob arguments reject bad input with hypernets_bad_input", {
  hg <- .blob_hg()
  hg$edge_data$kind <- c("a", "b", "c")
  expect_error(plot(hg, titles = "kind"), class = "hypernets_bad_input")
  expect_error(plot(hg, titles = "nope"), class = "hypernets_bad_input")
  expect_error(plot(hg, titles = NA), class = "hypernets_bad_input")
  expect_error(plot(hg, titles = FALSE, notes = c(Hint = "x")),
               class = "hypernets_bad_input")
  expect_error(plot(hg, titles = TRUE, notes = data.frame(group = "Hint")),
               class = "hypernets_bad_input")
  expect_error(plot(hg, titles = TRUE, notes = c("unnamed")),
               class = "hypernets_bad_input")
  expect_error(plot(hg, titles = TRUE, dismantled = TRUE),
               class = "hypernets_bad_input")
  expect_error(plot(hg, unit = 1), class = "hypernets_bad_input")
  expect_error(plot(hg, unit = c("a", "b")), class = "hypernets_bad_input")
  expect_error(plot(hg, title_prefix = 1), class = "hypernets_bad_input")
  expect_error(plot(hg, layout = "circle", pieces = "row"),
               class = "hypernets_bad_input")
  expect_error(plot(hg, pieces = "column"))
})
