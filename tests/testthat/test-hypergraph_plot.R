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
  discrete <- hypernets:::.thg_fill_colours(c("x", "y", "x"))
  expect_identical(discrete$colours, unname(discrete$palette[c("x", "y", "x")]))
  expect_identical(unname(discrete$palette), hypernets:::.thg_okabe_ito[1:2])
  numeric <- hypernets:::.thg_fill_colours(c(2, 12, 7))
  expect_null(numeric$palette)
  expect_identical(numeric$colours[1L], toupper("#F0E442"))
  expect_identical(numeric$colours[2L], toupper("#0072B2"))
  constant <- hypernets:::.thg_fill_colours(c(3, 3))
  expect_identical(constant$colours[1L], constant$colours[2L])
})

test_that("plot.net_hypergraph draws one hull per drawable hyperedge, largest first", {
  hg <- .plot_fixture()
  grDevices::pdf(NULL)
  before <- grDevices::dev.cur()
  p <- plot(hg)
  expect_identical(grDevices::dev.cur(), before)
  grDevices::dev.off()
  expect_s3_class(p, "ggplot")
  polygons <- Filter(function(l) inherits(l$geom, "GeomPolygon"), p$layers)
  expect_length(polygons, 1L)
  # e4 is a singleton and gets no hull; e1 and e2 (size 3) precede e3 (size 2)
  expect_identical(levels(polygons[[1L]]$data$hyperedge), c("e1", "e2", "e3"))
})

test_that("a hull is the member hull widened by exactly the padding", {
  xs <- c(0, 1, 0.5)
  ys <- c(0, 0, 1)
  ring <- hypernets:::.thg_hull(xs, ys, radius = 0.1, n_arc = 720L)
  # every member lies inside, at least `radius` from the outline
  inside <- vapply(seq_along(xs), function(i) {
    min(sqrt((ring$x - xs[i])^2 + (ring$y - ys[i])^2))
  }, numeric(1L))
  expect_true(all(inside >= 0.1 - 1e-9))
  # and no outline point is further than `radius` from its nearest member
  reach <- vapply(seq_len(nrow(ring)), function(j) {
    min(sqrt((ring$x[j] - xs)^2 + (ring$y[j] - ys)^2))
  }, numeric(1L))
  expect_equal(max(reach), 0.1, tolerance = 1e-9)
  # two members give a stadium: width across is twice the radius
  stadium <- hypernets:::.thg_hull(c(0, 1), c(0, 0), radius = 0.2, n_arc = 720L)
  expect_equal(diff(range(stadium$y)), 0.4, tolerance = 1e-6)
})

test_that("the bipartite layout places nodes and hyperedges in one frame", {
  hg <- .plot_fixture()
  pos <- hypernets:::.thg_positions(hg, "bipartite", seed = 2L)
  expect_identical(names(pos$nodes), c("node", "x", "y"))
  expect_identical(names(pos$edges), c("hyperedge", "x", "y"))
  expect_identical(pos$edges$hyperedge, colnames(hg$incidence))
  both <- rbind(pos$nodes[c("x", "y")], pos$edges[c("x", "y")])
  expect_true(all(both >= 0 & both <= 1))
  # one factor scales both axes, so the longer axis spans exactly [0, 1]
  expect_equal(max(diff(range(both$x)), diff(range(both$y))), 1)
  # a projection layout puts each hyperedge at the centroid of its members
  spring <- hypernets:::.thg_positions(hg, "spring", seed = 2L)
  e1 <- hg$hyperedges[[1L]]
  expect_equal(spring$edges$x[1L], mean(spring$nodes$x[e1]))
})

test_that("dismantled panels ink exactly the members of their hyperedge", {
  hg <- .plot_fixture()
  p <- plot(hg, dismantled = TRUE)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  text_layer <- which(vapply(p$layers, function(l) inherits(l$geom, "GeomText"),
                             logical(1L)))
  labelled <- built$data[[text_layer]]
  panel_names <- levels(p$layers[[1L]]$data$hyperedge)
  shown <- split(labelled$label, panel_names[as.integer(labelled$PANEL)])
  expected <- lapply(stats::setNames(hg$hyperedges, colnames(hg$incidence)),
                     function(idx) hg$nodes[idx])
  lapply(names(shown), function(e) expect_setequal(shown[[e]], expected[[e]]))
  expect_error(plot(hg, dismantled = TRUE, edge_labels = TRUE),
               class = "hypernets_bad_input")
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
  pos <- hypernets:::.thg_layout(hg, "spring", seed = 3L)
  expect_identical(names(pos), c("node", "x", "y"))
  expect_identical(pos$node, hg$nodes)
  expect_true(all(pos$x >= 0 & pos$x <= 1 & pos$y >= 0 & pos$y <= 1))
  again <- hypernets:::.thg_layout(hg, pos, seed = 1L)
  expect_equal(again$x, pos$x)
  layout_plot <- plot(hg, layout = pos)
  expect_s3_class(layout_plot, "ggplot")
  expect_error(plot(hg, layout = pos[-1, ]), class = "hypernets_bad_input")
})

test_that("plot.net_hypergraph rejects bad selectors", {
  hg <- .plot_fixture()
  expect_error(plot(hg, color_by = "nope"), class = "hypernets_bad_input")
  expect_error(plot(hg, color_by = 1:3), class = "hypernets_bad_input")
  expect_error(plot(hg, color_by = c(e1 = 1)), class = "hypernets_bad_input")
  expect_error(plot(hg, labels = c("A", "B")), class = "hypernets_bad_input")
  expect_error(plot(hg, alpha = 2))
  singletons <- group_hypergraph(
    data.frame(member = c("a", "b"), event = c("e1", "e2")), "member", "event"
  )
  expect_error(plot(singletons), class = "hypernets_bad_input")
})

test_that("one variable for colour and line type gives one merged legend", {
  hg <- .plot_fixture()
  merged <- plot(hg, color_by = "sector", linetype_by = "sector")
  # ggplot2 merges legends whose scales share a title and levels
  expect_identical(merged$scales$get_scales("fill")$name, "sector")
  expect_identical(merged$scales$get_scales("linetype")$name, "sector")
  apart <- plot(hg, color_by = "sector", linetype_by = c(e1 = "p", e2 = "q",
                                                         e3 = "p", e4 = "q"))
  expect_identical(apart$scales$get_scales("linetype")$name, "line type")
})

test_that("a pebble keeps every member inside with room, and stays convex", {
  set.seed(11)
  point_sets <- replicate(40, {
    k <- sample(2:7, 1L)
    list(x = stats::runif(k), y = stats::runif(k) * stats::runif(1L, 0.05, 1))
  }, simplify = FALSE)
  radius <- 0.045
  checks <- vapply(point_sets, function(pts) {
    vapply(c(2, 5, 8), function(detail) {
      shape <- hypernets:::.thg_pebble(pts$x, pts$y, radius, detail)
      # inside: each member lies on the inner side of every outline edge
      ex <- c(shape$x[-1L], shape$x[1L]) - shape$x
      ey <- c(shape$y[-1L], shape$y[1L]) - shape$y
      orientation <- sign(sum(shape$x * c(shape$y[-1L], shape$y[1L]) -
                                c(shape$x[-1L], shape$x[1L]) * shape$y))
      inside <- all(vapply(seq_along(pts$x), function(i) {
        all(orientation * (ex * (pts$y[i] - shape$y) -
                             ey * (pts$x[i] - shape$x)) >= -1e-9)
      }, logical(1L)))
      # distance to the outline's edges (point-to-segment), not its vertices
      room <- min(vapply(seq_along(pts$x), function(i) {
        along <- pmin(1, pmax(0, ((pts$x[i] - shape$x) * ex +
                                    (pts$y[i] - shape$y) * ey) / (ex^2 + ey^2)))
        min(sqrt((shape$x + along * ex - pts$x[i])^2 +
                   (shape$y + along * ey - pts$y[i])^2))
      }, numeric(1L)))
      # convex, hence no ripple: consecutive edges always turn the same way
      turns <- ex * c(ey[-1L], ey[1L]) - ey * c(ex[-1L], ex[1L])
      inside && room >= 0.6 * radius - 1e-6 &&
        all(orientation * turns >= -1e-12)
    }, logical(1L))
  }, logical(3L))
  expect_true(all(checks))
})

test_that("detail = Inf is the unsmoothed rounded hull", {
  xs <- c(0, 1, 0.4)
  ys <- c(0, 0.2, 0.9)
  expect_identical(hypernets:::.thg_pebble(xs, ys, 0.05, Inf),
                   hypernets:::.thg_hull(xs, ys, 0.05, n_arc = 48L))
})

test_that("outline = 'fill' maps the outline, a colour fixes it", {
  hg <- .plot_fixture()
  seams <- plot(hg, color_by = "sector")
  expect_identical(seams$layers[[1L]]$aes_params$colour, "white")
  expect_null(seams$scales$get_scales("colour"))
  own <- plot(hg, color_by = "sector", outline = "fill", alpha = 0.15)
  expect_identical(own$scales$get_scales("colour")$name, "sector")
  expect_error(plot(hg, detail = 0))
  expect_error(plot(hg, outline = c("white", "black")))
})

test_that("every hyperedge label lies inside its own pebble", {
  # two-member hyperedges put their bipartite vertex off the line between
  # the members; the label must fall back inside the shape
  memberships <- data.frame(
    researcher = c("a", "b", "c", "a", "d", "d", "e", "b", "e", "f", "f", "g"),
    project = c("p1", "p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4", "p5",
                "p6", "p6")
  )
  hg <- group_hypergraph(memberships, actor = "researcher", group = "project")
  inside <- vapply(1:6, function(seed) {
    p <- plot(hg, edge_labels = TRUE, seed = seed)
    shapes <- p$layers[[1L]]$data
    labels <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
    anchors <- labels[[1L]]$data
    drawn <- levels(shapes$hyperedge)
    all(vapply(seq_along(drawn), function(i) {
      ring <- shapes[shapes$hyperedge == drawn[i], , drop = FALSE]
      hypernets:::.thg_in_convex(ring$x, ring$y, anchors$x[i], anchors$y[i])
    }, logical(1L)))
  }, logical(1L))
  expect_true(all(inside))
  expect_true(hypernets:::.thg_in_convex(c(0, 1, 1, 0), c(0, 0, 1, 1), 0.5, 0.5))
  expect_false(hypernets:::.thg_in_convex(c(0, 1, 1, 0), c(0, 0, 1, 1), 1.5, 0.5))
})

test_that("a whole-number colour scale gets whole-number legend breaks", {
  hg <- .plot_fixture()
  p <- plot(hg, color_by = "size")
  # breaks exist only once the scale has been trained on the data
  breaks <- ggplot2::ggplot_build(p)$plot$scales$get_scales("fill")$get_breaks()
  breaks <- breaks[!is.na(breaks)]
  expect_true(length(breaks) >= 2L)
  expect_true(all(abs(breaks - round(breaks)) < 1e-12))
  fractional <- plot(hg, color_by = c(e1 = 0.1, e2 = 0.35, e3 = 0.6, e4 = 0.9))
  expect_s3_class(fractional$scales$get_scales("fill")$breaks, "waiver")
})

test_that("the force-directed layout is seeded, local, and pulls neighbours in", {
  ring <- matrix(0, 12, 12)
  ring[cbind(1:12, c(2:12, 1L))] <- 1
  ring <- ring + t(ring)
  set.seed(99)
  before <- stats::runif(1L)
  set.seed(99)
  first <- hypernets:::.thg_layout_fr(ring, seed = 4L)
  # the caller's stream is untouched: the next draw is the one it would be
  expect_identical(stats::runif(1L), before)
  expect_identical(hypernets:::.thg_layout_fr(ring, seed = 4L), first)
  expect_false(isTRUE(all.equal(hypernets:::.thg_layout_fr(ring, seed = 5L),
                                first)))
  gaps <- as.matrix(stats::dist(first))
  expect_lt(mean(gaps[ring == 1]), mean(gaps[ring == 0 & row(gaps) != col(gaps)]))
})

test_that("center pins the named nodes in the middle of every layout", {
  # a flower: every hyperedge shares the two outcome nodes
  memberships <- data.frame(
    event = c("out1", "out2", "a1", "a2", "out1", "out2", "b1", "b2",
              "out1", "out2", "c1", "c2", "out1", "out2", "d1", "d2", "d3"),
    step = rep(c("s1", "s2", "s3", "s4"), c(4, 4, 4, 5))
  )
  hg <- group_hypergraph(memberships, actor = "event", group = "step")
  # two leaves that would sit at the rim if left alone; the shared outcome
  # nodes land in the middle without help, so they would prove nothing
  middle <- c("a1", "d3")
  # "in the middle": every other node is further from the centred nodes than
  # they are from each other's centroid, and the other nodes surround them
  # (their centroid lies inside the other nodes' convex hull)
  surrounded <- function(layout, seed) {
    pos <- hypernets:::.thg_positions(hg, layout, seed, center = middle)$nodes
    inner <- pos$node %in% middle
    hub <- c(mean(pos$x[inner]), mean(pos$y[inner]))
    gap <- sqrt((pos$x - hub[1L])^2 + (pos$y - hub[2L])^2)
    outer <- pos[!inner, , drop = FALSE]
    hull <- grDevices::chull(outer$x, outer$y)
    max(gap[inner]) < min(gap[!inner]) &&
      hypernets:::.thg_in_convex(outer$x[hull], outer$y[hull], hub[1L], hub[2L])
  }
  closeness <- surrounded
  expect_true(all(vapply(1:5, function(s) closeness("bipartite", s), logical(1L))))
  expect_true(all(vapply(1:5, function(s) closeness("spring", s), logical(1L))))
  expect_true(closeness("circle", 1L))
  expect_s3_class(plot(hg, center = middle), "ggplot")
  expect_error(plot(hg, center = "nope"), class = "hypernets_bad_input")
  # a name missing from this hypergraph is skipped when another one matches
  expect_s3_class(plot(hg, center = c("a1", "nope")), "ggplot")
  expect_error(plot(hg, center = 1), class = "hypernets_bad_input")
  table_layout <- hypernets:::.thg_layout(hg, "spring", 1L)
  expect_error(plot(hg, layout = table_layout, center = middle),
               class = "hypernets_bad_input")
})

test_that("pinned vertices never move in the force-directed layout", {
  star <- matrix(0, 6, 6)
  star[1L, 2:6] <- 1
  star <- star + t(star)
  pinned <- c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE)
  xy <- hypernets:::.thg_layout_fr(star, seed = 3L, pinned = pinned)
  ring <- hypernets:::.thg_ring(pinned, radius = sqrt(6) / 4)
  expect_equal(xy[pinned, , drop = FALSE], ring[pinned, , drop = FALSE])
})
