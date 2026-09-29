testthat::skip_on_cran()

# Node overlays of plot.net_hypergraph(): circles sized by a value, direction
# triangles, transition arrows, haloed labels.

.nodes_fixture <- function() {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e", "f"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e4")
  )
  group_hypergraph(dat, "member", "event")
}

# Four nodes on the corners of a unit square, so every heading is known:
# a (0, 0), b (1, 0), c (0, 1), d (1, 1).
.square <- function() {
  hg <- group_hypergraph(
    data.frame(member = c("a", "b", "c", "b", "c", "d"),
               event = c("h1", "h1", "h1", "h2", "h2", "h2")),
    "member", "event")
  list(hg = hg,
       layout = data.frame(node = c("a", "b", "c", "d"), x = c(0, 1, 0, 1),
                           y = c(0, 0, 1, 1)))
}

.polygon_layer <- function(p, fill) {
  hits <- Filter(function(l) inherits(l$geom, "GeomPolygon") &&
                   identical(l$aes_params$fill, fill), p$layers)
  if (length(hits)) hits[[1L]]$data else NULL
}

# Radius of each drawn disc, in the hypergraph's node order.
.disc_radius <- function(p, fill = "#000000") {
  disc <- .polygon_layer(p, fill)
  unname(vapply(split(disc, disc$id), function(d) {
    max(sqrt((d$x - mean(d$x[-1L]))^2 + (d$y - mean(d$y[-1L]))^2))
  }, numeric(1L)))
}

test_that("default calls draw exactly what they drew before the overlays", {
  baseline <- readRDS(test_path("fixtures", "plot_hypergraph_baseline.rds"))
  skip_if_not(identical(baseline$ggplot2,
                        as.character(utils::packageVersion("ggplot2"))),
              "baseline was recorded under another ggplot2 version")
  hg <- .nodes_fixture()
  hg$edge_data <- data.frame(edge = c("e1", "e2", "e3", "e4"),
                             sector = c("x", "y", "x", "z"))
  calls <- list(
    default = function() plot(hg),
    size = function() plot(hg, color_by = "size", edge_labels = TRUE),
    sector = function() plot(hg, color_by = "sector", linetype_by = "sector",
                             outline = "fill", alpha = 0.15),
    nolabel = function() plot(hg, labels = FALSE, detail = Inf),
    circle = function() plot(hg, layout = "circle", center = c("b", "c"),
                             edge_labels = TRUE),
    spring = function() plot(hg, layout = "spring", labels = c(a = "A")),
    dismantled = function() plot(hg, dismantled = TRUE)
  )
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  now <- lapply(calls, function(f) {
    lapply(ggplot2::ggplot_build(f())$data, as.data.frame)
  })
  expect_identical(now, baseline$data)
})

# a direction table that draws circles; a -> b so every node has a heading
# or none
.moves_fixture <- function() {
  data.frame(from = c("a", "b", "c", "d", "d", "e"),
             to = c("b", "c", "a", "e", "b", "f"),
             weight = c(20, 9, 4, 8, 3, 1))
}

test_that("circle area follows the value, with a floor for rare nodes", {
  hg <- .nodes_fixture()
  counts <- c(a = 40, b = 25, c = 12, d = 30, e = 5, f = 0.01)
  moves <- .moves_fixture()
  p <- plot(hg, node_sizes = counts, direction = moves)
  expect_s3_class(p, "ggplot")
  radius <- .disc_radius(p)
  value <- unname(counts[hg$nodes])
  r_max <- max(radius)
  # monotone in the value
  expect_identical(order(radius), order(value))
  # area proportional to the value above the floor
  above <- radius > 0.3 * r_max + 1e-9
  expect_equal(radius[above]^2 / r_max^2, value[above] / max(value),
               tolerance = 1e-3)
  # the rare node sits on the floor
  expect_equal(radius[hg$nodes == "f"], 0.3 * r_max, tolerance = 1e-3)
  # a data.frame gives the same circles; extra names are skipped
  as_table <- plot(hg, node_sizes = data.frame(node = c(names(counts), "zz"),
                                               value = c(counts, 99)),
                   direction = moves)
  expect_equal(.disc_radius(as_table), radius)
  # a (state, trials) table is read by its one text and one numeric column
  as_pipeline <- plot(hg, node_sizes = data.frame(state = names(counts),
                                                  trials = unname(counts)),
                      direction = moves)
  expect_equal(.disc_radius(as_pipeline), radius)
  expect_match(p$labels$caption, "proportional to the event value",
               fixed = TRUE)
  titled <- plot(hg, node_sizes = counts, direction = moves,
                 size_title = "trials with the event")
  expect_identical(gsub("\n", " ", titled$labels$caption, fixed = TRUE),
                   "Circle area: trials with the event. Triangle: points to the event that most often follows it.")
  by_unit <- plot(hg, node_sizes = counts, direction = moves, unit = "trials")
  expect_identical(gsub("\n", " ", by_unit$labels$caption, fixed = TRUE),
                   "Circle area: trials with the event. Triangle: points to the event that most often follows it.")
  # a halo of eight white stamps and the bold label on top
  texts <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
  expect_length(texts, 9L)
  expect_identical(unname(vapply(texts, function(l) l$aes_params$colour,
                                 character(1L))), c(rep("white", 8L), "#2B2B2B"))
  expect_true(all(texts[[9L]]$data$vjust == -0.35))
})

test_that("a triangle points at the most frequent other successor", {
  sq <- .square()
  sizes <- c(a = 10, b = 10, c = 10, d = 10)
  moves <- data.frame(
    from = c("a", "a", "a", "b", "b", "c", "c"),
    to = c("b", "c", "a", "d", "a", "b", "b"),
    weight = c(5, 9, 100, 3, 3, 2, 2)
  )
  # a: c beats b, and its self-loop (100) is ignored -> heading pi / 2
  # b: a and d tie at 3 -> a by name -> heading pi
  # c: b twice, summed to 4 -> heading -pi / 4
  # d: a sink -> no triangle
  geometry <- hypernets:::.thg_node_geometry(
    data.frame(x = sq$layout$x, y = sq$layout$y), sq$hg$nodes,
    unname(sizes[sq$hg$nodes]),
    hypernets:::.thg_move_table(moves, "direction", sq$hg$nodes))
  expect_identical(geometry$target, c("c", "a", "b", NA))
  expect_equal(geometry$heading[1:3], c(pi / 2, pi, -pi / 4))

  p <- plot(sq$hg, layout = sq$layout, node_sizes = sizes, direction = moves)
  tri <- .polygon_layer(p, "#1F3A44")
  expect_setequal(unique(tri$id), 1:3)
  apex <- tri[!duplicated(tri$id), ]
  r <- geometry$radius[apex$id]
  expect_equal(apex$x, sq$layout$x[apex$id] + 0.95 * r * cos(geometry$heading[apex$id]))
  expect_equal(apex$y, sq$layout$y[apex$id] + 0.95 * r * sin(geometry$heading[apex$id]))
  expect_match(gsub("\n", " ", p$labels$caption, fixed = TRUE),
               "Triangle: points to the event that most often follows it.",
               fixed = TRUE)

  # ties are broken the same way whatever the row order
  shuffled <- plot(sq$hg, layout = sq$layout, node_sizes = sizes,
                   direction = moves[rev(seq_len(nrow(moves))), ])
  expect_identical(.polygon_layer(shuffled, "#1F3A44"), tri)
})

test_that("outside tips take the node colour and push an up-label below", {
  sq <- .square()
  sizes <- c(a = 10, b = 4, c = 10, d = 10)
  moves <- data.frame(from = c("a", "b"), to = c("c", "a"), weight = c(1, 1))
  p <- plot(sq$hg, layout = sq$layout, node_sizes = sizes, direction = moves,
            arrow_style = "outside", node_fill = "#0072B2")
  expect_false(is.null(.polygon_layer(p, "#0072B2")))
  labels <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)[[9L]]$data
  geometry <- hypernets:::.thg_node_geometry(
    data.frame(x = sq$layout$x, y = sq$layout$y), sq$hg$nodes,
    unname(sizes[sq$hg$nodes]),
    hypernets:::.thg_move_table(moves, "direction", sq$hg$nodes), "outside")
  # a points up at c: label below; the others above
  expect_identical(labels$vjust, c(1.35, -0.35, -0.35, -0.35))
  expect_equal(labels$y, sq$layout$y + c(-1, 1, 1, 1) * geometry$radius)
  custom <- plot(sq$hg, layout = sq$layout, node_sizes = sizes,
                 direction = moves, arrow_fill = "#D55E00")
  expect_false(is.null(.polygon_layer(custom, "#D55E00")))
})

test_that("transitions are curved arrows with width by weight and self-loops", {
  hg <- .nodes_fixture()
  counts <- c(a = 40, b = 25, c = 12, d = 30, e = 5, f = 2)
  moves <- data.frame(from = c("a", "b", "c", "a"), to = c("b", "c", "a", "a"),
                      weight = c(20, 9, 4, 5))
  p <- plot(hg, node_sizes = counts, direction = moves, transitions = moves)
  curves <- Filter(function(l) inherits(l$geom, "GeomCurve"), p$layers)
  expect_length(curves, 2L)
  expect_identical(curves[[1L]]$data$weight, c(20, 9, 4)[order(c("a", "b", "c"))])
  expect_identical(curves[[2L]]$data$from, "a")
  expect_true(p$scales$has_scale("linewidth"))
  expect_identical(p$scales$get_scales("linewidth")$name, "Transition weight")
  by_unit <- plot(hg, transitions = moves, unit = "trials")
  expect_identical(by_unit$scales$get_scales("linewidth")$name,
                   "Trials making the transition")
  # a weight column under another name is read as the weight
  renamed <- plot(hg, transitions = data.frame(from = moves$from, to = moves$to,
                                               trials = moves$weight))
  expect_identical(
    Filter(function(l) inherits(l$geom, "GeomCurve"), renamed$layers)[[1L]]$data$weight,
    curves[[1L]]$data$weight)
  # arrows sit under the circles
  kinds <- vapply(p$layers, function(l) class(l$geom)[1L], character(1L))
  expect_lt(max(which(kinds == "GeomCurve")),
            which(kinds == "GeomPolygon")[2L])
  # without circles the arrows still draw over plain points
  plain <- plot(hg, transitions = moves)
  expect_length(Filter(function(l) inherits(l$geom, "GeomCurve"),
                       plain$layers), 2L)
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("node overlays reject bad input with hypernets_bad_input", {
  hg <- .nodes_fixture()
  counts <- c(a = 40, b = 25, c = 12, d = 30, e = 5, f = 2)
  moves <- data.frame(from = "a", to = "b", weight = 1)
  expect_error(plot(hg, direction = moves), class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts[-1L]),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = replace(counts, 1L, -1)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts * 0),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = unname(counts)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = data.frame(name = "a", size = 1, n = 2)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = data.frame(node = "a", size = 1)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts,
                    direction = data.frame(from = "a", to = "b", count = 1,
                                           n = 2)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts,
                    direction = data.frame(from = "a", to = "b")),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts,
                    direction = data.frame(from = "a", to = "b", weight = -1)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, transitions = data.frame(from = "x", to = "y",
                                                 weight = 1)),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts, dismantled = TRUE),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts, node_fill = "nope"),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts, arrow_fill = c("red", "blue")),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts, size_title = 1),
               class = "hypernets_bad_input")
  expect_error(plot(hg, node_sizes = counts, arrow_style = "sideways"))
})
