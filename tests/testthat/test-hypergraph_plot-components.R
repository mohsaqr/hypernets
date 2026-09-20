testthat::skip_on_cran()

# A connected core of overlapping hyperedges plus hyperedges sharing no member
# with anything else -- the shape that made the reports unreadable.
.island_fixture <- function() {
  group_hypergraph(
    data.frame(
      member = c("a", "b", "c", "b", "c", "d", "c", "d", "e",
                 "i1", "i2", "i3", "j1", "j2", "j3", "k1", "k2"),
      event = c(rep(c("E1", "E2", "E3"), each = 3), rep("I1", 3),
                rep("I2", 3), rep("K1", 2)),
      stringsAsFactors = FALSE
    ),
    "member", "event"
  )
}

# A dense core with one hyperedge bridging out to a satellite cluster: the
# pendant that unbounded repulsion throws off the edge of the picture.
.bridge_fixture <- function() {
  group_hypergraph(
    data.frame(
      member = c("c1", "c2", "c3", "c4", "c2", "c3", "c4", "c5",
                 "c3", "c4", "c5", "c6", "c1", "c5", "c6", "c2",
                 "c4", "s1", "s1", "s2", "s3", "s4"),
      event = rep(c("C1", "C2", "C3", "C4", "BRIDGE", "SAT"),
                  c(4, 4, 4, 4, 2, 4)),
      stringsAsFactors = FALSE
    ),
    "member", "event"
  )
}

.frame_share <- function(xy, rows) {
  part <- xy[rows, , drop = FALSE]
  apply(part, 2L, function(z) diff(range(z))) /
    max(apply(xy, 2L, function(z) diff(range(z))))
}

test_that(".thg_components labels the connected components 1..K", {
  # two triangles and a lone vertex
  adjacency <- matrix(0, 7L, 7L)
  adjacency[cbind(c(1, 2, 3, 4, 5, 6), c(2, 3, 1, 5, 6, 4))] <- 1
  adjacency <- adjacency + t(adjacency)
  labels <- hypernets:::.thg_components(adjacency)
  expect_identical(labels, c(1L, 1L, 1L, 2L, 2L, 2L, 3L))
  # a path still settles, however long it takes to propagate
  path <- matrix(0, 6L, 6L)
  path[cbind(1:5, 2:6)] <- 1
  path <- path + t(path)
  expect_identical(hypernets:::.thg_components(path), rep(1L, 6L))
  expect_identical(hypernets:::.thg_components(matrix(0, 1L, 1L)), 1L)
})

test_that("shelf packing never overlaps two component boxes", {
  size <- cbind(c(4, 1, 3, 2, 0.5), c(3, 5, 1, 2, 0.5))
  gap <- 1.5
  offset <- hypernets:::.thg_shelf_offsets(size, gap)
  pairs <- utils::combn(nrow(size), 2L)
  clear <- apply(pairs, 2L, function(ij) {
    a <- ij[1L]
    b <- ij[2L]
    # boxes are disjoint when they are apart on either axis
    apart <- vapply(1:2, function(axis) {
      offset[a, axis] + size[a, axis] <= offset[b, axis] ||
        offset[b, axis] + size[b, axis] <= offset[a, axis]
    }, logical(1L))
    any(apart)
  })
  expect_true(all(clear))
})

test_that("a disconnected hypergraph spends its frame on the structure", {
  hg <- .island_fixture()
  core <- match(c("a", "b", "c", "d", "e"), hg$nodes)
  share <- vapply(1:20, function(seed) {
    pos <- hypernets:::.thg_positions(hg, "bipartite", seed = seed)
    .frame_share(cbind(pos$nodes$x, pos$nodes$y), core)
  }, numeric(2L))
  # Laying the whole graph out at once gave the connected core 0.05 of the
  # frame on average and 0.02 at worst, the islands having set the scale.
  # Packing the components gives it 0.43 on average and 0.24 at worst.
  expect_gt(mean(share), 0.25)
  expect_gt(min(share), 0.10)
})

test_that("packed components keep clear of one another", {
  hg <- .island_fixture()
  pos <- hypernets:::.thg_positions(hg, "bipartite", seed = 1L)
  xy <- cbind(pos$nodes$x, pos$nodes$y)
  parts <- list(match(c("a", "b", "c", "d", "e"), hg$nodes),
                match(c("i1", "i2", "i3"), hg$nodes),
                match(c("j1", "j2", "j3"), hg$nodes),
                match(c("k1", "k2"), hg$nodes))
  # no node of one component may fall inside another component's pebble: an
  # overlap between components would claim a membership that is not there
  intruders <- sum(vapply(seq_along(parts), function(k) {
    rows <- parts[[k]]
    ring <- hypernets:::.thg_pebble(xy[rows, 1L], xy[rows, 2L], 0.045,
                                   detail = 5)
    outside <- setdiff(seq_len(hg$n_nodes), rows)
    sum(vapply(outside, function(i)
      hypernets:::.thg_in_convex(ring$x, ring$y, xy[i, 1L], xy[i, 2L]),
      logical(1L)))
  }, numeric(1L)))
  expect_identical(intruders, 0)
})

test_that("components are packed into a frame worth drawing on", {
  hg <- .island_fixture()
  shapes <- vapply(1:10, function(seed) {
    pos <- hypernets:::.thg_positions(hg, "bipartite", seed = seed)
    xy <- cbind(pos$nodes$x, pos$nodes$y)
    span <- apply(xy, 2L, function(z) diff(range(z)))
    min(span) / max(span)
  }, numeric(1L))
  # a packed frame is roughly square, never a thin ribbon
  expect_true(all(shapes > 0.3))
})

test_that("bounded repulsion keeps a satellite cluster in the picture", {
  hg <- .bridge_fixture()
  bridge <- which(colnames(hg$incidence) == "BRIDGE")
  spans <- vapply(1:10, function(seed) {
    pos <- hypernets:::.thg_positions(hg, "bipartite", seed = seed)
    xy <- cbind(pos$nodes$x, pos$nodes$y)
    frame <- max(apply(xy, 2L, function(z) diff(range(z))))
    members <- xy[hg$hyperedges[[bridge]], , drop = FALSE]
    max(stats::dist(members)) / frame
  }, numeric(1L))
  # uncut, the bridging pebble averaged 0.43 of the frame and reached further
  expect_true(mean(spans) < 0.42)
  expect_lt(hypernets:::.thg_repulsion_cutoff, Inf)
})

test_that("more than nine components keep their members together", {
  # the packing maps rows back through split(), whose grouping vector becomes a
  # factor: if those levels ordered as character, component 10 would collect
  # component 2's coordinates and hyperedges would be torn apart
  hg <- group_hypergraph(
    do.call(rbind, lapply(1:12, function(k) data.frame(
      member = paste0("c", k, "_", 1:3), event = paste0("e", k),
      stringsAsFactors = FALSE))),
    "member", "event"
  )
  pos <- hypernets:::.thg_positions(hg, "bipartite", seed = 1L)
  xy <- cbind(pos$nodes$x, pos$nodes$y)
  spread <- vapply(seq_len(hg$n_hyperedges), function(k)
    max(stats::dist(xy[hg$hyperedges[[k]], , drop = FALSE])), numeric(1L))
  # each hyperedge is a component of its own, so it must stay a tight clump
  expect_lt(max(spread), 0.25)
})

test_that("packing leaves a connected hypergraph exactly as it was", {
  connected <- matrix(0, 5L, 5L)
  connected[cbind(c(1, 2, 3, 4), c(2, 3, 4, 5))] <- 1
  connected <- connected + t(connected)
  expect_identical(
    hypernets:::.thg_layout_packed(connected, seed = 3L),
    hypernets:::.thg_layout_fr(connected, seed = 3L)
  )
})

test_that("the packed layout is seeded and leaves the caller's stream alone", {
  hg <- .island_fixture()
  set.seed(99)
  before <- stats::runif(1L)
  set.seed(99)
  first <- hypernets:::.thg_positions(hg, "bipartite", seed = 7L)$nodes
  expect_identical(stats::runif(1L), before)
  expect_identical(hypernets:::.thg_positions(hg, "bipartite", seed = 7L)$nodes,
                   first)
  other <- hypernets:::.thg_positions(hg, "bipartite", seed = 8L)$nodes
  expect_false(isTRUE(all.equal(other$x, first$x)))
})

test_that("a wider padding pushes the components further apart", {
  hg <- .island_fixture()
  gap_at <- function(padding) {
    pos <- hypernets:::.thg_positions(hg, "bipartite", seed = 1L,
                                      padding = padding)
    xy <- cbind(pos$nodes$x, pos$nodes$y)
    islands <- list(match(c("i1", "i2", "i3"), hg$nodes),
                    match(c("j1", "j2", "j3"), hg$nodes),
                    match(c("k1", "k2"), hg$nodes),
                    match(c("a", "b", "c", "d", "e"), hg$nodes))
    centres <- t(vapply(islands, function(rows)
      colMeans(xy[rows, , drop = FALSE]), numeric(2L)))
    min(stats::dist(centres))
  }
  # the pebbles grow with `padding`, so the seam between components has to
  # grow with it too or the blobs would touch and imply a shared member
  expect_gt(gap_at(0.12), gap_at(0.045))
})

test_that("a disconnected hypergraph plots, dismantles and keeps its labels", {
  hg <- .island_fixture()
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  p <- plot(hg, color_by = "size", edge_labels = TRUE)
  expect_s3_class(p, "ggplot")
  expect_s3_class(ggplot2::ggplot_build(p), "ggplot_built")
  panels <- plot(hg, dismantled = TRUE, ncol = 2L)
  expect_s3_class(ggplot2::ggplot_build(panels), "ggplot_built")
  # a panel title is an identifier, so it must never be cut to the panel
  expect_identical(panels$theme$strip.clip, "off")
})

test_that("every hyperedge label still lands inside its own pebble", {
  hg <- .island_fixture()
  pos <- hypernets:::.thg_positions(hg, "bipartite", seed = 1L)
  drawable <- which(lengths(hg$hyperedges) >= 2L)
  inside <- vapply(drawable, function(k) {
    members <- hg$hyperedges[[k]]
    ring <- hypernets:::.thg_pebble(pos$nodes$x[members], pos$nodes$y[members],
                                   0.045, detail = 5)
    own <- c(pos$edges$x[k], pos$edges$y[k])
    centroid <- c(mean(pos$nodes$x[members]), mean(pos$nodes$y[members]))
    hypernets:::.thg_in_convex(ring$x, ring$y, own[1L], own[2L]) ||
      hypernets:::.thg_in_convex(ring$x, ring$y, centroid[1L], centroid[2L])
  }, logical(1L))
  expect_true(all(inside))
})
