testthat::skip_on_cran()

# The community figure fits the device: nothing is drawn within two pixels
# of any edge (a clipped legend, caption, title box or label runs into the
# edge).

.edge_ink <- function(p, width) {
  file <- tempfile(fileext = ".png")
  grDevices::png(file, width = width, height = 0.9 * width, units = "in",
                 res = 72)
  print(p)
  grDevices::dev.off()
  img <- png::readPNG(file)
  ink <- 1 - apply(img[, , 1:3, drop = FALSE], c(1L, 2L), min)
  rim <- c(1L, 2L)
  c(left = max(ink[, rim]), right = max(ink[, ncol(ink) + 1L - rim]),
    top = max(ink[rim, ]), bottom = max(ink[nrow(ink) + 1L - rim, ]))
}

test_that("the ring community figure fits at 6, 8 and 10 inches", {
  skip_if_not_installed("png")
  comm <- hg_communities(hon(ring_sequences, max_order = 2L),
                          trials = 3L, seed = 1L)
  p <- plot(comm)
  ink <- do.call(rbind, lapply(c(6, 8, 10), function(width) {
    data.frame(width = width, t(.edge_ink(p, width)))
  }))
  expect_true(all(ink[c("left", "right", "top", "bottom")] == 0),
              label = paste(capture.output(print(ink)), collapse = "\n"))
})

test_that("the community boxes clear the labels of the nodes at the rim", {
  comm <- hg_communities(hon(ring_sequences, max_order = 2L),
                          trials = 3L, seed = 1L)
  p <- plot(comm)
  boxes <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)[[1L]]$data
  default_gap <- plot(comm, title_gap = 0.06)
  near <- Filter(function(l) inherits(l$geom, "GeomLabel"),
                 default_gap$layers)[[1L]]$data
  # the same headings, further out
  expect_identical(boxes$hjust, near$hjust)
  expect_true(all(abs(boxes$y) + abs(boxes$x) != abs(near$y) + abs(near$x)))
})
