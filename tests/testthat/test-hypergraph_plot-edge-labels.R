skip_on_cran()
skip_if_not_installed("cograph")
skip_if_not_installed("ggplot2")

# A layout with one clearly outlying blob: that blob's label is the one the
# anchor rule pushes furthest out, so it is the one that used to be clipped.
petal_hg <- function() {
  memberships <- data.frame(
    actor = c("a1","a2","a3", "b1","b2","b3", "c1","c2","c3",
              "d1","d2","d3", "core","core","core","core"),
    group = c("alpha","alpha","alpha", "beta","beta","beta",
              "gamma","gamma","gamma", "delta","delta","delta",
              "alpha","beta","gamma","delta")
  )
  group_hypergraph(memberships, actor = "actor", group = "group")
}

edge_label_data <- function(p) {
  hits <- Filter(function(l) is.data.frame(l$data) &&
                   identical(names(l$data), c("x", "y", "label")), p$layers)
  if (length(hits) == 0L) NULL else hits[[length(hits)]]$data
}

test_that("edge_labels = FALSE is the default and writes nothing", {
  expect_null(edge_label_data(plot(petal_hg())))
  expect_null(edge_label_data(plot(petal_hg(), edge_labels = FALSE)))
})

test_that("edge_labels = TRUE writes one label per drawn hyperedge", {
  d <- edge_label_data(plot(petal_hg(), edge_labels = TRUE))
  expect_identical(nrow(d), 4L)
  expect_setequal(d$label, c("alpha", "beta", "gamma", "delta"))
})

test_that("every hyperedge label falls inside the panel", {
  # regression: the limits were fixed before the labels were placed, so the
  # outermost blob's label was pushed past them and silently dropped
  for (size in c(2.5, 4, 6)) {
    p <- plot(petal_hg(), edge_labels = TRUE, edge_label_size = size)
    d <- edge_label_data(p)
    built <- ggplot2::ggplot_build(p)
    xr <- built$layout$panel_params[[1L]]$x.range
    yr <- built$layout$panel_params[[1L]]$y.range
    expect_true(all(d$x >= xr[1L] & d$x <= xr[2L]),
                info = paste("x out of range at size", size))
    expect_true(all(d$y >= yr[1L] & d$y <= yr[2L]),
                info = paste("y out of range at size", size))
  }
})

test_that("a named vector renames the hyperedges", {
  d <- edge_label_data(plot(petal_hg(), edge_labels = c(alpha = "First",
                                                        beta = "Second")))
  expect_true(all(c("First", "Second") %in% d$label))
  # a hyperedge the vector does not name keeps its own name
  expect_true(all(c("gamma", "delta") %in% d$label))
})

test_that("edge_labels rejects an unnamed vector by class", {
  expect_error(plot(petal_hg(), edge_labels = c("one", "two")),
               class = "hypernets_bad_input")
})

test_that("node labels and hyperedge labels are independent", {
  p <- plot(petal_hg(), labels = FALSE, edge_labels = TRUE)
  expect_identical(nrow(edge_label_data(p)), 4L)
})
