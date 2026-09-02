test_that("association projection can count or collapse repeated hyperedges", {
  h <- group_hypergraph(
    data.frame(member = rep(c("a", "b", "c"), 2),
               event = rep(c("e1", "e2"), each = 3)),
    "member", "event"
  )
  counted <- hg_project(h, method = "association", what = "matrix")
  collapsed <- hg_project(h, method = "association", what = "matrix",
                          duplicate_edges = "collapse")
  expect_equal(counted["a", "b"], 1)
  expect_equal(collapsed["a", "b"], 0.5)
})

test_that("self-association contributes one unit per source", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c"),
               event = rep(c("e1", "e2"), each = 2)),
    "member", "event"
  )
  w <- hg_project(h, method = "association", what = "matrix",
                  self_association = TRUE,
                  edge_source = c(e1 = "s", e2 = "s"))
  expect_equal(w["a", "b"], 1)
  expect_equal(w["b", "c"], 1)
  expect_equal(unname(w["s", c("a", "b", "c")]), c(0.25, 0.5, 0.25))
  expect_equal(sum(w["s", ]), 1)
})

test_that("self-association supports sparse incidence", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "b", "c"),
               event = rep(c("e1", "e2"), each = 2)),
    "member", "event", sparse = TRUE
  )
  w <- hg_project(h, method = "association", what = "matrix",
                  self_association = TRUE,
                  edge_source = c(e1 = "s", e2 = "s"))
  expect_s4_class(w, "sparseMatrix")
  expect_equal(unname(as.matrix(w)["s", c("a", "b", "c")]),
               c(0.25, 0.5, 0.25))
})

test_that("hyperedge incidence threshold is configurable", {
  h <- group_hypergraph(
    data.frame(member = c("a", "b", "c", "b", "c", "d", "c", "e"),
               event = rep(c("e1", "e2", "e3"), c(3, 3, 2))),
    "member", "event"
  )
  expect_equal(hg_edges(h, s = 1)$n_incident_edges, c(2L, 2L, 2L))
  expect_equal(hg_edges(h, s = 2)$n_incident_edges, c(1L, 1L, 0L))
})
