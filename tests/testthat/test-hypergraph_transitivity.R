# hg_transitivity(): hand-computed values from the published definitions,
# the graph reduction (Klimm et al. 2021: extra overlap = ordinary
# clustering on a graph), dense/sparse parity, invariants and error classes.
# The XGI oracle lives in local_testing_and_equivalence/.

.tr_hg <- function(edges, sparse = FALSE) {
  long <- data.frame(
    member = unlist(edges),
    group = rep(names(edges), lengths(edges)),
    stringsAsFactors = FALSE
  )
  group_hypergraph(long, actor = "member", group = "group", sparse = sparse)
}

all_types <- c("projection", "extra_overlap", "two_node_union",
               "two_node_min", "two_node_max")

test_that("hand-computed values on a three-edge hypergraph", {
  # e1 = {a, b, c}, e2 = {a, d}, e3 = {b, d}
  hg <- .tr_hg(list(e1 = c("a", "b", "c"), e2 = c("a", "d"),
                    e3 = c("b", "d")))
  out <- hg_transitivity(hg, type = all_types)
  a <- subset(out, node == "a")
  # projection: neighbours b, c, d; linked pairs bc, bd of 3
  expect_equal(a$projection, 2 / 3, tolerance = 1e-14)
  # extra overlap of (e1, e2): D1 = {b, c}, D2 = {d};
  # N(D1) meets D2 in d, N(D2) = {a, b} meets D1 in b -> (1 + 1) / 3
  expect_equal(a$extra_overlap, 2 / 3, tolerance = 1e-14)
  # two-node: M(a) = {e1, e2}; b {e1, e3}: 1/3 | 1/2 | 1/2;
  # c {e1}: 1/2 | 1 | 1/2; d {e2, e3}: 1/3 | 1/2 | 1/2
  expect_equal(a$two_node_union, (1 / 3 + 1 / 2 + 1 / 3) / 3,
               tolerance = 1e-14)
  expect_equal(a$two_node_min, (1 / 2 + 1 + 1 / 2) / 3, tolerance = 1e-14)
  expect_equal(a$two_node_max, (1 / 2 + 1 / 2 + 1 / 2) / 3,
               tolerance = 1e-14)
  # c: hyperdegree 1 -> extra overlap 0 (Klimm et al.), neighbours a, b linked
  c_row <- subset(out, node == "c")
  expect_identical(c_row$extra_overlap, 0)
  expect_equal(c_row$projection, 1)
})

test_that("a single hyperedge closes every triangle but has no extra overlap", {
  hg <- .tr_hg(list(e1 = c("a", "b", "c", "d")))
  out <- hg_transitivity(hg, type = all_types)
  expect_true(all(out$projection == 1))
  expect_true(all(out$extra_overlap == 0))
  expect_true(all(out$two_node_union == 1))
})

test_that("on a graph, extra overlap and projection equal igraph local transitivity", {
  skip_if_not_installed("igraph")
  set.seed(7)
  g <- igraph::sample_gnm(25, 70)
  ends <- igraph::as_edgelist(g)
  edges <- lapply(seq_len(nrow(ends)), \(i) paste0("v", ends[i, ]))
  names(edges) <- sprintf("e%03d", seq_along(edges))
  hg <- .tr_hg(edges)
  out <- hg_transitivity(hg, type = c("projection", "extra_overlap"))
  local <- igraph::transitivity(g, type = "local", isolates = "NaN")
  names(local) <- paste0("v", seq_len(igraph::vcount(g)))
  reference <- unname(local[out$node])
  defined <- !is.na(reference)
  expect_equal(out$projection[defined], reference[defined], tolerance = 1e-12)
  expect_equal(out$extra_overlap[defined], reference[defined],
               tolerance = 1e-12)
  expect_identical(is.na(out$projection), !defined)
})

test_that("undefined values are NA, not 0", {
  # f has one neighbour (projection undefined); g sits only in a singleton
  hg <- .tr_hg(list(e1 = c("a", "b", "c"), e2 = c("c", "f"), e3 = "g"))
  out <- hg_transitivity(hg, type = all_types)
  expect_true(is.na(subset(out, node == "f")$projection))
  g_row <- subset(out, node == "g")
  expect_true(is.na(g_row$projection))
  expect_true(is.na(g_row$two_node_union))
  expect_identical(g_row$extra_overlap, 0)
  summary <- hg_transitivity(hg, type = all_types, what = "summary")
  expect_identical(summary$type, all_types)
  expect_identical(subset(summary, type == "projection")$n_undefined, 2L)
  expect_equal(subset(summary, type == "projection")$mean_local,
               mean(out$projection, na.rm = TRUE))
})

test_that("dense and sparse incidences give identical coefficients", {
  set.seed(3)
  edges <- lapply(seq_len(30), \(i) paste0("n", sample(15, sample(2:5, 1))))
  names(edges) <- sprintf("e%02d", seq_along(edges))
  dense <- hg_transitivity(.tr_hg(edges), type = all_types)
  sparse <- hg_transitivity(.tr_hg(edges, sparse = TRUE), type = all_types)
  expect_equal(sparse, dense, tolerance = 1e-14)
})

test_that("invariants: range, relabelling, duplicate hyperedges", {
  set.seed(11)
  edges <- lapply(seq_len(25), \(i) paste0("n", sample(12, sample(2:4, 1))))
  names(edges) <- sprintf("e%02d", seq_along(edges))
  base <- hg_transitivity(.tr_hg(edges), type = all_types)
  values <- unlist(base[all_types])
  expect_true(all(values[!is.na(values)] >= 0 & values[!is.na(values)] <= 1))

  # relabelling nodes and reordering hyperedges permutes the rows only
  perm <- sample(12)
  relabel <- setNames(paste0("m", perm), paste0("n", seq_len(12)))
  moved <- lapply(rev(edges), \(e) unname(relabel[e]))
  names(moved) <- rev(names(edges))
  shuffled <- hg_transitivity(.tr_hg(moved), type = all_types)
  back <- shuffled[match(relabel[base$node], shuffled$node), all_types]
  rownames(back) <- NULL
  expect_equal(back, base[all_types], tolerance = 1e-14)

  # the projection only sees who is linked, so a duplicate edge is invisible
  twice <- c(edges, list(dup = edges[[1]]))
  expect_equal(hg_transitivity(.tr_hg(twice))$projection, base$projection)
})

test_that("sort_by and n return the top rows", {
  hg <- .tr_hg(list(e1 = c("a", "b", "c"), e2 = c("a", "d"),
                    e3 = c("b", "d"), e4 = c("d", "e")))
  top <- hg_transitivity(hg, type = "two_node_union",
                         sort_by = "two_node_union", n = 2)
  full <- hg_transitivity(hg, type = "two_node_union")
  expect_identical(nrow(top), 2L)
  expect_identical(top$two_node_union,
                   sort(full$two_node_union, decreasing = TRUE)[1:2])
})

test_that("broken contracts raise classed errors", {
  expect_error(hg_transitivity(list(incidence = diag(2))),
               class = "hypernets_bad_input")
  hg <- .tr_hg(list(e1 = c("a", "b", "c")))
  expect_error(hg_transitivity(hg, n = 0), class = "hypernets_bad_input")
  expect_error(hg_transitivity(hg, type = "global"))
})
