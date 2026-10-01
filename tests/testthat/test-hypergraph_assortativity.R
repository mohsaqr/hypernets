# hg_assortativity() and hg_degree_correlation(): an enumeration reference
# written from Chodrow's (2020) definition (a different code path from the
# closed-form moments in the package), the graph reduction to Newman's
# coefficient, invariants and error classes. External oracles (XGI,
# hypergraphx) live in local_testing_and_equivalence/.

.as_hg <- function(edges, sparse = FALSE) {
  long <- data.frame(
    member = unlist(edges),
    group = rep(names(edges), lengths(edges)),
    stringsAsFactors = FALSE
  )
  group_hypergraph(long, actor = "member", group = "group", sparse = sparse)
}

# Enumerate the chosen pairs with their probability weights, then take the
# weighted Pearson correlation of the two coordinates.
.as_reference <- function(edges, type, scale) {
  nodes <- sort(unique(unlist(edges)))
  degree <- vapply(nodes, \(v) sum(vapply(edges, \(e) v %in% e, logical(1))),
                   numeric(1))
  score <- if (scale == "rank") rank(degree) else degree
  rows <- lapply(Filter(\(e) length(e) >= 2L, edges), \(e) {
    s <- score[e]
    m <- length(s)
    if (type == "uniform") {
      grid <- expand.grid(i = seq_len(m), j = seq_len(m))
      grid <- subset(grid, i != j)
      data.frame(x = s[grid$i], y = s[grid$j], w = 1 / (m * (m - 1)))
    } else if (type == "top_2") {
      top <- sort(s, decreasing = TRUE)[1:2]
      data.frame(x = top, y = rev(top), w = 0.5)
    } else {
      data.frame(x = max(s), y = min(s), w = 1)
    }
  })
  pairs <- do.call(rbind, rows)
  stats::cov.wt(pairs[c("x", "y")], wt = pairs$w, cor = TRUE)$cor[1, 2]
}

set.seed(19)
as_edges <- lapply(seq_len(40), \(i) paste0("n", sample(20, sample(2:5, 1))))
names(as_edges) <- sprintf("e%02d", seq_along(as_edges))

test_that("every type and scale matches the enumerated definition", {
  hg <- .as_hg(as_edges)
  types <- c("uniform", "top_2", "top_bottom")
  vapply(c("degree", "rank"), \(scale) {
    out <- hg_assortativity(hg, type = types, scale = scale)
    expect_identical(out$type, types)
    expect_identical(out$n_edges, rep(40L, 3L))
    reference <- vapply(types, \(t) .as_reference(as_edges, t, scale),
                        numeric(1))
    expect_equal(out$assortativity, unname(reference), tolerance = 1e-12)
    TRUE
  }, logical(1))
})

test_that("a star graph is perfectly disassortative", {
  edges <- lapply(paste0("leaf", 1:6), \(l) c("hub", l))
  names(edges) <- paste0("e", 1:6)
  out <- hg_assortativity(.as_hg(edges), type = c("uniform", "top_bottom"),
                          scale = "degree")
  expect_equal(out$assortativity[1], -1, tolerance = 1e-12)
  # top-bottom always picks (hub, leaf): both coordinates are constant
  expect_true(is.na(out$assortativity[2]))
})

test_that("on a graph, uniform equals Newman's degree assortativity", {
  skip_if_not_installed("igraph")
  set.seed(5)
  g <- igraph::sample_gnm(30, 60)
  ends <- igraph::as_edgelist(g)
  edges <- lapply(seq_len(nrow(ends)), \(i) paste0("v", ends[i, ]))
  names(edges) <- sprintf("e%03d", seq_along(edges))
  out <- hg_assortativity(.as_hg(edges), type = c("uniform", "top_2"),
                          scale = "degree")
  newman <- igraph::assortativity_degree(g, directed = FALSE)
  expect_equal(out$assortativity, c(newman, newman), tolerance = 1e-12)
})

test_that("undefined coefficients are NA (regular hypergraph)", {
  # every node has hyperdegree 2
  edges <- list(e1 = c("a", "b", "c"), e2 = c("d", "e", "f"),
                e3 = c("a", "d"), e4 = c("b", "e"), e5 = c("c", "f"))
  out <- hg_assortativity(.as_hg(edges), type = c("uniform", "top_bottom"))
  expect_true(all(is.na(out$assortativity)))
})

test_that("dense and sparse agree; relabelling does not matter", {
  types <- c("uniform", "top_2", "top_bottom")
  dense <- hg_assortativity(.as_hg(as_edges), type = types)
  sparse <- hg_assortativity(.as_hg(as_edges, sparse = TRUE), type = types)
  expect_equal(sparse, dense, tolerance = 1e-14)

  relabel <- setNames(paste0("z", sample(20)), paste0("n", seq_len(20)))
  moved <- lapply(rev(as_edges), \(e) unname(relabel[e]))
  names(moved) <- paste0("f", seq_along(moved))
  expect_equal(hg_assortativity(.as_hg(moved), type = types), dense,
               tolerance = 1e-12)
  values <- dense$assortativity
  expect_true(all(values >= -1 & values <= 1))
})

test_that("singleton hyperedges are ignored", {
  with_singletons <- c(as_edges, list(s1 = "n1", s2 = "n2"))
  a <- hg_assortativity(.as_hg(as_edges))
  b <- hg_assortativity(.as_hg(with_singletons), scale = "rank")
  expect_identical(b$n_edges, 40L)
  # degrees change (n1, n2 gain a hyperedge), so only the edge count is fixed
  expect_true(is.finite(b$assortativity))
  expect_true(is.finite(a$assortativity))
})

test_that("hg_degree_correlation: hand-computed Pearson across orders", {
  edges <- list(p1 = c("a", "b"), p2 = c("a", "c"), p3 = c("b", "c"),
                p4 = c("a", "d"), t1 = c("a", "b", "c"),
                t2 = c("b", "c", "d"), q1 = c("a", "b", "c", "d"))
  out <- hg_degree_correlation(.as_hg(edges))
  expect_identical(out$size_1, c(2L, 2L, 3L))
  expect_identical(out$size_2, c(3L, 4L, 4L))
  pair_degree <- c(a = 3, b = 2, c = 2, d = 1)
  triple_degree <- c(a = 1, b = 2, c = 2, d = 1)
  expect_equal(out$correlation[1], stats::cor(pair_degree, triple_degree),
               tolerance = 1e-14)
  # size 4: every node has degree 1 -> constant -> NA
  expect_true(all(is.na(out$correlation[2:3])))
  expect_identical(out$n_nodes, rep(4L, 3L))
})

test_that("broken contracts raise classed errors", {
  expect_error(hg_assortativity(list()), class = "hypernets_bad_input")
  expect_error(hg_assortativity(.as_hg(list(s1 = "a", s2 = "b"))),
               class = "hypernets_bad_input")
  expect_error(hg_degree_correlation(.as_hg(list(e1 = c("a", "b"),
                                                 e2 = c("b", "c")))),
               class = "hypernets_bad_input")
  expect_error(hg_degree_correlation(list()), class = "hypernets_bad_input")
})
