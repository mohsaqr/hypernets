# ---- build_hypergraph() tests --------------------------------------------

# Helpers ------------------------------------------------------------------

.hg_triangle_adj <- function() {
  # 3-node triangle (one 3-clique, three 2-edges)
  adj <- matrix(c(0, 1, 1,
                  1, 0, 1,
                  1, 1, 0), 3L, 3L)
  rownames(adj) <- colnames(adj) <- c("A", "B", "C")
  adj
}

.hg_two_triangles_adj <- function() {
  # Two disjoint triangles (A-B-C and D-E-F): 2 cliques, 6 pairwise edges
  adj <- matrix(0L, 6L, 6L)
  rownames(adj) <- colnames(adj) <- c("A", "B", "C", "D", "E", "F")
  for (e in list(c(1, 2), c(1, 3), c(2, 3),
                 c(4, 5), c(4, 6), c(5, 6))) {
    adj[e[1], e[2]] <- adj[e[2], e[1]] <- 1L
  }
  adj
}

.hg_chain_adj <- function() {
  # A-B-C (chain, no triangle)
  adj <- matrix(0L, 3L, 3L)
  rownames(adj) <- colnames(adj) <- c("A", "B", "C")
  adj[1, 2] <- adj[2, 1] <- 1L
  adj[2, 3] <- adj[3, 2] <- 1L
  adj
}

# Structure ----------------------------------------------------------------

test_that("build_hypergraph returns a net_hypergraph with required fields", {
  hg <- build_hypergraph(.hg_triangle_adj(), p = 1, max_size = 3L)
  expect_s3_class(hg, "net_hypergraph")
  expect_named(hg, c("hyperedges", "incidence", "nodes", "n_nodes",
                     "n_hyperedges", "size_distribution", "params"))
  expect_type(hg$hyperedges, "list")
  expect_true(is.matrix(hg$incidence))
  expect_equal(dim(hg$incidence), c(3L, hg$n_hyperedges))
  expect_identical(hg$nodes, c("A", "B", "C"))
})

# Triangle, p = 1, include_pairwise = TRUE ---------------------------------

test_that("single triangle p=1 yields 1 three-edge + 3 two-edges", {
  hg <- build_hypergraph(.hg_triangle_adj(), p = 1, include_pairwise = TRUE)
  expect_equal(hg$n_hyperedges, 4L)
  sizes <- vapply(hg$hyperedges, length, integer(1))
  expect_equal(sort(sizes), c(2L, 2L, 2L, 3L))
  # incidence column sums = hyperedge sizes
  expect_equal(unname(colSums(hg$incidence)), unname(sizes))
})

# Triangle, p = 0 ---------------------------------------------------------

test_that("p=0 drops all higher-order, keeps pairwise", {
  hg <- build_hypergraph(.hg_triangle_adj(), p = 0, include_pairwise = TRUE)
  expect_equal(hg$n_hyperedges, 3L)
  expect_true(all(vapply(hg$hyperedges, length, integer(1)) == 2L))
})

# Triangle, p = 1, include_pairwise = FALSE -------------------------------

test_that("include_pairwise=FALSE keeps only k>=3 hyperedges", {
  hg <- build_hypergraph(.hg_triangle_adj(), p = 1, include_pairwise = FALSE)
  expect_equal(hg$n_hyperedges, 1L)
  expect_equal(length(hg$hyperedges[[1]]), 3L)
  expect_equal(sort(hg$hyperedges[[1]]), 1:3)
})

# No triangles ------------------------------------------------------------

test_that("triangle-free network yields only 2-hyperedges", {
  hg <- build_hypergraph(.hg_chain_adj(), p = 1, max_size = 3L)
  expect_true(all(vapply(hg$hyperedges, length, integer(1)) == 2L))
  expect_equal(hg$n_hyperedges, 2L)
})

test_that("triangle-free + include_pairwise=FALSE yields empty hypergraph", {
  hg <- build_hypergraph(.hg_chain_adj(), p = 1,
                         include_pairwise = FALSE, max_size = 3L)
  expect_equal(hg$n_hyperedges, 0L)
  expect_equal(ncol(hg$incidence), 0L)
  expect_equal(nrow(hg$incidence), 3L)
  expect_length(hg$hyperedges, 0)
})

# Reproducibility ---------------------------------------------------------

test_that("seed gives reproducible sampling at intermediate p", {
  adj <- .hg_two_triangles_adj()
  hg1 <- build_hypergraph(adj, p = 0.5, seed = 42L)
  hg2 <- build_hypergraph(adj, p = 0.5, seed = 42L)
  expect_identical(hg1$hyperedges, hg2$hyperedges)
})

test_that("different seeds give (probably) different samples", {
  adj <- .hg_two_triangles_adj()
  # 2 triangles independently sampled at p=0.5: 4 outcomes,
  # only a 1/4 chance of matching by accident.
  set.seed(NULL)
  outcomes <- replicate(20, {
    hg1 <- build_hypergraph(adj, p = 0.5, seed = sample.int(1e6, 1))
    hg2 <- build_hypergraph(adj, p = 0.5, seed = sample.int(1e6, 1))
    identical(hg1$hyperedges, hg2$hyperedges)
  })
  # Out of 20 random pairs, at least one should differ.
  expect_false(all(outcomes))
})

# Incidence matrix structure ----------------------------------------------

test_that("incidence rows index nodes, columns index hyperedges", {
  hg <- build_hypergraph(.hg_two_triangles_adj(), p = 1)
  expect_equal(rownames(hg$incidence), hg$nodes)
  expect_equal(ncol(hg$incidence), hg$n_hyperedges)
  # node A (idx 1) belongs only to A-B, A-C edges and the A-B-C triangle
  a_membership <- which(hg$incidence["A", ] == 1L)
  for (j in a_membership) {
    expect_true(1L %in% hg$hyperedges[[j]])
  }
})

# Size distribution -------------------------------------------------------

test_that("size distribution matches hyperedge sizes", {
  hg <- build_hypergraph(.hg_two_triangles_adj(), p = 1, max_size = 3L)
  expect_equal(hg$size_distribution[["size_2"]], 6L)
  expect_equal(hg$size_distribution[["size_3"]], 2L)
})

# Higher max_size includes 4-cliques --------------------------------------

test_that("max_size = 4 includes 4-cliques", {
  # K4: complete graph on 4 nodes => one 4-clique, four 3-cliques, six 2-edges
  adj <- matrix(1L, 4L, 4L); diag(adj) <- 0L
  rownames(adj) <- colnames(adj) <- c("A", "B", "C", "D")
  hg <- build_hypergraph(adj, p = 1, max_size = 4L, include_pairwise = TRUE)
  sizes <- vapply(hg$hyperedges, length, integer(1))
  expect_equal(sum(sizes == 2L), 6L)
  expect_equal(sum(sizes == 3L), 4L)
  expect_equal(sum(sizes == 4L), 1L)
})

# Input flexibility -------------------------------------------------------

test_that("accepts a netobject", {
  adj <- .hg_two_triangles_adj()
  net <- .wrap_netobject(adj * 1.0, method = "manual", directed = FALSE)
  hg <- build_hypergraph(net, p = 1)
  expect_s3_class(hg, "net_hypergraph")
  expect_equal(hg$n_nodes, 6L)
})

test_that("accepts a simplicial_complex", {
  adj <- .hg_two_triangles_adj()
  sc <- build_simplicial(adj, type = "clique", max_dim = 2L)
  hg <- build_hypergraph(sc, p = 1)
  expect_s3_class(hg, "net_hypergraph")
  expect_equal(hg$n_nodes, 6L)
  expect_true(any(vapply(hg$hyperedges, length, integer(1)) == 3L))
})

# Input validation --------------------------------------------------------

test_that("rejects out-of-range p", {
  expect_error(build_hypergraph(.hg_triangle_adj(), p = -0.1))
  expect_error(build_hypergraph(.hg_triangle_adj(), p =  1.1))
})

test_that("rejects max_size < 2", {
  expect_error(build_hypergraph(.hg_triangle_adj(), max_size = 1L))
})

test_that("rejects bad input type", {
  expect_error(build_hypergraph("not a network"),
               "netobject, cograph_network")
})

# Print/summary methods don't error --------------------------------------

test_that("print and summary methods run without error", {
  hg <- build_hypergraph(.hg_two_triangles_adj(), p = 1)
  expect_invisible(print(hg))
  s <- summary(hg)
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s), c("node", "degree"))
})

test_that("print/summary handle empty hypergraph", {
  hg <- build_hypergraph(.hg_chain_adj(), p = 1,
                         include_pairwise = FALSE)
  expect_invisible(print(hg))
  s <- summary(hg)
  expect_s3_class(s, "data.frame")
})

test_that("as.data.frame(what = 'nodes') reports one row per node", {
  df <- data.frame(
    member = c("a", "b", "c", "b", "c", "d"),
    session = c("s1", "s1", "s1", "s2", "s2", "s2")
  )
  hg <- group_hypergraph(df, actor = "member", group = "session")
  nodes <- as.data.frame(hg, what = "nodes")
  expect_identical(names(nodes), c("node", "degree"))
  expect_identical(nodes$node, c("a", "b", "c", "d"))
  expect_identical(nodes$degree, c(1L, 2L, 2L, 1L))
  by_degree <- as.data.frame(hg, what = "nodes", sort_by = "degree", top = 2)
  expect_identical(by_degree$node, c("b", "c"))
  sparse <- group_hypergraph(df, actor = "member", group = "session",
                             sparse = TRUE)
  sparse_nodes <- as.data.frame(sparse, what = "nodes")
  expect_identical(sparse_nodes, nodes)
  # planted blocks ride along
  P <- matrix(c(.5, .05, .05, .5), 2, 2)
  sbm <- hg_sample_sbm(P = P, block_sizes = c(3, 3), d = 2, seed = 1)
  sbm_nodes <- as.data.frame(sbm, what = "nodes")
  expect_identical(names(sbm_nodes), c("node", "degree", "block"))
  expect_identical(sbm_nodes$block, rep(1:2, each = 3))
})

test_that("as.data.frame(what = 'memberships') lists every incidence cell", {
  dfw <- data.frame(person = c("A", "B", "A", "A", "B"),
                    grp = c("g1", "g1", "g1", "g2", "g2"),
                    n = c(2, 5, 3, 1, 4), stringsAsFactors = FALSE)
  hw <- group_hypergraph(dfw, actor = "person", group = "grp", weight = "n")
  cells <- as.data.frame(hw, what = "memberships")
  expect_identical(names(cells), c("node", "hyperedge", "weight"))
  expect_identical(cells$node, c("A", "B", "A", "B"))
  expect_identical(cells$hyperedge, c("g1", "g1", "g2", "g2"))
  expect_equal(cells$weight, c(5, 5, 1, 4))
  heaviest <- as.data.frame(hw, what = "memberships", sort_by = "weight", top = 1)
  expect_identical(heaviest$hyperedge, "g1")
  # INVARIANT: the table is the incidence matrix in long form, dense or sparse
  sparse <- group_hypergraph(dfw, actor = "person", group = "grp",
                             weight = "n", sparse = TRUE)
  sparse_cells <- as.data.frame(sparse, what = "memberships")
  expect_equal(sparse_cells, cells)
  expect_equal(sum(cells$weight), sum(hw$incidence))
  binary <- group_hypergraph(dfw, actor = "person", group = "grp")
  binary_cells <- as.data.frame(binary, what = "memberships")
  expect_equal(binary_cells$weight, rep(1, 4))
})
