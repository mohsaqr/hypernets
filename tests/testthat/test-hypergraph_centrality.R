# ---- hypergraph_centrality() tests ---------------------------------------

# Helpers ------------------------------------------------------------------

.hc_two_overlapping <- function() {
  d <- data.frame(
    member  = c("A", "B", "C",  "A", "B", "D"),
    session = c("S1", "S1", "S1", "S2", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  group_hypergraph(d, "member", "session")
}

# hypergraph_centrality() returns a tidy table (one row per node, one column
# per requested type). Several checks below are about a single centrality as
# a named-by-node vector; this recovers that view without changing what any
# test asserts mathematically.
.hc_vec <- function(cent, type) stats::setNames(cent[[type]], cent$node)

# the centrality columns, i.e. every column except the node id
.hc_types <- function(cent) setdiff(names(cent), "node")

# Structure ----------------------------------------------------------------

test_that("default returns a tidy table with all three types", {
  cent <- hypergraph_centrality(.hc_two_overlapping())
  expect_s3_class(cent, "data.frame")
  expect_named(cent, c("node", "clique", "Z", "H"))
  expect_identical(cent$node, c("A", "B", "C", "D"))
  for (nm in .hc_types(cent)) {
    expect_length(cent[[nm]], 4L)
    expect_true(all(is.finite(cent[[nm]])))
  }
})

test_that("single type returns a table with just that column", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), type = "clique")
  expect_named(cent, c("node", "clique"))
  expect_length(cent$clique, 4L)
})

test_that("sort_by ranks the table without changing the values", {
  cent <- hypergraph_centrality(.hc_two_overlapping())
  srt  <- hypergraph_centrality(.hc_two_overlapping(), sort_by = "clique")
  expect_identical(srt$clique, sort(cent$clique, decreasing = TRUE))
  expect_identical(.hc_vec(srt, "clique")[cent$node], .hc_vec(cent, "clique"))
  expect_error(hypergraph_centrality(.hc_two_overlapping(),
                                     type = "clique", sort_by = "Z"),
               "should be")
})

# CEC vs igraph validation ------------------------------------------------

test_that("clique centrality matches igraph::eigen_centrality on expansion", {
  skip_if_not_installed("igraph")
  hg  <- .hc_two_overlapping()
  cent <- hypergraph_centrality(hg, type = "clique")
  ours <- .hc_vec(cent, "clique")
  # Recompute via clique expansion + igraph
  net <- clique_expansion(hg)
  W   <- net$weights
  g   <- igraph::graph_from_adjacency_matrix(W, mode = "undirected",
                                              weighted = TRUE, diag = FALSE)
  igr <- igraph::eigen_centrality(g, directed = FALSE)$vector
  # Match up to sign and scale
  ours_n <- ours / sqrt(sum(ours^2))
  igr_n  <- igr  / sqrt(sum(igr^2))
  if (sum(ours_n * igr_n) < 0) ours_n <- -ours_n
  expect_equal(ours_n, igr_n[names(ours_n)],
               tolerance = 1e-5, ignore_attr = TRUE)
})

# Symmetry / structural validation ----------------------------------------

test_that("symmetric nodes get equal centrality (two overlapping triangles)", {
  # A & B are structurally identical (both in S1 and S2).
  # C & D are structurally identical (each in exactly one session).
  cent <- hypergraph_centrality(.hc_two_overlapping())
  for (nm in .hc_types(cent)) {
    v <- .hc_vec(cent, nm)
    expect_equal(v[["A"]], v[["B"]], tolerance = 1e-6)
    expect_equal(v[["C"]], v[["D"]], tolerance = 1e-6)
    # A/B should have higher centrality than C/D (in 2 edges vs 1)
    expect_gt(v[["A"]], v[["C"]])
  }
})

# ZEC manual validation on small uniform hypergraph -----------------------

test_that("ZEC power-iteration update matches manual formula on one triangle", {
  # Single triangle (A,B,C), uniform k=3.
  # ZEC eigen-equation: λ x_i = Π_{j≠i} x_j  ⇒ by symmetry x_A=x_B=x_C.
  # Thus the fixed point is (a, a, a), any positive a, normalized.
  d <- data.frame(member  = c("A", "B", "C"),
                  session = c("S1", "S1", "S1"),
                  stringsAsFactors = FALSE)
  hg <- group_hypergraph(d, "member", "session")
  zec <- hypergraph_centrality(hg, type = "Z")
  cent <- .hc_vec(zec, "Z")
  expected <- rep(1 / sqrt(3), 3)
  names(expected) <- c("A", "B", "C")
  expect_equal(cent, expected, tolerance = 1e-6)
})

test_that("ZEC distinguishes hub-like vs peripheral in a 4-node, 2-edge case", {
  # Edges (A,B,C) and (A,B,D). A and B are in both edges.
  zec <- hypergraph_centrality(.hc_two_overlapping(), type = "Z")
  cent <- .hc_vec(zec, "Z")
  expect_gt(cent[["A"]], cent[["C"]])
  expect_gt(cent[["B"]], cent[["D"]])
})

# HEC on uniform hypergraph ------------------------------------------------

test_that("HEC reproduces ZEC-like ranking on uniform hypergraph", {
  # For a uniform hypergraph the eigenvector direction is the same;
  # only the scale differs (H takes the k-1 root of the ZEC update).
  cent <- hypergraph_centrality(.hc_two_overlapping())
  # Both should rank A=B above C=D
  expect_equal(rank(cent$Z), rank(cent$H), ignore_attr = TRUE)
})

# Non-negativity, normalization --------------------------------------------

test_that("centralities are non-negative for non-negative hypergraphs", {
  cent <- hypergraph_centrality(.hc_two_overlapping())
  for (nm in .hc_types(cent)) {
    expect_true(all(cent[[nm]] >= -1e-10))
  }
})

test_that("normalize=TRUE gives unit-L2-norm vectors (clique at least)", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), normalize = TRUE)
  expect_equal(sqrt(sum(cent$clique^2)), 1, tolerance = 1e-6)
})

test_that("normalize=FALSE gives max-abs = 1", {
  cent <- hypergraph_centrality(.hc_two_overlapping(), normalize = FALSE)
  expect_equal(max(abs(cent$clique)), 1, tolerance = 1e-6)
})

# Empty & degenerate -----------------------------------------------------

test_that("empty hypergraph returns all-zero centralities", {
  inc <- matrix(0L, 3L, 0L, dimnames = list(c("A", "B", "C"), NULL))
  hg <- structure(
    list(hyperedges = list(), incidence = inc, nodes = c("A", "B", "C"),
         n_nodes = 3L, n_hyperedges = 0L,
         size_distribution = integer(0), params = list()),
    class = "net_hypergraph")
  cent <- hypergraph_centrality(hg)
  expect_identical(cent$node, c("A", "B", "C"))
  for (nm in .hc_types(cent)) {
    expect_true(all(cent[[nm]] == 0))
  }
})

test_that("single size-1 hyperedge contributes nothing (k < 2)", {
  # Hyperedge of size 1 shouldn't participate (no "other member" to pair with)
  d <- data.frame(member  = c("A"),
                  session = c("S1"),
                  stringsAsFactors = FALSE)
  hg <- group_hypergraph(d, "member", "session")
  expect_no_error(hypergraph_centrality(hg, type = c("Z", "H")))
})

# Input validation --------------------------------------------------------

test_that("rejects non-net_hypergraph input", {
  expect_error(hypergraph_centrality(matrix(0, 3, 3)), "net_hypergraph")
})

test_that("unknown type is rejected by match.arg", {
  expect_error(hypergraph_centrality(.hc_two_overlapping(), type = "bogus"))
})

# Reproducibility ---------------------------------------------------------

test_that("power iteration is deterministic given same input", {
  hg <- .hc_two_overlapping()
  c1 <- hypergraph_centrality(hg)
  c2 <- hypergraph_centrality(hg)
  expect_identical(c1, c2)
})

# Integration with group_hypergraph + bundled data ------------------------

test_that("runs on bundled human_long dataset without error", {
  data("human_long", package = "hypernets")
  # Use a small subset to keep the test fast
  sub <- head(human_long, 500)
  hg  <- group_hypergraph(sub, "code", "session_id")
  cent <- hypergraph_centrality(hg, type = "clique")
  expect_length(cent$clique, hg$n_nodes)
  expect_true(all(is.finite(cent$clique)))
})

# ---- pagerank (EDVW random walk, Chitra & Raphael 2019) -------------------

test_that("pagerank: symmetric two-node hypergraph is exactly uniform", {
  hg <- window_hypergraph(list(c("a", "b")), window = 2L)
  pr <- hypergraph_centrality(hg, type = "pagerank")
  expect_named(pr, c("node", "pagerank"))
  expect_equal(.hc_vec(pr, "pagerank"), c(a = 0.5, b = 0.5),
               tolerance = 1e-12)
})

test_that("pagerank matches a dense linear-system solve on EDVW walks", {
  set.seed(7)
  results <- lapply(1:10, function(i) {
    hg <- window_hypergraph(
      list(sample(letters[1:5], 60, replace = TRUE)), window = 3L)
    pr <- hypergraph_centrality(hg, type = "pagerank", damping = 0.85, tol = 1e-12)
    rw <- hypernets:::.hl_rw_transition(hg)
    n <- hg$n_nodes
    ref <- solve(diag(n) - 0.85 * t(rw$P), rep(0.15 / n, n))
    expect_equal(.hc_vec(pr, "pagerank"), ref, tolerance = 1e-7,
                 ignore_attr = TRUE)
    expect_equal(sum(pr$pagerank), 1, tolerance = 1e-12)
    expect_true(all(pr$pagerank > 0))
    hg
  })
  expect_length(results, 10L)
})

test_that("pagerank collapse theorem: binary incidence equals igraph on
           the weighted clique expansion", {
  skip_if_not_installed("igraph")
  df <- data.frame(
    p = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
    s = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(df, actor = "p", group = "s")
  pr <- hypergraph_centrality(hg, type = "pagerank")
  # Edge-independent case: W[u,v] = sum over shared edges of w(e)/delta(e),
  # self-loops included (the walk may stay in place)
  pat <- (hg$incidence > 0) * 1.0
  delta <- colSums(pat)
  W <- pat %*% (rep(1, hg$n_hyperedges) / delta * t(pat))
  g <- igraph::graph_from_adjacency_matrix(W, mode = "directed",
                                           weighted = TRUE, diag = TRUE)
  ref <- igraph::page_rank(g, damping = 0.85)$vector
  expect_equal(.hc_vec(pr, "pagerank"), ref, tolerance = 1e-6)
})

test_that("pagerank uses window counts by default and reacts to weights", {
  set.seed(11)
  traj <- c(rep(c("a", "b"), 6), "c", rep(c("a", "c"), 3))
  hg <- window_hypergraph(list(traj), window = 2L)
  pr_default <- hypergraph_centrality(hg, type = "pagerank")
  pr_counts <- hypergraph_centrality(
    hg, type = "pagerank",
    edge_weights = as.numeric(hg$window_counts))
  pr_unit <- hypergraph_centrality(
    hg, type = "pagerank", edge_weights = rep(1, hg$n_hyperedges))
  expect_identical(pr_default, pr_counts)
  expect_false(isTRUE(all.equal(pr_default, pr_unit)))
})

test_that("pagerank is invariant under state relabeling", {
  set.seed(41)
  traj <- sample(letters[1:4], 30, replace = TRUE)
  relabel <- c(a = "w", b = "x", c = "y", d = "z")
  hg1 <- window_hypergraph(list(traj), window = 3L)
  pr1 <- hypergraph_centrality(hg1, type = "pagerank")
  hg2 <- window_hypergraph(list(unname(relabel[traj])), window = 3L)
  pr2 <- hypergraph_centrality(hg2, type = "pagerank")
  expect_equal(unname(pr1$pagerank), unname(pr2$pagerank),
               tolerance = 1e-12)
})

test_that("pagerank handles nodes filtered out of every hyperedge", {
  hg <- window_hypergraph(list(c("a", "b", "a", "b", "a", "b", "z", "q")),
                          window = 2L, min_weight = 2L)
  pr <- hypergraph_centrality(hg, type = "pagerank")
  expect_equal(sum(pr$pagerank), 1, tolerance = 1e-12)
  expect_true(all(is.finite(pr$pagerank)))
  # dangling nodes keep at least the raw teleportation share
  expect_true(all(pr$pagerank >= (1 - 0.85) / hg$n_nodes))
})

test_that("pagerank degenerate and error paths", {
  adj <- matrix(0, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  hg0 <- build_hypergraph(adj, p = 1, include_pairwise = FALSE)
  pr0 <- hypergraph_centrality(hg0, type = "pagerank")
  expect_equal(unname(pr0$pagerank), rep(1 / 3, 3))
  hg <- window_hypergraph(list(c("a", "b", "c")), window = 2L)
  expect_error(hypergraph_centrality(hg, type = "pagerank", damping = 0),
               "damping")
  expect_error(hypergraph_centrality(hg, type = "pagerank", damping = 1),
               "damping")
  expect_error(hypergraph_centrality(hg, type = "pagerank",
                                     edge_weights = c(1, -1)),
               "edge_weights")
  expect_warning(
    hypergraph_centrality(hg, type = "pagerank", max_iter = 1L),
    class = "hypernets_no_converge"
  )
})

test_that("scalar edge_weights recycles across the Laplacian family too", {
  hg <- window_hypergraph(list(c("a", "b", "a", "b", "c", "a")), window = 2L)
  m <- hg$n_hyperedges
  pr_scalar <- hypergraph_centrality(hg, type = "pagerank", edge_weights = 1)
  pr_vector <- hypergraph_centrality(hg, type = "pagerank",
                                     edge_weights = rep(1, m))
  expect_identical(pr_scalar, pr_vector)
  laplacian_scalar <- hypergraph_laplacian(hg, edge_weights = 2)
  laplacian_vector <- hypergraph_laplacian(hg, edge_weights = rep(2, m))
  expect_identical(laplacian_scalar, laplacian_vector)
  expect_error(
    hypergraph_centrality(hg, type = "pagerank", edge_weights = -1),
    "edge_weights")
})
