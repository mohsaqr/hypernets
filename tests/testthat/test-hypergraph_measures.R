# ---- .hg_measures_fit() tests -----------------------------------------

# Helpers ------------------------------------------------------------------

.hm_two_overlapping <- function() {
  # Hyperedges: e1 = {A,B,C}, e2 = {A,B,D}.  Pair AB appears in both.
  d <- data.frame(
    member  = c("A", "B", "C",  "A", "B", "D"),
    session = c("S1", "S1", "S1", "S2", "S2", "S2"),
    stringsAsFactors = FALSE
  )
  group_hypergraph(d, "member", "session")
}

# Class + structure --------------------------------------------------------

test_that("returns a .hg_measures_fit object with expected fields", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_s3_class(m, "net_hg_measures")
  expect_named(m, c("hyperdegree", "node_strength", "max_edge_size",
                    "co_degree", "edge_sizes", "edge_pairwise_overlap",
                    "overlap_coefficient", "jaccard", "density",
                    "avg_edge_size", "size_distribution",
                    "intersection_profile", "pairwise_participation",
                    "n_nodes", "n_hyperedges"))
})

# Node-level measures ------------------------------------------------------

test_that("hyperdegree counts hyperedges containing each node", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_equal(m$hyperdegree[["A"]], 2L)  # in S1 and S2
  expect_equal(m$hyperdegree[["B"]], 2L)
  expect_equal(m$hyperdegree[["C"]], 1L)  # only S1
  expect_equal(m$hyperdegree[["D"]], 1L)  # only S2
})

test_that("node_strength = sum of edge sizes containing node", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  # Both S1 and S2 are size 3, so A in both => 6; C only in S1 (size 3) => 3
  expect_equal(m$node_strength[["A"]], 6)
  expect_equal(m$node_strength[["B"]], 6)
  expect_equal(m$node_strength[["C"]], 3)
  expect_equal(m$node_strength[["D"]], 3)
})

test_that("max_edge_size = largest hyperedge containing each node", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_true(all(m$max_edge_size == 3L))
})

# Co-degree (n x n) --------------------------------------------------------

test_that("co_degree[i,j] = #hyperedges containing both i and j", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_equal(m$co_degree["A", "B"], 2)  # AB in both S1, S2
  expect_equal(m$co_degree["A", "C"], 1)  # only S1
  expect_equal(m$co_degree["A", "D"], 1)  # only S2
  expect_equal(m$co_degree["C", "D"], 0)  # never together
  expect_true(all(diag(m$co_degree) == 0))
  expect_true(isSymmetric(m$co_degree))
})

# Hyperedge-level measures -------------------------------------------------

test_that("edge_sizes equals colSums(incidence)", {
  hg <- .hm_two_overlapping()
  m  <- .hg_measures_fit(hg)
  expect_equal(m$edge_sizes, c(3L, 3L))
})

test_that("edge_pairwise_overlap = |e_i and e_j|", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  # S1 and S2 share {A, B} => overlap 2
  expect_equal(m$edge_pairwise_overlap["S1", "S2"], 2)
  expect_equal(m$edge_pairwise_overlap["S2", "S1"], 2)
  expect_true(all(diag(m$edge_pairwise_overlap) == 0))
})

test_that("jaccard = overlap / union, overlap_coefficient = overlap / min", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  # |S1| = |S2| = 3, overlap = 2, union = 4
  expect_equal(m$jaccard["S1", "S2"], 2 / 4)
  expect_equal(m$overlap_coefficient["S1", "S2"], 2 / 3)
})

# Global measures ----------------------------------------------------------

test_that("density: uniform 3-uniform = m / choose(n, 3)", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_equal(m$density, 2 / choose(4L, 3L))
})

test_that("density: mixed sizes = sum(|e|) / (n * m)", {
  d <- data.frame(
    member = c("A", "B", "C", "A", "B"),
    grp    = c("g1", "g1", "g1", "g2", "g2"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  m <- .hg_measures_fit(hg)
  # n=3, m=2, sum(|e|) = 3 + 2 = 5
  expect_equal(m$density, 5 / (3 * 2))
})

test_that("avg_edge_size = mean of edge_sizes", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_equal(m$avg_edge_size, 3)
})

test_that("pairwise_participation = fraction of node pairs co-appearing", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  # 4 nodes => 6 pairs. AB AC AD BC BD = 5 pairs co-appear; CD does not.
  expect_equal(m$pairwise_participation, 5 / 6)
})

# Intersection profile -----------------------------------------------------

test_that("intersection_profile tabulates upper-tri overlap sizes", {
  d <- data.frame(
    member = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
    grp    = c("g1", "g1", "g1", "g2", "g2", "g3", "g3", "g3", "g3"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  m  <- .hg_measures_fit(hg)
  # g1 ∩ g2 = {A,B} (2);  g1 ∩ g3 = {C} (1);  g2 ∩ g3 = {} (0)
  expect_equal(m$intersection_profile[["overlap_0"]], 1L)
  expect_equal(m$intersection_profile[["overlap_1"]], 1L)
  expect_equal(m$intersection_profile[["overlap_2"]], 1L)
})

# Empty hypergraph ---------------------------------------------------------

test_that("empty hypergraph returns trivial zeros without error", {
  inc <- matrix(0L, 3L, 0L, dimnames = list(c("A", "B", "C"), NULL))
  hg <- structure(
    list(hyperedges = list(), incidence = inc, nodes = c("A", "B", "C"),
         n_nodes = 3L, n_hyperedges = 0L,
         size_distribution = integer(0), params = list()),
    class = "net_hg")
  m <- .hg_measures_fit(hg)
  expect_equal(m$n_hyperedges, 0L)
  expect_equal(m$density, 0)
  expect_equal(unname(m$hyperdegree), c(0L, 0L, 0L))
  expect_true(all(m$co_degree == 0))
  expect_equal(m$pairwise_participation, 0)
})

# Single hyperedge ---------------------------------------------------------

test_that("single 4-hyperedge: jaccard / overlap matrices are empty", {
  d <- data.frame(
    member = c("A", "B", "C", "D"),
    grp    = c("g1", "g1", "g1", "g1"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  m  <- .hg_measures_fit(hg)
  expect_equal(dim(m$jaccard), c(1L, 1L))
  expect_equal(m$jaccard[1, 1], 0)  # diag is zero
  expect_length(m$intersection_profile, 0L)  # need >=2 hyperedges
})

# Validation against manual matrix algebra --------------------------------

test_that("matches manual incidence-matrix computation", {
  hg <- .hm_two_overlapping()
  m  <- .hg_measures_fit(hg)
  B  <- (hg$incidence > 0) * 1
  expect_equal(m$hyperdegree, stats::setNames(as.integer(rowSums(B)), hg$nodes))
  expect_equal(m$edge_sizes, as.integer(colSums(B)))
  manual_co <- tcrossprod(B); diag(manual_co) <- 0
  expect_equal(unname(m$co_degree), unname(manual_co))
  manual_op <- crossprod(B); diag(manual_op) <- 0
  expect_equal(unname(m$edge_pairwise_overlap), unname(manual_op))
})

# Input validation --------------------------------------------------------

test_that("rejects non-net_hg input", {
  expect_error(.hg_measures_fit(matrix(0, 3, 3)), "net_hg")
  expect_error(.hg_measures_fit(list()), "net_hg")
})

# Print method ------------------------------------------------------------

test_that("print method runs without error and returns object invisibly", {
  m <- .hg_measures_fit(.hm_two_overlapping())
  expect_invisible(print(m))
})

test_that("print works on empty hypergraph", {
  inc <- matrix(0L, 2L, 0L, dimnames = list(c("A", "B"), NULL))
  hg <- structure(
    list(hyperedges = list(), incidence = inc, nodes = c("A", "B"),
         n_nodes = 2L, n_hyperedges = 0L,
         size_distribution = integer(0), params = list()),
    class = "net_hg")
  expect_invisible(print(.hg_measures_fit(hg)))
})
