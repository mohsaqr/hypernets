# ---- hg_clique_expansion() tests --------------------------------------------

# Helpers ------------------------------------------------------------------

.ce_hg_one_triangle <- function() {
  # Single 3-hyperedge (A,B,C). When include_pairwise=FALSE there is only
  # the triangle; pairs AB, AC, BC each get weight 1.
  inc <- matrix(c(1L, 1L, 1L), 3L, 1L,
                dimnames = list(c("A", "B", "C"), "h1"))
  structure(
    list(
      hyperedges        = list(1:3),
      incidence         = inc,
      nodes             = c("A", "B", "C"),
      n_nodes           = 3L,
      n_hyperedges      = 1L,
      size_distribution = setNames(1L, "size_3"),
      params            = list(source = "test")
    ),
    class = "net_hg"
  )
}

# Structure ----------------------------------------------------------------

test_that("returns a netobject + cograph_network", {
  net <- hg_clique_expansion(.ce_hg_one_triangle())
  expect_s3_class(net, "netobject")
  expect_s3_class(net, "cograph_network")
  expect_identical(net$method, "hg_clique_expansion")
  expect_false(net$directed)
})

# Weight semantics: single triangle ----------------------------------------

test_that("single triangle becomes K3 with all unit weights", {
  net <- hg_clique_expansion(.ce_hg_one_triangle())
  W <- net$weights
  expect_equal(W["A", "B"], 1)
  expect_equal(W["A", "C"], 1)
  expect_equal(W["B", "C"], 1)
  expect_true(isSymmetric(W))
  expect_true(all(diag(W) == 0))
})

# Weight semantics: overlapping hyperedges ---------------------------------

test_that("pair shared by k hyperedges has weight k", {
  # Two hyperedges: (A,B,C) and (A,B,D). AB shared by both => weight 2.
  inc <- matrix(0L, 4L, 2L,
                dimnames = list(c("A", "B", "C", "D"), c("h1", "h2")))
  inc[c("A", "B", "C"), "h1"] <- 1L
  inc[c("A", "B", "D"), "h2"] <- 1L
  hg <- structure(
    list(hyperedges = list(1:3, c(1L, 2L, 4L)),
         incidence = inc, nodes = c("A", "B", "C", "D"),
         n_nodes = 4L, n_hyperedges = 2L,
         size_distribution = setNames(2L, "size_3"),
         params = list()),
    class = "net_hg")
  W <- hg_clique_expansion(hg)$weights
  expect_equal(W["A", "B"], 2)  # shared by both
  expect_equal(W["A", "C"], 1)  # only h1
  expect_equal(W["B", "C"], 1)  # only h1
  expect_equal(W["A", "D"], 1)  # only h2
  expect_equal(W["B", "D"], 1)  # only h2
  expect_equal(W["C", "D"], 0)  # never co-occur
  expect_true(isSymmetric(W))
})

# Pairwise hyperedges --------------------------------------------------- ---

test_that("size-2 hyperedges become single edges", {
  inc <- matrix(0L, 3L, 2L,
                dimnames = list(c("A", "B", "C"), c("h1", "h2")))
  inc[c("A", "B"), "h1"] <- 1L
  inc[c("B", "C"), "h2"] <- 1L
  hg <- structure(
    list(hyperedges = list(1:2, 2:3), incidence = inc,
         nodes = c("A", "B", "C"), n_nodes = 3L, n_hyperedges = 2L,
         size_distribution = setNames(2L, "size_2"),
         params = list()),
    class = "net_hg")
  W <- hg_clique_expansion(hg)$weights
  expect_equal(W["A", "B"], 1)
  expect_equal(W["B", "C"], 1)
  expect_equal(W["A", "C"], 0)
})

# Empty hypergraph ---------------------------------------------------------

test_that("empty hypergraph -> all-zero adjacency", {
  inc <- matrix(0L, 3L, 0L, dimnames = list(c("A", "B", "C"), NULL))
  hg <- structure(
    list(hyperedges = list(), incidence = inc,
         nodes = c("A", "B", "C"), n_nodes = 3L, n_hyperedges = 0L,
         size_distribution = integer(0), params = list()),
    class = "net_hg")
  net <- hg_clique_expansion(hg)
  expect_equal(sum(net$weights), 0)
  expect_equal(dim(net$weights), c(3L, 3L))
  expect_equal(net$n_edges, 0L)
})

# Weighted vs unweighted incidence ----------------------------------------

test_that("weighted=FALSE binarises incidence before projecting", {
  inc <- matrix(0, 3L, 1L, dimnames = list(c("A", "B", "C"), "h1"))
  inc[, 1] <- c(2, 5, 3)  # weighted incidence
  hg <- structure(
    list(hyperedges = list(1:3), incidence = inc,
         nodes = c("A", "B", "C"), n_nodes = 3L, n_hyperedges = 1L,
         size_distribution = setNames(1L, "size_3"),
         params = list()),
    class = "net_hg")
  W_w  <- hg_clique_expansion(hg, weighted = TRUE)$weights
  W_uw <- hg_clique_expansion(hg, weighted = FALSE)$weights
  # Weighted: W[A,B] = 2*5 = 10
  expect_equal(W_w["A", "B"], 10)
  expect_equal(W_w["A", "C"], 6)
  expect_equal(W_w["B", "C"], 15)
  # Unweighted: just count of shared hyperedges
  expect_equal(W_uw["A", "B"], 1)
  expect_equal(W_uw["A", "C"], 1)
  expect_equal(W_uw["B", "C"], 1)
})

# Roundtrip with group_hypergraph -----------------------------------------

test_that("group_hypergraph -> hg_clique_expansion is consistent", {
  # 3 sessions: S1=(A,B,C), S2=(A,B), S3=(B,C,D)
  d <- data.frame(
    member  = c("A", "B", "C",  "A", "B",  "B", "C", "D"),
    session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "session")
  W  <- hg_clique_expansion(hg)$weights
  # AB shared by S1 and S2 => 2;  BC by S1 and S3 => 2;  CD by S3 => 1
  expect_equal(W["A", "B"], 2)
  expect_equal(W["B", "C"], 2)
  expect_equal(W["C", "D"], 1)
  expect_equal(W["A", "D"], 0)  # never together
})

# Roundtrip: network_hypergraph then expand ---------------------------------

test_that("network_hypergraph(p=0) -> hg_clique_expansion preserves the graph", {
  # Pure pairwise (p=0) hypergraph from a network. Expanding back must
  # match the binarised input adjacency.
  set.seed(7)
  n <- 6
  adj <- matrix(stats::rbinom(n * n, 1, 0.4), n, n)
  diag(adj) <- 0
  adj <- ((adj + t(adj)) > 0) * 1
  storage.mode(adj) <- "double"
  rownames(adj) <- colnames(adj) <- LETTERS[seq_len(n)]
  hg  <- network_hypergraph(adj, p = 0, include_pairwise = TRUE)
  W   <- hg_clique_expansion(hg)$weights
  expect_equal(W, adj, ignore_attr = TRUE)
})

# Input validation --------------------------------------------------------

test_that("rejects non-net_hg input", {
  expect_error(hg_clique_expansion(matrix(0, 3, 3)),
               "net_hg")
  expect_error(hg_clique_expansion(list()), "net_hg")
})

test_that("rejects non-logical weighted argument", {
  expect_error(hg_clique_expansion(.ce_hg_one_triangle(), weighted = "yes"))
})

# Records provenance in $params -------------------------------------------

test_that("params records source, weighted flag, hyperedge count + size dist", {
  net <- hg_clique_expansion(.ce_hg_one_triangle(), weighted = TRUE)
  expect_equal(net$params$source, "hg_clique_expansion")
  expect_true(net$params$weighted)
  expect_equal(net$params$n_hyperedges, 1L)
  expect_equal(unname(net$params$hypergraph_size_distribution), 1L)
})
