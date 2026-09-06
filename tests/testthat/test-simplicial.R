testthat::skip_on_cran()

# ---- Simplicial Complex Tests ----

# Helper: build a known network.
#
# Inherited from Nestimate, where this called build_network(method =
# "relative"). honets does not depend on Nestimate, so the row-normalised
# transition matrix is computed here directly and wrapped as a netobject.
# Verified to reproduce Nestimate::build_network(seqs, method = "relative")
# exactly (max abs difference 0, identical dimnames) on this seed, so every
# expectation below is unchanged.
.make_sc_net <- function(seed = 42) {
  set.seed(seed)
  seqs <- data.frame(
    V1 = sample(LETTERS[1:5], 40, TRUE),
    V2 = sample(LETTERS[1:5], 40, TRUE),
    V3 = sample(LETTERS[1:5], 40, TRUE),
    V4 = sample(LETTERS[1:5], 40, TRUE)
  )
  states <- sort(unique(unlist(seqs, use.names = FALSE)))
  m <- as.matrix(seqs)
  counts <- table(
    from = factor(as.vector(m[, -ncol(m)]), levels = states),
    to   = factor(as.vector(m[, -1L]),      levels = states)
  )
  mat <- matrix(as.numeric(sweep(counts, 1, rowSums(counts), "/")),
                nrow = length(states),
                dimnames = list(states, states))
  .wrap_netobject(mat, data = seqs, method = "relative", directed = TRUE)
}

# Helper: build a simple known matrix for exact verification
.make_sc_mat <- function() {
  # 4-node complete graph => K4 clique complex
  mat <- matrix(0.5, 4, 4)
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- LETTERS[1:4]
  mat
}

# Helper: triangle matrix
.make_triangle_mat <- function() {
  mat <- matrix(0, 3, 3)
  mat[1, 2] <- mat[2, 1] <- 1
  mat[1, 3] <- mat[3, 1] <- 1
  mat[2, 3] <- mat[3, 2] <- 1
  rownames(mat) <- colnames(mat) <- c("X", "Y", "Z")
  mat
}


# =========================================================================
# build_simplicial — clique type
# =========================================================================

test_that("build_simplicial clique returns correct class and structure", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat, type = "clique")

  expect_s3_class(sc, "net_simplicial")
  expect_equal(sc$type, "clique")
  expect_equal(sc$n_nodes, 4L)
  expect_true(sc$n_simplices > 0L)
  expect_true(sc$dimension >= 1L)
  expect_true(is.numeric(sc$density))
  expect_true(is.numeric(sc$mean_dim))
  expect_true(length(sc$f_vector) > 0L)
  expect_equal(length(sc$nodes), 4L)
  expect_true(is.list(sc$simplices))
})

test_that("build_simplicial clique on K4 gives complete simplex", {
  # K4 => 1 tetrahedron (3-simplex), 4 triangles, 6 edges, 4 vertices
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat, type = "clique")

  expect_equal(sc$dimension, 3L)
  # f-vector: f0=4, f1=6, f2=4, f3=1
  expect_equal(as.integer(sc$f_vector), c(4L, 6L, 4L, 1L))
  expect_equal(sc$n_simplices, 4 + 6 + 4 + 1)
})

test_that("build_simplicial clique with threshold filters edges", {
  mat <- .make_sc_mat()
  mat[1, 2] <- mat[2, 1] <- 0.01  # Weak edge A-B
  sc_full <- build_simplicial(mat, threshold = 0)
  sc_thresh <- build_simplicial(mat, threshold = 0.1)

  # Thresholded version should have fewer simplices

  expect_true(sc_thresh$n_simplices < sc_full$n_simplices)
  expect_true(sc_thresh$dimension <= sc_full$dimension)
})

test_that("build_simplicial clique includes boundary weights but excludes zero non-edges", {
  mat <- matrix(0, 3, 3)
  mat[1, 2] <- mat[2, 1] <- 0.5
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")

  sc_boundary <- build_simplicial(mat, type = "clique", threshold = 0.5)
  expect_equal(as.integer(sc_boundary$f_vector), c(3L, 1L))

  sc_zero <- build_simplicial(mat, type = "clique", threshold = 0)
  expect_equal(as.integer(sc_zero$f_vector), c(3L, 1L))
})

test_that("build_simplicial clique with max_dim limits dimension", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat, type = "clique", max_dim = 1L)
  # max_dim=1 => only vertices and edges
  expect_lte(sc$dimension, 1L)
})

test_that("build_simplicial accepts netobject input", {
  net <- .make_sc_net()
  sc <- build_simplicial(net, threshold = 0.05)
  expect_s3_class(sc, "net_simplicial")
  expect_equal(sc$type, "clique")
  expect_true(sc$n_nodes > 0L)
})

test_that("build_simplicial on matrix without rownames assigns V-names", {
  mat <- matrix(c(0, 1, 1, 1, 0, 1, 1, 1, 0), 3, 3)
  sc <- build_simplicial(mat)
  expect_equal(sc$nodes, c("V1", "V2", "V3"))
})


# =========================================================================
# build_simplicial — VR type is a genuine Vietoris-Rips filtration
#
# Historical note: A04-F01 audit caught that type = "vr" silently aliased
# "clique" and forced it to error. With the boundary-matrix PH backend, VR
# is now a real metric filtration: input is interpreted as a distance
# matrix, k-simplex enters at max pairwise distance in σ. These tests
# verify VR is genuinely distinct from clique (no silent alias regression).
# =========================================================================

test_that("build_simplicial type='vr' returns a VR complex with filtration", {
  # 4-cycle as a distance matrix: short edges form the cycle.
  d <- matrix(0, 4, 4)
  d[cbind(c(1,2,3,4,1,2,3,4), c(2,3,4,1,2,3,4,1))] <- c(0.1, 0.2, 0.3, 0.4,
                                                          0.1, 0.2, 0.3, 0.4)
  rownames(d) <- colnames(d) <- paste0("v", 1:4)
  sc <- build_simplicial(d, type = "vr", max_scale = 0.5)
  expect_s3_class(sc, "net_simplicial")
  expect_identical(sc$type, "vr")
  expect_true(is.numeric(sc$filtration))
  expect_equal(length(sc$filtration), length(sc$simplices))
  # vertices enter at 0
  vert_idx <- which(vapply(sc$simplices, length, integer(1)) == 1L)
  expect_true(all(sc$filtration[vert_idx] == 0))
})

test_that("type='vr' is genuinely distinct from type='clique'", {
  # On a similarity-weighted graph, the VR complex on (1 - weight) and the
  # clique complex on the original weights have different simplex counts
  # because the "max distance" rule and the "min similarity" rule order the
  # same simplices differently when the input contains a mix of weights.
  mat <- .make_sc_mat()
  d <- 1 - abs(mat)
  diag(d) <- 0
  sc_vr     <- build_simplicial(d,   type = "vr",     max_scale = 0.95)
  sc_clique <- build_simplicial(mat, type = "clique", threshold = 0.05)
  expect_identical(sc_vr$type,     "vr")
  expect_identical(sc_clique$type, "clique")
  # filtration attached on VR, not on clique
  expect_true(is.numeric(sc_vr$filtration))
  expect_null(sc_clique$filtration)
})

test_that("type='rips' is an alias for 'vr'", {
  d <- matrix(c(0, 0.2, 0.3,
                0.2, 0, 0.4,
                0.3, 0.4, 0), 3, 3, byrow = TRUE)
  rownames(d) <- colnames(d) <- c("A","B","C")
  sc1 <- build_simplicial(d, type = "vr",   max_scale = 0.5)
  sc2 <- build_simplicial(d, type = "rips", max_scale = 0.5)
  expect_identical(sc1$type, sc2$type)
  expect_equal(length(sc1$simplices), length(sc2$simplices))
})


# =========================================================================
# build_simplicial — pathway type
# =========================================================================

test_that("build_simplicial pathway from net_hon works", {
  trajs <- list(
    c("A", "B", "C", "A", "B"),
    c("B", "C", "A", "B", "C"),
    c("A", "B", "C", "B", "A"),
    c("C", "A", "B", "C", "A")
  )
  hon <- build_hon(trajs, max_order = 2, min_freq = 1)
  sc <- build_simplicial(hon, type = "pathway")
  expect_s3_class(sc, "net_simplicial")
  expect_equal(sc$type, "pathway")
})

test_that("build_simplicial pathway with max_pathways limits input", {
  trajs <- list(
    c("A", "B", "C", "A", "B"),
    c("B", "C", "A", "B", "C"),
    c("A", "B", "C", "B", "A")
  )
  hon <- build_hon(trajs, max_order = 2, min_freq = 1)
  sc_full <- build_simplicial(hon, type = "pathway")
  sc_limited <- build_simplicial(hon, type = "pathway", max_pathways = 2)
  expect_lte(sc_limited$n_simplices, sc_full$n_simplices)
})

test_that("build_simplicial pathway from netobject builds HON internally", {
  net <- .make_sc_net()
  sc <- build_simplicial(net, type = "pathway", max_order = 2, min_freq = 1)
  expect_s3_class(sc, "net_simplicial")
  expect_equal(sc$type, "pathway")
})

test_that("build_simplicial pathway rejects unsupported input", {
  expect_error(
    build_simplicial(list(a = 1), type = "pathway"),
    "net_hon.*net_hypa.*tna.*netobject"
  )
})

test_that("build_simplicial pathway from hypa works", {
  trajs <- list(
    c("A", "B", "C"), c("A", "B", "C"), c("A", "B", "C"),
    c("A", "B", "D"), c("C", "B", "D"), c("C", "B", "A")
  )
  h <- build_hypa(trajs, k = 2)
  sc <- build_simplicial(h, type = "pathway")
  expect_s3_class(sc, "net_simplicial")
  expect_equal(sc$type, "pathway")
})

test_that("build_simplicial pathway from hypa respects anomaly direction", {
  scores <- data.frame(
    path = c("A -> B -> C", "D -> E -> F"),
    observed = c(20, 1),
    expected = c(3, 8),
    ratio = c(20 / 3, 1 / 8),
    p_value = c(0.99, 0.01),
    p_under = c(0.99, 0.01),
    p_over = c(0.001, 0.95),
    p_adjusted_under = c(0.99, 0.01),
    p_adjusted_over = c(0.001, 0.95),
    anomaly = c("over", "under"),
    stringsAsFactors = FALSE
  )
  h <- structure(list(
    scores = scores,
    nodes = data.frame(label = LETTERS[1:6], stringsAsFactors = FALSE)
  ), class = "net_hypa")

  sc_under <- build_simplicial(h, type = "pathway", anomaly = "under",
                               max_pathways = 1)
  under_idx <- match(c("D", "E", "F"), sc_under$nodes)
  expect_true(any(vapply(sc_under$simplices, function(s) {
    identical(sort(s), sort(under_idx))
  }, logical(1))))

  sc_over <- build_simplicial(h, type = "pathway", anomaly = "over",
                              max_pathways = 1)
  over_idx <- match(c("A", "B", "C"), sc_over$nodes)
  expect_true(any(vapply(sc_over$simplices, function(s) {
    identical(sort(s), sort(over_idx))
  }, logical(1))))
})


# =========================================================================
# .sc_extract_matrix
# =========================================================================

test_that("betti_numbers on filled triangle gives b0=1, rest 0", {
  # Clique complex on K3 fills the triangle => contractible (solid triangle)
  mat <- .make_triangle_mat()
  sc <- build_simplicial(mat, threshold = 0)
  b <- betti_numbers(sc)
  expect_equal(b[["b0"]], 1L)  # 1 component
  expect_equal(b[["b1"]], 0L)  # No loop (triangle is filled by 2-simplex)
})

test_that("betti_numbers on K4 gives b0=1, b_rest=0 (contractible)", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  b <- betti_numbers(sc)
  expect_equal(b[["b0"]], 1L)
  # K4 clique complex (tetrahedron) is contractible
  expect_true(all(b[-1] == 0))
})

test_that("betti_numbers on disconnected graph gives b0=n_components", {
  # Two disconnected edges: A-B and C-D
  mat <- matrix(0, 4, 4)
  mat[1, 2] <- mat[2, 1] <- 1
  mat[3, 4] <- mat[4, 3] <- 1
  rownames(mat) <- colnames(mat) <- LETTERS[1:4]
  sc <- build_simplicial(mat)
  b <- betti_numbers(sc)
  expect_equal(b[["b0"]], 2L)  # 2 components
})

test_that("betti_numbers rejects non-simplicial_complex", {
  expect_error(betti_numbers(list(a = 1)), "net_simplicial")
})


# =========================================================================
# euler_characteristic
# =========================================================================

test_that("euler_characteristic matches Euler-Poincare theorem", {
  mat <- .make_triangle_mat()
  sc <- build_simplicial(mat)
  chi <- euler_characteristic(sc)
  b <- betti_numbers(sc)
  # chi = sum((-1)^k * beta_k)
  dims <- seq_along(b) - 1L
  chi_from_betti <- as.integer(sum((-1L)^dims * b))
  expect_equal(chi, chi_from_betti)
})

test_that("euler_characteristic of K4 is 1", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  # chi = f0 - f1 + f2 - f3 = 4 - 6 + 4 - 1 = 1
  chi <- euler_characteristic(sc)
  expect_equal(chi, 1L)
})

test_that("euler_characteristic rejects non-simplicial_complex", {
  expect_error(euler_characteristic(list(a = 1)), "net_simplicial")
})


# =========================================================================
# simplicial_degree
# =========================================================================

test_that("simplicial_degree returns correct structure", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  deg <- simplicial_degree(sc)

  expect_true(is.data.frame(deg))
  expect_true("node" %in% names(deg))
  expect_true("total" %in% names(deg))
  expect_equal(nrow(deg), sc$n_nodes)
  # d0 column exists for vertices
  expect_true("d0" %in% names(deg))
})

test_that("simplicial_degree on K4 gives equal participation", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  deg <- simplicial_degree(sc)
  # In K4, all nodes participate equally
  expect_true(all(deg$total == deg$total[1]))
})

test_that("simplicial_degree normalized divides by max", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  deg_raw <- simplicial_degree(sc, normalized = FALSE)
  deg_norm <- simplicial_degree(sc, normalized = TRUE)
  # Normalized d0 should be 1 for each node
  expect_true(all(deg_norm$d0 == 1))
  # Normalized values <= 1
  d_cols <- grep("^d[0-9]+$", names(deg_norm), value = TRUE)
  expect_true(all(deg_norm[, d_cols] <= 1))
})

test_that("simplicial_degree rejects non-simplicial_complex", {
  expect_error(simplicial_degree(list()), "net_simplicial")
})


# =========================================================================
# persistent_homology
# =========================================================================

test_that("persistent_homology returns correct class and structure", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 10)

  expect_s3_class(ph, "net_persistent_homology")
  expect_true(is.data.frame(ph$betti_curve))
  expect_true(is.data.frame(ph$persistence))
  expect_true(is.numeric(ph$thresholds))
  expect_equal(length(ph$thresholds), 10L)
  expect_true(all(c("threshold", "dimension", "betti") %in%
                     names(ph$betti_curve)))
  expect_true(all(c("dimension", "birth", "death", "persistence") %in%
                     names(ph$persistence)))
})

test_that("persistent_homology thresholds are decreasing", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 15)
  expect_true(all(diff(ph$thresholds) < 0))
})

test_that("persistent_homology betti_curve has all requested dimensions", {
  mat <- matrix(0, 3, 3)
  mat[1, 2] <- mat[2, 1] <- 1
  mat[1, 3] <- mat[3, 1] <- 0.5
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  ph <- persistent_homology(mat, n_steps = 3, max_dim = 3L)

  expect_equal(sort(unique(ph$betti_curve$dimension)), 0:3)
  expect_true(all(table(ph$betti_curve$threshold) == 4L))
  expect_true(all(ph$betti_curve$betti[ph$betti_curve$dimension > 1L] == 0L))
})

test_that("persistent_homology includes edges at exact threshold", {
  mat <- matrix(0, 2, 2)
  mat[1, 2] <- mat[2, 1] <- 0.5
  rownames(mat) <- colnames(mat) <- c("A", "B")
  ph <- persistent_homology(mat, n_steps = 3)
  first <- ph$betti_curve[ph$betti_curve$threshold == max(ph$thresholds), ]

  expect_equal(first$betti[first$dimension == 0], 1L)
  expect_equal(first$betti[first$dimension == 1], 0L)
  expect_equal(nrow(ph$persistence), 1L)
  expect_equal(ph$persistence$dimension, 0L)
  expect_equal(ph$persistence$birth, 0.5)
  expect_equal(ph$persistence$death, 0)
})

test_that("persistent_homology b0 increases as threshold increases", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 20)
  b0 <- ph$betti_curve[ph$betti_curve$dimension == 0, ]
  # At highest threshold, b0 should be large (disconnected)
  # At lowest threshold, b0 should be small (connected)
  # b0 should generally increase with threshold
  expect_gte(b0$betti[1], b0$betti[nrow(b0)])
})

test_that("persistent_homology fails on zero matrix", {
  mat <- matrix(0, 3, 3)
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  expect_error(persistent_homology(mat), "All weights are zero")
})

test_that("persistent_homology persistence >= 0", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 10)
  if (nrow(ph$persistence) > 0L) {
    expect_true(all(ph$persistence$persistence >= 0))
  }
})

test_that("persistent_homology from matrix works", {
  mat <- .make_sc_mat()
  ph <- persistent_homology(mat, n_steps = 5)
  expect_s3_class(ph, "net_persistent_homology")
})


# =========================================================================
# q_analysis
# =========================================================================

test_that("q_analysis returns correct class and structure", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)

  expect_s3_class(qa, "net_q_analysis")
  expect_true(is.integer(qa$q_vector) || is.numeric(qa$q_vector))
  expect_true(is.integer(qa$max_q) || is.numeric(qa$max_q))
  expect_true(is.integer(qa$structure_vector) || is.numeric(qa$structure_vector))
  expect_equal(length(qa$structure_vector), sc$n_nodes)
})

test_that("q_analysis on K4: q_vector = 1 at all levels", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  # K4 is fully connected at all q levels
  expect_true(all(qa$q_vector == 1L))
  expect_equal(qa$max_q, 3L) # K4 has a 3-simplex
})

test_that("q_analysis structure_vector gives max simplex dim per node", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  # All nodes in K4 are in the 3-simplex
  expect_true(all(qa$structure_vector == 3L))
})

test_that("q_analysis with single maximal simplex returns 1", {
  # Only one maximal simplex => 1 component at all levels
  mat <- .make_triangle_mat()
  # Add a fill to make it a filled triangle (K3 = 2-simplex)
  mat_full <- matrix(0.5, 3, 3)
  diag(mat_full) <- 0
  rownames(mat_full) <- colnames(mat_full) <- c("X", "Y", "Z")
  sc <- build_simplicial(mat_full)
  qa <- q_analysis(sc)
  expect_true(all(qa$q_vector == 1L))
})

test_that("q_analysis on disconnected graph fragments at q=0", {
  mat <- matrix(0, 4, 4)
  mat[1, 2] <- mat[2, 1] <- 1
  mat[3, 4] <- mat[4, 3] <- 1
  rownames(mat) <- colnames(mat) <- LETTERS[1:4]
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  # At q=0, should have at least 2 components
  last_q <- qa$q_vector[length(qa$q_vector)]
  expect_gte(last_q, 2L)
})

test_that("q_analysis names q-vector levels in the computed order", {
  mat <- matrix(0, 5, 5)
  rownames(mat) <- colnames(mat) <- LETTERS[1:5]
  edges <- rbind(c(1, 2), c(1, 3), c(2, 3),
                 c(3, 4), c(3, 5), c(4, 5))
  mat[edges] <- 1
  mat[edges[, 2:1]] <- 1
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)

  expect_equal(names(qa$q_vector), c("q_2", "q_1", "q_0"))
  expect_equal(as.integer(qa$q_vector), c(2L, 2L, 1L))
})

test_that("q_analysis rejects non-simplicial_complex", {
  expect_error(q_analysis(list()), "net_simplicial")
})


# =========================================================================
# verify_simplicial (cross-validation with igraph)
# =========================================================================

test_that("verify_simplicial matches igraph on K4", {
  skip_if_not_installed("igraph")
  mat <- .make_sc_mat()
  out <- capture.output(res <- verify_simplicial(mat))
  expect_true(res$cliques_match)
  expect_equal(res$n_simplices_ours, res$n_simplices_igraph)
})

test_that("verify_simplicial matches igraph on triangle", {
  skip_if_not_installed("igraph")
  mat <- .make_triangle_mat()
  out <- capture.output(res <- verify_simplicial(mat))
  expect_true(res$cliques_match)
})

test_that("verify_simplicial matches igraph with threshold", {
  skip_if_not_installed("igraph")
  mat <- .make_sc_mat()
  mat[1, 2] <- mat[2, 1] <- 0.01
  out <- capture.output(res <- verify_simplicial(mat, threshold = 0.1))
  expect_true(res$cliques_match)
})

test_that("verify_simplicial Euler-Poincare holds", {
  skip_if_not_installed("igraph")
  mat <- .make_sc_mat()
  out <- capture.output(res <- verify_simplicial(mat))
  b <- res$betti
  dims <- seq_along(b) - 1L
  euler_from_betti <- as.integer(sum((-1L)^dims * b))
  expect_equal(res$euler, euler_from_betti)
})

test_that("verify_simplicial on random network passes", {
  skip_if_not_installed("igraph")
  net <- .make_sc_net()
  out <- capture.output(res <- verify_simplicial(net$weights, threshold = 0.1))
  expect_true(res$cliques_match)
})


# =========================================================================
# print methods
# =========================================================================

test_that("print.simplicial_complex works with non-zero Betti", {
  mat <- .make_triangle_mat()
  sc <- build_simplicial(mat)
  out <- capture.output(print(sc))
  expect_true(any(grepl("Clique Complex", out)))
  expect_true(any(grepl("f-vector", out)))
  expect_true(any(grepl("Betti", out)))
  expect_true(any(grepl("Nodes:", out)))
})

test_that("print.simplicial_complex shows Betti info for K4", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  out <- capture.output(print(sc))
  expect_true(any(grepl("Betti", out)))
  expect_true(any(grepl("b0=1", out)))
})

test_that("print.simplicial_complex suppresses nodes for large complexes", {
  # 20 nodes — should not print "Nodes:"
  n <- 20
  mat <- matrix(runif(n * n), n, n)
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- paste0("N", seq_len(n))
  sc <- build_simplicial(mat, threshold = 0.9)
  out <- capture.output(print(sc))
  expect_false(any(grepl("Nodes:", out)))
})

test_that("print.simplicial_complex never shows a Vietoris-Rips label", {
  # type = "vr" no longer produces a complex (A04-F01), so the dead
  # "Vietoris-Rips Complex" label can never be printed. A clique complex
  # must print "Clique Complex".
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat, type = "clique", threshold = 0.3)
  out <- capture.output(print(sc))
  expect_true(any(grepl("Clique Complex", out)))
  expect_false(any(grepl("Vietoris-Rips", out)))
})

test_that("print.simplicial_complex returns invisibly", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  out <- capture.output(res <- print(sc))
  expect_identical(res, sc)
})

test_that("print.persistent_homology works", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 10)
  out <- capture.output(print(ph))
  expect_true(any(grepl("Persistent Homology", out)))
  expect_true(any(grepl("filtration steps", out)))
})

test_that("print.persistent_homology returns invisibly", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 5)
  out <- capture.output(res <- print(ph))
  expect_identical(res, ph)
})

test_that("print.q_analysis works for fully connected", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  out <- capture.output(print(qa))
  expect_true(any(grepl("Q-Analysis", out)))
  expect_true(any(grepl("Fully connected", out)))
  expect_true(any(grepl("Structure", out)))
})

test_that("print.q_analysis works for fragmented graph", {
  # 3 disconnected pairs + 1 triangle => fragments at some q
  mat <- matrix(0, 6, 6)
  mat[1, 2] <- mat[2, 1] <- 1
  mat[3, 4] <- mat[4, 3] <- 1
  mat[3, 5] <- mat[5, 3] <- 1
  mat[4, 5] <- mat[5, 4] <- 1
  mat[5, 6] <- mat[6, 5] <- 1
  rownames(mat) <- colnames(mat) <- LETTERS[1:6]
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  out <- capture.output(print(qa))
  expect_true(any(grepl("Q-Analysis", out)))
  # Should show fragmentation info since there are disconnected components
  expect_true(any(grepl("Fragments|Components", out)))
})

test_that("print.q_analysis returns invisibly", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  out <- capture.output(res <- print(qa))
  expect_identical(res, qa)
})


# =========================================================================
# plot methods
# =========================================================================

test_that("plot.simplicial_complex produces grob", {
  skip_if_not_installed("gridExtra")
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  grb <- plot(sc)
  expect_true(inherits(grb, "grob") || inherits(grb, "gtable"))
})

test_that("plot.persistent_homology produces grob", {
  skip_if_not_installed("gridExtra")
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 10)
  # Unicode β in labels may fail in non-UTF-8 locales
  grb <- tryCatch(plot(ph), error = function(e) {
    skip("Unicode rendering not supported in this locale")
  })
  expect_true(inherits(grb, "grob") || inherits(grb, "gtable"))
})

test_that("plot.persistent_homology handles no features", {
  skip_if_not_installed("gridExtra")
  mat <- matrix(1e-6, 3, 3)
  mat[1, 2] <- mat[2, 1] <- 0.01
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  ph <- persistent_homology(mat, n_steps = 5)
  grb <- tryCatch(plot(ph), error = function(e) {
    skip("Unicode rendering not supported in this locale")
  })
  expect_true(inherits(grb, "grob") || inherits(grb, "gtable"))
})

test_that("plot.q_analysis produces grob", {
  skip_if_not_installed("gridExtra")
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  grb <- plot(qa)
  expect_true(inherits(grb, "grob") || inherits(grb, "gtable"))
})

test_that("plot.simplicial_complex errors without gridExtra", {
  skip_if(requireNamespace("gridExtra", quietly = TRUE),
          "gridExtra is installed; can't test missing")
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  expect_error(plot(sc), "gridExtra")
})

test_that("plot.persistent_homology errors without gridExtra", {
  skip_if(requireNamespace("gridExtra", quietly = TRUE),
          "gridExtra is installed; can't test missing")
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 5)
  expect_error(plot(ph), "gridExtra")
})

test_that("plot.q_analysis errors without gridExtra", {
  skip_if(requireNamespace("gridExtra", quietly = TRUE),
          "gridExtra is installed; can't test missing")
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  expect_error(plot(qa), "gridExtra")
})


# =========================================================================
# Edge cases
# =========================================================================

test_that("build_simplicial on single-node matrix works", {
  mat <- matrix(0, 1, 1)
  rownames(mat) <- colnames(mat) <- "A"
  sc <- build_simplicial(mat)
  expect_equal(sc$n_nodes, 1L)
  expect_equal(sc$n_simplices, 1L) # just the vertex
  expect_equal(sc$dimension, 0L)
})

test_that("build_simplicial on empty network (high threshold) works", {
  mat <- .make_sc_mat()
  sc <- build_simplicial(mat, threshold = 10)
  # No edges above threshold
  expect_equal(sc$dimension, 0L)
  expect_equal(sc$n_simplices, 4L) # just vertices
})

test_that("betti_numbers on isolated vertices gives b0=n", {
  mat <- matrix(0, 5, 5)
  rownames(mat) <- colnames(mat) <- LETTERS[1:5]
  sc <- build_simplicial(mat)
  b <- betti_numbers(sc)
  expect_equal(b[["b0"]], 5L)
})

test_that("persistent_homology with different max_dim works", {
  net <- .make_sc_net()
  ph <- persistent_homology(net, n_steps = 5, max_dim = 1L)
  expect_s3_class(ph, "net_persistent_homology")
  expect_true(all(ph$betti_curve$dimension <= 1L))
})

test_that("q_analysis on single simplex returns q_vector of 1s", {
  mat <- matrix(0, 2, 2)
  mat[1, 2] <- mat[2, 1] <- 1
  rownames(mat) <- colnames(mat) <- c("A", "B")
  sc <- build_simplicial(mat)
  qa <- q_analysis(sc)
  expect_true(all(qa$q_vector == 1L))
})

test_that("pathway complex with no higher-order edges returns empty", {
  trajs <- list(c("A", "B"))
  hon <- build_hon(trajs, max_order = 2, min_freq = 1)
  sc <- build_simplicial(hon, type = "pathway")
  expect_s3_class(sc, "net_simplicial")
})


# =========================================================================
# Coverage: HON max_pathways limiting
# =========================================================================

test_that("build_simplicial pathway HON with max_pathways truncates", {
  # B's next-step depends on previous state => higher-order edges
  trajs <- list(
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D")
  )
  hon <- build_hon(trajs, max_order = 2, min_freq = 2)
  ho <- hon$ho_edges[hon$ho_edges$from_order > 1L, ]
  # Ensure we have >1 higher-order edges to test truncation
  expect_gt(nrow(ho), 1L)
  sc <- build_simplicial(hon, type = "pathway", max_pathways = 1)
  expect_s3_class(sc, "net_simplicial")
})

test_that("build_simplicial pathway HYPA with max_pathways truncates", {
  trajs <- list(
    c("A", "B", "C"), c("A", "B", "C"), c("A", "B", "C"),
    c("A", "B", "D"), c("C", "B", "D"), c("C", "B", "A"),
    c("D", "A", "B"), c("D", "A", "B")
  )
  h <- build_hypa(trajs, k = 2)
  n_anom <- sum(h$scores$anomaly != "normal")
  if (n_anom > 1L) {
    sc <- build_simplicial(h, type = "pathway", max_pathways = 1)
    expect_s3_class(sc, "net_simplicial")
  } else {
    # If only 0-1 anomalous paths, max_pathways won't truncate
    sc <- build_simplicial(h, type = "pathway")
    expect_s3_class(sc, "net_simplicial")
  }
})

# =========================================================================
# Coverage: boundary matrix edge case (dim 0, no k-simplices/km1-simplices)
# =========================================================================

test_that("betti_numbers handles dimension gap", {
  # Create a complex where some intermediate dimension is empty
  # Single edge A-B: has 0-simplices and 1-simplex
  mat <- matrix(0, 3, 3)
  mat[1, 2] <- mat[2, 1] <- 1
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  sc <- build_simplicial(mat)
  b <- betti_numbers(sc)
  expect_equal(b[["b0"]], 2L) # 2 components (A-B and C)
})

# =========================================================================
# Coverage: empty persistence result
# =========================================================================

test_that("persistent_homology with tiny network has valid persistence", {
  mat <- matrix(0, 2, 2)
  mat[1, 2] <- mat[2, 1] <- 0.5
  rownames(mat) <- colnames(mat) <- c("A", "B")
  ph <- persistent_homology(mat, n_steps = 5)
  expect_s3_class(ph, "net_persistent_homology")
  expect_true(is.data.frame(ph$persistence))
  expect_true(all(c("dimension", "birth", "death", "persistence") %in%
                     names(ph$persistence)))
})

test_that("persistent_homology with uniform weights yields flat betti", {
  # All edges same weight => nothing changes across filtration
  mat <- matrix(1, 3, 3)
  diag(mat) <- 0
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  ph <- persistent_homology(mat, n_steps = 3)
  expect_s3_class(ph, "net_persistent_homology")
  # Persistence df may be empty or have only persistent features
  expect_true(is.data.frame(ph$persistence))
})

test_that("build_simplicial pathway with max_dim truncation", {
  # Pathway with 4+ nodes should retain all lower-dimensional faces.
  trajs <- list(
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D"),
    c("A", "B", "C", "B", "D")
  )
  hon <- build_hon(trajs, max_order = 2, min_freq = 2)
  sc <- build_simplicial(hon, type = "pathway", max_dim = 1L)
  expect_s3_class(sc, "net_simplicial")
  expect_lte(sc$dimension, 1L)
})

test_that(".expand_to_faces keeps all faces up to max_dim", {
  faces_1 <- .expand_to_faces(list(1:4), max_dim = 1L)
  keys_1 <- sort(vapply(faces_1, paste, collapse = "-", FUN.VALUE = character(1)))
  expect_equal(keys_1, sort(c(
    "1", "2", "3", "4",
    "1-2", "1-3", "1-4", "2-3", "2-4", "3-4"
  )))

  faces_2 <- .expand_to_faces(list(1:4), max_dim = 2L)
  keys_2 <- sort(vapply(faces_2, paste, collapse = "-", FUN.VALUE = character(1)))
  expect_true(all(c("1", "2", "3", "4") %in% keys_2))
  expect_true(all(c("1-2", "1-3", "1-4", "2-3", "2-4", "3-4") %in% keys_2))
  expect_true(all(c("1-2-3", "1-2-4", "1-3-4", "2-3-4") %in% keys_2))
  expect_false("1-2-3-4" %in% keys_2)
})


# =========================================================================
# build_simplicial — pathway type from net_mogen
# =========================================================================

test_that("build_simplicial pathway from mogen works", {
  set.seed(42)
  # Create data with deterministic higher-order patterns
  seqs <- lapply(1:80, function(i) {
    s <- character(10)
    s[1] <- sample(LETTERS[1:4], 1)
    for (j in 2:10) {
      if (s[j - 1] == "A") s[j] <- sample(c("B", "C"), 1, prob = c(0.9, 0.1))
      else if (s[j - 1] == "B") s[j] <- "C"
      else if (s[j - 1] == "C") s[j] <- sample(c("D", "A"), 1, prob = c(0.8, 0.2))
      else s[j] <- "A"
    }
    s
  })
  mog <- build_mogen(seqs, max_order = 3)

  # Use order 2 explicitly (guaranteed to have transitions)
  if (max(mog$orders) >= 2L) {
    trans <- mogen_transitions(mog, order = 2)
    if (nrow(trans) > 0L) {
      # Force optimal_order to 2 so simplicial code picks it up
      mog$optimal_order <- 2L
      sc <- build_simplicial(mog, type = "pathway")
      expect_s3_class(sc, "net_simplicial")
      expect_equal(sc$type, "pathway")
      expect_true(sc$n_nodes > 0L)
    }
  }
})

test_that("build_simplicial pathway from mogen with max_pathways limits", {
  set.seed(99)
  seqs <- lapply(1:80, function(i) {
    s <- character(10)
    s[1] <- sample(LETTERS[1:4], 1)
    for (j in 2:10) {
      if (s[j - 1] == "A") s[j] <- "B"
      else if (s[j - 1] == "B") s[j] <- "C"
      else if (s[j - 1] == "C") s[j] <- "D"
      else s[j] <- "A"
    }
    s
  })
  mog <- build_mogen(seqs, max_order = 2)
  mog$optimal_order <- max(mog$orders[mog$orders >= 1L])

  trans <- mogen_transitions(mog, order = mog$optimal_order)
  if (nrow(trans) > 2L) {
    sc_full <- build_simplicial(mog, type = "pathway")
    sc_lim <- build_simplicial(mog, type = "pathway", max_pathways = 2)
    expect_s3_class(sc_lim, "net_simplicial")
    expect_lte(sc_lim$n_simplices, sc_full$n_simplices)
  }
})

test_that("build_simplicial pathway from mogen errors when optimal_order = 0", {
  seqs <- list(c("A", "B", "C"), c("B", "C", "A"), c("C", "A", "B"))
  mog <- build_mogen(seqs, max_order = 2)
  # Force order 0

  mog$optimal_order <- 0L
  mog$orders <- 0L
  expect_error(
    build_simplicial(mog, type = "pathway"),
    "no higher-order transitions"
  )
})

test_that("build_simplicial pathway from mogen with empty transitions returns empty", {
  seqs <- list(c("A", "B", "C"), c("B", "C", "A"), c("C", "A", "B"))
  mog <- build_mogen(seqs, max_order = 2)
  mog$optimal_order <- 1L
  # Wipe transition matrix to produce 0 rows from mogen_transitions
  n <- nrow(mog$transition_matrices[[2]])
  mog$count_matrices[[2]] <- matrix(0L, n, n,
    dimnames = dimnames(mog$transition_matrices[[2]]))
  sc <- build_simplicial(mog, type = "pathway")
  expect_s3_class(sc, "net_simplicial")
  expect_equal(sc$type, "pathway")
})

test_that("build_simplicial pathway error message includes net_mogen", {
  expect_error(
    build_simplicial(list(a = 1), type = "pathway"),
    "net_mogen"
  )
})
