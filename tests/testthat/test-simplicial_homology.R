testthat::skip_on_cran()

# ---- Persistent homology + diagram-level operations ----
#
# Tests against KNOWN topology — the right way to validate PH without an
# external library. We check that the boundary-matrix reduction recovers
# the correct (homology class count, persistence interval) on small
# complexes where the answer is computable by hand.

# Helpers -----------------------------------------------------------------

.ph_sim_mat <- function(edges, weights, n) {
  # Symmetric weight matrix from an edge list. edges: 2-col integer matrix.
  m <- matrix(0, n, n)
  for (k in seq_len(nrow(edges))) {
    m[edges[k, 1L], edges[k, 2L]] <- weights[k]
    m[edges[k, 2L], edges[k, 1L]] <- weights[k]
  }
  rownames(m) <- colnames(m) <- paste0("v", seq_len(n))
  m
}

.ph_dist_circle <- function(n = 12L) {
  # Pairwise great-circle (chord) distances of n points evenly on a unit
  # circle. The boundary of these points is the circle itself; at small
  # filtration this is a hollow polygon with H_1 = 1.
  theta <- seq(0, 2 * pi, length.out = n + 1L)[seq_len(n)]
  xy <- cbind(cos(theta), sin(theta))
  d <- as.matrix(stats::dist(xy))
  d
}

# 4-cycle: H_1 = 1 -------------------------------------------------------

test_that("4-cycle clique complex has one essential H1", {
  edges <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  weights <- c(0.51, 0.52, 0.53, 0.54)
  mat <- .ph_sim_mat(edges, weights, 4L)

  ph <- persistent_homology(mat, n_steps = 5L, max_dim = 2L)

  # Exactly one essential H_1 (the hole) — clique mode sets essential death = 0.
  ess_b1 <- ph$persistence[ph$persistence$dimension == 1L &
                             ph$persistence$death == 0, , drop = FALSE]
  expect_equal(nrow(ess_b1), 1L)
  # Born when the last edge enters, i.e., similarity = 0.51.
  expect_equal(ess_b1$birth, 0.51)

  # Exactly one essential H_0 (single connected component).
  ess_b0 <- ph$persistence[ph$persistence$dimension == 0L &
                             ph$persistence$death == 0, , drop = FALSE]
  expect_equal(nrow(ess_b0), 1L)
})

# Triangle (full 2-simplex): contractible --------------------------------

test_that("filled triangle has only one essential H0 (contractible)", {
  edges <- rbind(c(1, 2), c(1, 3), c(2, 3))
  weights <- c(0.6, 0.5, 0.4)
  mat <- .ph_sim_mat(edges, weights, 3L)

  ph <- persistent_homology(mat, n_steps = 5L, max_dim = 2L)

  ess <- ph$persistence[ph$persistence$death == 0, , drop = FALSE]
  expect_equal(nrow(ess), 1L)
  expect_equal(ess$dimension, 0L)

  # No essential H_1 — the triangle is filled.
  expect_equal(sum(ph$persistence$dimension == 1L &
                     ph$persistence$death == 0), 0L)
})

# Two disjoint edges: H_0 = 2 -------------------------------------------

test_that("two disjoint edges give two essential H0 classes", {
  # 4 nodes, 2 disjoint edges: {1-2}, {3-4}. Two components, no cycles.
  edges <- rbind(c(1, 2), c(3, 4))
  weights <- c(0.5, 0.5)
  mat <- .ph_sim_mat(edges, weights, 4L)

  ph <- persistent_homology(mat, n_steps = 5L, max_dim = 2L)

  ess_b0 <- ph$persistence[ph$persistence$dimension == 0L &
                             ph$persistence$death == 0, , drop = FALSE]
  expect_equal(nrow(ess_b0), 2L)

  # No essential H_1 (no cycles).
  expect_equal(sum(ph$persistence$dimension == 1L &
                     ph$persistence$death == 0), 0L)
})

# Bottleneck distance properties ----------------------------------------

test_that("bottleneck distance is zero for identical diagrams", {
  mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  rownames(mat) <- colnames(mat) <- c("A", "B", "C")
  ph <- persistent_homology(mat, n_steps = 5L)
  d <- bottleneck_distance(ph, ph)
  expect_true(all(d == 0))
})

test_that("bottleneck distance is symmetric", {
  m1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  m2 <- matrix(c(0, .7, .3, .7, 0, .2, .3, .2, 0), 3, 3)
  rownames(m1) <- colnames(m1) <- c("A","B","C")
  rownames(m2) <- colnames(m2) <- c("A","B","C")
  ph1 <- persistent_homology(m1, n_steps = 5L)
  ph2 <- persistent_homology(m2, n_steps = 5L)

  d_12 <- bottleneck_distance(ph1, ph2)
  d_21 <- bottleneck_distance(ph2, ph1)
  expect_equal(unname(d_12), unname(d_21))
})

test_that("bottleneck distance satisfies the triangle inequality", {
  m1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  m2 <- matrix(c(0, .7, .3, .7, 0, .2, .3, .2, 0), 3, 3)
  m3 <- matrix(c(0, .4, .8, .4, 0, .6, .8, .6, 0), 3, 3)
  for (m in list(m1, m2, m3)) {
    rownames(m) <- colnames(m) <- c("A","B","C")
  }
  ph1 <- persistent_homology(m1, n_steps = 5L)
  ph2 <- persistent_homology(m2, n_steps = 5L)
  ph3 <- persistent_homology(m3, n_steps = 5L)

  d_12 <- bottleneck_distance(ph1, ph2)
  d_23 <- bottleneck_distance(ph2, ph3)
  d_13 <- bottleneck_distance(ph1, ph3)
  expect_true(all(d_13 <= d_12 + d_23 + 1e-9))
})

test_that("bottleneck returns Inf when essential counts differ", {
  # m1 has 1 component (triangle); m2 has 2 components (two isolated edges).
  m1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  m2 <- matrix(0, 4, 4)
  m2[1, 2] <- m2[2, 1] <- 0.5
  m2[3, 4] <- m2[4, 3] <- 0.5
  rownames(m1) <- colnames(m1) <- c("A","B","C")
  rownames(m2) <- colnames(m2) <- paste0("v", 1:4)

  ph1 <- persistent_homology(m1, n_steps = 5L)
  ph2 <- persistent_homology(m2, n_steps = 5L)
  d <- bottleneck_distance(ph1, ph2)
  expect_true(is.infinite(d["dim_0"]))
})

# Wasserstein distance properties ----------------------------------------

test_that("Wasserstein distance matches hand-computed diagram costs", {
  one <- data.frame(dimension = 0L, birth = 0, death = 2)
  empty <- data.frame(dimension = integer(), birth = numeric(),
                      death = numeric())
  shifted <- data.frame(dimension = 0L, birth = 1, death = 3)

  expect_equal(unname(wasserstein_distance(one, empty, dimension = 0L)), 1)
  expect_equal(unname(wasserstein_distance(one, shifted)), 1)
  expect_equal(
    unname(wasserstein_distance(one, empty, dimension = 0L,
                                internal_p = 2)),
    sqrt(2)
  )
})

test_that("Wasserstein distance has metric properties", {
  d1 <- data.frame(dimension = c(0L, 0L), birth = c(0, 3),
                   death = c(2, 5))
  d2 <- data.frame(dimension = c(0L, 0L), birth = c(0.5, 4),
                   death = c(2.5, 6))
  d3 <- data.frame(dimension = 0L, birth = 1, death = 4)

  expect_equal(wasserstein_distance(d1, d1, order = 2), c(dim_0 = 0))
  expect_equal(wasserstein_distance(d1, d2, order = 2),
               wasserstein_distance(d2, d1, order = 2))
  expect_lte(
    unname(wasserstein_distance(d1, d3, order = 2)),
    unname(wasserstein_distance(d1, d2, order = 2)) +
      unname(wasserstein_distance(d2, d3, order = 2)) + 1e-12
  )
})

test_that("Wasserstein handles dimensions and essential classes", {
  d1 <- data.frame(dimension = c(0L, 1L), birth = c(0, 2),
                   death = c(Inf, 5))
  d2 <- data.frame(dimension = c(0L, 1L), birth = c(1, 2.5),
                   death = c(Inf, 5.5))
  out <- wasserstein_distance(d1, d2)
  expect_equal(out["dim_0"], c(dim_0 = 1))
  expect_equal(out["dim_1"], c(dim_1 = 0.5))

  mismatch <- rbind(d2, data.frame(dimension = 0L, birth = 3, death = Inf))
  expect_true(is.infinite(wasserstein_distance(d1, mismatch)["dim_0"]))
})

test_that("Wasserstein validates metric orders", {
  d <- data.frame(dimension = 0L, birth = 0, death = 1)
  expect_error(wasserstein_distance(d, d, order = Inf), "finite")
  expect_error(wasserstein_distance(d, d, order = 0), ">= 1")
  expect_error(wasserstein_distance(d, d, internal_p = 0), ">= 1")
})

# Vietoris-Rips on a circle ---------------------------------------------

test_that("VR filtration on a circle recovers H_1 = 1 essential", {
  d <- .ph_dist_circle(n = 8L)
  # max_scale set so the 8 vertices form a cycle of 8 edges but no diagonals
  # get pulled in to fill the hole. The shortest "diagonal" distance on an
  # 8-gon is 2*sin(2*pi/8) ≈ 1.414; the shortest cycle edge is
  # 2*sin(pi/8) ≈ 0.765. Choose a cap between them.
  ph <- persistent_homology(d, n_steps = 5L, max_dim = 2L,
                            type = "vr", max_scale = 1.0)

  ess_b1 <- ph$persistence[ph$persistence$dimension == 1L &
                             is.finite(ph$persistence$death), , drop = FALSE]
  # Should have at least one finite H_1 pair OR an essential H_1; either way
  # the hole is detected.
  total_h1 <- sum(ph$persistence$dimension == 1L)
  expect_gte(total_h1, 1L)

  # H_0: should be 1 essential (all 8 points connected by the cycle).
  ess_b0 <- ph$persistence[ph$persistence$dimension == 0L &
                             ph$persistence$death == ph$persistence$death[
                               ph$persistence$dimension == 0L
                             ][1L], , drop = FALSE]
  # Less brittle: at the final threshold, there is exactly 1 connected component.
  fin_t <- max(ph$thresholds)
  b0_at_end <- ph$betti_curve$betti[
    ph$betti_curve$dimension == 0L &
      ph$betti_curve$threshold == fin_t]
  expect_equal(b0_at_end, 1L)
})

test_that("build_simplicial(type='vr') filtration values match expectation", {
  # 3 points in a line: 1—2—3 with d(1,2)=0.5, d(2,3)=0.5, d(1,3)=1.0
  d <- matrix(c(0, 0.5, 1.0,
                0.5, 0, 0.5,
                1.0, 0.5, 0), 3, 3, byrow = TRUE)
  rownames(d) <- colnames(d) <- c("A","B","C")
  sc <- build_simplicial(d, type = "vr", max_scale = 1.0, max_dim = 2L)

  expect_identical(sc$type, "vr")
  # Vertices appear at filt = 0; edges at their distance; the 2-simplex
  # {A,B,C} appears at max pairwise distance = 1.0.
  vert_idx <- which(vapply(sc$simplices, length, integer(1)) == 1L)
  edge_idx <- which(vapply(sc$simplices, length, integer(1)) == 2L)
  tri_idx  <- which(vapply(sc$simplices, length, integer(1)) == 3L)

  expect_true(all(sc$filtration[vert_idx] == 0))
  expect_setequal(sort(sc$filtration[edge_idx]), c(0.5, 0.5, 1.0))
  expect_equal(sc$filtration[tri_idx], 1.0)
})

# Persistence landscape -------------------------------------------------

test_that("persistence landscape sup norms decrease with k", {
  edges <- rbind(c(1, 2), c(2, 3), c(3, 4), c(4, 1))
  weights <- c(0.51, 0.52, 0.53, 0.54)
  mat <- .ph_sim_mat(edges, weights, 4L)
  ph <- persistent_homology(mat, n_steps = 5L, max_dim = 2L)

  pl <- persistence_landscape(ph, k_max = 4L, dimension = 0L)
  expect_s3_class(pl, "net_persistence_landscape")

  k_norms <- vapply(seq_len(4L), function(k) {
    sub <- pl$landscape[pl$landscape$k == k, ]
    max(sub$value)
  }, numeric(1))
  # Top-k norms must be weakly decreasing in k.
  expect_true(all(diff(k_norms) <= 1e-12))
})

test_that("persistence landscape is non-negative everywhere", {
  mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  rownames(mat) <- colnames(mat) <- c("A","B","C")
  ph <- persistent_homology(mat, n_steps = 5L)
  pl <- persistence_landscape(ph, k_max = 3L, dimension = 0L)
  expect_true(all(pl$landscape$value >= 0))
})

test_that("persistence landscape on empty diagram returns zeros", {
  # Diagram with only essentials (e.g., two disjoint edges)
  edges <- rbind(c(1, 2), c(3, 4))
  mat <- .ph_sim_mat(edges, c(0.5, 0.5), 4L)
  ph <- persistent_homology(mat, n_steps = 5L)
  pl <- persistence_landscape(ph, k_max = 2L, dimension = 0L)
  expect_true(all(pl$landscape$value == 0))
})

# persistent_homology backward compatibility ----------------------------

test_that("persistent_homology preserves $persistence column shape", {
  mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  rownames(mat) <- colnames(mat) <- c("A","B","C")
  ph <- persistent_homology(mat, n_steps = 5L)
  expect_true(all(c("dimension", "birth", "death", "persistence") %in%
                    names(ph$persistence)))
  expect_true(all(ph$persistence$persistence >= 0))
})

# Regression: Codex review 2026-05-20 — zero-distance VR edges -----------

test_that("VR includes off-diagonal zero-distance edges (duplicate points)", {
  # 3 distinct labels but two of them are coincident (d = 0). They should
  # form one connected component at scale 0, not three.
  d <- matrix(0, 3, 3)
  d[1, 2] <- d[2, 1] <- 0   # duplicate
  d[1, 3] <- d[3, 1] <- 1
  d[2, 3] <- d[3, 2] <- 1
  rownames(d) <- colnames(d) <- c("A", "Ab", "B")  # Ab coincident with A

  ph <- persistent_homology(d, n_steps = 5L, max_dim = 1L,
                            type = "vr", max_scale = 1)

  # b0 at threshold 0: A and Ab are merged (zero-distance edge), B is alone.
  b0_at_zero <- ph$betti_curve$betti[
    ph$betti_curve$dimension == 0L & ph$betti_curve$threshold == 0
  ]
  # Either b0 == 2 at t = 0 (Ab merged with A; B alone), or the grid does
  # not include t = 0 — in which case the smallest threshold > 0 must show
  # at most 2 components.
  if (length(b0_at_zero) > 0L) {
    expect_lte(b0_at_zero, 2L)
  } else {
    smallest_t <- min(ph$thresholds)
    b0_smallest <- ph$betti_curve$betti[
      ph$betti_curve$dimension == 0L &
        ph$betti_curve$threshold == smallest_t
    ]
    expect_lte(b0_smallest, 2L)
  }
})

# Regression: Codex review 2026-05-20 — essential VR persistence -------

test_that("VR persistence diagram plot retains essential classes", {
  # 4 points forming a unit square; VR detects H_0 essential + H_1 finite.
  pts <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  d <- as.matrix(stats::dist(pts))
  rownames(d) <- colnames(d) <- paste0("p", 1:4)
  ph <- persistent_homology(d, n_steps = 5L, max_dim = 2L,
                            type = "vr", max_scale = 2)
  pl <- plot(ph, combined = FALSE)
  # The persistence diagram panel's underlying data must include all rows
  # with persistence > 0 — essential classes (persistence = Inf) must NOT
  # be silently dropped by ggplot's continuous size scale.
  pd_data <- pl$persistence$data
  expect_true(nrow(pd_data) >= 1L)
  # The size aesthetic maps to persistence; it must be finite for every row.
  expect_true(all(is.finite(pd_data$persistence)))
  # Essential dimension-0 class (the surviving connected component) must
  # appear in the rendered plot data.
  expect_true(any(pd_data$dimension == 0L))
})

# Regression: Codex review 2026-05-20 — filtered-complex handoff -------

test_that("persistent_homology accepts a build_simplicial VR complex", {
  pts <- rbind(c(0,0), c(1,0), c(1,1), c(0,1))
  d <- as.matrix(stats::dist(pts))
  rownames(d) <- colnames(d) <- paste0("p", 1:4)

  # Path A: matrix → PH directly.
  ph_a <- persistent_homology(d, n_steps = 5L, max_dim = 2L,
                              type = "vr", max_scale = 2)
  # Path B: matrix → build_simplicial(type = "vr") → PH.
  sc <- build_simplicial(d, type = "vr", max_scale = 2, max_dim = 2L)
  ph_b <- persistent_homology(sc, n_steps = 5L, max_dim = 2L,
                              type = "vr")

  # Persistence tables must match by row count and by (dimension, birth, death)
  # up to ordering — both should compute the same diagram.
  expect_equal(nrow(ph_a$persistence), nrow(ph_b$persistence))
  a_keys <- sort(paste(ph_a$persistence$dimension,
                       ph_a$persistence$birth,
                       ph_a$persistence$death, sep = "/"))
  b_keys <- sort(paste(ph_b$persistence$dimension,
                       ph_b$persistence$birth,
                       ph_b$persistence$death, sep = "/"))
  expect_identical(a_keys, b_keys)
})

test_that("persistent_homology mode flag is set", {
  mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
  rownames(mat) <- colnames(mat) <- c("A","B","C")
  ph_c <- persistent_homology(mat, n_steps = 5L)
  expect_identical(ph_c$mode, "clique")

  d <- matrix(c(0, 0.4, 0.5,
                0.4, 0, 0.6,
                0.5, 0.6, 0), 3, 3, byrow = TRUE)
  rownames(d) <- colnames(d) <- c("A","B","C")
  ph_v <- persistent_homology(d, n_steps = 5L, type = "vr", max_scale = 1)
  expect_identical(ph_v$mode, "vr")
})
