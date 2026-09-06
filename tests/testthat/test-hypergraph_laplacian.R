# ---- Tests for hypergraph_laplacian / hypergraph_cluster /
# ---- hypergraph_transduction

# Two-community hypergraph with a bridge: {a,b,c} x2, {d,e,f} x2, {c,d}
.hl_planted <- function() {
  events <- data.frame(
    person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
               "d", "e", "f", "c", "d"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
                "m4", "m4", "m4", "m5", "m5"),
    stringsAsFactors = FALSE
  )
  group_hypergraph(events, actor = "person", group = "meeting")
}

# Small weighted hypergraph (EDVW cells) - connected
.hl_weighted <- function() {
  events <- data.frame(
    person = c("a", "b", "c", "a", "b", "d", "c", "d", "e", "e", "a"),
    meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
                "m4", "m4"),
    hours = c(2, 1, 1, 3, 2, 1, 2, 2, 4, 1, 1),
    stringsAsFactors = FALSE
  )
  group_hypergraph(events, actor = "person", group = "meeting",
                   weight = "hours")
}

# ---- formula: Zhou Theta against independent brute-force ---------------

test_that("zhou Laplacian matches the paper formula computed brute-force", {
  hg <- .hl_weighted()
  H <- (hg$incidence > 0) * 1.0
  n <- nrow(H); m <- ncol(H)
  set.seed(7)
  w <- runif(m, 0.5, 2)
  d <- as.vector(H %*% w)
  delta <- colSums(H)
  Theta_bf <- outer(seq_len(n), seq_len(n), Vectorize(function(u, v) {
    sum(w * H[u, ] * H[v, ] / delta) / sqrt(d[u] * d[v])
  }))
  L <- hypergraph_laplacian(hg, type = "zhou", edge_weights = w)
  expect_equal(unclass(L), diag(n) - Theta_bf,
               ignore_attr = TRUE, tolerance = 1e-12)
})

test_that("random_walk Laplacian matches brute-force EDVW construction", {
  hg <- .hl_weighted()
  gamma <- hg$incidence * 1.0
  H <- (gamma > 0) * 1.0
  n <- nrow(gamma); m <- ncol(gamma)
  we <- apply(gamma, 2, function(col) {
    x <- col[col > 0]
    sqrt(mean((x - mean(x))^2)) + 1
  })
  delta_g <- colSums(gamma)
  d_v <- as.vector(H %*% we)
  P_bf <- outer(seq_len(n), seq_len(n), Vectorize(function(v, u) {
    sum(H[v, ] * we / d_v[v] * gamma[u, ] / delta_g)
  }))
  expect_equal(rowSums(P_bf), rep(1, n), tolerance = 1e-12)
  pi_bf <- Re(eigen(t(P_bf))$vectors[, 1])
  pi_bf <- pi_bf / sum(pi_bf)
  G <- diag(sqrt(pi_bf)) %*% P_bf %*% diag(1 / sqrt(pi_bf))
  L_bf <- diag(n) - (G + t(G)) / 2
  L <- hypergraph_laplacian(hg, type = "random_walk")
  expect_equal(unclass(L), L_bf, ignore_attr = TRUE, tolerance = 1e-10)
  expect_equal(attr(L, "edge_weights"), we, ignore_attr = TRUE)
})

# ---- invariants --------------------------------------------------------

test_that("both Laplacians are symmetric PSD with spectrum in [0, 2]", {
  for (hg in list(.hl_planted(), .hl_weighted())) {
    for (type in c("zhou", "random_walk")) {
      L <- hypergraph_laplacian(hg, type = type)
      expect_true(isSymmetric(unclass(L)))
      ev <- eigen(unclass(L), symmetric = TRUE, only.values = TRUE)$values
      expect_gte(min(ev), -1e-10)
      expect_lte(max(ev), 2 + 1e-10)
      expect_lt(abs(min(ev)), 1e-10)  # connected: smallest eigenvalue 0
      pi_v <- attr(L, "pi")
      expect_equal(sum(pi_v), 1, tolerance = 1e-10)
      expect_true(all(pi_v > 0))
    }
  }
})

test_that("binary incidence + unit weights: zhou == random_walk", {
  hg <- .hl_planted()
  Lz <- hypergraph_laplacian(hg, type = "zhou")
  Lrw <- hypergraph_laplacian(hg, type = "random_walk")
  expect_equal(unclass(Lz), unclass(Lrw), ignore_attr = TRUE,
               tolerance = 1e-10)
})

test_that("Laplacian is permutation-equivariant in the node order", {
  hg <- .hl_weighted()
  set.seed(3)
  perm <- sample(hg$n_nodes)
  hg2 <- hg
  hg2$incidence <- hg$incidence[perm, , drop = FALSE]
  hg2$nodes <- hg$nodes[perm]
  laplacian_original <- hypergraph_laplacian(hg, type = "random_walk")
  L1 <- unclass(laplacian_original)
  laplacian_permuted <- hypergraph_laplacian(hg2, type = "random_walk")
  L2 <- unclass(laplacian_permuted)
  expect_equal(L2, L1[perm, perm], ignore_attr = TRUE, tolerance = 1e-12)
})

# ---- clustering --------------------------------------------------------

test_that("hypergraph_cluster recovers the planted two communities", {
  hg <- .hl_planted()
  cl <- hypergraph_cluster(hg, k = 2, seed = 1)
  expect_s3_class(cl, "net_hypergraph_cluster")
  a <- as.data.frame(cl)
  expect_equal(names(a), c("node", "cluster", "pi", "dim1", "dim2"))
  expect_equal(nrow(a), 6L)
  expect_equal(sum(a$pi), 1, tolerance = 1e-10)
  expect_equal(dim(cl$embedding), c(6L, 2L))
  # {a,b,c} together, {d,e,f} together
  expect_length(unique(a$cluster[a$node %in% c("a", "b", "c")]), 1L)
  expect_length(unique(a$cluster[a$node %in% c("d", "e", "f")]), 1L)
  expect_false(a$cluster[a$node == "a"] == a$cluster[a$node == "d"])
  # first-appearance labelling: node "a" is in "Cluster 1"
  expect_identical(a$cluster[a$node == "a"], "Cluster 1")
})

test_that("clustering is stable across seeds on separated structure", {
  hg <- .hl_planted()
  parts <- lapply(c(1L, 42L, 99L), function(s) {
    as.data.frame(hypergraph_cluster(hg, k = 2, seed = s))$cluster
  })
  expect_identical(parts[[1]], parts[[2]])
  expect_identical(parts[[1]], parts[[3]])
})

test_that("RDC-SymNMF follows Algorithm 2 and decreases Eq. 16", {
  hg <- .hl_planted()
  cl <- hypergraph_cluster(hg, k = 2, type = "random_walk",
                           algorithm = "symnmf", nstart = 3, seed = 7,
                           max_iter = 1000, tol = 1e-8)
  expect_s3_class(cl, "net_hypergraph_cluster")
  expect_identical(cl$algorithm, "symnmf")
  expect_true(all(cl$embedding >= 0))
  expect_equal(dim(cl$embedding), c(6L, 2L))
  expect_true(all(diff(cl$params$objective_history) <= 1e-8))
  rw_laplacian <- hypergraph_laplacian(hg, "random_walk")
  expect_equal(cl$params$objective,
               sum((diag(hg$n_nodes) -
                      unclass(rw_laplacian) -
                      tcrossprod(cl$embedding))^2),
               tolerance = 1e-8)
  assignment <- max.col(cl$embedding, ties.method = "first")
  expect_identical(cl$clusters$cluster,
                   paste("Cluster", match(assignment, unique(assignment))))
  expect_length(unique(cl$clusters$cluster[1:3]), 1L)
  expect_length(unique(cl$clusters$cluster[4:6]), 1L)
  expect_false(cl$clusters$cluster[1] == cl$clusters$cluster[4])
})

test_that("RDC-SymNMF is deterministic for a fixed seed", {
  args <- list(hg = .hl_planted(), k = 2, algorithm = "symnmf",
               nstart = 2, seed = 11, max_iter = 100)
  a <- do.call(hypergraph_cluster, args)
  b <- do.call(hypergraph_cluster, args)
  expect_identical(a$clusters, b$clusters)
  expect_equal(a$embedding, b$embedding)
  expect_equal(a$params$objective_history, b$params$objective_history)
})

test_that("J-NMF and JS-NMF implement Eqs. 18 and 19", {
  hg <- .hl_planted()
  S <- matrix(0, hg$n_nodes, hg$n_nodes,
              dimnames = list(hg$nodes, hg$nodes))
  S[1:3, 1:3] <- 1
  S[4:6, 4:6] <- 1
  diag(S) <- 0
  for (method in c("joint", "joint_symmetric")) {
    fit <- hypergraph_joint_cluster(
      hg, S, k = 2, method = method, nstart = 2, seed = 3,
      max_iter = 1000, tol = 1e-8
    )
    expect_s3_class(fit, "net_hypergraph_cluster")
    expect_identical(fit$algorithm, method)
    expect_true(all(fit$embedding >= 0))
    expect_true(tail(fit$params$objective_history, 1) <=
                  fit$params$objective_history[1])
    M <- fit$params$M
    Mt <- fit$params$Mtilde
    if (method == "joint") {
      X <- t(hg$incidence * 1.0)
      expected_objective <- sum((X - fit$params$Z %*% t(M))^2) +
        sum((S - M %*% t(Mt))^2) + sum((M - Mt)^2)
      expect_error(plot(fit, what = "spectrum"), "no Laplacian spectrum")
    } else {
      rw_laplacian <- hypergraph_laplacian(hg, type = "random_walk")
      C <- diag(hg$n_nodes) - unclass(rw_laplacian)
      expected_objective <- sum((C - M %*% t(fit$params$Mhat))^2) +
        sum((M - fit$params$Mhat)^2) +
        sum((S - M %*% t(Mt))^2) + sum((M - Mt)^2)
    }
    expect_equal(fit$params$objective, expected_objective,
                 tolerance = 1e-8)
    expect_length(unique(fit$clusters$cluster[1:3]), 1L)
    expect_length(unique(fit$clusters$cluster[4:6]), 1L)
    expect_false(fit$clusters$cluster[1] == fit$clusters$cluster[4])
  }
})

test_that("J-NMF Eq. 18 reduces exactly to its ordinary-NMF objective", {
  hg <- .hl_planted()
  S <- matrix(0, hg$n_nodes, hg$n_nodes)
  joint <- .hl_joint_nmf_fit(hg, S, k = 2, method = "joint",
                             alpha = 0, beta = 0, gamma = 0,
                             type = "random_walk", edge_weights = NULL,
                             seed = 17, max_iter = 25, tol = 1e-20)
  X <- t(hg$incidence * 1.0)
  expect_equal(joint$objective, sum((X - joint$Z %*% t(joint$M))^2),
               tolerance = 1e-12)
  # Mtilde is absent from the reduced objective, as stated after Eq. 18.
  shifted <- joint$Mtilde + 100
  expect_equal(joint$objective,
               sum((X - joint$Z %*% t(joint$M))^2) +
                 0 * sum((S - joint$M %*% t(shifted))^2) +
                 0 * sum((joint$M - shifted)^2),
               tolerance = 1e-12)
})

test_that("joint clustering aligns named relations and validates inputs", {
  hg <- .hl_planted()
  S <- diag(hg$n_nodes)
  dimnames(S) <- list(rev(hg$nodes), rev(hg$nodes))
  aligned_fit <- hypergraph_joint_cluster(
    hg, S, 2, nstart = 1, seed = 1, max_iter = 2
  )
  expect_s3_class(aligned_fit, "net_hypergraph_cluster")
  expect_error(hypergraph_joint_cluster(hg, matrix(-1, 6, 6), 2),
               "non-negative")
  expect_error(hypergraph_joint_cluster(hg, diag(5), 2), "n_nodes")
})

test_that("summary.net_hypergraph_cluster returns tidy shares", {
  cl <- hypergraph_cluster(.hl_planted(), k = 2, seed = 1)
  s <- summary(cl)
  expect_true(is.data.frame(s))
  expect_equal(names(s), c("cluster", "size", "share"))
  expect_equal(sum(s$share), 1, tolerance = 1e-12)
  expect_invisible(print(cl))
})

test_that("plot.net_hypergraph_cluster panels are ggplots", {
  cl <- hypergraph_cluster(.hl_planted(), k = 2, seed = 1)
  spectrum_plot <- plot(cl, what = "spectrum")
  expect_s3_class(spectrum_plot, "ggplot")
  embedding_plot <- plot(cl, what = "embedding")
  expect_s3_class(embedding_plot, "ggplot")
})

# ---- transduction ------------------------------------------------------

test_that("transduction closed form equals the Neumann series", {
  hg <- .hl_planted()
  zhou_laplacian <- hypergraph_laplacian(hg, type = "zhou")
  L <- unclass(zhou_laplacian)
  n <- nrow(L)
  S <- diag(n) - L
  xi <- 0.9
  y <- as.numeric(seq_len(n) == 1L)
  f_closed <- unname((1 - xi) * solve(diag(n) - xi * S, y))
  f_series <- numeric(n)
  term <- y
  # geometric-series reference; sequential by construction
  for (t in 0:3000) {
    f_series <- f_series + (1 - xi) * xi^t * term
    term <- as.vector(S %*% term)
  }
  expect_equal(f_closed, f_series, tolerance = 1e-10)
})

test_that("balanced seeds classify the planted communities correctly", {
  hg <- .hl_planted()
  tr <- hypergraph_transduction(
    hg, labels = c(a = "x", b = "x", e = "y", f = "y"), xi = 0.9
  )
  p <- as.data.frame(tr)
  expect_identical(p$predicted[p$node %in% c("a", "b", "c")],
                   rep("x", 3L))
  expect_identical(p$predicted[p$node %in% c("e", "f")], rep("y", 2L))
  expect_true(all(p$margin >= 0))
  expect_invisible(print(tr))
})

test_that("transduction score accessor is tidy long format", {
  hg <- .hl_planted()
  tr <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"))
  sc <- as.data.frame(tr, what = "scores")
  expect_equal(names(sc), c("node", "class", "score"))
  expect_equal(nrow(sc), hg$n_nodes * 2L)
  s <- summary(tr)
  expect_equal(names(s),
               c("class", "n_labeled", "n_predicted", "mean_margin"))
  p <- plot(tr)
  expect_s3_class(p, "ggplot")
})

test_that("transduction works with the random_walk Laplacian and weights", {
  hg <- .hl_weighted()
  tr <- hypergraph_transduction(hg, labels = c(a = "x", e = "y"),
                                type = "random_walk")
  expect_s3_class(tr, "net_hypergraph_transduction")
  predictions <- as.data.frame(tr)
  expect_equal(nrow(predictions), hg$n_nodes)
})

# ---- error paths -------------------------------------------------------

test_that("disconnected hypergraphs raise a classed condition", {
  events <- data.frame(person = c("a", "b", "c", "d"),
                       meeting = c("m1", "m1", "m2", "m2"),
                       stringsAsFactors = FALSE)
  hg <- group_hypergraph(events, actor = "person", group = "meeting")
  expect_error(hypergraph_laplacian(hg),
               class = "honets_hypergraph_disconnected")
  expect_error(hypergraph_cluster(hg, k = 2),
               class = "honets_hypergraph_disconnected")
})

test_that("argument contracts are enforced", {
  hg <- .hl_planted()
  expect_error(hypergraph_cluster(hg, k = 1), "between 2 and")
  expect_error(hypergraph_cluster(hg, k = 6), "between 2 and")
  expect_error(hypergraph_laplacian(hg, edge_weights = c(1, 2)),
               "one per hyperedge")
  expect_error(hypergraph_laplacian(hg, edge_weights = rep(-1, 5)),
               "positive")
  expect_error(hypergraph_transduction(hg, labels = c(a = "x")),
               "two distinct classes")
  expect_error(hypergraph_transduction(hg, labels = c(zz = "x", a = "y")),
               "Unknown node names")
  expect_error(hypergraph_transduction(hg, labels = c(a = "x", d = "y"),
                                       xi = 1.5),
               "xi")
})

# ---- class-mass normalization ------------------------------------------

test_that("class_mass decision rule matches a hand-computed fixture", {
  # class A has 10x the spread mass; raw argmax picks A everywhere, the
  # normalized rule recovers row 2 for B: decision = F / colSums(F)
  f <- matrix(c(0.60, 0.50, 0.40,
                0.05, 0.09, 0.02), nrow = 3,
              dimnames = list(c("n1", "n2", "n3"), c("A", "B")))
  lab <- c("A", NA, NA)
  raw <- .hl_score_predictions(f, lab, "none")
  cmn <- .hl_score_predictions(f, lab, "class_mass")
  expect_identical(raw$predicted, c("A", "A", "A"))
  # hand: colSums = (1.5, 0.16); F/mass = A: .4 .3333 .2667 / B: .3125 .5625 .125
  expect_identical(cmn$predicted, c("A", "B", "A"))
  expect_equal(cmn$score, c(0.4, 0.5625, 4 / 15), tolerance = 1e-12)
  expect_equal(cmn$margin, c(0.4 - 0.3125, 0.5625 - 1 / 3, 4 / 15 - 0.125),
               tolerance = 1e-12)
})

test_that("class_mass predictions are invariant to per-class score scaling", {
  set.seed(7)
  f <- matrix(runif(20, 0.1, 1), nrow = 10,
              dimnames = list(paste0("n", seq_len(10)), c("A", "B")))
  scaled <- f
  scaled[, "A"] <- scaled[, "A"] * 7
  lab <- rep(NA_character_, 10)
  expect_identical(.hl_score_predictions(f, lab, "class_mass")$predicted,
                   .hl_score_predictions(scaled, lab, "class_mass")$predicted)
})

test_that("transduction engines accept normalization end to end", {
  hg <- .hl_planted()
  raw <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"))
  cmn <- hypergraph_transduction(hg, labels = c(a = "x", d = "y"),
                                 normalization = "class_mass")
  # class_mass recovers the planted partition {a,b,c} / {d,e,f}; the raw
  # rule does not (y's larger spread mass pulls b and c across the bridge)
  expect_identical(cmn$predictions$predicted,
                   c("x", "x", "x", "y", "y", "y"))
  expect_identical(raw$predictions$predicted,
                   c("x", "y", "y", "y", "y", "y"))
  expect_identical(cmn$normalization, "class_mass")
  expect_identical(raw$normalization, "none")
  # scores stay the raw spread scores under either rule
  expect_identical(cmn$scores, raw$scores)
  expect_error(hypergraph_transduction(hg, labels = c(a = "x", d = "y"),
                                       normalization = "bogus"))
})

test_that("zero-mass classes are refused with a classed condition", {
  f <- matrix(c(0.5, 0.4, 0, 0), nrow = 2,
              dimnames = list(c("n1", "n2"), c("A", "B")))
  expect_error(.hl_score_predictions(f, c("A", NA), "class_mass"),
               class = "honets_bad_input")
})
