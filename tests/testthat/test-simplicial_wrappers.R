testthat::skip_on_cran()

# ---- Simplicial wrappers: identity with the Nestimate estimators ----
#
# simplicial(), hg_homology(), hg_landscape(), hg_bottleneck(), hg_betti(),
# hg_euler(), hg_qanalysis() and hg_degree() return the estimator's result
# unchanged; identical() against the direct call on several inputs.

.sw_mat <- function(seed = 1L, n = 7L) {
  set.seed(seed)
  m <- matrix(stats::runif(n * n), n, n)
  m <- (m + t(m)) / 2
  diag(m) <- 0
  dimnames(m) <- list(letters[seq_len(n)], letters[seq_len(n)])
  m
}

test_that("simplicial() is the complex builder for every type", {
  m <- .sw_mat()
  expect_identical(.unclass_sc(simplicial(m, threshold = 0.5)),
                   Nestimate::build_simplicial(m, threshold = 0.5))
  expect_identical(.unclass_sc(simplicial(m, threshold = 0.3, max_dim = 2L)),
                   Nestimate::build_simplicial(m, threshold = 0.3,
                                               max_dim = 2L))
  d <- 1 - m
  diag(d) <- 0
  expect_identical(.unclass_sc(simplicial(d, type = "vr", max_scale = 0.6)),
                   Nestimate::build_simplicial(d, type = "vr",
                                               max_scale = 0.6))
  seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
               c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
  h <- hon(seqs, max_order = 2L)
  expect_identical(.unclass_sc(simplicial(h, type = "pathway")),
                   Nestimate::build_simplicial(h, type = "pathway"))
  net <- Nestimate::build_network(.ho_wide_sequences(seqs),
                                  method = "relative")
  expect_identical(.unclass_sc(simplicial(net, threshold = 0.2)),
                   Nestimate::build_simplicial(net, threshold = 0.2))
})

test_that("the measures of a complex are the estimator's", {
  sc <- simplicial(.sw_mat(), threshold = 0.4)
  expect_identical(hg_betti(sc), Nestimate::betti_numbers(sc))
  expect_identical(hg_euler(sc), Nestimate::euler_characteristic(sc))
  expect_identical(.ho_unresult(hg_qanalysis(sc)), Nestimate::q_analysis(sc))
  ref <- Nestimate::simplicial_degree(sc, normalized = TRUE)
  rownames(ref) <- NULL
  expect_identical(hg_degree(sc, normalized = TRUE), ref)
  expect_identical(hg_get(sc, what = "degree", normalized = TRUE), ref)
  expect_identical(hg_get(sc, what = "degree", top = 2L),
                   utils::head(hg_degree(sc), 2L))
  expect_error(hg_degree(list()), class = "hypernets_bad_input")
})

test_that("hg_homology, hg_landscape and hg_bottleneck are the estimator's", {
  m <- .sw_mat()
  ph <- hg_homology(m, n_steps = 8L, max_dim = 2L)
  expect_identical(.ho_unresult(ph),
                   Nestimate::persistent_homology(m, n_steps = 8L,
                                                  max_dim = 2L))
  d <- 1 - m
  diag(d) <- 0
  expect_identical(.ho_unresult(hg_homology(d, type = "vr", max_scale = 0.7)),
                   Nestimate::persistent_homology(d, type = "vr",
                                                  max_scale = 0.7))
  expect_identical(.ho_unresult(hg_landscape(ph, k_max = 3L, dimension = 0L)),
                   Nestimate::persistence_landscape(.ho_unresult(ph),
                                                    k_max = 3L,
                                                    dimension = 0L))
  ph2 <- hg_homology(.sw_mat(2L), n_steps = 8L, max_dim = 2L)
  expect_identical(hg_bottleneck(ph, ph2),
                   Nestimate::bottleneck_distance(ph, ph2))
  expect_identical(hg_bottleneck(ph, ph2, dimension = 0L),
                   Nestimate::bottleneck_distance(ph, ph2, dimension = 0L))
  # invariant: a diagram is at distance 0 from itself
  self <- hg_bottleneck(ph, ph)
  expect_equal(unname(self), rep(0, length(self)))
})

test_that("simplicial(type = \"window\") is the co-occurrence complex of windows", {
  seqs <- list(c("a", "b", "c", "a", "b", "c"), c("a", "b", "d", "a", "b", "d"),
               c("c", "d", "e", "c"))
  sc <- simplicial(seqs, type = "window", window = 3L)
  expect_s3_class(sc, "simplicial_complex")
  expect_identical(sc$type, "window")
  keys <- vapply(sc$simplices, \(s) paste(sort(sc$nodes[s]), collapse = ","),
                 character(1L))
  # closed under subsets: every face of every simplex is present
  faces <- unique(unlist(lapply(sc$simplices, \(s) {
    unlist(lapply(seq_along(s), \(k) combn(sort(sc$nodes[s]), k,
                                           \(x) paste(x, collapse = ","))))
  })))
  expect_setequal(keys, faces)
  # the window sets: {a,b,c}, {a,b,d}, {c,d,e} and their faces
  expect_true(all(c("a,b,c", "a,b,d", "c,d,e") %in% keys))
  expect_false("a,c,d" %in% keys)
  expect_identical(max(lengths(sc$simplices)), 3L)       # window - 1 = dim 2
  # min_count keeps only sets seen often enough; monotone
  strict <- simplicial(seqs, type = "window", window = 3L, min_count = 3L)
  expect_lte(strict$n_simplices, sc$n_simplices)
  expect_false("c,d,e" %in% vapply(strict$simplices, \(s)
    paste(sort(strict$nodes[s]), collapse = ","), character(1L)))
  # Euler-Poincare holds on the result
  expect_equal(hg_euler(sc), sum((-1)^(seq_along(hg_betti(sc)) - 1L) *
                                       hg_betti(sc)))
  # long input gives the same complex
  long <- data.frame(id = rep(c("s1", "s2", "s3"), c(6, 6, 4)),
                     t = c(1:6, 1:6, 1:4), code = unlist(seqs))
  expect_identical(simplicial(long, type = "window", window = 3L, actor = "id",
                              action = "code", time = "t"), sc)
  expect_error(simplicial(seqs, type = "window", window = 1L))
  expect_error(simplicial(seqs, type = "window", min_count = 0))
})

test_that("a window complex prints its construction and counts", {
  sc <- simplicial(list(c("a", "b", "c", "a", "b")), type = "window")
  expect_output(print(sc), "Window co-occurrence complex \\(window = 3, min_count = 1\\)")
  expect_output(print(sc), "Betti")
  expect_invisible(print(sc))
})

test_that("direction = 'both' joins a pair only when both weights reach threshold", {
  w <- matrix(c(0,   0.5, 0.5, 0,
                0.5, 0,   0.1, 0.5,
                0.5, 0.5, 0,   0.5,
                0,   0.5, 0.5, 0), 4, byrow = TRUE,
              dimnames = list(letters[1:4], letters[1:4]))
  keys <- function(sc) vapply(sc$simplices, \(s) paste(sort(sc$nodes[s]), collapse = ","),
                              character(1L))
  either <- simplicial(w, threshold = 0.3)
  both <- simplicial(w, threshold = 0.3, direction = "both")
  # b -> c is 0.1, c -> b is 0.5: joined by "either", not by "both"
  expect_true("b,c" %in% keys(either))
  expect_false("b,c" %in% keys(both))
  expect_true("a,b,c" %in% keys(either))
  expect_false("a,b,c" %in% keys(both))
  expect_true(all(keys(both) %in% keys(either)))               # "both" is a subcomplex
  expect_identical(simplicial(w, direction = "either"), simplicial(w))  # default unchanged
  expect_error(simplicial(w, direction = "sideways"))
  expect_error(simplicial(1:3, threshold = 0.3, direction = "both"),
               class = "hypernets_bad_input")
})

test_that("hg_homology() of a window complex follows the count filtration", {
  set.seed(4)
  seqs <- lapply(1:25, \(i) sample(c("a", "b", "c", "d", "e"), 12, replace = TRUE,
                                   prob = c(.3, .25, .2, .15, .1)))
  sc <- simplicial(seqs, type = "window", window = 3L)
  expect_length(sc$count, sc$n_simplices)
  # counts are monotone: a face is never rarer than a simplex containing it
  keys <- vapply(sc$simplices, \(s) paste(sort(s), collapse = ","), "")
  ok <- vapply(seq_along(sc$simplices), \(i) {
    s <- sc$simplices[[i]]
    if (length(s) < 2L) return(TRUE)
    faces <- lapply(seq_along(s), \(j) paste(sort(s[-j]), collapse = ","))
    all(sc$count[match(unlist(faces), keys)] >= sc$count[i])
  }, logical(1L))
  expect_true(all(ok))
  ph <- hg_homology(sc, max_dim = 2L)
  curve <- hg_get(ph, what = "betti")
  # at every count t the curve equals the Betti numbers of the complex built
  # with min_count = t (an independent code path)
  lapply(unique(curve$threshold), function(t) {
    direct <- as.integer(hg_betti(simplicial(seqs, type = "window", window = 3L,
                                              min_count = t)))
    direct <- c(direct, rep(0L, 3L - length(direct)))[1:3]
    at_t <- curve[curve$threshold == t, ]
    expect_identical(at_t$betti[order(at_t$dimension)], direct,
                     label = sprintf("Betti at count %s", t))
  })
  pers <- hg_get(ph)
  expect_true(all(pers$birth >= pers$death))
})

test_that("a clique complex of sequences is the complex of their transition network", {
  seqs <- ring_sequences[1:60]
  net <- Nestimate::build_network(.ho_wide_sequences(seqs), method = "relative")
  expect_identical(simplicial(seqs, threshold = 0.2),
                   simplicial(net, threshold = 0.2))
  expect_identical(simplicial(seqs, threshold = 0.2, direction = "both"),
                   simplicial(net, threshold = 0.2, direction = "both"))
  # long input with the shared sequence arguments gives the same complex
  long <- data.frame(id = rep(seq_along(seqs), lengths(seqs)),
                     step = unlist(lapply(lengths(seqs), seq_len)),
                     code = unlist(seqs))
  expect_identical(
    simplicial(long, threshold = 0.2, actor = "id", action = "code",
               time = "step"),
    simplicial(net, threshold = 0.2))
})

test_that("a Vietoris-Rips complex takes points as well as distances", {
  set.seed(4)
  points <- matrix(stats::rnorm(24), 8, 3,
                   dimnames = list(sprintf("p%d", 1:8), c("x", "y", "z")))
  distances <- as.matrix(stats::dist(points))
  expect_identical(simplicial(points, type = "vr", max_scale = 1.5),
                   simplicial(distances, type = "vr", max_scale = 1.5))
  expect_identical(simplicial(stats::dist(points), type = "vr",
                              max_scale = 1.5),
                   simplicial(distances, type = "vr", max_scale = 1.5))
  frame <- data.frame(id = rownames(points), points)
  expect_identical(simplicial(frame, type = "vr", max_scale = 1.5),
                   simplicial(distances, type = "vr", max_scale = 1.5))
  expect_identical(hg_homology(points, type = "vr", max_dim = 2L),
                   hg_homology(distances, type = "vr", max_dim = 2L))
  expect_error(simplicial(data.frame(a = letters[1:3], b = letters[3:1]),
                          type = "vr"), class = "hypernets_bad_input")
})

test_that("hg_get(what = \"betti\") tabulates the Betti numbers", {
  sc <- simplicial(.sw_mat(), threshold = 0.4)
  betti <- hg_get(sc, what = "betti")
  expect_named(betti, c("dimension", "betti"))
  expect_identical(betti$betti, as.integer(hg_betti(sc)))
})

test_that("plot() draws the maximal simplices in order", {
  sc <- simplicial(ring_sequences, type = "window", window = 3L,
                   min_count = 200L)
  facets <- .sc_facets(sc)
  # no drawn simplex is a face of another
  expect_false(any(vapply(seq_along(facets), \(i) {
    any(vapply(facets[-i], \(g) all(facets[[i]] %in% g), logical(1L)))
  }, logical(1L))))
  # the most frequent window set comes first
  counts <- sc$count[match(vapply(facets, paste, character(1L), collapse = ","),
                           vapply(sc$simplices, \(s) paste(sc$nodes[sort(s)],
                                                            collapse = ","),
                                  character(1L)))]
  expect_false(is.unsorted(rev(counts[lengths(facets) == max(lengths(facets))])))
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_s3_class(plot(sc), "ggplot")
  expect_no_error(plot(sc, dismantled = TRUE, top = 3))
  expect_error(plot(sc, top = 0))
})
