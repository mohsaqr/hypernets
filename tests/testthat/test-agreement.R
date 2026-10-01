# hg_agreement(), hg_stability(), hg_seeds(), and the eigenvalue gap /
# keyword collapse accessors.

toy_hg <- function() {
  text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
}

test_that("hg_agreement matches the hand-computed adjusted Rand index", {
  # a = (1,1,1,2,2), b = (1,1,2,2,2). Contingency: a1b1=2, a1b2=1, a2b2=2.
  # sum_ij = choose(2,2)+choose(1,2)+choose(2,2) = 2; sum_a = choose(3,2)+
  # choose(2,2) = 4; sum_b = 4; expected = 4*4/choose(5,2) = 1.6;
  # ARI = (2 - 1.6) / ((4+4)/2 - 1.6) = 0.4/2.4 = 1/6.
  x <- data.frame(node = letters[1:5],
                  cluster = c("g1", "g1", "g1", "g2", "g2"))
  y <- data.frame(node = letters[1:5],
                  cluster = c("g1", "g1", "g2", "g2", "g2"))
  out <- hg_agreement(x, y)
  expect_s3_class(out, "data.frame")
  expect_identical(nrow(out), 1L)
  expect_identical(out$n, 5L)
  expect_equal(out$agreement, 4 / 5)
  expect_equal(out$ari, 1 / 6)
})

test_that("hg_agreement ari is label-permutation invariant, agreement is not", {
  x <- data.frame(node = letters[1:6],
                  cluster = c("A", "A", "A", "B", "B", "B"))
  y <- data.frame(node = letters[1:6],
                  cluster = c("B", "B", "B", "A", "A", "A"))
  out <- hg_agreement(x, y)
  expect_equal(out$ari, 1)
  expect_equal(out$agreement, 0)
  self <- hg_agreement(x, x)
  expect_equal(self$agreement, 1)
  expect_equal(self$ari, 1)
})

test_that("hg_agreement prefers `predicted` and joins on shared nodes", {
  pred <- data.frame(node = c("a", "b", "c", "zzz"),
                     cluster = c("wrong", "wrong", "wrong", "wrong"),
                     predicted = c("A", "A", "B", "B"))
  ref <- data.frame(node = c("a", "b", "c"),
                    cluster = c("A", "B", "B"))
  out <- hg_agreement(pred, ref)
  expect_identical(out$n, 3L)
  expect_equal(out$agreement, 2 / 3)
})

test_that("hg_agreement what = 'table' reproduces the contingency counts", {
  x <- data.frame(node = letters[1:5],
                  cluster = c("g1", "g1", "g1", "g2", "g2"))
  y <- data.frame(node = letters[1:5],
                  cluster = c("h1", "h1", "h2", "h2", "h2"))
  tab <- hg_agreement(x, y, what = "table")
  expect_identical(names(tab), c("label_x", "label_y", "n"))
  expect_identical(tab$n, c(2L, 1L, 2L))
  expect_identical(tab$label_x, c("g1", "g1", "g2"))
  expect_identical(tab$label_y, c("h1", "h2", "h2"))
})

test_that("hg_agreement raises classed errors on bad input", {
  good <- data.frame(node = c("a", "b"), cluster = c("A", "B"))
  expect_error(hg_agreement(good, data.frame(node = c("a", "b"))),
               class = "hypernets_bad_input")
  expect_error(hg_agreement(good, data.frame(node = c("x", "y"),
                                             cluster = c("A", "B"))),
               class = "hypernets_bad_input")
})

test_that("hg_stability(resample = 'seeds') keeps the two-seed solver check", {
  hg <- toy_hg()
  out <- hg_stability(hg, k = 2, type = "random_walk", resample = "seeds")
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("k", "identical_partition", "ari"))
  expect_identical(nrow(out), 1L)
  expect_true(out$identical_partition)
  expect_equal(out$ari, 1)
  two <- hg_stability(hg, k = c(2, 3), type = "random_walk",
                      resample = "seeds")
  expect_identical(two$k, c(2, 3))
})

test_that("hg_stability validates its contract", {
  hg <- toy_hg()
  expect_error(hg_stability(hg, k = 1), "at least 2")
  expect_error(hg_stability(hg, k = 2, resample = "seeds", seeds = c(1, 1)),
               "distinct")
  expect_error(hg_stability(hg, k = 2, resample = "seeds", what = "clusters"),
               class = "hypernets_bad_input")
})

# base kmeans() (Hartigan-Wong) warns "did not converge in 100 iterations"
# or "Quick-TRANSfer stage steps exceeded maximum" on the tied embedding rows
# of planted blocks; hg_cluster() passes those on. Muffle exactly those two
# kmeans warnings, nothing else.
quiet_kmeans <- function(expr) {
  withCallingHandlers(expr, warning = function(w) {
    if (grepl("did not converge in|Quick-TRANSfer stage steps exceeded",
              conditionMessage(w))) {
      invokeRestart("muffleWarning")
    }
  })
}

# three planted topics of 12 documents each, tied by one shared word
planted_hg <- function() {
  old <- if (exists(".Random.seed", envir = globalenv())) .Random.seed
  on.exit(if (!is.null(old)) assign(".Random.seed", old, envir = globalenv()))
  set.seed(42L)
  vocab <- list(a = paste0("a", 1:15), b = paste0("b", 1:15),
                c = paste0("c", 1:15))
  docs <- unlist(lapply(names(vocab), \(t) vapply(seq_len(12), \(i)
    paste(c(sample(vocab[[t]], 6), "common"), collapse = " "), "")))
  names(docs) <- paste0(rep(names(vocab), each = 12), seq_len(12))
  text_hypergraph(docs)
}

test_that("subsampling stability recovers planted topics and the eigengap", {
  hg <- planted_hg()
  out <- quiet_kmeans(hg_stability(hg, k = 2:5, n_boot = 20))
  expect_named(out, c("k", "n_runs", "n_failed", "mean_jaccard",
                      "min_jaccard", "n_stable", "n_dissolved", "eigengap"))
  expect_identical(out$k, 2:5)
  expect_identical(out$n_runs, rep(20L, 4))
  expect_identical(out$n_failed, rep(0L, 4))
  at3 <- subset(out, k == 3)
  expect_equal(at3$mean_jaccard, 1)
  expect_identical(at3$n_stable, 3L)
  # the planted k is the most stable resolution and the largest eigengap
  expect_identical(out$k[which.max(out$mean_jaccard)], 3L)
  expect_identical(out$k[which.max(out$eigengap)], 3L)
  # the eigengap is the gap column of hg_cluster(what = "eigenvalues")
  spectrum <- quiet_kmeans(hg_cluster(hg, k = 5, seed = 1,
                                     what = "eigenvalues"))
  expect_equal(out$eigengap, spectrum$gap[2:5])
  clusters <- quiet_kmeans(hg_stability(hg, k = 3, n_boot = 20,
                                        what = "clusters"))
  expect_named(clusters, c("k", "cluster", "size", "jaccard", "n_dissolved",
                           "n_recovered", "n_runs"))
  expect_identical(clusters$size, rep(12L, 3))
  expect_equal(clusters$jaccard, rep(1, 3))
  expect_identical(clusters$n_recovered, rep(20L, 3))
})

test_that("planted stability holds across seeds; unstructured text is less stable", {
  hg <- planted_hg()
  across <- vapply(c(1L, 7L, 19L, 101L), \(s)
    quiet_kmeans(hg_stability(hg, k = 3, n_boot = 15, seed = s))$mean_jaccard,
    numeric(1))
  expect_equal(across, rep(1, 4))
  # words drawn from one shared pool: no topic structure to recover
  old <- if (exists(".Random.seed", envir = globalenv())) .Random.seed
  set.seed(3)
  pool <- paste0("w", seq_len(40))
  noise <- vapply(seq_len(36), \(i) paste(sample(pool, 7), collapse = " "), "")
  if (!is.null(old)) assign(".Random.seed", old, envir = globalenv())
  names(noise) <- paste0("d", seq_len(36))
  noise_hg <- text_hypergraph(noise)
  noisy <- vapply(c(1L, 7L, 19L), \(s)
    quiet_kmeans(hg_stability(noise_hg, k = 3, n_boot = 15,
                              seed = s))$mean_jaccard,
    numeric(1))
  expect_true(all(noisy < 0.75))
})

test_that("hg_stability is reproducible under a seed and restores the caller's RNG", {
  hg <- planted_hg()
  set.seed(123)
  before <- stats::runif(1)
  set.seed(123)
  a <- quiet_kmeans(hg_stability(hg, k = 2:4, n_boot = 10, seed = 5))
  expect_identical(stats::runif(1), before)
  b <- quiet_kmeans(hg_stability(hg, k = 2:4, n_boot = 10, seed = 5))
  expect_identical(a, b)
})

test_that("a subsample that disconnects the hypergraph is counted and warned", {
  # two topics joined only through the bridge document's two words
  docs <- c(
    a1 = "apple pear plum", a2 = "apple pear fig", a3 = "pear plum fig",
    a4 = "apple plum fig", bridge = "fig kiwi",
    b1 = "kiwi lime lemon", b2 = "kiwi lemon melon", b3 = "lime lemon melon",
    b4 = "kiwi lime melon"
  )
  hg <- text_hypergraph(docs)
  expect_warning(
    out <- hg_stability(hg, k = 2, n_boot = 30, fraction = 0.6),
    class = "hypernets_hypergraph_disconnected"
  )
  expect_gt(out$n_failed, 0L)
  expect_identical(out$n_runs + out$n_failed, 30L)
})

test_that("hg_stability raises classed errors for a bad resampling design", {
  hg <- planted_hg()
  expect_error(hg_stability(hg, k = 2, fraction = 1),
               class = "hypernets_bad_input")
  expect_error(hg_stability(hg, k = 2, fraction = 0),
               class = "hypernets_bad_input")
  expect_error(hg_stability(hg, k = 2, n_boot = 0),
               class = "hypernets_bad_input")
  # floor(0.1 * 36) = 3 nodes cannot hold 3 clusters
  expect_error(hg_stability(hg, k = 3, fraction = 0.1),
               class = "hypernets_bad_input")
})

test_that("the best-match Jaccard follows Hennig's definition", {
  full <- c(x = "A", y = "A", z = "B", w = "B")
  sub <- c(x = "1", y = "2", z = "2", w = "2")
  # A = {x, y}: with 1 = {x} 1/2, with 2 = {y, z, w} 1/4 -> 0.5
  # B = {z, w}: with 2 -> 2/3
  expect_equal(.thg_best_jaccard(full, sub, c("A", "B")),
               c(A = 0.5, B = 2 / 3))
  # a cluster with no member in the subsample scores 0
  expect_equal(.thg_best_jaccard(full, sub, c("A", "B", "C")),
               c(A = 0.5, B = 2 / 3, C = 0))
})

test_that("hg_seeds picks the top-pi nodes per cluster, deterministically", {
  pool <- data.frame(
    node = c("d1", "d2", "d3", "d4", "d5", "d6"),
    cluster = c("A", "A", "A", "B", "B", "B"),
    pi = c(0.3, 0.1, 0.2, 0.05, 0.15, 0.15)
  )
  seeds <- hg_seeds(pool, n = 2)
  # A: d1 (0.3), d3 (0.2). B: pi tie 0.15/0.15 broken alphabetically ->
  # d5 then d6.
  expect_identical(seeds, c(d1 = "A", d3 = "A", d5 = "B", d6 = "B"))
  # n larger than the cluster returns the whole cluster
  all_seeds <- hg_seeds(pool, n = 10)
  expect_identical(length(all_seeds), 6L)
})

test_that("hg_seeds feeds hg_classify directly", {
  hg <- toy_hg()
  pool <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk",
                     what = "embedding")
  seeds <- hg_seeds(pool, n = 1)
  expect_identical(length(seeds), 2L)
  fit <- hg_classify(hg, labels = seeds, type = "random_walk")
  full <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk")
  fit_vs_full <- hg_agreement(fit, full)
  expect_equal(fit_vs_full$agreement, 1)
})

test_that("hg_seeds validates its contract", {
  expect_error(hg_seeds(data.frame(node = "a", cluster = "A")), "embedding")
  expect_error(hg_seeds(data.frame(node = "a", cluster = "A", pi = 1),
                        n = 0), "positive")
})

test_that("classifiers accept labels as a tidy data.frame", {
  hg <- toy_hg()
  as_vector <- c(cooking_1 = "cooking", space_1 = "space")
  as_frame <- data.frame(node = c("cooking_1", "space_1"),
                         label = c("cooking", "space"))
  from_frame <- hg_classify(hg, labels = as_frame)
  from_vector <- hg_classify(hg, labels = as_vector)
  expect_identical(from_frame, from_vector)
  # a hg_cluster() result is a valid labels input directly
  topics <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk")
  one_per <- subset(topics, !duplicated(cluster))
  from_cluster <- hg_classify(hg, labels = one_per)
  from_coerced <- hg_classify(hg, labels = .thg_labels_input(one_per))
  expect_identical(from_cluster$predicted, from_coerced$predicted)
  expect_error(hg_classify(hg, labels = data.frame(node = "cooking_1")),
               "labels")
})

test_that("eigenvalue table carries the gap column", {
  hg <- toy_hg()
  eig <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk",
                    what = "eigenvalues")
  expect_identical(names(eig), c("index", "value", "gap"))
  expect_equal(head(eig$gap, -1), diff(eig$value))
  expect_true(is.na(tail(eig$gap, 1)))
})

test_that("hg_keywords collapse = TRUE matches the long form", {
  hg <- toy_hg()
  topics <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk")
  long <- hg_keywords(hg, topics, n = 3)
  wide <- hg_keywords(hg, topics, n = 3, collapse = TRUE)
  expect_identical(names(wide), c("type", "cluster", "size", "words"))
  expect_identical(wide$size, as.integer(table(topics$cluster)[wide$cluster]))
  expect_identical(nrow(wide), length(unique(long$cluster)))
  rebuilt <- aggregate(word ~ cluster, data = long, paste, collapse = ", ")
  expect_identical(wide$words, rebuilt$word)
  expect_true(all(wide$type == "mass"))
})

test_that("hg_agreement(what = 'mapping') names the best-matching label per row", {
  x <- data.frame(node = c("a", "b", "c", "d", "e"),
                  cluster = c("Cluster 1", "Cluster 1", "Cluster 1",
                              "Cluster 2", "Cluster 2"))
  y <- data.frame(node = c("a", "b", "c", "d", "e"),
                  cluster = c("p", "p", "q", "q", "q"))
  m <- hg_agreement(x, y, what = "mapping")
  expect_named(m, c("label_x", "n", "label_y", "overlap", "share"))
  expect_identical(m$label_x, c("Cluster 1", "Cluster 2"))
  expect_identical(m$n, c(3L, 2L))
  expect_identical(m$label_y, c("p", "q"))
  expect_identical(m$overlap, c(2L, 2L))
  expect_equal(m$share, c(2 / 3, 1))
  # natural order of label_x
  x$cluster <- sub("Cluster 2", "Cluster 10", x$cluster)
  natural <- hg_agreement(x, y, what = "mapping")
  expect_identical(natural$label_x, c("Cluster 1", "Cluster 10"))
})

test_that("hg_agreement summary counts majority-aligned nodes", {
  x <- data.frame(node = letters[1:5],
                  cluster = c("g1", "g1", "g1", "g2", "g2"))
  y <- data.frame(node = letters[1:5],
                  cluster = c("h1", "h1", "h2", "h2", "h2"))
  out <- hg_agreement(x, y)
  mapping <- hg_agreement(x, y, what = "mapping")
  expect_named(out, c("n", "agreement", "aligned", "ari"))
  expect_identical(out$aligned, 4L)
  expect_identical(out$aligned, sum(mapping$overlap))
  # INVARIANT: aligned is bounded by n and equals n for identical labelings
  same <- hg_agreement(x, x)
  expect_identical(same$aligned, 5L)
  expect_true(out$aligned <= out$n)
})

test_that("hg_agreement node and label selectors read any two tables", {
  predictions <- data.frame(node = c("d1", "d2", "d3", "d4"),
                            predicted = c("early", "early", "late", "late"))
  corpus <- data.frame(doc = c("d1", "d2", "d3", "d4"),
                       year = c(2020L, 2020L, 2020L, 2024L))
  tab <- hg_agreement(predictions, corpus, node = c("node", "doc"),
                      label = c("predicted", "year"), what = "table")
  expect_identical(tab$label_x, c("early", "late", "late"))
  expect_identical(tab$label_y, c("2020", "2020", "2024"))
  expect_identical(tab$n, c(2L, 1L, 1L))
  # one name applies to both sides
  renamed <- data.frame(id = c("d1", "d2", "d3", "d4"),
                        predicted = c("early", "early", "late", "late"))
  names(corpus) <- c("id", "predicted")
  corpus$predicted <- c("early", "early", "early", "late")
  one_name <- hg_agreement(renamed, corpus, node = "id")
  expect_identical(one_name$n, 4L)
  expect_equal(one_name$agreement, 3 / 4)
  expect_error(hg_agreement(predictions, corpus, node = "id"),
               class = "hypernets_bad_input")
  expect_error(hg_agreement(predictions, corpus, node = c("node", "id"),
                            label = c("predicted", "nope")),
               class = "hypernets_bad_input")
  expect_error(hg_agreement(predictions, corpus, node = c("a", "b", "c")),
               class = "hypernets_bad_input")
})
