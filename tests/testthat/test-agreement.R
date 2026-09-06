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
  expect_equal(hg_agreement(x, x)$agreement, 1)
  expect_equal(hg_agreement(x, x)$ari, 1)
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
               class = "honets_bad_input")
  expect_error(hg_agreement(good, data.frame(node = c("x", "y"),
                                             cluster = c("A", "B"))),
               class = "honets_bad_input")
})

test_that("hg_stability reports one tidy row per k and finds planted stability", {
  hg <- toy_hg()
  out <- hg_stability(hg, k = 2, type = "random_walk")
  expect_s3_class(out, "data.frame")
  expect_identical(names(out), c("k", "identical_partition", "ari"))
  expect_identical(nrow(out), 1L)
  expect_true(out$identical_partition)
  expect_equal(out$ari, 1)
  two <- hg_stability(hg, k = c(2, 3), type = "random_walk")
  expect_identical(two$k, c(2, 3))
})

test_that("hg_stability validates its contract", {
  hg <- toy_hg()
  expect_error(hg_stability(hg, k = 1), "at least 2")
  expect_error(hg_stability(hg, k = 2, seeds = c(1, 1)), "distinct")
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
  expect_identical(length(hg_seeds(pool, n = 10)), 6L)
})

test_that("hg_seeds feeds hg_classify directly", {
  hg <- toy_hg()
  pool <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk",
                     what = "embedding")
  seeds <- hg_seeds(pool, n = 1)
  expect_identical(length(seeds), 2L)
  fit <- hg_classify(hg, labels = seeds, type = "random_walk")
  full <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk")
  expect_equal(hg_agreement(fit, full)$agreement, 1)
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
  expect_identical(hg_classify(hg, labels = as_frame),
                   hg_classify(hg, labels = as_vector))
  # a hg_cluster() result is a valid labels input directly
  topics <- hg_cluster(hg, k = 2, seed = 1, type = "random_walk")
  one_per <- subset(topics, !duplicated(cluster))
  expect_identical(hg_classify(hg, labels = one_per)$predicted,
                   hg_classify(hg,
                               labels = .thg_labels_input(one_per))$predicted)
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
  expect_identical(hg_agreement(x, y, what = "mapping")$label_x,
                   c("Cluster 1", "Cluster 10"))
})
