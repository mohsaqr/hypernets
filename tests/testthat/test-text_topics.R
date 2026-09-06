skip_on_cran()

.topic_hg <- function() {
  text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night",
    space_3 = "the stars and the telescope"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
}
.topic_clusters <- data.frame(
  node = c("cooking_1", "cooking_2", "space_1", "space_2", "space_3"),
  cluster = c("food", "food", "sky", "sky", "sky")
)

test_that("hg_topic_sizes counts documents and weights them", {
  hg <- .topic_hg()
  sizes <- hg_topic_sizes(hg, .topic_clusters)
  expect_s3_class(sizes, "honets_topic_sizes")
  expect_named(sizes, c("topic", "n", "share"))
  expect_identical(sizes$topic, c("food", "sky"))
  expect_identical(sizes$n, c(2L, 3L))
  expect_equal(sizes$share, c(0.4, 0.6))
  w <- c(cooking_1 = 10, cooking_2 = 10, space_1 = 1, space_2 = 1, space_3 = 1)
  weighted <- hg_topic_sizes(hg, .topic_clusters, weights = w)
  expect_named(weighted, c("topic", "n", "share", "weighted_n",
                           "weighted_share"))
  expect_equal(weighted$weighted_n, c(20, 3))
  expect_equal(weighted$weighted_share, c(20, 3) / 23)
  # natural order and the plot
  many <- stats::setNames(c("Cluster 10", "Cluster 2", "Cluster 1",
                            "Cluster 10", "Cluster 2"), .topic_clusters$node)
  expect_identical(hg_topic_sizes(hg, many)$topic,
                   c("Cluster 1", "Cluster 2", "Cluster 10"))
  expect_s3_class(plot(sizes), "ggplot")
  expect_s3_class(plot(weighted), "ggplot")
  expect_error(hg_topic_sizes(hg, .topic_clusters, weights = c(cooking_1 = 1)),
               class = "honets_bad_input")
  expect_error(hg_topic_sizes(hg, c(zz = "a", cooking_1 = "b")),
               class = "honets_bad_input")
  expect_identical(hypergraph_topic_sizes, hg_topic_sizes)
})

test_that("hg_topic_quality scores coherence by NPMI and exclusivity by share", {
  hg <- .topic_hg()
  q <- hg_topic_quality(hg, .topic_clusters, n = 2)
  expect_s3_class(q, "honets_topic_quality")
  expect_named(q, c("topic", "size", "n_words", "coherence", "exclusivity"))
  expect_identical(q$topic, c("food", "sky"))
  expect_identical(q$size, c(2L, 3L))
  # sky's two most owned words with support are "stars" and "telescope",
  # both in all three sky documents: p_ij = p_i = p_j = 1, NPMI 1 by
  # convention, and both owned outright
  sky <- subset(q, topic == "sky")
  expect_equal(sky$coherence, 1)
  expect_equal(sky$exclusivity, 1)
  # hand NPMI: two words each in one of two documents, never together
  present <- rbind(a = c(1, 0), b = c(0, 1))
  expect_equal(.thg_npmi(present), -1)
  # two words, one document each, always together
  expect_equal(.thg_npmi(rbind(a = c(1, 1), b = c(0, 0))), 1)
  # independence: word 1 in docs 1,2; word 2 in docs 1,3 of 4 -> p_ij = 1/4
  # = p_i p_j -> NPMI 0
  indep <- rbind(c(1, 1), c(1, 0), c(0, 1), c(0, 0))
  expect_equal(.thg_npmi(indep), 0)
  expect_true(is.na(.thg_npmi(matrix(1, 3, 1))))
  expect_s3_class(plot(q), "ggplot")
  expect_error(hg_topic_quality(hg, c(zz = "a", cooking_1 = "b")),
               class = "honets_bad_input")
  expect_identical(hypergraph_topic_quality, hg_topic_quality)
})

test_that("hg_membership gives fuzzy weights that sum to one and favour the hard label", {
  hg <- .topic_hg()
  topics <- hg_cluster(hg, k = 2, seed = 1)
  m <- hg_membership(hg, topics)
  expect_s3_class(m, "honets_membership")
  expect_named(m, c("node", "cluster", "topic", "membership"))
  expect_identical(nrow(m), 5L * 2L)
  sums <- tapply(m$membership, m$node, sum)
  expect_equal(as.numeric(sums), rep(1, 5))
  expect_true(all(m$membership >= 0 & m$membership <= 1))
  # the hard label is the topic of maximal membership for every document
  best <- vapply(split(m, m$node), \(d) d$topic[which.max(d$membership)],
                 character(1))
  hard <- stats::setNames(topics$cluster, topics$node)
  expect_identical(best[names(hard)], hard)
  # a document exactly at a centre has membership 1 there: a topic with one
  # document is its own centre
  alone <- c(cooking_1 = "A", cooking_2 = "B", space_1 = "B", space_2 = "B",
             space_3 = "B")
  m1 <- hg_membership(hg, alone)
  expect_equal(subset(m1, node == "cooking_1" & topic == "A")$membership, 1)
  expect_s3_class(plot(m), "ggplot")
  expect_error(hg_membership(hg, c(zz = "a", cooking_1 = "b")),
               class = "honets_bad_input")
  expect_identical(hypergraph_membership, hg_membership)
})
