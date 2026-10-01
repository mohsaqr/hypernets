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
  expect_s3_class(sizes, "hypernets_topic_sizes")
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
  many_sizes <- hg_topic_sizes(hg, many)
  expect_identical(many_sizes$topic,
                   c("Cluster 1", "Cluster 2", "Cluster 10"))
  sizes_plot <- plot(sizes)
  expect_s3_class(sizes_plot, "ggplot")
  weighted_plot <- plot(weighted)
  expect_s3_class(weighted_plot, "ggplot")
  expect_error(hg_topic_sizes(hg, .topic_clusters, weights = c(cooking_1 = 1)),
               class = "hypernets_bad_input")
  expect_error(hg_topic_sizes(hg, c(zz = "a", cooking_1 = "b")),
               class = "hypernets_bad_input")
  expect_identical(hg_topic_sizes, hg_topic_sizes)
})

test_that("hg_topic_quality keeps the pre-0.6.0 measures behind explicit arguments", {
  hg <- .topic_hg()
  q <- hg_topic_quality(hg, .topic_clusters, n = 2, coherence = "npmi_cluster",
                        exclusivity = "share", sort_by = "share")
  expect_s3_class(q, "hypernets_topic_quality")
  expect_named(q, c("topic", "size", "n_words", "coherence", "exclusivity",
                    "coherence_type", "exclusivity_type"))
  expect_identical(q$topic, c("food", "sky"))
  expect_identical(q$size, c(2L, 3L))
  # sky's two most owned words with support are "stars" and "telescope",
  # both in all three sky documents: p_ij = p_i = p_j = 1, NPMI 1 by
  # convention, and both owned outright
  sky <- subset(q, topic == "sky")
  expect_equal(sky$coherence, 1)
  expect_equal(sky$exclusivity, 1)
  expect_identical(unique(q$coherence_type), "npmi_cluster")
  expect_identical(unique(q$exclusivity_type), "share")
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
  quality_plot <- plot(q)
  expect_s3_class(quality_plot, "ggplot")
})

test_that("UMass and corpus NPMI match hand computation on the toy corpus", {
  hg <- .topic_hg()
  # soup: cooking_1, cooking_2 (D = 2); salt: cooking_2 (D = 1); together in
  # one document, N = 5. UMass conditions on the higher-ranked word (soup).
  words <- data.frame(topic = "t", word = c("soup", "salt"))
  umass <- hg_topic_quality(hg, words = words, exclusivity = "none")
  expect_equal(umass$coherence, log((1 + 0.01) / (2 + 0.01)))
  expect_true(is.na(umass$size))
  expect_true(is.na(umass$exclusivity))
  npmi <- hg_topic_quality(hg, words = words, coherence = "npmi",
                           exclusivity = "none")
  expect_equal(npmi$coherence, log(0.2 / (0.4 * 0.2)) / -log(0.2))
  # rank order matters for UMass: salt first conditions on salt (D = 1)
  flipped <- data.frame(topic = "t", word = c("soup", "salt"), rank = 2:1)
  expect_equal(hg_topic_quality(hg, words = flipped,
                                exclusivity = "none")$coherence,
               log((1 + 0.01) / (1 + 0.01)))
  # a single word has no pair
  one <- hg_topic_quality(hg, words = data.frame(topic = "t", word = "soup"),
                          exclusivity = "none")
  expect_true(is.na(one$coherence))
})

test_that("FREX matches the stm formula on a hand-built distribution", {
  # A = (4, 2, 1), B = (1, 2, 4) over words x, y, z. For A: beta ranks and
  # exclusivity ranks are both (3, 2, 1) / 3, so FREX(x) = 1, FREX(y) = 2/3;
  # the sum over the top two words is 5/3 whatever the weight.
  counts <- rbind(A = c(x = 4, y = 2, z = 1), B = c(x = 1, y = 2, z = 4))
  expect_equal(.thg_frex(counts, "A", c("x", "y"), 0.7), 5 / 3)
  expect_equal(.thg_frex(counts, "A", c("x", "y"), 0.2), 5 / 3)
  # word share: mean of 4/5 and 2/4
  expect_equal(.thg_word_share(counts, "A", c("x", "y")), mean(c(0.8, 0.5)))
})

test_that("hg_topic_quality defaults to UMass + FREX on cluster word distributions", {
  hg <- .topic_hg()
  q <- hg_topic_quality(hg, .topic_clusters, n = 3)
  expect_identical(unique(q$coherence_type), "umass")
  expect_identical(unique(q$exclusivity_type), "frex")
  expect_identical(q$n_words, c(3L, 3L))
  # sky's three most frequent words: stars (3 tokens), telescope (3), then
  # the first one-token word in vocabulary order
  expect_true(all(q$exclusivity > 0 & q$exclusivity <= 3))
  expect_true(all(q$coherence <= 0))
})

test_that("topic quality invariants: relabelling, bounds, owned vocabulary", {
  hg <- .topic_hg()
  base <- hg_topic_quality(hg, .topic_clusters, n = 3, coherence = "npmi")
  relabelled <- transform(.topic_clusters,
                          cluster = ifelse(cluster == "food", "z_food", "a_sky"))
  moved <- hg_topic_quality(hg, relabelled, n = 3, coherence = "npmi")
  expect_equal(moved$coherence[match(c("z_food", "a_sky"), moved$topic)],
               base$coherence[match(c("food", "sky"), base$topic)])
  expect_equal(moved$exclusivity[match(c("z_food", "a_sky"), moved$topic)],
               base$exclusivity[match(c("food", "sky"), base$topic)])
  expect_true(all(base$coherence >= -1 & base$coherence <= 1))
  # UMass is never positive (D(m, l) <= D(l))
  umass <- hg_topic_quality(hg, .topic_clusters, n = 4)
  expect_true(all(umass$coherence <= 0))
  # topics that share no words own their vocabulary: share exclusivity 1
  owned <- hg_topic_quality(hg, .topic_clusters, n = 2,
                            exclusivity = "share", sort_by = "share")
  expect_equal(owned$exclusivity, c(1, 1))
})

test_that("hg_topic_quality raises classed errors on bad input", {
  hg <- .topic_hg()
  expect_error(hg_topic_quality(hg, c(zz = "a", cooking_1 = "b")),
               class = "hypernets_bad_input")
  expect_error(hg_topic_quality(hg), class = "hypernets_bad_input")
  words <- data.frame(topic = "t", word = c("soup", "salt"))
  # exclusivity needs the cluster distributions
  expect_error(hg_topic_quality(hg, words = words),
               class = "hypernets_bad_input")
  expect_error(hg_topic_quality(hg, words = words, coherence = "npmi_cluster",
                                exclusivity = "none"),
               class = "hypernets_bad_input")
  expect_error(hg_topic_quality(hg, words = data.frame(topic = "t",
                                                       word = "zebra"),
                                exclusivity = "none"),
               class = "hypernets_bad_input")
  expect_error(hg_topic_quality(hg, words = data.frame(topic = "t")),
               class = "hypernets_bad_input")
  expect_error(hg_topic_quality(hg, .topic_clusters, words = words),
               class = "hypernets_bad_input")
  expect_error(hg_topic_quality(hg, .topic_clusters, frexw = 1),
               class = "hypernets_bad_input")
  word_hg <- text_hypergraph(c(
    cooking_1 = "simmer the soup with onions",
    space_1 = "the telescope and the soup"
  ), nodes = "word")
  expect_error(hg_topic_quality(word_hg, words = words,
                                exclusivity = "none"),
               class = "hypernets_bad_input")
})

test_that("hg_membership gives fuzzy weights that sum to one and favour the hard label", {
  hg <- .topic_hg()
  topics <- hg_cluster(hg, k = 2, seed = 1)
  m <- hg_membership(hg, topics)
  expect_s3_class(m, "hypernets_membership")
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
  membership_plot <- plot(m)
  expect_s3_class(membership_plot, "ggplot")
  expect_error(hg_membership(hg, c(zz = "a", cooking_1 = "b")),
               class = "hypernets_bad_input")
  expect_identical(hg_membership, hg_membership)
})

test_that("words = takes hg_keywords() output and a terms()-style matrix as they are", {
  hg <- text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots at night",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars at night",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
  topics <- hg_cluster(hg, k = 2, seed = 1)
  keywords <- hg_keywords(hg, topics, n = 3L)
  long <- data.frame(topic = keywords$cluster, word = keywords$word,
                     rank = keywords$rank, stringsAsFactors = FALSE)
  expect_identical(hg_topic_quality(hg, topics, words = keywords),
                   hg_topic_quality(hg, topics, words = long))
  mat <- matrix(c("soup", "night", "telescope", "stars"), nrow = 2L,
                dimnames = list(NULL, c("Topic 1", "Topic 2")))
  as_long <- data.frame(topic = c("Topic 1", "Topic 1", "Topic 2", "Topic 2"),
                        word = c("soup", "night", "telescope", "stars"),
                        rank = c(1L, 2L, 1L, 2L), stringsAsFactors = FALSE)
  expect_identical(hg_topic_quality(hg, words = mat, exclusivity = "none"),
                   hg_topic_quality(hg, words = as_long, exclusivity = "none"))
  expect_error(hg_topic_quality(hg, words = matrix(1:4, 2L), exclusivity = "none"),
               class = "hypernets_bad_input")
})
