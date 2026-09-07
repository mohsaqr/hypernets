# Delegating verbs: parity with the hypergraph-family engines, tidy shapes, conditions.

verb_corpus <- c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
)
verb_stop <- c("the", "with", "and", "a", "this", "at", "on", "all")

test_that("hg_measures matches the engine and returns tidy tables", {
  hg <- text_hypergraph(verb_corpus, stop_words = verb_stop)
  m <- hypergraph_measures(hg)

  nodes <- hg_measures(hg, what = "nodes")
  expect_identical(nodes$node, names(m$hyperdegree))
  expect_identical(nodes$hyperdegree, unname(as.integer(m$hyperdegree)))
  expect_identical(nodes$strength, unname(as.numeric(m$node_strength)))
  expect_identical(nrow(nodes), hg$n_nodes)

  edges <- hg_measures(hg, what = "edges")
  expect_identical(edges$edge, colnames(hg$incidence))
  expect_identical(edges$size, unname(as.integer(m$edge_sizes)))

  overlap <- hg_measures(hg, what = "overlap")
  expect_identical(
    nrow(overlap),
    (hg$n_hyperedges * (hg$n_hyperedges - 1L)) %/% 2L
  )
  expect_true(all(overlap$jaccard >= 0 & overlap$jaccard <= 1))

  summary_tab <- hg_measures(hg, what = "summary")
  expect_identical(names(summary_tab), c("measure", "value"))
  expect_identical(
    summary_tab$value[summary_tab$measure == "n_nodes"],
    as.numeric(hg$n_nodes)
  )
})

test_that("hg_centrality matches the engine value for value", {
  hg <- text_hypergraph(verb_corpus, stop_words = verb_stop)
  direct <- hypergraph_centrality(hg, type = c("clique", "Z", "H"))

  tab <- hg_centrality(hg)
  expect_identical(names(tab), c("node", "clique", "Z", "H"))
  expect_identical(tab$node, direct$node)
  expect_identical(tab$clique, direct$clique)
  expect_identical(tab$Z, direct$Z)
  expect_identical(tab$H, direct$H)

  one <- hg_centrality(hg, type = "clique")
  expect_identical(names(one), c("node", "clique"))
})

test_that("hg_cluster matches the seeded engine partition and is reproducible", {
  hg <- text_hypergraph(verb_corpus, stop_words = verb_stop)
  direct <- hypergraph_cluster(hg, k = 2, type = "random_walk",
                                          seed = 7)
  expected <- direct$clusters
  rownames(expected) <- NULL

  tab <- hg_cluster(hg, k = 2, type = "random_walk", seed = 7)
  expect_identical(tab, expected)
  again <- hg_cluster(hg, k = 2, type = "random_walk", seed = 7)
  expect_identical(tab, again)
  expect_identical(nrow(tab), hg$n_nodes)
  expect_identical(length(unique(tab$cluster)), 2L)
})

test_that("hg_classify matches the engine and preserves given labels", {
  hg <- text_hypergraph(verb_corpus, stop_words = verb_stop)
  labels <- c(cooking_1 = "cooking", space_1 = "space")
  direct <- hypergraph_transduction(hg, labels = labels)
  expected <- direct$predictions
  rownames(expected) <- NULL

  tab <- hg_classify(hg, labels = labels)
  expect_identical(tab, expected)
  expect_identical(names(tab),
                   c("node", "label", "predicted", "score", "margin"))
  expect_identical(
    tab$label[tab$node == "cooking_1"],
    "cooking"
  )
  expect_true(is.na(tab$label[tab$node == "cooking_2"]))
})

test_that("every verb rejects a non-hypergraph with a classed error", {
  expect_error(hg_measures(42), class = "hypernets_bad_input")
  expect_error(hg_centrality(list()), class = "hypernets_bad_input")
  expect_error(hg_cluster("x", k = 2), class = "hypernets_bad_input")
  expect_error(hg_classify(NULL, labels = c(a = "x")),
               class = "hypernets_bad_input")
})

test_that("hg_centrality sort_by and n select without user-side subsetting", {
  hg <- text_hypergraph(verb_corpus, stop_words = verb_stop)
  full <- hg_centrality(hg, type = "clique")
  top <- hg_centrality(hg, type = "clique", sort_by = "clique", n = 3)
  expect_identical(nrow(top), 3L)
  expect_identical(top$clique, sort(full$clique, decreasing = TRUE)[1:3])
  expect_true(all(diff(top$clique) <= 0))
  expect_error(
    hg_centrality(hg, type = "clique", sort_by = "Z"),
    "clique"
  )
})

# ---- hg_cluster what= and hg_keywords ----------------------------------

test_that("hg_cluster returns embedding and eigenvalues via `what`", {
  hg <- text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
  parts <- hg_cluster(hg, k = 2, seed = 1)
  emb <- hg_cluster(hg, k = 2, seed = 1, what = "embedding")
  expect_named(emb, c("node", "cluster", "pi", "dim1", "dim2"))
  expect_identical(emb$cluster, parts$cluster)
  expect_identical(emb$node, parts$node)
  eig <- hg_cluster(hg, k = 2, seed = 1, what = "eigenvalues")
  expect_named(eig, c("index", "value", "gap"))
  expect_identical(nrow(eig), hg$n_nodes)
  expect_true(all(diff(eig$value) >= -1e-12))  # ascending spectrum
})

test_that("hg_keywords matches a hand-computed fixture", {
  hg <- text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"),
  weight = "n")
  clusters <- data.frame(node = c("cooking_1", "cooking_2",
                                  "space_1", "space_2"),
                         cluster = c("food", "food", "sky", "sky"))
  kw <- hg_keywords(hg, clusters, n = Inf)
  # hand: "soup" count 1+1 in food, 0 in sky -> score 2, share 1
  food_soup <- subset(kw, cluster == "food" & word == "soup")
  expect_identical(food_soup$rank, 1L)
  expect_identical(food_soup$score, 2)
  expect_identical(food_soup$share, 1)
  # hand: "night" once in each cluster -> share 0.5 both sides
  expect_identical(subset(kw, cluster == "food" & word == "night")$share, 0.5)
  expect_identical(subset(kw, cluster == "sky" & word == "night")$share, 0.5)
  # "telescope" is sky-only, count 2
  sky_tel <- subset(kw, cluster == "sky" & word == "telescope")
  expect_identical(sky_tel$score, 2)
  expect_identical(sky_tel$share, 1)
  # invariant: with every node assigned and n = Inf, per-word scores sum
  # to the word's total corpus mass
  total_by_word <- tapply(kw$score, kw$word, sum)
  expect_equal(as.numeric(total_by_word[colnames(hg$incidence)]),
               as.numeric(colSums(hg$incidence)), tolerance = 1e-12)
  expect_true(all(kw$share > 0 & kw$share <= 1))
})

test_that("hg_keywords works on sparse hypergraphs and vector input", {
  hg <- text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"),
  weight = "n", sparse = TRUE)
  labels <- c(cooking_1 = "food", cooking_2 = "food",
              space_1 = "sky", space_2 = "sky")
  kw <- hg_keywords(hg, labels, n = 2)
  expect_identical(nrow(kw), 4L)
  expect_true(all(kw$rank %in% c(1L, 2L)))
  expect_error(hg_keywords(hg, c(zz = "a", cooking_1 = "b")),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(hg, labels, n = 0), "positive")
})

# ---- hg_keywords type = ------------------------------------------------

.kw_fixture <- function(...) {
  text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"),
  ...)
}
.kw_clusters <- data.frame(node = c("cooking_1", "cooking_2",
                                    "space_1", "space_2"),
                           cluster = c("food", "food", "sky", "sky"))

test_that("type = 'frequency' is the raw count whatever the weighting", {
  hg_n <- .kw_fixture(weight = "n")
  hg_tfidf <- .kw_fixture(weight = "tfidf")
  kw_n <- hg_keywords(hg_n, .kw_clusters, n = Inf)
  kw_freq_n <- hg_keywords(hg_n, .kw_clusters, n = Inf, type = "frequency")
  kw_freq_tfidf <- hg_keywords(hg_tfidf, .kw_clusters, n = Inf,
                               type = "frequency")
  values <- c("cluster", "rank", "word", "score", "share")
  # under weight = "n" the incidence mass IS the count
  expect_equal(kw_freq_n[values], kw_n[values])
  # and the count does not move when the incidence weighting changes
  expect_equal(kw_freq_tfidf[values], kw_freq_n[values])
  # invariant: per-word counts sum to the corpus token count of the word
  hg <- .kw_fixture()
  totals <- tapply(kw_freq_tfidf$score, kw_freq_tfidf$word, sum)
  layer <- as.data.frame(hg, what = "weights")
  corpus <- tapply(layer$n, layer$word, sum)
  expect_equal(as.numeric(totals[names(corpus)]), as.numeric(corpus))
})

test_that("type = 'ctfidf' reproduces BERTopic's ClassTfidfTransformer", {
  # oracle: bertopic 0.17.4 ClassTfidfTransformer().fit_transform() on the
  # 2 x 17 class-term count matrix of this fixture (python, 2026-09-05)
  hg <- .kw_fixture()
  kw <- hg_keywords(hg, .kw_clusters, n = Inf, type = "ctfidf")
  pick <- function(cl, w) subset(kw, cluster == cl & word == w)$score
  expect_equal(pick("food", "soup"), 0.358351893845611, tolerance = 1e-12)
  expect_equal(pick("food", "carrots"), 0.23978952727983707,
               tolerance = 1e-12)
  expect_equal(pick("food", "night"), 0.1791759469228055, tolerance = 1e-12)
  expect_equal(pick("sky", "telescope"), 0.358351893845611,
               tolerance = 1e-12)
  expect_equal(pick("sky", "night"), 0.1791759469228055, tolerance = 1e-12)
  # the shared word is split evenly, cluster-only words are fully owned
  expect_equal(subset(kw, word == "night")$share, c(0.5, 0.5))
  expect_true(all(subset(kw, word != "night")$share == 1))
  # hand formula: tf (L1 per class) * log(A / f_w + 1); both clusters hold
  # 10 tokens after stop words, so A = trunc(mean(c(10, 10))) = 10
  expect_equal(pick("food", "soup"), (2 / 10) * log(10 / 2 + 1),
               tolerance = 1e-12)
})

test_that("type = 'centrality' equals hypergraph_centrality on the cluster", {
  hg <- .kw_fixture()
  kw <- hg_keywords(hg, .kw_clusters, n = Inf, type = "centrality",
                    centrality = "clique")
  layer <- as.data.frame(hg, what = "weights")
  sky <- subset(layer, doc %in% c("space_1", "space_2"))
  sky_hg <- group_hypergraph(sky, actor = "word", group = "doc",
                             weight = "weight")
  direct <- hypergraph_centrality(sky_hg, type = "clique")
  got <- subset(kw, cluster == "sky")
  expect_equal(stats::setNames(got$score, got$word)[direct$node],
               stats::setNames(direct$clique, direct$node))
  # only the cluster's own vocabulary is ranked
  expect_false("soup" %in% got$word)
  pr <- hg_keywords(hg, .kw_clusters, n = 2, type = "centrality",
                    centrality = "pagerank", collapse = TRUE)
  expect_named(pr, c("type", "cluster", "size", "words"))
  expect_identical(nrow(pr), 2L)
})

test_that("hg_keywords type = contracts are enforced by class", {
  hg <- .kw_fixture()
  # a window hypergraph has words as nodes and no per-document token layer
  window <- text_hypergraph(c(a = "one two three four", b = "two three five"),
                            construction = "window", window = 2)
  word_groups <- stats::setNames(rep(c("x", "y"), length.out = window$n_nodes),
                                 window$nodes)
  expect_error(hg_keywords(window, word_groups, type = "frequency"),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(window, word_groups, type = "ctfidf"),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(window, word_groups, type = "centrality"),
               class = "hypernets_bad_input")
  # mass still works on any hypergraph
  mass <- hg_keywords(window, word_groups)
  expect_s3_class(mass, "data.frame")
  expect_error(hg_keywords(hg, .kw_clusters, type = "attention"),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(hg, .kw_clusters,
                           scores = data.frame(node = "cooking_1")),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(hg, .kw_clusters,
                           scores = data.frame(node = "zz", word = "w",
                                               attention = 1)),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(hg, .kw_clusters, type = "nope"),
               class = "hypernets_bad_input")
})

test_that("`scores` aggregates a hand-built external table as its own block", {
  hg <- .kw_fixture()
  att <- data.frame(
    node = c("cooking_1", "cooking_2", "space_1", "space_1"),
    word = c("soup", "soup", "stars", "soup"),
    attention = c(0.5, 0.25, 0.9, 0.1)
  )
  kw <- hg_keywords(hg, .kw_clusters, n = Inf, scores = att)
  expect_identical(unique(kw$type), "attention")
  both <- hg_keywords(hg, .kw_clusters, n = Inf, type = "frequency",
                      scores = att)
  expect_identical(unique(both$type), c("frequency", "attention"))
  expect_equal(subset(kw, cluster == "food" & word == "soup")$score, 0.75)
  expect_equal(subset(kw, cluster == "sky" & word == "soup")$score, 0.1)
  expect_equal(subset(kw, cluster == "food" & word == "soup")$share,
               0.75 / 0.85)
  expect_equal(subset(kw, cluster == "sky")$word, c("stars", "soup"))
})

test_that("several `type`s stack into one table and plot", {
  hg <- .kw_fixture()
  one <- lapply(c("frequency", "ctfidf", "centrality"),
                \(t) hg_keywords(hg, .kw_clusters, n = 3, type = t))
  many <- hg_keywords(hg, .kw_clusters, n = 3,
                      type = c("frequency", "ctfidf", "centrality"))
  expect_s3_class(many, "hypernets_keywords")
  expect_named(many, c("type", "cluster", "size", "rank", "word", "score",
                       "share", "n_docs"))
  many_table <- as.data.frame(many)
  stacked <- do.call(rbind, one)
  stacked_table <- as.data.frame(stacked)
  expect_equal(many_table, stacked_table, ignore_attr = TRUE)
  expect_identical(unique(many$type), c("frequency", "ctfidf", "centrality"))
  wide <- hg_keywords(hg, .kw_clusters, n = 3, collapse = TRUE,
                      type = c("ctfidf", "frequency"))
  expect_named(wide, c("type", "cluster", "size", "words"))
  expect_identical(wide$size, rep(2L, 4L))
  expect_identical(wide$type, rep(c("ctfidf", "frequency"), each = 2L))
  p <- plot(many)
  expect_s3_class(p, "ggplot")
  expect_error(plot(wide), class = "hypernets_bad_input")
  expect_error(hg_keywords(hg, .kw_clusters, type = c("mass", "nope")),
               class = "hypernets_bad_input")
  expect_error(hg_keywords(hg, .kw_clusters, type = c("mass", "mass")),
               class = "hypernets_bad_input")
})

test_that("sort_by = 'share' and min_docs rank and filter as documented", {
  hg <- .kw_fixture(weight = "n")
  by_score <- hg_keywords(hg, .kw_clusters, n = Inf)
  by_share <- hg_keywords(hg, .kw_clusters, n = Inf, sort_by = "share")
  # same rows, different order: share descending within cluster
  expect_setequal(paste(by_score$cluster, by_score$word),
                  paste(by_share$cluster, by_share$word))
  food <- subset(by_share, cluster == "food")
  expect_true(all(diff(food$share) <= 0))
  # "night" is shared (share 0.5) so it ranks last under share, and is in
  # one document per cluster
  expect_identical(tail(food$word, 1), "night")
  expect_identical(subset(by_score, word == "night")$n_docs, c(1, 1))
  expect_identical(subset(by_score, cluster == "food" & word == "soup")$n_docs, 2)
  # support floor: only "soup" (food) and "telescope"/"stars" (sky) are in
  # two documents of their cluster
  two <- hg_keywords(hg, .kw_clusters, n = Inf, min_docs = 2)
  expect_setequal(two$word, c("soup", "telescope", "stars"))
  expect_error(hg_keywords(hg, .kw_clusters, min_docs = 0), "min_docs")
  expect_error(hg_keywords(hg, .kw_clusters, sort_by = "nope"))
  p <- plot(by_share, value = "share")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$x, "share")
  two_col <- plot(by_share, ncol = 2)
  expect_s3_class(two_col, "ggplot")
  expect_error(plot(by_share, ncol = 0), "ncol")
})

test_that("clusters come out in natural order, not string order", {
  hg <- .kw_fixture()
  many <- stats::setNames(c("Cluster 10", "Cluster 2", "Cluster 1", "Cluster 10"),
                          c("cooking_1", "cooking_2", "space_1", "space_2"))
  kw <- hg_keywords(hg, many, n = 2, collapse = TRUE)
  expect_identical(kw$cluster, c("Cluster 1", "Cluster 2", "Cluster 10"))
  named <- stats::setNames(c("b", "a", "b", "a"), names(many))
  named_kw <- hg_keywords(hg, named, n = 1, collapse = TRUE)
  expect_identical(named_kw$cluster, c("a", "b"))
  expect_identical(.thg_kw_natural(c("x2", "x10", "x1")), c("x1", "x2", "x10"))
})

# ---- construction = "sentence" and type = "sentence_centrality" -----------

.sent_docs <- c(
  cooking_1 = "Simmer the soup slowly. Add onions and carrots to the soup.",
  cooking_2 = "This soup recipe needs salt; serve on a cold night.",
  space_1 = "The telescope revealed a distant galaxy. Stars everywhere!",
  space_2 = "Astronomers aimed the telescope at the stars all night"
)
.sent_sw <- c("the", "with", "and", "a", "this", "at", "on", "all", "to")

test_that("construction = 'sentence' binds the words of each sentence", {
  sh <- text_hypergraph(.sent_docs, construction = "sentence",
                        stop_words = .sent_sw)
  expect_s3_class(sh, "text_hypergraph")
  expect_identical(sh$text$construction, "sentence")
  expect_identical(sh$text$nodes, "word")
  sents <- as.data.frame(sh, what = "sentences")
  expect_named(sents, c("edge", "doc", "sentence", "n_tokens"))
  # 2 + 2 + 2 + 1 sentences, one hyperedge each, named doc#index
  expect_identical(nrow(sents), 7L)
  expect_identical(sh$n_hyperedges, 7L)
  expect_identical(sents$edge[1:2], c("cooking_1#1", "cooking_1#2"))
  expect_identical(table(sents$doc)[["space_2"]], 1L)
  # "soup" sits in three sentences across two documents; "onions" in one
  expect_identical(sum(sh$incidence["soup", ] > 0), 3L)
  expect_identical(sum(sh$incidence["onions", ] > 0), 1L)
  # incidence weight is the in-sentence count
  two <- text_hypergraph(c(a = "salt salt soup. soup"),
                         construction = "sentence")
  expect_identical(as.numeric(two$incidence["salt", "a#1"]), 2)
  expect_identical(as.numeric(two$incidence["soup", "a#2"]), 1)
  # the vocabulary equals the bag construction's under the same filters
  bag <- text_hypergraph(.sent_docs, stop_words = .sent_sw)
  sh_vocab <- as.data.frame(sh, what = "vocabulary")
  bag_vocab <- as.data.frame(bag, what = "vocabulary")
  expect_identical(sh_vocab$word, bag_vocab$word)
  # sparse gives the same incidence
  sparse <- text_hypergraph(.sent_docs, construction = "sentence",
                            stop_words = .sent_sw, sparse = TRUE)
  expect_equal(as.matrix(sparse$incidence), sh$incidence)
  # min_count filters inside sentences too
  mc <- text_hypergraph(.sent_docs, construction = "sentence",
                        stop_words = .sent_sw, min_count = 2L)
  expect_setequal(mc$nodes, c("soup", "telescope", "stars", "night"))
  expect_error(text_hypergraph(.sent_docs, construction = "sentence",
                               weight = "tfidf"), class = "hypernets_bad_input")
  expect_error(as.data.frame(bag, what = "sentences"),
               class = "hypernets_bad_input")
  expect_output(print(sh), "sentence hyperedges: 7 sentences")
})

test_that("a sentence hypergraph scopes centrality to the cluster's sentences", {
  hg <- text_hypergraph(.sent_docs, stop_words = .sent_sw)
  sh <- text_hypergraph(.sent_docs, construction = "sentence",
                        stop_words = .sent_sw)
  clusters <- data.frame(node = names(.sent_docs),
                         cluster = c("food", "food", "sky", "sky"))
  kw <- hg_keywords(sh, clusters, n = Inf, type = "centrality",
                    centrality = "clique")
  # direct: the sky sentences only
  sents <- as.data.frame(sh, what = "sentences")
  sky_edges <- sents$edge[sents$doc %in% c("space_1", "space_2")]
  sh_table <- as.data.frame(sh)
  sky <- subset(sh_table, edge %in% sky_edges)
  sky_hg <- group_hypergraph(sky, actor = "word", group = "edge",
                             weight = "weight")
  direct <- hypergraph_centrality(sky_hg, type = "clique")
  got <- subset(kw, cluster == "sky")
  expect_equal(stats::setNames(got$score, got$word)[direct$node],
               stats::setNames(direct$clique, direct$node))
  # support counts documents, not sentences: "soup" is in 3 sentences of
  # 2 food documents
  expect_identical(subset(kw, cluster == "food" & word == "soup")$n_docs, 2)
  # scope changes the ranking: "add"/"carrots" share a sentence with soup
  doc_scope <- hg_keywords(hg, clusters, n = 3, type = "centrality",
                           centrality = "clique")
  expect_false(identical(subset(kw, cluster == "food")$word[1:3],
                         subset(doc_scope, cluster == "food")$word[1:3]))
  # mass under sentence scope sums the in-sentence counts per document
  mass <- hg_keywords(sh, clusters, n = Inf)
  expect_identical(subset(mass, cluster == "food" & word == "soup")$score, 3)
  other <- text_hypergraph(c(zz = "unrelated words here."),
                           construction = "sentence")
  expect_error(hg_keywords(other, clusters), class = "hypernets_bad_input")
})

test_that("hg_cluster(edge_weights) reaches both engines and 'idf' reads the vocabulary", {
  docs <- c(a = "soup salt onions soup", b = "soup salt carrots",
            c = "stars galaxy telescope", d = "stars telescope night",
            e = "night soup stars")
  dense <- text_hypergraph(docs, weight = "tfidf")
  sparse <- text_hypergraph(docs, weight = "tfidf", sparse = TRUE)
  # the default cut ignores the cells: tfidf and counts agree exactly
  counts <- text_hypergraph(docs, weight = "n")
  e_default <- hg_cluster(dense, k = 2, seed = 1, what = "eigenvalues")
  e_counts <- hg_cluster(counts, k = 2, seed = 1, what = "eigenvalues")
  expect_equal(e_default$value, e_counts$value)
  # idf weights change the spectrum, identically on both engines
  e_dense <- hg_cluster(dense, k = 2, seed = 1, edge_weights = "idf",
                        what = "eigenvalues")
  e_sparse <- hg_cluster(sparse, k = 2, seed = 1, edge_weights = "idf",
                         what = "eigenvalues")
  # the sparse solver returns the leading eigenvalues only
  n_common <- min(nrow(e_dense), nrow(e_sparse))
  expect_equal(head(e_dense$value, n_common), head(e_sparse$value, n_common),
               tolerance = 1e-8)
  e_unweighted <- hg_cluster(dense, k = 2, seed = 1, what = "eigenvalues")
  expect_false(isTRUE(all.equal(e_dense$value, e_unweighted$value)))
  # explicit numeric weights equal to the idf give the same answer
  vocab <- as.data.frame(dense, what = "vocabulary")
  idf <- vocab$idf
  names(idf) <- vocab$word
  e_num <- hg_cluster(dense, k = 2, seed = 1,
                      edge_weights = as.numeric(idf[colnames(dense$incidence)]),
                      what = "eigenvalues")
  expect_equal(e_num$value, e_dense$value)
  # a scalar recycles to unit weights, i.e. the default
  e_scalar <- hg_cluster(dense, k = 2, seed = 1, edge_weights = 2,
                         what = "eigenvalues")
  e_unit <- hg_cluster(dense, k = 2, seed = 1, what = "eigenvalues")
  expect_equal(e_scalar$value, e_unit$value)
  expect_error(hg_cluster(counts, k = 2, edge_weights = "idf"),
               class = "hypernets_bad_input")
  expect_error(hg_cluster(dense, k = 2, edge_weights = c(1, 2)),
               class = "hypernets_bad_input")
  expect_error(hg_cluster(dense, k = 2, edge_weights = -1),
               class = "hypernets_bad_input")
})

test_that("the keyword print method is compact and the default centrality is pagerank", {
  hg <- .kw_fixture()
  kw <- hg_keywords(hg, .kw_clusters, n = 3, type = c("frequency", "ctfidf"))
  out <- capture.output(print(kw))
  expect_true(any(grepl("^ *type +cluster +size +words", out)))
  expect_identical(sum(grepl("Cluster|food|sky", out)), 4L)
  expect_true(any(grepl("12 rows in the long form", out)))
  long <- as.data.frame(kw)
  expect_identical(nrow(long), 12L)
  pr <- hg_keywords(hg, .kw_clusters, n = Inf, type = "centrality")
  explicit <- hg_keywords(hg, .kw_clusters, n = Inf, type = "centrality",
                          centrality = "pagerank")
  expect_equal(pr$score, explicit$score)
})

test_that("hg_relations is the bibliometric co-occurrence of topics through words", {
  hg <- .kw_fixture(weight = "n")
  rel <- hg_relations(hg, .kw_clusters)
  expect_named(rel, c("source", "target", "weight"))
  expect_identical(nrow(rel), 1L)
  expect_identical(rel$source, "food")
  expect_identical(rel$target, "sky")
  # hand: topic x word document counts, then the sum over words of the
  # products; only "night" is in both topics, one document each
  counts <- rbind(food = colSums((hg$incidence[c("cooking_1", "cooking_2"), ] > 0) * 1),
                  sky = colSums((hg$incidence[c("space_1", "space_2"), ] > 0) * 1))
  co <- tcrossprod(counts)
  expect_equal(rel$weight, co["food", "sky"])
  expect_equal(rel$weight, 1)
  d <- diag(co)
  association <- hg_relations(hg, .kw_clusters, similarity = "association")
  expect_equal(association$weight,
               co["food", "sky"] / (d[["food"]] * d[["sky"]]))
  cosine <- hg_relations(hg, .kw_clusters, similarity = "cosine")
  expect_equal(cosine$weight,
               co["food", "sky"] / sqrt(d[["food"]] * d[["sky"]]))
  jaccard <- hg_relations(hg, .kw_clusters, similarity = "jaccard")
  expect_equal(jaccard$weight,
               co["food", "sky"] / (d[["food"]] + d[["sky"]] - co["food", "sky"]))
  inclusion <- hg_relations(hg, .kw_clusters, similarity = "inclusion")
  expect_equal(inclusion$weight,
               co["food", "sky"] / min(d[["food"]], d[["sky"]]))
  equivalence <- hg_relations(hg, .kw_clusters, similarity = "equivalence")
  expect_equal(equivalence$weight,
               co["food", "sky"]^2 / (d[["food"]] * d[["sky"]]))
  # three topics: pairs in natural order, zero-weight pairs dropped
  three <- stats::setNames(c("Cluster 1", "Cluster 2", "Cluster 10", "Cluster 10"),
                           .kw_clusters$node)
  rel3 <- hg_relations(hg, three)
  expect_true(all(paste(rel3$source, rel3$target) %in%
                    c("Cluster 1 Cluster 2", "Cluster 1 Cluster 10",
                      "Cluster 2 Cluster 10")))
  expect_true(all(rel3$weight > 0))
  net <- hg_relations(hg, .kw_clusters, what = "network")
  expect_s3_class(net, "cograph_network")
  expect_identical(nrow(net$edges), 1L)
  expect_setequal(net$nodes$name, c("food", "sky"))
  expect_identical(net$nodes$size, c(2L, 2L))
  expect_error(hg_relations(hg, c(zz = "a", cooking_1 = "b")),
               class = "hypernets_bad_input")
  expect_identical(hypergraph_relations, hg_relations)
})
