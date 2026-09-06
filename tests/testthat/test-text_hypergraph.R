# Construction: tokenization, weighting, both node orientations, conditions.

tiny_corpus <- c(a = "salt and soup", b = "soup and stars")

test_that("counts and the weights table are exact for a hand-built corpus", {
  hg <- text_hypergraph(tiny_corpus)
  expect_s3_class(hg, "text_hypergraph")
  expect_s3_class(hg, "net_hypergraph")
  expect_identical(
    as.data.frame(hg),
    data.frame(
      doc = c("a", "a", "a", "b", "b", "b"),
      word = c("and", "salt", "soup", "and", "soup", "stars"),
      n = 1L,
      weight = 1
    )
  )
})

test_that("repeated tokens are counted, not deduplicated", {
  hg <- text_hypergraph(c(d = "soup soup onion"))
  tab <- as.data.frame(hg)
  expect_identical(tab$n, c(1L, 2L))
  expect_identical(tab$word, c("onion", "soup"))
})

test_that("tf-idf matches the hand-computed smoothed formula", {
  hg <- text_hypergraph(tiny_corpus, weight = "tfidf")
  tab <- as.data.frame(hg)
  # N = 2; df(and) = df(soup) = 2 -> idf = log(3/3) + 1 = 1
  # df(salt) = df(stars) = 1 -> idf = log(3/2) + 1
  rare <- log(3 / 2) + 1
  expect_equal(
    tab$weight,
    c(1, rare, 1, 1, 1, rare)
  )
  vocab <- as.data.frame(hg, what = "vocabulary")
  expect_identical(vocab$word, c("and", "salt", "soup", "stars"))
  expect_identical(vocab$doc_freq, c(2L, 1L, 2L, 1L))
  expect_equal(vocab$idf, c(1, rare, 1, rare))
})

test_that("node orientation controls which entity is the vertex set", {
  doc_hg <- text_hypergraph(tiny_corpus, nodes = "doc")
  word_hg <- text_hypergraph(tiny_corpus, nodes = "word")
  expect_identical(doc_hg$n_nodes, 2L)
  expect_identical(doc_hg$n_hyperedges, 4L)
  expect_identical(word_hg$n_nodes, 4L)
  expect_identical(word_hg$n_hyperedges, 2L)
  expect_identical(sort(doc_hg$nodes), c("a", "b"))
  expect_identical(sort(word_hg$nodes), c("and", "salt", "soup", "stars"))
})

test_that("total incidence weight equals the weights table in both modes", {
  doc_hg <- text_hypergraph(tiny_corpus, nodes = "doc", weight = "tfidf")
  word_hg <- text_hypergraph(tiny_corpus, nodes = "word", weight = "tfidf")
  expect_equal(sum(doc_hg$incidence), sum(as.data.frame(doc_hg)$weight))
  expect_equal(sum(word_hg$incidence), sum(as.data.frame(word_hg)$weight))
  expect_true(all(doc_hg$incidence >= 0))
})

test_that("stop words and min_count filter the vocabulary", {
  hg <- text_hypergraph(tiny_corpus, stop_words = "and")
  expect_false("and" %in% as.data.frame(hg, what = "vocabulary")$word)

  hg2 <- text_hypergraph(tiny_corpus, min_count = 2L)
  expect_identical(
    as.data.frame(hg2, what = "vocabulary")$word,
    c("and", "soup")
  )
})

test_that("data.frame input keeps IDs and carries metadata into documents", {
  articles <- data.frame(
    key = c("x1", "x2"),
    abstract = c("salt and soup", "soup and stars"),
    year = c(2020L, 2021L)
  )
  hg <- text_hypergraph(articles, column = "abstract", id = "key")
  docs <- as.data.frame(hg, what = "documents")
  expect_identical(docs$doc, c("x1", "x2"))
  expect_identical(docs$year, c(2020L, 2021L))
  expect_identical(docs$n_tokens, c(3L, 3L))
  expect_identical(docs$n_types, c(3L, 3L))
})

test_that("documents emptied by filtering are dropped with a classed warning", {
  expect_warning(
    hg <- text_hypergraph(c(a = "salt and soup", b = "and", c = "soup"),
                          stop_words = "and"),
    class = "honets_dropped_documents"
  )
  expect_identical(as.data.frame(hg, what = "documents")$doc, c("a", "c"))
})

test_that("contract violations raise classed errors", {
  expect_error(
    text_hypergraph(c(a = "and", b = "and"), stop_words = "and"),
    class = "honets_empty_corpus"
  )
  expect_error(
    text_hypergraph(data.frame(txt = "salt")),
    class = "honets_bad_input"
  )
  expect_error(
    text_hypergraph(data.frame(id = c("a", "a"), txt = c("x", "y")),
                    column = "txt", id = "id"),
    class = "honets_bad_input"
  )
})

test_that("construction is deterministic", {
  expect_identical(
    text_hypergraph(tiny_corpus, weight = "tfidf"),
    text_hypergraph(tiny_corpus, weight = "tfidf")
  )
})

test_that("print announces the corpus and delegates to the engine", {
  hg <- text_hypergraph(tiny_corpus)
  expect_output(print(hg), "Text hypergraph: 2 documents, 4 words")
  expect_invisible(print(hg))
})

test_that("covid_abstracts dataset is intact", {
  expect_identical(dim(covid_abstracts), c(165L, 4L))
  expect_identical(names(covid_abstracts),
                   c("doc", "title", "abstract", "year"))
  expect_identical(anyDuplicated(covid_abstracts$doc), 0L)
  expect_identical(sort(unique(covid_abstracts$year)), 2020:2024)
})

test_that("curly apostrophes keep possessives as one token", {
  hg <- text_hypergraph(c(d = "the children\u2019s teacher's plan"))
  expect_identical(
    as.data.frame(hg, what = "vocabulary")$word,
    c("children's", "plan", "teacher's", "the")
  )
})

test_that("stop_words_en is a clean, deterministic function-word list", {
  s <- stop_words_en()
  expect_identical(s, sort(s))
  expect_identical(anyDuplicated(s), 0L)
  expect_identical(s, tolower(s))
  expect_true(all(c("the", "and", "of") %in% s))
  expect_false(any(c("covid", "study", "learning") %in% s))
})

test_that("sparse storage is chosen automatically above a million cells", {
  small <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
  expect_false(inherits(small$incidence, "dgCMatrix"))
  expect_true(inherits(
    text_hypergraph(c(a = "salt and soup", b = "soup and stars"), sparse = TRUE)$incidence,
    "dgCMatrix"))
  # 1,200 documents x 900 words crosses the threshold (alphabetic words:
  # the tokeniser drops digits)
  set.seed(1)
  vocab <- head(apply(expand.grid(letters, letters, letters), 1, paste,
                      collapse = ""), 900)
  docs <- vapply(seq_len(1200), \(i) paste(sample(vocab, 12), collapse = " "),
                 character(1))
  big <- text_hypergraph(docs)
  expect_true(inherits(big$incidence, "dgCMatrix"))
  expect_false(inherits(text_hypergraph(docs, sparse = FALSE)$incidence,
                        "dgCMatrix"))
  expect_error(text_hypergraph(docs, sparse = NA), "sparse")
})
