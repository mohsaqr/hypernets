# Construction: tokenization, weighting, both node orientations, conditions.

tiny_corpus <- c(a = "salt and soup", b = "soup and stars")

test_that("counts and the weights table are exact for a hand-built corpus", {
  hg <- text_hypergraph(tiny_corpus)
  expect_s3_class(hg, "text_hypergraph")
  expect_s3_class(hg, "net_hypergraph")
  tab <- as.data.frame(hg)
  expect_identical(
    tab,
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
  doc_tab <- as.data.frame(doc_hg)
  word_tab <- as.data.frame(word_hg)
  expect_equal(sum(doc_hg$incidence), sum(doc_tab$weight))
  expect_equal(sum(word_hg$incidence), sum(word_tab$weight))
  expect_true(all(doc_hg$incidence >= 0))
})

test_that("stop words and min_count filter the vocabulary", {
  hg <- text_hypergraph(tiny_corpus, stop_words = "and")
  vocab <- as.data.frame(hg, what = "vocabulary")
  expect_false("and" %in% vocab$word)

  hg2 <- text_hypergraph(tiny_corpus, min_count = 2L)
  vocab2 <- as.data.frame(hg2, what = "vocabulary")
  expect_identical(vocab2$word, c("and", "soup"))
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
    class = "hypernets_dropped_documents"
  )
  docs <- as.data.frame(hg, what = "documents")
  expect_identical(docs$doc, c("a", "c"))
})

test_that("contract violations raise classed errors", {
  expect_error(
    text_hypergraph(c(a = "and", b = "and"), stop_words = "and"),
    class = "hypernets_empty_corpus"
  )
  expect_error(
    text_hypergraph(data.frame(txt = "salt")),
    class = "hypernets_bad_input"
  )
  expect_error(
    text_hypergraph(data.frame(id = c("a", "a"), txt = c("x", "y")),
                    column = "txt", id = "id"),
    class = "hypernets_bad_input"
  )
})

test_that("construction is deterministic", {
  first <- text_hypergraph(tiny_corpus, weight = "tfidf")
  second <- text_hypergraph(tiny_corpus, weight = "tfidf")
  expect_identical(first, second)
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
  vocab <- as.data.frame(hg, what = "vocabulary")
  expect_identical(vocab$word, c("children's", "plan", "teacher's", "the"))
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
  small_sparse <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"),
                                  sparse = TRUE)
  expect_true(inherits(small_sparse$incidence, "dgCMatrix"))
  # 1,200 documents x 900 words crosses the threshold (alphabetic words:
  # the tokeniser drops digits)
  set.seed(1)
  vocab <- head(apply(expand.grid(letters, letters, letters), 1, paste,
                      collapse = ""), 900)
  docs <- vapply(seq_len(1200), \(i) paste(sample(vocab, 12), collapse = " "),
                 character(1))
  big <- text_hypergraph(docs)
  expect_true(inherits(big$incidence, "dgCMatrix"))
  big_dense <- text_hypergraph(docs, sparse = FALSE)
  expect_false(inherits(big_dense$incidence, "dgCMatrix"))
  expect_error(text_hypergraph(docs, sparse = NA), "sparse")
})

# --- min_chars -------------------------------------------------------------

test_that("text_hypergraph(min_chars) gates the vocabulary", {
  x <- c("M. Wathelet and J. N. Cunha wrote on energy",
         "a coal fired plant emits carbon dioxide")
  vocab <- function(k) {
    hg <- text_hypergraph(x, min_count = 1L, min_chars = k)
    sort(as.data.frame(hg, what = "vocabulary")$word)
  }
  expect_true(all(c("a", "j", "m", "n") %in% vocab(1)))
  expect_false(any(nchar(vocab(3)) < 3))
  expect_true(all(c("energy", "carbon", "wathelet") %in% vocab(3)))
  # the default keeps everything
  expect_identical(vocab(1), sort(as.data.frame(
    text_hypergraph(x, min_count = 1L), what = "vocabulary")$word))
})

test_that("text_hypergraph(min_chars) shrinks the vocabulary monotonically", {
  x <- c("a bo cat food plates schools", "cat food and plates for schools")
  sizes <- vapply(1:6, \(k) nrow(as.data.frame(
    text_hypergraph(x, min_count = 1L, min_chars = k), what = "vocabulary")),
    integer(1))
  expect_false(is.unsorted(rev(sizes)))
})

test_that("text_hypergraph rejects a bad min_chars", {
  expect_error(text_hypergraph("some text here", min_chars = 0))
  expect_error(text_hypergraph("some text here", min_chars = c(2, 3)))
})

test_that("min_chars is refused for construction = 'knn'", {
  expect_error(
    text_hypergraph(c("a b", "c d"), construction = "knn", min_chars = 3L),
    class = "hypernets_bad_input"
  )
})

test_that("the storage rule does not overflow on a large corpus", {
  skip_on_cran()
  # 60000 docs x 70000 words = 4.2e9 cells: the integer product overflows to
  # NA, which `isTRUE()` would read as "dense" and try to allocate ~34 Gb.
  expect_true(hypernets:::.thg_choose_sparse(60000L, 70000L, "bag"))
  expect_true(hypernets:::.thg_choose_sparse(46341L, 46341L, "sentence"))
  expect_false(hypernets:::.thg_choose_sparse(10L, 10L, "bag"))
  expect_false(hypernets:::.thg_choose_sparse(60000L, 70000L, "window"))
  # the decision is never NA, whatever the scale
  sizes <- c(1L, 1000L, 46341L, 60000L, .Machine$integer.max)
  decisions <- vapply(sizes,
                      \(n) hypernets:::.thg_choose_sparse(n, n, "bag"),
                      logical(1))
  expect_false(anyNA(decisions))
})

test_that("the storage rule turns sparse exactly at the threshold", {
  skip_on_cran()
  expect_false(hypernets:::.thg_choose_sparse(1000L, 999L, "bag"))
  expect_true(hypernets:::.thg_choose_sparse(1000L, 1000L, "bag"))
})

test_that("the vocabulary filter is one deterministic ranking", {
  skip_on_cran()
  total <- c(rare = 1L, common = 10L, mid = 5L, tie_b = 5L)
  # min_count alone reproduces the pre-0.4.6 rule exactly
  expect_identical(
    hypernets:::.thg_keep_vocabulary(total, min_count = 5L),
    sort(names(total)[total >= 5L])
  )
  # max_words caps the head; ties (mid, tie_b at 5) break alphabetically
  expect_identical(hypernets:::.thg_keep_vocabulary(total, max_words = 1L),
                   "common")
  expect_identical(hypernets:::.thg_keep_vocabulary(total, max_words = 2L),
                   sort(c("common", "mid")))
  # the whole vocabulary survives the defaults
  expect_identical(hypernets:::.thg_keep_vocabulary(total), sort(names(total)))
  expect_identical(hypernets:::.thg_keep_vocabulary(total, coverage = 1),
                   sort(names(total)))
})

test_that("coverage retains at least the share of tokens it promises", {
  skip_on_cran()
  set.seed(11)
  total <- stats::setNames(as.integer(stats::rpois(400, 8) + 1L),
                           sprintf("w%03d", seq_len(400)))
  shares <- c(0.5, 0.75, 0.9, 0.95, 0.99, 1)
  kept <- lapply(shares, \(s) hypernets:::.thg_keep_vocabulary(total,
                                                               coverage = s))
  retained <- vapply(kept,
                     \(k) sum(total[k]) / sum(total),
                     numeric(1))
  expect_true(all(retained >= shares - sqrt(.Machine$double.eps)))
  # and it is the FEWEST such words: drop the rarest kept word and the share
  # falls below the target (the full vocabulary at coverage = 1 excepted)
  minimal <- vapply(seq_along(shares), \(i) {
    k <- kept[[i]]
    # coverage = 1 must keep everything, so there is no word to drop
    if (isTRUE(all.equal(shares[[i]], 1))) return(length(k) == length(total))
    rarest <- k[[which.min(total[k])]]
    sum(total[setdiff(k, rarest)]) / sum(total) < shares[[i]]
  }, logical(1))
  expect_true(all(minimal))
  # smaller coverage never keeps more words
  expect_false(is.unsorted(vapply(kept, length, integer(1))))
  # and the filters compose: the strictest wins
  both <- hypernets:::.thg_keep_vocabulary(total, coverage = 0.9,
                                           max_words = 10L)
  expect_length(both, 10L)
})

test_that("max_words and coverage prune the constructor's vocabulary", {
  skip_on_cran()
  set.seed(3)
  words <- sprintf("%s%s", rep(letters, each = 26), rep(letters, 26))
  p <- 1 / seq_along(words)
  docs <- vapply(seq_len(120L),
                 \(i) paste(sample(words, 60, replace = TRUE, prob = p),
                            collapse = " "),
                 character(1))
  corpus <- data.frame(id = sprintf("d%03d", seq_len(120L)), text = docs)

  full <- suppressMessages(
    text_hypergraph(corpus, column = "text", id = "id", nodes = "doc"))
  capped <- suppressMessages(
    text_hypergraph(corpus, column = "text", id = "id", nodes = "doc",
                    max_words = 50L))
  covered <- suppressMessages(
    text_hypergraph(corpus, column = "text", id = "id", nodes = "doc",
                    coverage = 0.9))

  expect_identical(capped$n_hyperedges, 50L)
  expect_lt(covered$n_hyperedges, full$n_hyperedges)
  expect_gte(covered$text$token_share, 0.9)
  expect_identical(full$text$token_share, 1)
  # the kept vocabulary is exactly the surviving hyperedges
  expect_setequal(as.data.frame(capped, what = "vocabulary")$word,
                  colnames(capped$incidence))
  # pruning announces itself
  expect_message(
    text_hypergraph(corpus, column = "text", id = "id", nodes = "doc",
                    max_words = 50L),
    "vocabulary pruned to 50"
  )
})

test_that("the vocabulary filters reject bad input", {
  skip_on_cran()
  corpus <- data.frame(id = c("a", "b"),
                       text = c("night owl sings", "night crow calls"))
  expect_error(
    text_hypergraph(corpus, column = "text", id = "id", coverage = 0),
    "`coverage` must be a single share"
  )
  expect_error(
    text_hypergraph(corpus, column = "text", id = "id", coverage = 1.5),
    "`coverage` must be a single share"
  )
  expect_error(
    text_hypergraph(corpus, column = "text", id = "id", max_words = 0),
    "`max_words` must be a single count"
  )
  expect_error(
    text_hypergraph(corpus, column = "text", id = "id", min_count = 99L),
    class = "hypernets_empty_corpus"
  )
})
