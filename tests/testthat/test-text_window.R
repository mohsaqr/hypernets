# Windowed construction: hand-computed windows, conservation invariant, and
# the reduction identity to window_hypergraph().

test_that("sliding windows are hand-computed correctly", {
  hg <- text_hypergraph(c(d = "a b c a b"), construction = "window",
                        window = 2)
  expect_identical(
    as.data.frame(hg),
    data.frame(
      edge = c("a+b", "a+b", "a+c", "a+c", "b+c", "b+c"),
      word = c("a", "b", "a", "c", "b", "c"),
      weight = c(2, 2, 1, 1, 1, 1)
    )
  )
  expect_identical(hg$text$n_windows, 4L)
})

test_that("tumbling windows include the trailing partial chunk", {
  hg <- text_hypergraph(c(d = "a b c a b"), construction = "window",
                        window = 2, window_mode = "tumbling")
  expect_identical(
    as.data.frame(hg),
    data.frame(
      edge = c("a+b", "a+b", "a+c", "a+c", "b"),
      word = c("a", "b", "a", "c", "b"),
      weight = c(1, 1, 1, 1, 1)
    )
  )
  expect_identical(hg$text$n_windows, 3L)
})

test_that("repeated tokens inside a window collapse to a singleton edge", {
  hg <- text_hypergraph(c(d = "a a b"), construction = "window", window = 2)
  expect_identical(
    as.data.frame(hg),
    data.frame(
      edge = c("a", "a+b", "a+b"),
      word = c("a", "a", "b"),
      weight = c(1, 1, 1)
    )
  )
})

test_that("a document shorter than the window forms one whole-document window", {
  hg <- text_hypergraph(c(d = "a b"), construction = "window", window = 5)
  expect_identical(
    as.data.frame(hg),
    data.frame(edge = c("a+b", "a+b"), word = c("a", "b"), weight = c(1, 1))
  )
  expect_identical(hg$text$n_windows, 1L)
})

test_that("windows never cross document boundaries", {
  hg <- text_hypergraph(c(x = "a b", y = "c d"), construction = "window",
                        window = 2)
  expect_identical(sort(unique(as.data.frame(hg)$edge)), c("a+b", "c+d"))
})

test_that("window counts are conserved", {
  corpus <- c(d1 = "a b c a b c a", d2 = "b c b c b")
  hg <- text_hypergraph(corpus, construction = "window", window = 3)
  # sliding, full windows: (7 - 3 + 1) + (5 - 3 + 1) = 8
  tab <- as.data.frame(hg)
  per_edge <- aggregate(weight ~ edge, data = tab, FUN = max)
  expect_identical(sum(per_edge$weight), 8)
  expect_identical(hg$text$n_windows, 8L)
  # every member of an edge carries the same weight (the window count)
  spread <- aggregate(weight ~ edge, data = tab, FUN = \(w) diff(range(w)))
  expect_true(all(spread$weight == 0))
})

test_that("min_count filtering closes the gap before windowing", {
  hg <- text_hypergraph(c(d = "a q b a b"), construction = "window",
                        window = 2, min_count = 2L)
  tab <- as.data.frame(hg)
  expect_false("q" %in% tab$word)
  # filtered sequence a b a b -> windows (a,b),(b,a),(a,b) -> a+b weight 3
  expect_identical(
    tab,
    data.frame(edge = c("a+b", "a+b"), word = c("a", "b"),
               weight = c(3, 3))
  )
})

# The reduction identity: text_hypergraph(construction = "window") is the
# text front end of window_hypergraph(). On token sequences with no repeated
# word inside a window, and with every document long enough to fill a window
# (tumbling: lengths that are multiples of `window`), the two constructions
# produce the same incidence matrix. The hyperedge NAMES differ by design:
# window_hypergraph() numbers them (h1, h2, ...), text_hypergraph() names them
# by their word content.
.unnamed_incidence <- function(hg) {
  m <- as.matrix(hg$incidence)
  m <- m[order(rownames(m)), , drop = FALSE]
  unname(m[, order(apply(m, 2, paste, collapse = "\x01")), drop = FALSE])
}

test_that("window construction equals window_hypergraph() on its shared domain", {
  corpus <- c(d1 = "salt soup onion carrot pepper",
              d2 = "soup stars sky night",
              d3 = "stars galaxy telescope night sky")
  tokens <- strsplit(corpus, " ")
  for (w in 2:3) {
    text <- text_hypergraph(corpus, construction = "window", window = w)
    engine <- window_hypergraph(tokens, window = w, step = 1L)
    expect_identical(text$nodes, engine$nodes)
    expect_identical(text$n_hyperedges, engine$n_hyperedges)
    expect_identical(.unnamed_incidence(text), .unnamed_incidence(engine))
  }
  tumbling_corpus <- c(d1 = "salt soup onion carrot", d2 = "soup stars sky night")
  text <- text_hypergraph(tumbling_corpus, construction = "window",
                          window = 2, window_mode = "tumbling")
  engine <- window_hypergraph(strsplit(tumbling_corpus, " "), window = 2,
                              step = 2L)
  expect_identical(.unnamed_incidence(text), .unnamed_incidence(engine))
  # and the two constructions differ exactly where the definitions do: a
  # repeated word inside a window is collapsed here, counted by the engine
  rep_corpus <- c(d = "a a b")
  text <- text_hypergraph(rep_corpus, construction = "window", window = 2)
  engine <- window_hypergraph(list(c("a", "a", "b")), window = 2, step = 1L)
  expect_identical(text$n_hyperedges, 2L)
  expect_identical(engine$n_hyperedges, 2L)
  expect_identical(max(text$incidence), 1)
  expect_identical(max(engine$incidence), 2)
})

test_that("windowed hypergraphs feed the verbs and print distinctly", {
  hg <- text_hypergraph(c(d1 = "a b c a b", d2 = "b c d"),
                        construction = "window", window = 2)
  nodes <- hg_measures(hg, what = "nodes")
  expect_identical(nrow(nodes), hg$n_nodes)
  expect_output(print(hg), "windowed hyperedges: w = 2, sliding")
})

test_that("window construction rejects tf-idf and window < 2", {
  expect_error(
    text_hypergraph(c(d = "a b c"), construction = "window",
                    weight = "tfidf"),
    class = "hypernets_bad_input"
  )
  expect_error(
    text_hypergraph(c(d = "a b c"), construction = "window", window = 1)
  )
})
