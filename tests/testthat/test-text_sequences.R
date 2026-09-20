testthat::skip_on_cran()

# ===========================================================================
# hg_sequences(): the text family -> memory family bridge.
#
# The reference for the round trip is the join a caller had to write by hand
# before this verb existed. It is kept inside a function body on purpose:
# merge / order / split on a result is exactly the ritual the verb removes
# from the public surface.
# ===========================================================================

seq_posts <- data.frame(
  student = rep(c("ana", "bo", "cy"), each = 4L),
  turn = rep(1:4, times = 3L),
  phase = rep(c("early", "early", "late", "late"), times = 3L),
  text = c(
    "simmer the soup with onions and carrots",
    "this soup recipe needs salt on a cold night",
    "the telescope revealed a distant galaxy and stars",
    "astronomers aimed the telescope at the stars all night",
    "a galaxy of stars seen through the telescope",
    "salt and onions make the soup",
    "the soup simmered all night with carrots",
    "astronomers watched the distant galaxy",
    "onions carrots and salt in the soup",
    "the distant stars and the galaxy at night",
    "soup with salt on a cold night",
    "astronomers and the telescope watched stars"
  ),
  stringsAsFactors = FALSE
)
seq_stop <- c("the", "with", "and", "a", "this", "at", "on", "all", "of",
              "through", "in", "make")

seq_hg <- function() {
  text_hypergraph(seq_posts, column = "text", stop_words = seq_stop)
}

# The hand-written join the verb replaces: merge the documents table onto the
# cluster assignments, sort by actor then turn, split the state by actor.
hand_written_join <- function(hg, clusters, actor, order_by) {
  docs <- as.data.frame(hg, what = "documents")
  joined <- merge(docs, clusters, by.x = "doc", by.y = "node")
  joined <- joined[order(joined[[actor]], joined[[order_by]]), ]
  lapply(split(joined$cluster, joined[[actor]]), as.character)
}

# The same reshape applied to the verb's output, so the two can be compared.
as_sequence_list <- function(x) {
  lapply(split(x$action, x$actor), as.character)
}

test_that("hg_sequences() reproduces the hand-written join exactly", {
  hg <- seq_hg()
  topics <- hg_cluster(hg, k = 2L, seed = 1L)

  seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
  expect_identical(
    as_sequence_list(seqs),
    hand_written_join(hg, topics, "student", "turn")
  )
})

test_that("hg_sequences() returns the canonical long sequence table", {
  hg <- seq_hg()
  topics <- hg_cluster(hg, k = 2L, seed = 1L)
  seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")

  expect_s3_class(seqs, "data.frame")
  expect_identical(class(seqs), "data.frame")
  expect_identical(names(seqs), c("actor", "time", "action"))
  expect_type(seqs$actor, "character")
  expect_type(seqs$action, "character")
  # `time` keeps the type of the ordering column.
  expect_identical(seqs$time, rep(1:4, times = 3L))
  expect_identical(rownames(seqs), as.character(seq_len(nrow(seqs))))
})

test_that("hg_sequences() preserves every document and orders within actor", {
  hg <- seq_hg()
  topics <- hg_cluster(hg, k = 2L, seed = 1L)
  docs <- as.data.frame(hg, what = "documents")

  # Invariant 1: one row per document, none added, none dropped.
  seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
  expect_identical(nrow(seqs), nrow(docs))
  expect_identical(as.vector(table(seqs$actor)), rep(4L, 3L))
  expect_identical(sort(seqs$action), sort(topics$cluster))

  # Invariant 2: strictly increasing time within every actor, and the actors
  # come out in blocks rather than interleaved.
  expect_true(all(vapply(
    split(seqs$time, seqs$actor),
    function(t) !is.unsorted(t, strictly = TRUE),
    logical(1L)
  )))
  expect_identical(seqs$actor, sort(seqs$actor))

  # Invariant 3: the order of the input rows does not matter.
  shuffled <- seq_posts[c(7L, 2L, 11L, 4L, 9L, 1L, 6L, 12L, 3L, 8L, 5L, 10L), ]
  hg_shuffled <- text_hypergraph(shuffled, column = "text",
                                 stop_words = seq_stop)
  states <- data.frame(
    node = as.data.frame(hg_shuffled, what = "documents")$doc,
    cluster = as.data.frame(hg_shuffled, what = "documents")$phase,
    stringsAsFactors = FALSE
  )
  expect_identical(
    hg_sequences(hg_shuffled, states, actor = "student", order_by = "turn"),
    hg_sequences(hg, actor = "student", order_by = "turn", state = "phase")
  )
})

test_that("hg_sequences() takes a state already in the documents table", {
  hg <- seq_hg()
  seqs <- hg_sequences(hg, actor = "student", order_by = "turn",
                       state = "phase")

  expect_identical(names(seqs), c("actor", "time", "action"))
  expect_identical(seqs$action, rep(c("early", "early", "late", "late"),
                                    times = 3L))
})

test_that("hg_sequences() feeds build_hon() with no further coercion", {
  hg <- seq_hg()
  topics <- hg_cluster(hg, k = 2L, seed = 1L)
  seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")

  hon <- build_hon(seqs, action = "action", actor = "actor", time = "time",
                   max_order = 2L)
  expect_s3_class(hon, "net_hon")
  expect_identical(hon$n_trajectories, 3L)

  # The network is the one the hand-assembled list of sequences would give.
  reference <- build_hon(hand_written_join(hg, topics, "student", "turn"),
                         max_order = 2L)
  expect_identical(hon$matrix, reference$matrix)
  expect_identical(hon$edges, reference$edges)

  # The other verbs that take the long form read the same table: the
  # bootstrap, and the sequence -> hypergraph bridge.
  boot <- bootstrap_hon(seqs, action = "action", actor = "actor",
                        time = "time", n_boot = 10L, max_order = 2L,
                        seed = 1L)
  expect_s3_class(boot, "net_hon_boot")
  wh <- window_hypergraph(seqs, window = 2L, action = "action",
                          actor = "actor", time = "time")
  expect_s3_class(wh, "net_hypergraph")
})

test_that("hg_sequences() raises classed conditions on bad input", {
  hg <- seq_hg()
  topics <- hg_cluster(hg, k = 2L, seed = 1L)

  expect_error(
    hg_sequences(hg, topics, actor = "teacher", order_by = "turn"),
    class = "hypernets_bad_input"
  )
  expect_error(
    hg_sequences(hg, topics, actor = "student", order_by = "when"),
    class = "hypernets_bad_input"
  )
  # No clusters and no state: nothing says what the states are.
  expect_error(
    hg_sequences(hg, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )
  expect_error(
    hg_sequences(hg, actor = "student", order_by = "turn", state = "topic"),
    class = "hypernets_bad_input"
  )
  # A partition of something other than the documents (here, the words).
  words <- hg_cluster(text_hypergraph(seq_posts, column = "text",
                                      nodes = "word", stop_words = seq_stop),
                      k = 2L, seed = 1L)
  expect_error(
    hg_sequences(hg, words, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )
  # Not a hypergraph at all, and a hypergraph with no documents table.
  expect_error(
    hg_sequences(seq_posts, topics, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )
  expect_error(
    hg_sequences(
      group_hypergraph(data.frame(member = c("a", "b"), group = c("g", "g"))),
      topics, actor = "student", order_by = "turn"
    ),
    class = "hypernets_bad_input"
  )
})

test_that("hg_sequences() refuses missing and ambiguous input", {
  hg <- seq_hg()
  docs <- as.data.frame(hg, what = "documents")
  topics <- hg_cluster(hg, k = 2L, seed = 1L)

  gappy <- seq_posts
  gappy$student[2L] <- NA_character_
  hg_gappy <- text_hypergraph(gappy, column = "text", stop_words = seq_stop)
  expect_error(
    hg_sequences(hg_gappy, topics, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )

  duplicated_nodes <- rbind(topics, topics)
  expect_error(
    hg_sequences(hg, duplicated_nodes, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )

  partial <- data.frame(node = docs$doc, cluster = docs$phase,
                        stringsAsFactors = FALSE)
  partial$cluster[3L] <- NA_character_
  expect_error(
    hg_sequences(hg, partial, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )

  expect_error(
    hg_sequences(hg, topics, actor = c("student", "turn"), order_by = "turn"),
    "single column name"
  )
})

test_that("hg_sequences() names the state column with `state`", {
  hg <- seq_hg()
  docs <- as.data.frame(hg, what = "documents")
  states <- data.frame(node = docs$doc, topic = docs$phase,
                       stringsAsFactors = FALSE)

  expect_identical(
    hg_sequences(hg, states, actor = "student", order_by = "turn",
                 state = "topic"),
    hg_sequences(hg, actor = "student", order_by = "turn", state = "phase")
  )
  expect_error(
    hg_sequences(hg, states, actor = "student", order_by = "turn"),
    class = "hypernets_bad_input"
  )
})

test_that("hg_sequences() sorts calendar times and breaks ties by document", {
  hg <- seq_hg()
  docs <- as.data.frame(hg, what = "documents")
  dated <- seq_posts
  dated$when <- as.Date("2026-01-01") + rep(c(3L, 1L, 4L, 2L), times = 3L)
  hg_dated <- text_hypergraph(dated, column = "text", stop_words = seq_stop)
  states <- data.frame(node = docs$doc, cluster = docs$phase,
                       stringsAsFactors = FALSE)

  # Turns 1..4 are phases early, early, late, late and fall on days 3, 1, 4, 2:
  # sorted by date the turns come out 2, 4, 1, 3 -- early, late, early, late.
  seqs <- hg_sequences(hg_dated, states, actor = "student", order_by = "when")
  expect_s3_class(seqs$time, "Date")
  expect_identical(seqs$action,
                   rep(c("early", "late", "early", "late"), times = 3L))

  # Every document at the same time: the tie-break is the document id, which
  # does not depend on the order the rows arrived in.
  flat <- seq_posts
  flat$turn <- 1L
  hg_flat <- text_hypergraph(flat, column = "text", stop_words = seq_stop)
  flat_states <- data.frame(
    node = as.data.frame(hg_flat, what = "documents")$doc,
    cluster = as.data.frame(hg_flat, what = "documents")$phase,
    stringsAsFactors = FALSE
  )
  expect_identical(
    hg_sequences(hg_flat, flat_states, actor = "student", order_by = "turn"),
    hg_sequences(hg_flat, flat_states, actor = "student", order_by = "turn")
  )
})
