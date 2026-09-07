# ---- window_hypergraph() tests -------------------------------------------

test_that("hand-computed case is exact (sets, occurrence incidence, counts)", {
  # windows of c(a b a c), w = 2, step = 1: (a,b) (b,a) (a,c)
  # -> {a,b} from 2 windows, {a,c} from 1
  hg <- window_hypergraph(list(s1 = c("a", "b", "a", "c")), window = 2L)
  expect_s3_class(hg, "net_hypergraph")
  expect_identical(hg$nodes, c("a", "b", "c"))
  expect_identical(hg$hyperedges, list(c(1L, 2L), c(1L, 3L)))
  expect_identical(hg$window_counts, c(2L, 1L))
  expected_inc <- matrix(c(2, 2, 0, 1, 0, 1), nrow = 3,
                         dimnames = list(c("a", "b", "c"), c("h1", "h2")))
  expect_identical(hg$incidence, expected_inc)
  expect_identical(hg$n_hyperedges, 2L)
  expect_identical(hg$size_distribution, c(size_2 = 2L))
})

test_that("within-window repeats count occurrences, not presence", {
  # window (a,a) -> set {a} with incidence cell 2
  hg <- window_hypergraph(list(c("a", "a")), window = 2L)
  expect_identical(hg$hyperedges, list(1L))
  expect_identical(as.vector(hg$incidence), 2)
  expect_identical(hg$window_counts, 1L)
})

test_that("window-count and occurrence conservation hold across random data", {
  set.seed(303)
  results <- lapply(1:20, function(i) {
    n <- sample(4:30, 1)
    w <- sample(2:4, 1)
    s <- sample(1:3, 1)
    traj <- sample(letters[1:5], n, replace = TRUE)
    hg <- window_hypergraph(list(traj), window = w, step = s)
    n_windows <- length(seq.int(1L, n - w + 1L, by = s))
    # No NAs: every full window is non-empty and no hyperedge is dropped
    expect_identical(sum(hg$window_counts), n_windows)
    expect_identical(hg$params$n_windows, n_windows)
    expect_identical(sum(hg$incidence), as.numeric(n_windows * w))
    expect_identical(hg$params$n_empty_windows, 0L)
    expect_identical(hg$params$n_dropped, 0L)
    hg
  })
  expect_length(results, 20L)
})

test_that("w = 2 sliding reduces to adjacent pairwise co-occurrence counts", {
  set.seed(99)
  traj <- sample(c("a", "b", "c", "d"), 40, replace = TRUE)
  hg <- window_hypergraph(list(traj), window = 2L)
  # Independent count of unordered adjacent pairs
  pairs <- paste(pmin(traj[-length(traj)], traj[-1L]),
                 pmax(traj[-length(traj)], traj[-1L]), sep = "|")
  expected <- table(pairs)
  got <- stats::setNames(
    hg$window_counts,
    vapply(hg$hyperedges,
           function(idx) paste(hg$nodes[idx][c(1L, length(idx))],
                               collapse = "|"),
           character(1L))
  )
  expect_identical(as.integer(expected[names(got)]), as.integer(got))
  expect_identical(sum(hg$window_counts), length(traj) - 1L)
})

test_that("step = window gives tumbling windows", {
  hg <- window_hypergraph(list(c("a", "b", "c", "d", "e", "f")),
                          window = 2L, step = 2L)
  # windows (a,b) (c,d) (e,f): three disjoint size-2 hyperedges, weight 1
  expect_identical(hg$n_hyperedges, 3L)
  expect_identical(hg$window_counts, c(1L, 1L, 1L))
  expect_identical(sort(unlist(hg$hyperedges)), 1:6)
  # trailing partial window is never formed
  hg2 <- window_hypergraph(list(c("a", "b", "c", "d", "e")),
                           window = 2L, step = 2L)
  expect_identical(hg2$params$n_windows, 2L)
})

test_that("long and wide input give the same hypergraph", {
  set.seed(17)
  seqs <- lapply(1:6, function(i) sample(letters[1:4], 8, replace = TRUE))
  wide <- as.data.frame(do.call(rbind, seqs), stringsAsFactors = FALSE)
  long <- data.frame(
    act = unlist(seqs),
    id  = rep(paste0("s", 1:6), each = 8L),
    t   = rep(1:8, times = 6L),
    stringsAsFactors = FALSE
  )
  hg_w <- window_hypergraph(wide, window = 3L)
  hg_l <- window_hypergraph(long, action = "act", actor = "id", time = "t",
                            window = 3L)
  expect_identical(hg_w$hyperedges,    hg_l$hyperedges)
  expect_identical(hg_w$incidence,     hg_l$incidence)
  expect_identical(hg_w$window_counts, hg_l$window_counts)
  expect_identical(hg_w$nodes,         hg_l$nodes)
})

test_that("long format with time = NULL preserves row order within actor", {
  long <- data.frame(act = c("a", "b", "c", "x", "y", "z"),
                     id  = c("s1", "s1", "s1", "s2", "s2", "s2"),
                     stringsAsFactors = FALSE)
  hg <- window_hypergraph(long, action = "act", actor = "id", window = 3L)
  expect_identical(
    sort(vapply(hg$hyperedges,
                function(idx) paste(hg$nodes[idx], collapse = ""),
                character(1L))),
    c("abc", "xyz")
  )
})

test_that("state relabeling permutes but preserves structure", {
  set.seed(41)
  traj <- sample(letters[1:4], 25, replace = TRUE)
  relabel <- c(a = "w", b = "x", c = "y", d = "z")
  hg1 <- window_hypergraph(list(traj), window = 3L)
  hg2 <- window_hypergraph(list(unname(relabel[traj])), window = 3L)
  expect_identical(sort(lengths(hg1$hyperedges)), sort(lengths(hg2$hyperedges)))
  expect_identical(sort(hg1$window_counts), sort(hg2$window_counts))
  expect_identical(unname(relabel[hg1$nodes]), hg2$nodes)
  # relabeling here is order-preserving (a<b<c<d -> w<x<y<z), so the full
  # object must match up to dimnames
  expect_identical(unname(hg1$incidence), unname(hg2$incidence))
})

test_that("NA states are excluded; all-NA windows are counted, not silent", {
  hg <- window_hypergraph(list(c("a", NA, NA, "b")), window = 2L)
  # windows (a,NA) (NA,NA) (NA,b) -> {a}, empty, {b}
  expect_identical(hg$params$n_windows, 3L)
  expect_identical(hg$params$n_empty_windows, 1L)
  expect_identical(sum(hg$window_counts), 2L)
  expect_identical(hg$nodes, c("a", "b"))
})

test_that("wide input strips trailing NAs but keeps internal ones", {
  wide <- data.frame(t1 = c("a", "x"), t2 = c(NA, "y"),
                     t3 = c("b", NA), t4 = c("c", NA),
                     stringsAsFactors = FALSE)
  hg <- window_hypergraph(wide, window = 2L)
  # row 1: a NA b c -> windows {a} {b} {b,c}; row 2: x y -> {x,y}
  expect_identical(hg$params$n_windows, 4L)
  expect_identical(hg$params$n_short_sequences, 0L)
  expect_identical(sum(hg$window_counts), 4L)
})

test_that("sequences shorter than the window are counted, not silent", {
  hg <- window_hypergraph(list(c("a", "b", "c"), c("z")), window = 2L)
  expect_identical(hg$params$n_sequences, 2L)
  expect_identical(hg$params$n_short_sequences, 1L)
  expect_false("z" %in% hg$nodes)
})

test_that("min_size drops small hyperedges and records the count", {
  hg <- window_hypergraph(list(c("a", "a", "b", "c")), window = 2L,
                          min_size = 2L)
  # windows (a,a) (a,b) (b,c): {a} dropped
  expect_identical(hg$params$n_dropped, 1L)
  expect_identical(hg$n_hyperedges, 2L)
  expect_true(all(lengths(hg$hyperedges) >= 2L))
  expect_identical(colnames(hg$incidence), c("h1", "h2"))
})

test_that("error paths: invalid arguments and empty results", {
  seqs <- list(c("a", "b", "c"))
  expect_error(window_hypergraph(seqs, window = 1L), "window")
  expect_error(window_hypergraph(seqs, window = 2.5), "window")
  expect_error(window_hypergraph(seqs, step = 0L), "step")
  expect_error(window_hypergraph(seqs, min_size = 0L), "min_size")
  expect_error(window_hypergraph(seqs, actor = "id"), "actor")
  expect_error(window_hypergraph(data.frame(a = "x"), action = "missing"),
               "action")
  expect_error(window_hypergraph(list("a"), window = 2L), "shorter")
  expect_error(window_hypergraph(list(c(NA_character_, NA_character_)),
                                 window = 2L), "NA")
  expect_error(window_hypergraph(list(c("a", "a", "a")), window = 2L,
                                 min_size = 2L), "min_size")
  expect_error(window_hypergraph(1:5), "data")
})

test_that("as.data.frame accessor is tidy for all constructors", {
  hg <- window_hypergraph(list(c("a", "b", "a", "c")), window = 2L)
  df <- as.data.frame(hg)
  expect_identical(names(df), c("hyperedge", "size", "states", "weight"))
  expect_identical(nrow(df), hg$n_hyperedges)
  expect_identical(df$weight, c(2, 1))
  expect_identical(df$states, c("a, b", "a, c"))
  # clique constructor: unweighted hyperedges -> NA weights
  adj <- matrix(1, 3, 3, dimnames = list(letters[1:3], letters[1:3]))
  diag(adj) <- 0
  hg_c <- build_hypergraph(adj, p = 1, max_size = 3L)
  df_c <- as.data.frame(hg_c)
  expect_identical(nrow(df_c), hg_c$n_hyperedges)
  expect_true(all(is.na(df_c$weight)))
})

test_that("print and summary run and report the windowed source", {
  hg <- window_hypergraph(list(c("a", "b", "a", "c")), window = 2L)
  expect_invisible(print(hg))
  expect_output(print(hg), "windowed sequences, window = 2, step = 1")
  expect_output(summary(hg), "Hypergraph summary")
  # bipartite print line no longer vanishes either
  bp <- group_hypergraph(
    data.frame(p = c("a", "b", "a"), g = c("g1", "g1", "g2")),
    actor = "p", group = "g")
  expect_output(print(bp), "group membership")
})

test_that("Laplacian family defaults to window counts as hyperedge weights", {
  set.seed(11)
  traj <- c(rep(c("a", "b"), 6), "c", rep(c("a", "c"), 3))
  hg <- window_hypergraph(list(traj), window = 2L)
  expect_true(length(unique(hg$window_counts)) > 1L)  # non-trivial weights
  for (ty in c("zhou", "random_walk")) {
    # for loop kept: two assertions over a 2-level argument, no data to grow
    default_laplacian <- hypergraph_laplacian(hg, type = ty)
    counts_laplacian <- hypergraph_laplacian(
      hg, type = ty, edge_weights = as.numeric(hg$window_counts))
    expect_identical(default_laplacian, counts_laplacian)
  }
  zhou_default <- hypergraph_laplacian(hg)
  zhou_unit <- hypergraph_laplacian(hg, edge_weights = rep(1, hg$n_hyperedges))
  expect_false(isTRUE(all.equal(
    zhou_default,
    zhou_unit,
    check.attributes = FALSE
  )))
  cl <- hypergraph_cluster(hg, k = 2L, seed = 5)
  expect_s3_class(cl, "net_hypergraph_cluster")
})

test_that("downstream verbs consume windowed hypergraphs", {
  set.seed(23)
  seqs <- lapply(1:8, function(i) sample(letters[1:5], 10, replace = TRUE))
  hg <- window_hypergraph(seqs, window = 3L)
  m <- hypergraph_measures(hg)
  expect_s3_class(m, "net_hypergraph_measures")
  ce <- hypergraph_centrality(hg, type = "clique")
  expect_s3_class(ce, "data.frame")
  expect_identical(nrow(ce), hg$n_nodes)
  expect_identical(ce$node, hg$nodes)
  net <- clique_expansion(hg)
  expect_s3_class(net, "netobject")
})

test_that("works on the bundled human_long dataset (long format)", {
  data("human_long", package = "honets")
  hg <- window_hypergraph(human_long, action = "code", actor = "session_id",
                          time = "timestamp", window = 3L)
  expect_s3_class(hg, "net_hypergraph")
  expect_identical(hg$params$n_sequences, length(unique(human_long$session_id)))
  expect_identical(sum(hg$window_counts) + hg$params$n_empty_windows,
                   hg$params$n_windows)
  expect_true(all(hg$nodes %in% unique(human_long$code)))
})

test_that("as.data.frame sort_by orders deterministically, largest first", {
  hg <- window_hypergraph(list(c("a", "b", "a", "b", "a", "c")), window = 2L)
  df <- as.data.frame(hg, sort_by = "weight")
  expect_true(all(diff(df$weight) <= 0))
  expect_identical(df$states[1], "a, b")
  expect_identical(nrow(df), hg$n_hyperedges)
  expect_error(as.data.frame(hg, sort_by = "nope"), "arg")
})

test_that("min_weight keeps only recurrent hyperedges", {
  # sets: {a,b} x3, {a,c} x1
  hg <- window_hypergraph(list(c("a", "b", "a", "b", "a", "c")), window = 2L,
                          min_weight = 2L)
  expect_identical(hg$n_hyperedges, 1L)
  hyperedges <- as.data.frame(hg)
  expect_identical(hyperedges$states, "a, b")
  expect_identical(hg$params$n_dropped, 1L)
  expect_error(window_hypergraph(list(c("a", "b")), min_weight = 0L),
               "min_weight")
  expect_error(window_hypergraph(list(c("a", "b", "c")), window = 2L,
                                 min_weight = 5L), "min_weight")
})
