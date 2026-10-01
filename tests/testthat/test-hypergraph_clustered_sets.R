# ---- group_hypergraph() on Nestimate clusterings ----------------------------
# The objects are built by structure (as Nestimate returns them), so these
# tests need no Nestimate. Loading Nestimate re-registers S3 methods of
# hypernets classes for the rest of the session, so the test on real
# Nestimate objects lives in local_testing_and_equivalence/
# (test-equiv-clustered-sets-eventdata.R).

.cs_data <- function() {
  data.frame(
    V1 = c("a", "a", "b", "a", "c", "b", "a"),
    V2 = c("b", "b", "c", "c", "a", "c", "b"),
    V3 = c("a", NA, "a", "b", "b", "b", "b"),
    stringsAsFactors = FALSE
  )
}
.cs_assign <- c(1L, 1L, 2L, 1L, 2L, 2L, 1L)

.cs_mmm <- function() {
  structure(list(data = .cs_data(), assignments = .cs_assign, k = 2L,
                 posterior = cbind(as.numeric(.cs_assign == 1L),
                                   as.numeric(.cs_assign == 2L))),
            class = "net_mmm")
}
.cs_clustering <- function() {
  structure(list(data = .cs_data(), k = 2L,
                 assignments = stats::setNames(.cs_assign, seq_along(.cs_assign))),
            class = "net_clustering")
}
.cs_group <- function(names = c("Cluster 1", "Cluster 2")) {
  d <- .cs_data()
  states <- c("a", "b", "c")
  nets <- lapply(1:2, function(g) {
    structure(list(data = d[.cs_assign == g, , drop = FALSE],
                   weights = matrix(0, 3, 3, dimnames = list(states, states))),
              class = c("netobject", "cograph_network"))
  })
  names(nets) <- names
  structure(nets, class = "netobject_group")
}

test_that("each clustering class is coerced to sequences per group", {
  for (x in list(.cs_mmm(), .cs_clustering(), .cs_group())) {
    grouped <- .coerce_grouped_sequences(x)
    expect_identical(names(grouped), c("Cluster 1", "Cluster 2"))
    expect_identical(lengths(grouped), c(`Cluster 1` = 4L, `Cluster 2` = 3L))
    expect_identical(grouped[["Cluster 1"]][[2L]], c("a", "b"))  # NA dropped
    expect_identical(grouped[["Cluster 2"]][[1L]], c("b", "c", "a"))
  }
})

test_that("a netobject_group keeps its (renamed) group names", {
  grouped <- .coerce_grouped_sequences(.cs_group(c("Short", "Long")))
  expect_identical(names(grouped), c("Short", "Long"))
  hg <- group_hypergraph(.cs_group(c("Short", "Long")))
  expect_setequal(unique(hg_get(hg, what = "sets")$group), c("Short", "Long"))
})

test_that("integer-coded states are decoded to labels", {
  coded <- .cs_mmm()
  coded$data[] <- lapply(coded$data, function(col) match(col, c("a", "b", "c")))
  coded$models <- list(list(weights = matrix(0, 3, 3,
    dimnames = list(c("a", "b", "c"), c("a", "b", "c")))))
  expect_identical(.coerce_grouped_sequences(coded), .coerce_grouped_sequences(.cs_mmm()))
})

test_that("top sets and counts match a hand count", {
  # Cluster 1: {a,b} {a,b} {a,b,c} {a,b}  -> a+b: 3, a+b+c: 1
  # Cluster 2: {a,b,c} {a,b,c} {b,c}      -> a+b+c: 2, b+c: 1
  sets <- hg_get(group_hypergraph(.cs_mmm()), what = "sets")
  expect_identical(sets$group, c("Cluster 1", "Cluster 1", "Cluster 2", "Cluster 2"))
  expect_identical(sets$set, c("a + b", "a + b + c", "a + b + c", "b + c"))
  expect_identical(sets$count, c(3L, 1L, 2L, 1L))
  expect_identical(sets$size, c(2L, 3L, 3L, 2L))
  expect_equal(sets$share, c(3 / 4, 1 / 4, 2 / 3, 1 / 3))
  states <- hg_get(group_hypergraph(.cs_mmm()), what = "state_counts")
  expect_identical(states$count[states$group == "Cluster 1"], c(4L, 4L, 1L))
  expect_identical(states$node[states$group == "Cluster 1"], c("a", "b", "c"))
})

test_that("top keeps the most frequent sets; ties go by set name", {
  sets <- hg_get(group_hypergraph(.cs_mmm(), top = 1), what = "sets")
  expect_identical(sets$set, c("a + b", "a + b + c"))
  tie <- structure(list(data = data.frame(V1 = c("c", "a", "b"), V2 = c("d", "b", "c")),
                        assignments = c(1L, 1L, 1L), k = 1L), class = "net_mmm")
  expect_identical(hg_get(group_hypergraph(tie, top = 2), what = "sets")$set,
                   c("a + b", "b + c"))
})

test_that("states keeps a subset of states before counting", {
  sets <- hg_get(group_hypergraph(.cs_mmm(), states = c("a", "c")), what = "sets")
  # Cluster 1: {a} x3, {a,c}; Cluster 2: {a,c} {a,c} {c}
  expect_identical(sets$set, c("a", "a + c", "a + c", "c"))
  expect_identical(sets$count, c(3L, 1L, 2L, 1L))
})

test_that("the three classes give the same hypergraph", {
  a <- group_hypergraph(.cs_mmm())
  b <- group_hypergraph(.cs_clustering())
  c <- group_hypergraph(.cs_group())
  expect_identical(a$incidence, b$incidence)
  expect_identical(a$incidence, c$incidence)
  expect_identical(a$edge_data, c$edge_data)
  expect_identical(hg_get(a, what = "sets"), hg_get(c, what = "sets"))
})

test_that("invariant: counts never exceed group sizes, shares sum to at most 1", {
  set.seed(7)
  d <- as.data.frame(matrix(sample(letters[1:4], 200, TRUE), 40, 5))
  x <- structure(list(data = d, assignments = sample(1:3, 40, TRUE), k = 3L),
                 class = "net_mmm")
  sets <- hg_get(group_hypergraph(x, top = 50), what = "sets")
  totals <- tapply(sets$count, sets$group, sum)
  expect_identical(as.integer(totals), as.integer(tabulate(x$assignments, 3L)))
  expect_true(all(abs(tapply(sets$share, sets$group, sum) - 1) < 1e-12))
  # permuting the sequences leaves the result unchanged
  perm <- sample(40)
  y <- x
  y$data <- d[perm, ]
  y$assignments <- x$assignments[perm]
  expect_identical(hg_get(group_hypergraph(y, top = 50), what = "sets"), sets)
})

test_that("plot(group =) draws one group's sets, sized and titled", {
  hg <- group_hypergraph(.cs_mmm())
  p <- plot(hg, group = "Cluster 1")
  expect_s3_class(p, "ggplot")
  expect_identical(p$labels$title, "Cluster 1 (4 sequences)")
  expect_identical(p$labels$subtitle, "Shared by every set: a, b")
  # one group: no title boxes by default, legends stacked
  geoms <- vapply(p$layers, \(l) class(l$geom)[1L], character(1L))
  expect_false("GeomLabel" %in% geoms)
  expect_true("GeomLabel" %in% vapply(plot(hg, group = "Cluster 1", titles = "count")$layers,
                                      \(l) class(l$geom)[1L], character(1L)))
  expect_identical(p$theme$legend.box, "vertical")
  # the group's hyperedges are named by what they add, in name order
  expect_identical(sort(colnames(hg_subset(hg, where = list(group = "Cluster 1"))$incidence)),
                   c("Cluster 1: a + b", "Cluster 1: a + b + c"))
  expect_s3_class(plot(hg), "ggplot")
  expect_error(plot(hg, group = "Cluster 9"), class = "hypernets_bad_input")
  expect_error(plot(hg, group = c("Cluster 1", "Cluster 2")), class = "hypernets_bad_input")
  plain <- group_hypergraph(data.frame(a = c("x", "y"), g = c("1", "1")), "a", "g")
  expect_error(plot(plain, group = "1"), class = "hypernets_bad_input")
})

test_that("malformed input raises hypernets_bad_input", {
  bad <- .cs_mmm()
  bad$assignments <- bad$assignments[-1L]
  expect_error(group_hypergraph(bad), class = "hypernets_bad_input")
  bad <- .cs_mmm()
  bad$assignments[1L] <- 5L
  expect_error(group_hypergraph(bad), class = "hypernets_bad_input")
  bad <- .cs_mmm()
  bad$data <- NULL
  expect_error(group_hypergraph(bad), class = "hypernets_bad_input")
  expect_error(group_hypergraph(.cs_group(c("A", "A"))), class = "hypernets_bad_input")
  unnamed <- .cs_group()
  names(unnamed) <- NULL
  expect_error(group_hypergraph(unnamed), class = "hypernets_bad_input")
  expect_error(group_hypergraph(.cs_mmm(), top = 0), class = "hypernets_bad_input")
  expect_error(group_hypergraph(.cs_mmm(), states = 1), class = "hypernets_bad_input")
  expect_error(group_hypergraph(.cs_mmm(), states = "zzz"), class = "hypernets_bad_input")
  expect_error(group_hypergraph(.cs_data(), "V1", "V2", states = "a"),
               class = "hypernets_bad_input")
  plain <- group_hypergraph(data.frame(a = c("x", "y"), g = c("1", "1")), "a", "g")
  expect_error(hg_get(plain, what = "sets"), class = "hypernets_bad_input")
  expect_error(hg_get(plain, what = "state_counts"), class = "hypernets_bad_input")
})

test_that("data.frame input is unchanged (frozen before the clustering branch)", {
  frozen <- readRDS(test_path("fixtures", "group_hypergraph_frozen.rds"))
  # Frozen when the class was net_hypergraph; the rename to net_hg (0.6.0) is
  # names only, so relabel the class and compare everything else as is.
  frozen <- lapply(frozen, \(x) structure(x, class = "net_hg"))
  df <-data.frame(person = c("Alice", "Bob", "Carol", "Alice", "Bob", "Dave", "Carol", "Dave", "Eve"),
                   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3"),
                   w = c(1, 2, 3, 1, 1, 2, 2, 1, 5), day = c(1, 1, 1, 2, 2, 3, 3, 3, 3))
  now <- list(
    a = group_hypergraph(df, actor = "person", group = "session"),
    b = group_hypergraph(df, actor = "person", group = "session", weight = "w"),
    c = group_hypergraph(df, actor = "person", group = "session", sparse = TRUE),
    d = group_hypergraph(data.frame(from = c("a", "b"), to = c("b", "c")), from = "from", to = "to"),
    e = group_hypergraph(data.frame(ref = c("a; b", "b;c", ""), blk = c("x", "y", "z")),
                         actor = "ref", group = "blk", separator = ";"),
    f = group_hypergraph(df, "person", "session",
                         nodes = c("Alice", "Bob", "Carol", "Dave", "Eve", "Zed")),
    g = group_hypergraph(data.frame(from = c("a", "b"), to = c("b", "c")))
  )
  for (n in names(frozen)) expect_identical(now[[n]], frozen[[n]], label = n)
})
