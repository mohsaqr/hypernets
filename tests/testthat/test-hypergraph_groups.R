# ---- group_hypergraph() tests --------------------------------------------

# Helpers ------------------------------------------------------------------

.bg_sample_data <- function() {
  data.frame(
    member = c("Alice", "Bob", "Carol", "Alice", "Bob",
               "Dave", "Carol", "Dave", "Eve"),
    session = c("S1", "S1", "S1", "S2", "S2",
                "S3", "S3", "S3", "S3"),
    stringsAsFactors = FALSE
  )
}

# Structure ----------------------------------------------------------------

test_that("returns a net_hypergraph with required fields", {
  hg <- group_hypergraph(.bg_sample_data(), actor = "member", group = "session")
  expect_s3_class(hg, "net_hypergraph")
  expect_named(hg, c("hyperedges", "incidence", "nodes", "n_nodes",
                     "n_hyperedges", "size_distribution", "params"))
  expect_equal(hg$n_nodes, 5L)
  expect_equal(hg$n_hyperedges, 3L)
  expect_setequal(hg$nodes, c("Alice", "Bob", "Carol", "Dave", "Eve"))
})

# Hyperedge content --------------------------------------------------------

test_that("each group becomes a hyperedge spanning its members", {
  hg <- group_hypergraph(.bg_sample_data(), actor = "member", group = "session")
  # Map back: hyperedges are integer indices into hg$nodes
  members_by_session <- lapply(hg$hyperedges, function(idx) sort(hg$nodes[idx]))
  names(members_by_session) <- colnames(hg$incidence)
  expect_equal(members_by_session[["S1"]], c("Alice", "Bob", "Carol"))
  expect_equal(members_by_session[["S2"]], c("Alice", "Bob"))
  expect_equal(members_by_session[["S3"]], c("Carol", "Dave", "Eve"))
})

# Incidence matrix ---------------------------------------------------------

test_that("incidence is binary by default with correct dimensions and sums", {
  hg <- group_hypergraph(.bg_sample_data(), actor = "member", group = "session")
  expect_equal(dim(hg$incidence), c(5L, 3L))
  expect_true(all(hg$incidence %in% c(0L, 1L)))
  # column sums = group sizes
  expect_equal(unname(colSums(hg$incidence)), c(3L, 2L, 3L))
  # row sums = member hyperdegree (num groups they appeared in)
  expect_equal(hg$incidence["Alice", , drop = TRUE], c(S1 = 1L, S2 = 1L, S3 = 0L))
  expect_equal(rowSums(hg$incidence)[["Carol"]], 2L)
})

# Size distribution --------------------------------------------------------

test_that("size_distribution counts hyperedges by size", {
  hg <- group_hypergraph(.bg_sample_data(), actor = "member", group = "session")
  expect_equal(hg$size_distribution[["size_2"]], 1L)
  expect_equal(hg$size_distribution[["size_3"]], 2L)
})

# Weighted incidence -------------------------------------------------------

test_that("weight column produces weighted incidence (sum per cell)", {
  d <- data.frame(
    member = c("A", "B", "A", "A", "B"),
    grp    = c("g1", "g1", "g1", "g2", "g2"),
    n      = c(2, 5, 3, 1, 4),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, actor = "member", group = "grp", weight = "n")
  # A in g1: 2 + 3 = 5; B in g1: 5; A in g2: 1; B in g2: 4
  expect_equal(hg$incidence["A", "g1"], 5)
  expect_equal(hg$incidence["B", "g1"], 5)
  expect_equal(hg$incidence["A", "g2"], 1)
  expect_equal(hg$incidence["B", "g2"], 4)
  # hyperedges still based on membership (any nonzero), not weight
  expect_setequal(hg$nodes[hg$hyperedges[[1]]], c("A", "B"))
})

# NA handling --------------------------------------------------------------

test_that("rows with NA in member or group are dropped silently", {
  d <- data.frame(
    member = c("A", "B", NA,  "C", "D"),
    grp    = c("g1", "g1", "g1", NA,  "g2"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, actor = "member", group = "grp")
  expect_equal(hg$n_nodes, 3L)         # A, B, D (C dropped because group NA)
  expect_setequal(hg$nodes, c("A", "B", "D"))
  expect_equal(hg$n_hyperedges, 2L)    # g1, g2
  # n_observations records post-drop count
  expect_equal(hg$params$n_observations, 3L)
})

test_that("all-NA data raises an error", {
  d <- data.frame(member = NA, grp = NA, stringsAsFactors = FALSE)
  expect_error(group_hypergraph(d, "member", "grp"),
               "No complete observations")
})

# Repeated rows are deduplicated in binary mode ----------------------------

test_that("duplicate rows do not duplicate membership in binary mode", {
  d <- data.frame(
    member = c("A", "A", "A", "B"),
    grp    = c("g1", "g1", "g1", "g1"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  expect_equal(hg$incidence["A", "g1"], 1L)
  expect_equal(hg$incidence["B", "g1"], 1L)
  expect_equal(hg$n_hyperedges, 1L)
  expect_equal(length(hg$hyperedges[[1]]), 2L)
})

# Single-member groups -----------------------------------------------------

test_that("singleton groups become size-1 hyperedges", {
  d <- data.frame(
    member = c("A", "B", "C"),
    grp    = c("g1", "g2", "g3"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(d, "member", "grp")
  expect_equal(hg$n_hyperedges, 3L)
  expect_true(all(vapply(hg$hyperedges, length, integer(1)) == 1L))
})

test_that("an explicit node universe preserves isolates", {
  d <- data.frame(member = c("A", "B"), grp = c("g1", "g1"))
  hg <- group_hypergraph(d, "member", "grp", nodes = c("A", "B", "C"))
  expect_equal(hg$nodes, c("A", "B", "C"))
  expect_equal(unname(rowSums(hg$incidence)), c(1, 1, 0))
  expect_equal(hg$n_nodes, 3L)
  expect_error(group_hypergraph(d, "member", "grp", nodes = c("A", "C")),
               "Every observed member")
  expect_error(group_hypergraph(d, "member", "grp", nodes = c("A", "B", "B")),
               "unique")
})

test_that("sparse group incidence is identical to dense incidence", {
  d <- rbind(.bg_sample_data(), .bg_sample_data()[1, , drop = FALSE])
  dense <- group_hypergraph(d, "member", "session")
  sparse <- group_hypergraph(d, "member", "session", sparse = TRUE)
  expect_s4_class(sparse$incidence, "sparseMatrix")
  expect_equal(as.matrix(sparse$incidence), dense$incidence)
  expect_equal(sparse$hyperedges, dense$hyperedges)

  weighted <- transform(d, n = seq_len(nrow(d)))
  dense_w <- group_hypergraph(weighted, "member", "session", "n")
  sparse_w <- group_hypergraph(weighted, "member", "session", "n",
                               sparse = TRUE)
  expect_equal(as.matrix(sparse_w$incidence), dense_w$incidence)
})

# Validation ---------------------------------------------------------------

test_that("missing member column raises error", {
  d <- data.frame(member = c("A"), grp = c("g1"), stringsAsFactors = FALSE)
  expect_error(group_hypergraph(d, actor = "playor", group = "grp"))
})

test_that("missing group column raises error", {
  d <- data.frame(member = c("A"), grp = c("g1"), stringsAsFactors = FALSE)
  expect_error(group_hypergraph(d, actor = "member", group = "guppy"))
})

test_that("non-data.frame input rejected", {
  expect_error(group_hypergraph(list(member = "A", grp = "g"), "member", "grp"))
})

# print and summary work --------------------------------------------------

test_that("print and summary work via shared net_hypergraph methods", {
  hg <- group_hypergraph(.bg_sample_data(), actor = "member", group = "session")
  expect_invisible(print(hg))
  # summary now returns a tidy node-degree data.frame (visible)
  s <- summary(hg)
  expect_s3_class(s, "data.frame")
  expect_setequal(names(s), c("node", "degree"))
})

# Integration with bundled dataset ----------------------------------------

test_that("works on bundled human_long dataset (long-format event data)", {
  data("human_long", package = "hypernets")
  hg <- group_hypergraph(human_long, actor = "code", group = "session_id")
  expect_s3_class(hg, "net_hypergraph")
  expect_gt(hg$n_nodes, 0L)
  expect_gt(hg$n_hyperedges, 0L)
  # Each session is a hyperedge; members = codes appearing in that session
  expect_equal(ncol(hg$incidence), hg$n_hyperedges)
  expect_equal(nrow(hg$incidence), hg$n_nodes)
})

# Hyperedge attributes -----------------------------------------------------

test_that("columns constant within a group become hyperedge attributes", {
  dat <- data.frame(
    member = c("a", "b", "c", "b", "c", "d", "d", "e"),
    event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
    kind = c("x", "x", "x", "x", "x", "x", "y", "y"),
    seat = c("p", "q", "r", "p", "q", "r", "p", "q"),
    stringsAsFactors = FALSE
  )
  hg <- group_hypergraph(dat, actor = "member", group = "event")
  expect_identical(names(hg$edge_data), c("edge", "kind"))
  expect_identical(hg$edge_data$edge, c("e1", "e2", "e3"))
  expect_identical(hg$edge_data$kind, c("x", "x", "y"))
  kept <- hg_subset(hg, where = c(kind = "y"))
  expect_identical(kept$nodes, c("d", "e"))
  # an NA inside a group does not make the attribute vary
  dat$kind[2L] <- NA
  with_na <- group_hypergraph(dat, actor = "member", group = "event")
  expect_identical(with_na$edge_data$kind, c("x", "x", "y"))
  # sparse and dense agree
  sparse <- group_hypergraph(dat, actor = "member", group = "event",
                             sparse = TRUE)
  expect_identical(sparse$edge_data, with_na$edge_data)
})

test_that("INVARIANT: a plain member/group table keeps the original layout", {
  hg <- group_hypergraph(.bg_sample_data(), actor = "member", group = "session")
  expect_false("edge_data" %in% names(hg))
  weighted <- data.frame(member = c("a", "b"), session = c("s", "s"),
                         w = c(1, 2))
  hg_w <- group_hypergraph(weighted, actor = "member", group = "session",
                           weight = "w")
  expect_false("edge_data" %in% names(hg_w))
})

test_that("an edge list carries its attributes onto the size-2 hyperedges", {
  contacts <- data.frame(from = c("a", "b"), to = c("b", "c"),
                         kind = c("call", "mail"), stringsAsFactors = FALSE)
  hg <- group_hypergraph(contacts, from = "from", to = "to")
  expect_identical(hg$edge_data,
                   data.frame(edge = c("e1", "e2"), kind = c("call", "mail"),
                              stringsAsFactors = FALSE))
})

test_that("an unaddressable dense incidence is refused, not attempted", {
  skip_on_cran()
  # 50000 members x 50000 groups = 2.5e9 cells: past the integer cell index,
  # and ~20 Gb if it were allocated.
  n <- 50000L
  d <- data.frame(member = sprintf("m%d", seq_len(n)),
                  group = sprintf("g%d", seq_len(n)))
  expect_error(group_hypergraph(d, actor = "member", group = "group"),
               class = "hypernets_dense_too_large")
  expect_s3_class(group_hypergraph(d, actor = "member", group = "group",
                                   sparse = TRUE),
                  "net_hypergraph")
})
