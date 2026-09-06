testthat::skip_on_cran()

# ---- Tidy accessors: every net_* result class ----------------------------
#
# Taxonomy contract (honets 0.2.0): every result object in every structure
# family is reachable with as.data.frame(), a secondary table is selected
# with what =, and the return is always a base data.frame with >= 1 column
# and no row names. Reaching into a result with $ is never required.

.ac_mat <- function(seed = 1L, n = 8L) {
  set.seed(seed)
  m <- matrix(stats::runif(n * n), n, n)
  m <- (m + t(m)) / 2
  diag(m) <- 0
  dimnames(m) <- list(LETTERS[seq_len(n)], LETTERS[seq_len(n)])
  m
}

.ac_seqs <- function() split(human_long$code, human_long$session_id)

.ac_wide <- function() {
  s <- .ac_seqs()
  L <- max(lengths(s))
  as.data.frame(do.call(rbind, lapply(s, function(x) c(x, rep(NA, L - length(x))))),
                stringsAsFactors = FALSE)
}

# A tidy table is a base data.frame with rows, columns and default row names.
expect_tidy <- function(d, min_rows = 1L) {
  expect_s3_class(d, "data.frame")
  expect_identical(class(d), "data.frame")
  expect_gte(nrow(d), min_rows)
  expect_gte(ncol(d), 1L)
  expect_identical(attr(d, "row.names"), seq_len(nrow(d)))
  invisible(d)
}

# ---- memory family -------------------------------------------------------

test_that("as.data.frame.net_hon returns rules and nodes", {
  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  rules <- as.data.frame(h)
  expect_tidy(rules)
  expect_named(rules, c("path", "from", "to", "count", "probability",
                        "from_order", "to_order"))
  expect_identical(rules, h$ho_edges, ignore_attr = "row.names")

  nodes <- as.data.frame(h, what = "nodes")
  expect_tidy(nodes)
  expect_named(nodes, c("id", "label", "name"))
  expect_identical(nrow(nodes), h$n_nodes)
})

test_that("as.data.frame.net_hon filters by order and sorts", {
  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  ho <- as.data.frame(h, order_min = 2L)
  expect_true(all(ho$from_order >= 2L))
  all_rules <- as.data.frame(h)
  expect_lt(nrow(ho), nrow(all_rules))

  sorted <- as.data.frame(h, sort_by = "count")
  expect_identical(sorted$count, sort(sorted$count, decreasing = TRUE))
  # sorting is a permutation, not a filter
  expect_identical(sum(sorted$count), sum(all_rules$count))
})

test_that("as.data.frame.net_honem returns embeddings and variance", {
  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  em <- build_honem(h, dim = 4L)
  d <- as.data.frame(em)
  expect_tidy(d)
  expect_named(d, c("node", "dim1", "dim2", "dim3", "dim4"))
  expect_identical(d$node, em$nodes)

  v <- as.data.frame(em, what = "variance")
  expect_tidy(v)
  expect_named(v, c("dim", "singular_value", "proportion"))
  # proportions are a distribution over the retained dimensions
  expect_equal(sum(v$proportion), 1)
  expect_true(all(v$proportion >= 0))
})

test_that("as.data.frame.net_hypa returns scores, over and under", {
  hp <- build_hypa(.ac_seqs(), order = 2L)
  all_s <- as.data.frame(hp)
  expect_tidy(all_s)
  over <- as.data.frame(hp, what = "over")
  expect_tidy(over, min_rows = 0L)
  under <- as.data.frame(hp, what = "under")
  expect_tidy(under, min_rows = 0L)
  expect_lte(nrow(over) + nrow(under), nrow(all_s))
  expect_identical(names(over), names(all_s))

  # sort_by = "p" orders by the two-sided p_value, smallest first
  by_p <- as.data.frame(hp, sort_by = "p")
  expect_false(is.unsorted(by_p$p_value))
  by_ratio <- as.data.frame(hp, sort_by = "ratio")
  expect_identical(by_ratio$ratio, sort(by_ratio$ratio, decreasing = TRUE))
})

test_that("as.data.frame.net_mogen returns the order table and transitions", {
  mo <- build_mogen(.ac_seqs(), max_order = 2L)
  ord <- as.data.frame(mo)
  expect_tidy(ord)
  expect_named(ord, c("order", "log_likelihood", "aic", "bic", "dof",
                      "layer_dof", "optimal"))
  expect_identical(sum(ord$optimal), 1L)
  expect_identical(ord$order[ord$optimal], mo$optimal_order)

  tr <- as.data.frame(mo, what = "transitions")
  expect_tidy(tr)
  tr_direct <- mogen_transitions(mo, order = mo$optimal_order)
  expect_identical(tr, tr_direct)
})

test_that("as.data.frame.net_markov_order returns orders and the null", {
  mk <- markov_order_test(.ac_wide(), max_order = 2L, n_perm = 20L, seed = 1L)
  ord <- as.data.frame(mk)
  expect_tidy(ord)
  expect_identical(ord, mk$test_table, ignore_attr = "row.names")

  nul <- as.data.frame(mk, what = "null")
  expect_tidy(nul)
  expect_named(nul, c("order", "replicate", "g2"))
  # one row per (tested order, permutation replicate)
  expect_identical(nrow(nul), length(unlist(mk$permutation_null)))
  expect_true(all(nul$replicate >= 1L))
})

test_that("as.data.frame.net_path_dependence filters and sorts", {
  pd <- path_dependence(.ac_wide(), order = 2L)
  d <- as.data.frame(pd)
  expect_tidy(d)
  expect_true(all(c("context", "n", "KL") %in% names(d)))

  filt <- as.data.frame(pd, min_count = 50L)
  expect_true(all(filt$n >= 50L))
  expect_lte(nrow(filt), nrow(d))

  srt <- as.data.frame(pd, sort_by = "KL")
  expect_identical(srt$KL, sort(srt$KL, decreasing = TRUE))
})

# ---- simplicial family ---------------------------------------------------

test_that("as.data.frame.net_simplicial returns simplices and the f-vector", {
  sc <- build_simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  d <- as.data.frame(sc)
  expect_tidy(d)
  expect_named(d, c("id", "dim", "size", "members"))
  expect_identical(nrow(d), sc$n_simplices)
  expect_identical(d$size, d$dim + 1L)

  fv <- as.data.frame(sc, what = "f_vector")
  expect_tidy(fv)
  expect_named(fv, c("dim", "count"))
  # the f-vector must account for every simplex, and match the per-dim counts
  expect_identical(sum(fv$count), sc$n_simplices)
  expect_identical(as.integer(table(d$dim)), fv$count[fv$count > 0L])
})

test_that("as.data.frame.net_simplicial filters by dimension", {
  sc <- build_simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  d0 <- as.data.frame(sc, dim = 0L)
  expect_true(all(d0$dim == 0L))
  # every node contributes exactly one 0-simplex
  expect_identical(nrow(d0), sc$n_nodes)
})

test_that("as.data.frame.net_q_analysis returns q levels and node max-q", {
  sc <- build_simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  qa <- q_analysis(sc)
  lv <- as.data.frame(qa)
  expect_tidy(lv)
  expect_named(lv, c("q", "components"))
  expect_identical(max(lv$q), qa$max_q)
  expect_true(all(lv$components >= 1L))

  nd <- as.data.frame(qa, what = "nodes")
  expect_tidy(nd)
  expect_named(nd, c("node", "max_q"))
  expect_identical(nd$node, sc$nodes)
})

test_that("as.data.frame.net_persistent_homology returns diagram and curves", {
  ph <- persistent_homology(.ac_mat(), n_steps = 6L, max_dim = 2L)
  pd <- as.data.frame(ph)
  expect_tidy(pd)
  expect_named(pd, c("dimension", "birth", "death", "persistence"))
  # persistence is death - birth, and never negative
  expect_true(all(pd$persistence >= 0))

  bc <- as.data.frame(ph, what = "betti")
  expect_tidy(bc)
  expect_named(bc, c("threshold", "dimension", "betti"))

  srt <- as.data.frame(ph, sort_by = "persistence")
  expect_identical(srt$persistence, sort(srt$persistence, decreasing = TRUE))
  dim1 <- as.data.frame(ph, dimension = 1L)
  expect_true(all(dim1$dimension == 1L))
})

test_that("as.data.frame.net_persistence_landscape returns the grid", {
  ph <- persistent_homology(.ac_mat(), n_steps = 6L, max_dim = 2L)
  pl <- persistence_landscape(ph, k_max = 3L)
  d <- as.data.frame(pl)
  expect_tidy(d)
  expect_named(d, c("k", "t", "value"))
  expect_identical(sort(unique(d$k)), 1:3)
  # landscape functions are non-negative by construction
  expect_true(all(d$value >= 0))
  k2 <- as.data.frame(pl, k = 2L)
  expect_true(all(k2$k == 2L))
})

# ---- hypergraph family ---------------------------------------------------

test_that("as.data.frame.net_hypergraph_measures returns all three tables", {
  hg <- build_hypergraph(.ac_mat(), threshold = 0.5)
  hm <- hypergraph_measures(hg)

  nd <- as.data.frame(hm)
  expect_tidy(nd)
  expect_named(nd, c("node", "hyperdegree", "node_strength", "max_edge_size"))
  expect_identical(nrow(nd), hg$n_nodes)

  ed <- as.data.frame(hm, what = "edges")
  expect_tidy(ed)
  expect_named(ed, c("hyperedge", "size"))
  expect_identical(nrow(ed), hg$n_hyperedges)
  # summed hyperedge sizes must equal summed hyperdegrees (both count
  # (node, hyperedge) memberships)
  expect_identical(sum(ed$size), sum(nd$hyperdegree))

  gl <- as.data.frame(hm, what = "global")
  expect_tidy(gl)
  expect_named(gl, c("measure", "value"))
  expect_identical(gl$value[gl$measure == "n_nodes"], as.numeric(hg$n_nodes))

  srt <- as.data.frame(hm, sort_by = "hyperdegree")
  expect_identical(srt$hyperdegree, sort(srt$hyperdegree, decreasing = TRUE))
})

# ---- error paths ---------------------------------------------------------

test_that("accessors reject unknown what = and bad filters", {
  sc <- build_simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  expect_error(as.data.frame(sc, what = "nope"), "should be one of")
  expect_error(as.data.frame(sc, dim = -1L), "`dim` must be")

  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  expect_error(as.data.frame(h, what = "nope"), "should be one of")
  expect_error(as.data.frame(h, order_min = 0L), "`order_min` must be")
  expect_error(as.data.frame(h, sort_by = "nope"), "should be one of")

  ph <- persistent_homology(.ac_mat(), n_steps = 6L, max_dim = 2L)
  expect_error(
    as.data.frame(ph, what = "betti", sort_by = "persistence"),
    class = "honets_bad_input"
  )

  pd <- path_dependence(.ac_wide(), order = 2L)
  expect_error(as.data.frame(pd, min_count = 0L), "`min_count` must be")
})

test_that("every exported net_* result class has an as.data.frame method", {
  # The taxonomy contract: no result object forces the user to reach in
  # with $. If a new net_* class appears without an accessor, this fails.
  methods <- as.character(utils::methods("as.data.frame"))
  covered <- sub("^as\\.data\\.frame\\.", "", methods)
  expected <- c(
    "net_hon", "net_hon_boot", "net_hon_compare", "net_honem", "net_hypa",
    "net_markov_order", "net_mogen", "net_path_dependence",
    "net_simplicial", "net_q_analysis", "net_persistent_homology",
    "net_persistence_landscape",
    "net_hypergraph", "net_hypergraph_cluster", "net_hypergraph_measures",
    "net_hypergraph_transduction"
  )
  expect_setequal(intersect(expected, covered), expected)
})

# ---- top = ---------------------------------------------------------------
#
# `top` is the argument that removes head() from the public surface. Its
# contract: applied LAST, after `what`, after every filter and after
# `sort_by`, so `sort_by` and `top` compose. A test that only checked
# nrow() would pass on an accessor that truncated before sorting, which is
# the bug this argument exists to prevent -- so every case below compares
# against head() of the fully-ordered table.

# Every accessor that can return many rows, as (label, full, topped) thunks.
.ac_top_cases <- function() {
  seqs <- .ac_seqs()
  h  <- build_hon(seqs, max_order = 2L, min_freq = 20L)
  em <- build_honem(h, dim = 4L)
  hp <- build_hypa(seqs, order = 2L, min_count = 20L)
  mo <- build_mogen(seqs, max_order = 2L)
  mk <- markov_order_test(seqs, max_order = 2L, n_perm = 10L, seed = 1L)
  pd <- path_dependence(.ac_wide(), order = 2L, min_count = 20L)
  bs <- bootstrap_hon(seqs, n_boot = 10L, max_order = 2L, min_freq = 20L,
                      seed = 1L)
  grp <- rep(c("a", "b"), length.out = length(seqs))
  cm <- compare_hon(seqs[grp == "a"], seqs[grp == "b"], n_perm = 10L,
                    max_order = 2L, min_freq = 10L, seed = 1L)
  m  <- .ac_mat()
  sc <- build_simplicial(m, type = "clique", threshold = 0.5)
  ph <- persistent_homology(m, n_steps = 6L, max_dim = 2L)
  pl <- persistence_landscape(ph, dimension = 1L)
  qa <- q_analysis(sc)
  hg <- window_hypergraph(seqs, window = 3L)
  hm <- hypergraph_measures(hg)
  cl <- hypergraph_cluster(hg, k = 3L, seed = 1L)
  lb <- stats::setNames(rep(NA_character_, hg$n_nodes), hg$nodes)
  lb[c(1L, 2L)] <- c("x", "y")
  tr <- hypergraph_transduction(hg, labels = lb, xi = 0.9)

  list(
    list("net_hon rules",     function(n) as.data.frame(h, sort_by = "count", top = n),
                              as.data.frame(h, sort_by = "count")),
    list("net_hon nodes",     function(n) as.data.frame(h, what = "nodes", top = n),
                              as.data.frame(h, what = "nodes")),
    list("net_honem emb",     function(n) as.data.frame(em, top = n),
                              as.data.frame(em)),
    list("net_honem var",     function(n) as.data.frame(em, what = "variance", top = n),
                              as.data.frame(em, what = "variance")),
    list("net_hypa",          function(n) as.data.frame(hp, sort_by = "ratio", top = n),
                              as.data.frame(hp, sort_by = "ratio")),
    list("mogen_transitions", function(n) mogen_transitions(mo, order = 2L, top = n),
                              mogen_transitions(mo, order = 2L)),
    list("net_mogen orders",  function(n) as.data.frame(mo, top = n),
                              as.data.frame(mo)),
    list("net_mogen trans",   function(n) as.data.frame(mo, what = "transitions",
                                                        order = 2L, top = n),
                              as.data.frame(mo, what = "transitions", order = 2L)),
    list("path_counts",       function(n) path_counts(seqs, k = 3L, top = n),
                              path_counts(seqs, k = 3L)),
    list("net_markov_order",  function(n) as.data.frame(mk, top = n),
                              as.data.frame(mk)),
    list("net_markov null",   function(n) as.data.frame(mk, what = "null", top = n),
                              as.data.frame(mk, what = "null")),
    list("net_path_dep",      function(n) as.data.frame(pd, sort_by = "KL", top = n),
                              as.data.frame(pd, sort_by = "KL")),
    list("net_hon_boot",      function(n) as.data.frame(bs, sort_by = "support", top = n),
                              as.data.frame(bs, sort_by = "support")),
    list("net_hon_compare",   function(n) as.data.frame(cm, sort_by = "p_adj", top = n),
                              as.data.frame(cm, sort_by = "p_adj")),
    list("hon_centrality",    function(n) hon_centrality(h, sort_by = "pagerank", top = n),
                              hon_centrality(h, sort_by = "pagerank")),
    list("net_simplicial",    function(n) as.data.frame(sc, top = n),
                              as.data.frame(sc)),
    list("net_simplicial fv", function(n) as.data.frame(sc, what = "f_vector", top = n),
                              as.data.frame(sc, what = "f_vector")),
    # simplicial_degree's default carries the stale row names it inherited
    # from Nestimate 0.9.0 (pinned by the cross-package identity test); the
    # truncated return resets them, exactly as .ho_top() does.
    list("simplicial_degree", function(n) simplicial_degree(sc, top = n),
                              simplicial_degree(sc)),
    list("net_q_analysis",    function(n) as.data.frame(qa, top = n),
                              as.data.frame(qa)),
    list("net_q_analysis nd", function(n) as.data.frame(qa, what = "nodes", top = n),
                              as.data.frame(qa, what = "nodes")),
    list("net_persist_homol", function(n) as.data.frame(ph, sort_by = "persistence",
                                                        top = n),
                              as.data.frame(ph, sort_by = "persistence")),
    list("net_pers_landscape",function(n) as.data.frame(pl, top = n),
                              as.data.frame(pl)),
    list("net_hypergraph",    function(n) as.data.frame(hg, sort_by = "weight", top = n),
                              as.data.frame(hg, sort_by = "weight")),
    list("hypergraph_centr",  function(n) hypergraph_centrality(hg, sort_by = "clique",
                                                                top = n),
                              hypergraph_centrality(hg, sort_by = "clique")),
    list("net_hg_measures",   function(n) as.data.frame(hm, sort_by = "hyperdegree",
                                                        top = n),
                              as.data.frame(hm, sort_by = "hyperdegree")),
    list("net_hg_cluster",    function(n) as.data.frame(cl, top = n),
                              as.data.frame(cl)),
    list("net_hg_transd",     function(n) as.data.frame(tr, top = n),
                              as.data.frame(tr)),
    list("net_hg_transd sc",  function(n) as.data.frame(tr, what = "scores", top = n),
                              as.data.frame(tr, what = "scores"))
  )
}

test_that("top = returns the first n rows of the sorted table, everywhere", {
  for (case in .ac_top_cases()) {
    label <- case[[1L]]
    topped <- case[[2L]]
    full <- case[[3L]]
    n <- min(3L, nrow(full))
    expect_gte(nrow(full), n)
    # the contract: top-n OF THE ORDERED TABLE, not of an arbitrary order
    expect_identical(topped(n), `rownames<-`(utils::head(full, n), NULL),
                     info = label)
    expect_identical(nrow(topped(n)), n, info = label)
  }
})

test_that("top = NULL is the default and changes nothing", {
  for (case in .ac_top_cases()) {
    expect_identical(case[[2L]](NULL), case[[3L]], info = case[[1L]])
  }
})

test_that("top = larger than the table returns the whole table", {
  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  full <- as.data.frame(h, sort_by = "count")
  over_topped <- as.data.frame(h, sort_by = "count", top = nrow(full) + 100L)
  expect_identical(over_topped, full)
})

test_that("top = rejects a non-whole, zero, negative or vector value", {
  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  expect_error(as.data.frame(h, top = 0L),      "`top` must be")
  expect_error(as.data.frame(h, top = -3L),     "`top` must be")
  expect_error(as.data.frame(h, top = 2.5),     "`top` must be")
  expect_error(as.data.frame(h, top = c(2L, 3L)), "`top` must be")
  expect_error(as.data.frame(h, top = "3"),     "`top` must be")
  expect_error(as.data.frame(h, top = Inf),     "`top` must be")
  # the same guard is shared by every accessor, so one non-as.data.frame
  # verb is checked too
  expect_error(hon_centrality(h, top = 0L), "`top` must be")
  expect_error(path_counts(.ac_seqs(), k = 2L, top = 0L), "`top` must be")
})

test_that("sort_by and top compose rather than fighting", {
  # The regression this exists for: truncating BEFORE sorting would give the
  # top n of construction order, re-sorted -- a different, wrong answer.
  h <- build_hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  by_count <- as.data.frame(h, sort_by = "count", top = 5L)
  expect_identical(by_count$count, sort(by_count$count, decreasing = TRUE))
  all_rules <- as.data.frame(h)
  expect_identical(max(by_count$count), max(all_rules$count))

  unsorted_then_topped <- as.data.frame(h, top = 5L)
  expect_false(identical(unsorted_then_topped, by_count))

  # top interacts with a filter too: order_min narrows, then top truncates
  ho <- as.data.frame(h, order_min = 2L, sort_by = "count", top = 4L)
  expect_true(all(ho$from_order >= 2L))
  expect_identical(nrow(ho), 4L)
  ho_full <- as.data.frame(h, order_min = 2L, sort_by = "count")
  expect_identical(ho, `rownames<-`(utils::head(ho_full, 4L), NULL))
})
