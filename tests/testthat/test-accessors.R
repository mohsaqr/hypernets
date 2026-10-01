testthat::skip_on_cran()

# ---- Tidy accessors: every net_* result class ----------------------------
#
# Taxonomy contract (hypernets 0.2.0): every result object in every structure
# family is reachable with hg_get(), a secondary table is selected
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

test_that("hg_get.net_hon returns rules and nodes", {
  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  rules <- hg_get(h)
  expect_tidy(rules)
  expect_named(rules, c("path", "from", "to", "count", "probability",
                        "from_order", "to_order"))
  expect_identical(rules, h$ho_edges, ignore_attr = "row.names")

  nodes <- hg_get(h, what = "nodes")
  expect_tidy(nodes)
  expect_named(nodes, c("id", "node"))
  expect_identical(nrow(nodes), h$n_nodes)
})

test_that("hg_get.net_hon filters by order and sorts", {
  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  ho <- hg_get(h, order_min = 2L)
  expect_true(all(ho$from_order >= 2L))
  all_rules <- hg_get(h)
  expect_lt(nrow(ho), nrow(all_rules))

  sorted <- hg_get(h, sort_by = "count")
  expect_identical(sorted$count, sort(sorted$count, decreasing = TRUE))
  # sorting is a permutation, not a filter
  expect_identical(sum(sorted$count), sum(all_rules$count))
})

test_that("hg_get.net_honem returns embeddings and variance", {
  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  em <- honem(h, dim = 4L)
  d <- hg_get(em)
  expect_tidy(d)
  expect_named(d, c("node", "dim1", "dim2", "dim3", "dim4"))
  expect_identical(d$node, em$nodes)

  v <- hg_get(em, what = "variance")
  expect_tidy(v)
  expect_named(v, c("dimension", "singular_value", "proportion"))
  # proportions are a distribution over the retained dimensions
  expect_equal(sum(v$proportion), 1)
  expect_true(all(v$proportion >= 0))
})

test_that("hg_get.net_hypa returns scores, over and under", {
  hp <- hypa(.ac_seqs(), order = 2L)
  all_s <- hg_get(hp, what = "scores")
  expect_tidy(all_s)
  anomalies <- hg_get(hp)
  expect_tidy(anomalies, min_rows = 0L)
  over <- hg_get(hp, what = "over")
  expect_tidy(over, min_rows = 0L)
  under <- hg_get(hp, what = "under")
  expect_tidy(under, min_rows = 0L)
  expect_lte(nrow(over) + nrow(under), nrow(all_s))
  expect_identical(names(over), names(all_s))

  # order_by = "sig" orders scores by the two-sided p_value, smallest first;
  # "ratio" by the absolute log deviation, largest first
  by_p <- hg_get(hp, what = "scores", order_by = "sig")
  expect_false(is.unsorted(by_p$p_value))
  by_ratio <- hg_get(hp, what = "scores", order_by = "ratio")
  expect_false(is.unsorted(rev(abs(log(by_ratio$ratio)))))
})

test_that("hg_get.net_mogen returns the order table and transitions", {
  mo <- mogen(.ac_seqs(), max_order = 2L)
  ord <- hg_get(mo)
  expect_tidy(ord)
  expect_named(ord, c("order", "log_likelihood", "aic", "bic", "df",
                      "layer_df", "optimal"))
  expect_identical(sum(ord$optimal), 1L)
  expect_identical(ord$order[ord$optimal], mo$optimal_order)

  tr <- hg_get(mo, what = "transitions")
  expect_tidy(tr)
  tr_direct <- Nestimate::mogen_transitions(mo, order = mo$optimal_order)
  expect_identical(tr, tr_direct)
})

test_that("hg_get.net_markov_order returns orders and the null", {
  mk <- markov_order(.ac_wide(), max_order = 2L, n_perm = 20L, seed = 1L)
  ord <- hg_get(mk)
  expect_tidy(ord)
  # the estimator's test table in the package vocabulary
  expect_identical(ord, .ho_rename(mk$test_table,
                                   c(loglik = "log_likelihood", AIC = "aic",
                                     BIC = "bic", p_permutation = "p_value")),
                   ignore_attr = "row.names")
  expect_named(ord, c("order", "log_likelihood", "aic", "bic", "df", "g2",
                      "p_value", "p_asymptotic", "significant"))

  nul <- hg_get(mk, what = "null")
  expect_tidy(nul)
  expect_named(nul, c("order", "replicate", "g2"))
  # one row per (tested order, permutation replicate)
  expect_identical(nrow(nul), length(unlist(mk$permutation_null)))
  expect_true(all(nul$replicate >= 1L))
})

test_that("hg_get.net_path_dependence filters and sorts", {
  pd <- memory(.ac_wide(), order = 2L)
  d <- hg_get(pd)
  expect_tidy(d)
  expect_true(all(c("context", "count", "kl") %in% names(d)))

  filt <- hg_get(pd, min_count = 50L)
  expect_true(all(filt$count >= 50L))
  expect_lte(nrow(filt), nrow(d))

  srt <- hg_get(pd, sort_by = "kl")
  expect_identical(srt$kl, sort(srt$kl, decreasing = TRUE))
})

# ---- simplicial family ---------------------------------------------------

test_that("hg_get.simplicial_complex returns simplices and the f-vector", {
  sc <- simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  d <- hg_get(sc)
  expect_tidy(d)
  expect_named(d, c("simplex", "dimension", "size", "members"))
  expect_identical(nrow(d), sc$n_simplices)
  expect_identical(d$size, d$dim + 1L)

  fv <- hg_get(sc, what = "f_vector")
  expect_tidy(fv)
  expect_named(fv, c("dimension", "count"))
  # the f-vector must account for every simplex, and match the per-dim counts
  expect_identical(sum(fv$count), sc$n_simplices)
  expect_identical(as.integer(table(d$dim)), fv$count[fv$count > 0L])
})

test_that("hg_get.simplicial_complex filters by dimension", {
  sc <- simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  d0 <- hg_get(sc, dimension = 0L)
  expect_true(all(d0$dimension == 0L))
  # every node contributes exactly one 0-simplex
  expect_identical(nrow(d0), sc$n_nodes)
})

test_that("hg_get.q_analysis returns q levels and node max-q", {
  sc <- simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  qa <- hg_qanalysis(sc)
  lv <- hg_get(qa)
  expect_tidy(lv)
  expect_named(lv, c("q", "components"))
  expect_identical(max(lv$q), qa$max_q)
  expect_true(all(lv$components >= 1L))

  nd <- hg_get(qa, what = "nodes")
  expect_tidy(nd)
  expect_named(nd, c("node", "max_q"))
  expect_identical(nd$node, sc$nodes)
})

test_that("hg_get.persistent_homology returns diagram and curves", {
  ph <- hg_homology(.ac_mat(), n_steps = 6L, max_dim = 2L)
  pd <- hg_get(ph)
  expect_tidy(pd)
  expect_named(pd, c("dimension", "birth", "death", "persistence"))
  # persistence is death - birth, and never negative
  expect_true(all(pd$persistence >= 0))

  bc <- hg_get(ph, what = "betti")
  expect_tidy(bc)
  expect_named(bc, c("threshold", "dimension", "betti"))

  srt <- hg_get(ph, sort_by = "persistence")
  expect_identical(srt$persistence, sort(srt$persistence, decreasing = TRUE))
  dim1 <- hg_get(ph, dimension = 1L)
  expect_true(all(dim1$dimension == 1L))
})

test_that("hg_get.persistence_landscape returns the grid", {
  ph <- hg_homology(.ac_mat(), n_steps = 6L, max_dim = 2L)
  pl <- hg_landscape(ph, k_max = 3L)
  d <- hg_get(pl)
  expect_tidy(d)
  expect_named(d, c("k", "t", "value"))
  expect_identical(sort(unique(d$k)), 1:3)
  # landscape functions are non-negative by construction
  expect_true(all(d$value >= 0))
  k2 <- hg_get(pl, k = 2L)
  expect_true(all(k2$k == 2L))
})

# ---- hypergraph family ---------------------------------------------------

test_that("hg_get.net_hg_measures returns all three tables", {
  hg <- network_hypergraph(.ac_mat(), threshold = 0.5)
  hm <- .hg_measures_fit(hg)

  nd <- hg_get(hm)
  expect_tidy(nd)
  expect_named(nd, c("node", "hyperdegree", "node_strength", "max_edge_size"))
  expect_identical(nrow(nd), hg$n_nodes)

  ed <- hg_get(hm, what = "edges")
  expect_tidy(ed)
  expect_named(ed, c("hyperedge", "size"))
  expect_identical(nrow(ed), hg$n_hyperedges)
  # summed hyperedge sizes must equal summed hyperdegrees (both count
  # (node, hyperedge) memberships)
  expect_identical(sum(ed$size), sum(nd$hyperdegree))

  gl <- hg_get(hm, what = "global")
  expect_tidy(gl)
  expect_named(gl, c("measure", "value"))
  expect_identical(gl$value[gl$measure == "n_nodes"], as.numeric(hg$n_nodes))

  srt <- hg_get(hm, sort_by = "hyperdegree")
  expect_identical(srt$hyperdegree, sort(srt$hyperdegree, decreasing = TRUE))
})

# ---- error paths ---------------------------------------------------------

test_that("accessors reject unknown what = and bad filters", {
  sc <- simplicial(.ac_mat(), type = "clique", threshold = 0.5)
  expect_error(hg_get(sc, what = "nope"), "should be one of")
  expect_error(hg_get(sc, dimension = -1L), "`dimension` must be")

  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  expect_error(hg_get(h, what = "nope"), "should be one of")
  expect_error(hg_get(h, order_min = 0L), "`order_min` must be")
  expect_error(hg_get(h, sort_by = "nope"), "should be one of")

  ph <- hg_homology(.ac_mat(), n_steps = 6L, max_dim = 2L)
  expect_error(
    hg_get(ph, what = "betti", sort_by = "persistence"),
    class = "hypernets_bad_input"
  )

  pd <- memory(.ac_wide(), order = 2L)
  expect_error(hg_get(pd, min_count = 0L), "`min_count` must be")
})

test_that("every hypernets result class has an hg_get method", {
  # The taxonomy contract: no result object forces the user to reach in with $.
  #
  # The class list is DISCOVERED from the package itself: every hypernets
  # result class has a print method by the same taxonomy, so hypernets' own
  # print methods enumerate its classes. The namespace's S3 registry is the
  # parsed NAMESPACE, so it reads the same under devtools::load_all() and an
  # installed build.
  registry <- getNamespaceInfo("hypernets", "S3methods")
  declared <- function(generic) registry[registry[, 1L] == generic, 2L]
  classes <- setdiff(sort(unique(declared("print"))),
                     # print-only views of data.frame results that are
                     # returned as the table itself
                     c("hypernets_keywords"))
  covered <- declared("hg_get")
  expect_true(length(classes) >= 8L)   # guard against discovering nothing
  expect_setequal(intersect(classes, covered), classes)
})

test_that("no as.data.frame method is registered", {
  registry <- getNamespaceInfo("hypernets", "S3methods")
  expect_length(registry[registry[, 1L] == "as.data.frame", 2L], 0L)
  exports <- getNamespaceExports("hypernets")
  expect_length(grep("^as\\.data\\.frame", exports, value = TRUE), 0L)
})

test_that("hg_get has a method for every imported result class, and nothing else", {
  # The memory and simplicial estimators return their own classes. hypernets
  # adds hg_get() for them and registers NO print/summary/plot method on
  # them, which would overwrite the estimator's own whenever both are loaded.
  imported_classes <- c(
    "net_hon", "net_honem", "net_hypa", "net_mogen", "net_markov_order",
    "net_markov_order_group", "net_path_dependence", "net_markov_stability",
    "simplicial_complex", "persistent_homology", "q_analysis",
    "persistence_landscape"
  )
  ours <- getNamespaceInfo("hypernets", "S3methods")
  theirs <- getNamespaceInfo("Nestimate", "S3methods")
  ours_on <- function(generic) ours[ours[, 1L] == generic, 2L]
  expect_setequal(intersect(ours_on("hg_get"), imported_classes),
                  imported_classes)
  shared <- merge(as.data.frame(ours[, 1:2, drop = FALSE]),
                  as.data.frame(theirs[, 1:2, drop = FALSE]))
  expect_identical(nrow(shared), 0L)
})

test_that("hg_get.default names the class it cannot read", {
  expect_error(hg_get(data.frame(a = 1)), class = "hypernets_bad_input")
  expect_error(hg_get(list(1)), "`list`")
  expect_error(hg_get(1:3), "`integer`")
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
  h  <- hon(seqs, max_order = 2L, min_freq = 20L)
  em <- honem(h, dim = 4L)
  hp <- hypa(seqs, order = 2L, min_count = 20L)
  mo <- mogen(seqs, max_order = 2L)
  mk <- markov_order(seqs, max_order = 2L, n_perm = 10L, seed = 1L)
  pd <- memory(.ac_wide(), order = 2L, min_count = 20L)
  bs <- hg_bootstrap(seqs, n_boot = 10L, max_order = 2L, min_freq = 20L,
                      seed = 1L)
  grp <- rep(c("a", "b"), length.out = length(seqs))
  cm <- hg_compare(hon(seqs, group = grp, max_order = 2L, min_freq = 10L),
                   n_perm = 10L, seed = 1L)
  m  <- .ac_mat()
  sc <- simplicial(m, type = "clique", threshold = 0.5)
  ph <- hg_homology(m, n_steps = 6L, max_dim = 2L)
  pl <- hg_landscape(ph, dimension = 1L)
  qa <- hg_qanalysis(sc)
  hg <- window_hypergraph(seqs, window = 3L)
  hm <- .hg_measures_fit(hg)
  cl <- .hg_cluster_fit(hg, k = 3L, seed = 1L)
  lb <- stats::setNames(rep(NA_character_, hg$n_nodes), hg$nodes)
  lb[c(1L, 2L)] <- c("x", "y")
  tr <- .hg_transduction_fit(hg, labels = lb, xi = 0.9)

  list(
    list("net_hon rules",     function(n) hg_get(h, sort_by = "count", top = n),
                              hg_get(h, sort_by = "count")),
    list("net_hon nodes",     function(n) hg_get(h, what = "nodes", top = n),
                              hg_get(h, what = "nodes")),
    list("net_honem emb",     function(n) hg_get(em, top = n),
                              hg_get(em)),
    list("net_honem var",     function(n) hg_get(em, what = "variance", top = n),
                              hg_get(em, what = "variance")),
    list("net_hypa",          function(n) hg_get(hp, what = "scores", order_by = "ratio", top = n),
                              hg_get(hp, what = "scores", order_by = "ratio")),
    list("net_mogen orders",  function(n) hg_get(mo, top = n),
                              hg_get(mo)),
    list("net_mogen trans",   function(n) hg_get(mo, what = "transitions",
                                                        order = 2L, top = n),
                              hg_get(mo, what = "transitions", order = 2L)),
    list("net_mogen paths",   function(n) hg_get(mo, what = "paths", k = 3L,
                                                        top = n),
                              hg_get(mo, what = "paths", k = 3L)),
    list("net_markov_order",  function(n) hg_get(mk, top = n),
                              hg_get(mk)),
    list("net_markov null",   function(n) hg_get(mk, what = "null", top = n),
                              hg_get(mk, what = "null")),
    list("net_path_dep",      function(n) hg_get(pd, sort_by = "kl", top = n),
                              hg_get(pd, sort_by = "kl")),
    list("net_hon_boot",      function(n) hg_get(bs, sort_by = "support", top = n),
                              hg_get(bs, sort_by = "support")),
    list("net_hon_compare",   function(n) hg_get(cm, sort_by = "p_adj", top = n),
                              hg_get(cm, sort_by = "p_adj")),
    list("hg_centrality",    function(n) hg_centrality(h, sort_by = "pagerank", top = n),
                              hg_centrality(h, sort_by = "pagerank")),
    list("simplicial_complex",    function(n) hg_get(sc, top = n),
                              hg_get(sc)),
    list("simplicial_complex fv", function(n) hg_get(sc, what = "f_vector", top = n),
                              hg_get(sc, what = "f_vector")),
    list("q_analysis",    function(n) hg_get(qa, top = n),
                              hg_get(qa)),
    list("q_analysis nd", function(n) hg_get(qa, what = "nodes", top = n),
                              hg_get(qa, what = "nodes")),
    list("net_persist_homol", function(n) hg_get(ph, sort_by = "persistence",
                                                        top = n),
                              hg_get(ph, sort_by = "persistence")),
    list("net_pers_landscape",function(n) hg_get(pl, top = n),
                              hg_get(pl)),
    list("net_hg",    function(n) hg_get(hg, sort_by = "weight", top = n),
                              hg_get(hg, sort_by = "weight")),
    list("hypergraph_centr",  function(n) .hg_centrality_fit(hg, sort_by = "clique",
                                                                top = n),
                              .hg_centrality_fit(hg, sort_by = "clique")),
    list("net_hg_measures",   function(n) hg_get(hm, sort_by = "hyperdegree",
                                                        top = n),
                              hg_get(hm, sort_by = "hyperdegree")),
    list("net_hg_cluster",    function(n) hg_get(cl, top = n),
                              hg_get(cl)),
    list("net_hg_transd",     function(n) hg_get(tr, top = n),
                              hg_get(tr)),
    list("net_hg_transd sc",  function(n) hg_get(tr, what = "scores", top = n),
                              hg_get(tr, what = "scores"))
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
  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  full <- hg_get(h, sort_by = "count")
  over_topped <- hg_get(h, sort_by = "count", top = nrow(full) + 100L)
  expect_identical(over_topped, full)
})

test_that("top = rejects a non-whole, zero, negative or vector value", {
  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  expect_error(hg_get(h, top = 0L),      "`top` must be")
  expect_error(hg_get(h, top = -3L),     "`top` must be")
  expect_error(hg_get(h, top = 2.5),     "`top` must be")
  expect_error(hg_get(h, top = c(2L, 3L)), "`top` must be")
  expect_error(hg_get(h, top = "3"),     "`top` must be")
  expect_error(hg_get(h, top = Inf),     "`top` must be")
  # the same guard is shared by every accessor, so one non-hg_get
  # verb is checked too
  expect_error(hg_centrality(h, top = 0L), "`top` must be")
})

test_that("sort_by and top compose rather than fighting", {
  # The regression this exists for: truncating BEFORE sorting would give the
  # top n of construction order, re-sorted -- a different, wrong answer.
  h <- hon(.ac_seqs(), max_order = 2L, min_freq = 20L)
  by_count <- hg_get(h, sort_by = "count", top = 5L)
  expect_identical(by_count$count, sort(by_count$count, decreasing = TRUE))
  all_rules <- hg_get(h)
  expect_identical(max(by_count$count), max(all_rules$count))

  unsorted_then_topped <- hg_get(h, top = 5L)
  expect_false(identical(unsorted_then_topped, by_count))

  # top interacts with a filter too: order_min narrows, then top truncates
  ho <- hg_get(h, order_min = 2L, sort_by = "count", top = 4L)
  expect_true(all(ho$from_order >= 2L))
  expect_identical(nrow(ho), 4L)
  ho_full <- hg_get(h, order_min = 2L, sort_by = "count")
  expect_identical(ho, `rownames<-`(utils::head(ho_full, 4L), NULL))
})
