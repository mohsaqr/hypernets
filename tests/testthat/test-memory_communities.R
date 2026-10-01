# Tests for hg_communities(): the map equation for memory networks
# (Rosvall et al. 2014; Edler, Bohlin & Rosvall 2017).

# ---- fixtures ---------------------------------------------------------------

# A HON-shaped object carrying only the state network, for hand-built
# networks whose flows are known in closed form.
.cm_fake_hon <- function(W) structure(list(weights = W), class = "net_hon")

# Sparse memory network with the flows of Edler et al. (2017), Fig. 4: five
# physical nodes, a red module (a, b and the centre i) and a blue module
# (c, d and i); in the centre the walker switches module at rate r / 2.
# Every state node has flow 1/6 and each module exit flow r / 12.
.cm_edler <- function(r) {
  s <- c("a", "b", "b -> i", "c", "d", "d -> i")
  W <- matrix(0, 6, 6, dimnames = list(s, s))
  W["a", "b"] <- 1
  W["b", "b -> i"] <- 1
  W["b -> i", "a"] <- 1 - r / 2
  W["b -> i", "c"] <- r / 2
  W["c", "d"] <- 1
  W["d", "d -> i"] <- 1
  W["d -> i", "c"] <- 1 - r / 2
  W["d -> i", "a"] <- r / 2
  W
}

.cm_H <- function(p) {
  p <- p / sum(p)
  -sum(p * log2(p))
}

.cm_toy_seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
                     c("c", "h", "d", "c", "h", "d", "c"),
                     c("a", "h", "b", "a", "h", "b"),
                     c("c", "h", "d", "c", "h", "d"))

# Planted memory modules: a ring of K cliques, neighbouring cliques share one
# node. Inside a clique the walker moves uniformly; on reaching a shared node
# it switches clique with probability eps. A second-order state "u -> v"
# knows its clique, a first-order view of a shared node does not.
.cm_sim_ring <- function(n_seq = 200L, len = 40L, K = 4L, m = 4L,
                         eps = 0.05, seed = 1L) {
  set.seed(seed)
  shared <- sprintf("s%d", seq_len(K))
  cliques <- lapply(seq_len(K), function(k) {
    c(shared[k], sprintf("p%d_%d", k, seq_len(m - 2L)), shared[k %% K + 1L])
  })
  lapply(seq_len(n_seq), function(i) {
    cm <- sample.int(K, 1L)
    start <- list(v = sample(cliques[[cm]], 1L), cm = cm)
    # a simulated walk is sequential: the steps are carried by Reduce()
    steps <- Reduce(function(state, t) {
      v <- sample(setdiff(cliques[[state$cm]], state$v), 1L)
      other <- setdiff(which(vapply(cliques, function(cl) v %in% cl,
                                    logical(1L))), state$cm)
      switch_now <- length(other) > 0L && stats::runif(1L) < eps
      list(v = v, cm = if (switch_now) other else state$cm)
    }, seq.int(2L, len), start, accumulate = TRUE)
    vapply(steps, function(state) state$v, character(1L))
  })
}

.cm_ring_truth <- function(states, K = 4L, m = 4L) {
  shared <- sprintf("s%d", seq_len(K))
  cliques <- lapply(seq_len(K), function(k) {
    c(shared[k], sprintf("p%d_%d", k, seq_len(m - 2L)), shared[k %% K + 1L])
  })
  vapply(strsplit(states, " -> ", fixed = TRUE), function(p) {
    last2 <- utils::tail(p, 2L)
    ks <- which(vapply(cliques, function(cl) all(last2 %in% cl), logical(1L)))
    if (length(ks) == 1L) ks else NA_integer_
  }, integer(1L))
}

# All set partitions of n items as restricted growth strings.
.cm_set_partitions <- function(n) {
  grow <- function(prefix) {
    if (length(prefix) == n) return(list(prefix))
    unlist(lapply(seq_len(max(prefix) + 1L), function(k) grow(c(prefix, k))),
           recursive = FALSE)
  }
  grow(1L)
}

# ---- codelength: hand-computed values from the paper ------------------------

test_that("codelength reproduces Edler et al. (2017) Eq. 31 and Eq. 32", {
  for (r in c(1, 0.5, 0.1)) {
    hon <- .cm_fake_hon(.cm_edler(r))
    one <- hg_communities(hon, partition = c(a = 1, b = 1, "b -> i" = 1,
                                              c = 1, d = 1, "d -> i" = 1))
    two <- hg_communities(hon, partition = c(a = 1, b = 1, "b -> i" = 1,
                                              c = 2, d = 2, "d -> i" = 2))
    cl1 <- hg_get(one, what = "codelength")
    cl2 <- hg_get(two, what = "codelength")
    # Eq. 31: shared codewords for the centre give H(1/6,1/6,1/6,1/6,2/6)
    eq31 <- .cm_H(c(1, 1, 1, 1, 2) / 6)
    # Eq. 32: (r/6) H(r/12, r/12) + 2 (6 + r)/12 H(1/6, 1/6, 1/6, r/12)
    eq32 <- (r / 6) * .cm_H(c(r, r) / 12) +
      2 * (6 + r) / 12 * .cm_H(c(1 / 6, 1 / 6, 1 / 6, r / 12))
    expect_equal(cl1$codelength[1L], eq31, tolerance = 1e-12)
    expect_equal(cl2$codelength[1L], eq32, tolerance = 1e-12)
  }
  # the paper's printed approximations
  expect_equal(round(.cm_H(c(1, 1, 1, 1, 2) / 6), 2), 2.25)
  two_r1 <- hg_communities(.cm_fake_hon(.cm_edler(1)),
                            partition = c(a = 1, b = 1, "b -> i" = 1,
                                          c = 2, d = 2, "d -> i" = 2))
  expect_equal(round(hg_get(two_r1, what = "codelength")$codelength[1L],
                     2), 2.44)
})

test_that("state nodes of one physical node share a codeword (Eq. 29 vs 31)", {
  W <- .cm_edler(1)
  fl <- hypernets:::.hcm_flow(unname(W), 0.15)
  expect_equal(fl$node_flow, rep(1 / 6, 6), tolerance = 1e-12)
  # six virtual physical nodes: Eq. 29, log2(6) ~ 2.58
  virt <- hypernets:::.hcm_codelength(fl, seq_len(6), rep(1L, 6))
  expect_equal(virt$codelength, log2(6), tolerance = 1e-12)
  # proper physical nodes (b -> i and d -> i are both i): Eq. 31
  phys <- c(1L, 2L, 3L, 4L, 5L, 3L)
  real <- hypernets:::.hcm_codelength(fl, phys, rep(1L, 6))
  expect_equal(real$codelength, .cm_H(c(1, 1, 1, 1, 2) / 6), tolerance = 1e-12)
  expect_lt(real$codelength, virt$codelength)
  # two modules: Eq. 30 equals Eq. 32 because no module holds two state
  # nodes of one physical node
  two_virt <- hypernets:::.hcm_codelength(fl, seq_len(6), rep(1:2, each = 3))
  two_real <- hypernets:::.hcm_codelength(fl, phys, rep(1:2, each = 3))
  expect_equal(two_virt$codelength, two_real$codelength, tolerance = 1e-12)
})

test_that("first-order map equation matches Rosvall & Bergstrom (2008) by hand", {
  # three directed 3-cycles joined asymmetrically, so that under unrecorded
  # teleportation the module enter flows are not a permutation of the exit
  # flows (they are with two modules, where exit_1 = enter_2)
  W <- matrix(0, 9, 9)
  W[cbind(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 3, 6, 9, 1, 4),
          c(2, 3, 1, 5, 6, 4, 8, 9, 7, 4, 7, 1, 7, 1))] <-
    c(rep(1, 9), 1, 1, 1, 3, 2)
  fl <- hypernets:::.hcm_flow(W, 0.15)
  m <- rep(1:3, each = 3)
  cl <- hypernets:::.hcm_codelength(fl, seq_len(9), m)
  # expanded map equation, recomputed independently from the link flows
  lf <- matrix(0, 9, 9)
  lf[cbind(fl$from, fl$to)] <- fl$flow
  q_exit <- vapply(1:3, function(k) sum(lf[m == k, m != k]), numeric(1))
  q_enter <- vapply(1:3, function(k) sum(lf[m != k, m == k]), numeric(1))
  plogp <- function(x) ifelse(x > 0, x * log2(x), 0)
  expect_gt(abs(sum(plogp(q_enter)) - sum(plogp(q_exit))), 1e-4)
  p <- fl$node_flow
  p_mod <- vapply(1:3, function(k) sum(p[m == k]), numeric(1))
  hand <- plogp(sum(q_enter)) - sum(plogp(q_enter)) +
    sum(plogp(q_exit + p_mod)) - sum(plogp(q_exit)) - sum(plogp(p))
  expect_equal(cl$codelength, hand, tolerance = 1e-12)
})

# ---- flow and invariants ----------------------------------------------------

test_that("flows are distributions and the one-module code is H(physical flow)", {
  hon <- hon(.cm_sim_ring(n_seq = 60L, seed = 2L), max_order = 2L)
  fl <- hypernets:::.hcm_flow(unname(hon$weights), 0.15)
  expect_equal(sum(fl$flow), 1, tolerance = 1e-12)
  expect_equal(sum(fl$node_flow), 1, tolerance = 1e-12)
  expect_true(all(fl$flow >= 0))
  comm <- hg_communities(hon, partition = stats::setNames(
    rep(1L, nrow(hon$weights)), rownames(hon$weights)))
  phys <- hg_get(comm, what = "physical")
  cl <- hg_get(comm, what = "codelength")
  expect_equal(cl$codelength[1L], .cm_H(phys$flow), tolerance = 1e-12)
  expect_equal(cl$one_level_codelength[1L], .cm_H(phys$flow),
               tolerance = 1e-12)
})

test_that("relabelling modules or reordering states leaves the codelength", {
  hon <- hon(.cm_sim_ring(n_seq = 60L, seed = 3L), max_order = 2L)
  fl <- hypernets:::.hcm_flow(unname(hon$weights), 0.15)
  states <- rownames(hon$weights)
  phys <- match(sub("^.* -> ", "", states), unique(sub("^.* -> ", "", states)))
  set.seed(9)
  m <- sample.int(4L, length(states), replace = TRUE)
  base <- hypernets:::.hcm_codelength(fl, phys, m)$codelength
  relab <- hypernets:::.hcm_codelength(fl, phys, c(7, 3, 11, 5)[m])$codelength
  expect_equal(relab, base, tolerance = 1e-12)
  perm <- sample.int(length(states))
  Wp <- unname(hon$weights)[perm, perm]
  flp <- hypernets:::.hcm_flow(Wp, 0.15)
  expect_equal(flp$node_flow, fl$node_flow[perm], tolerance = 1e-10)
  expect_equal(hypernets:::.hcm_codelength(flp, phys[perm], m[perm])$codelength,
               base, tolerance = 1e-10)
})

# ---- search -----------------------------------------------------------------

test_that("the search finds the brute-force optimum on a small memory network", {
  hon <- hon(.cm_toy_seqs, max_order = 2L)
  states <- rownames(hon$weights)
  phys <- match(sub("^.* -> ", "", states), unique(sub("^.* -> ", "", states)))
  fl <- hypernets:::.hcm_flow(unname(hon$weights), 0.15)
  parts <- .cm_set_partitions(length(states))
  expect_length(parts, 877L)  # Bell(7)
  brute <- min(vapply(parts, function(p) {
    hypernets:::.hcm_codelength(fl, phys, p)$codelength
  }, numeric(1L)))
  comm <- hg_communities(hon, trials = 5L)
  found <- hg_get(comm, what = "codelength")$codelength[1L]
  expect_equal(found, brute, tolerance = 1e-10)
})

test_that("the bundled ring data are the planted-module simulation", {
  expect_identical(ring_sequences, .cm_sim_ring(seed = 1L))
  # the bundled truth agrees with the rule it states, on every state it lists
  expect_identical(ring_communities$community,
                   .cm_ring_truth(ring_communities$node))
  expect_false(anyNA(ring_communities$community))
})

test_that("the planted groups of ring_sequences are recovered", {
  skip_on_cran()
  comm <- hg_communities(hon(ring_sequences, max_order = 2L),
                          trials = 10L, seed = 1L)
  agree <- hg_agreement(hg_get(comm), ring_communities, node = "node",
                        label = "community")
  expect_gt(agree$ari, 0.95)
  codes <- hg_get(comm, what = "codelength")
  expect_lt(codes$codelength[1L], codes$codelength[2L] - 1)
})

test_that("planted memory modules are recovered where first order fails", {
  skip_on_cran()
  # five simulated data sets, each searched from its own seed
  res <- lapply(1:5, function(s) {
    hon <- hon(.cm_sim_ring(seed = s), max_order = 2L)
    comm <- hg_communities(hon, trials = 5L, seed = s)
    st <- hg_get(comm)
    truth <- .cm_ring_truth(st$node)
    known <- !is.na(truth)
    cl <- hg_get(comm, what = "codelength")
    data.frame(memory = cl$codelength[1L], first = cl$codelength[2L],
               k_first = cl$n_communities[2L],
               ari = hypernets:::.thg_ari(st$community[known], truth[known]))
  })
  res <- do.call(rbind, res)
  expect_true(all(res$memory < res$first - 1))
  expect_true(all(res$k_first == 1L))
  expect_true(all(res$ari > 0.95))
})

test_that("overlap: shared physical nodes sit in two modules", {
  skip_on_cran()
  hon <- hon(.cm_sim_ring(seed = 1L), max_order = 2L)
  comm <- hg_communities(hon, trials = 3L)
  over <- hg_get(comm, what = "physical", overlapping = TRUE)
  expect_setequal(unique(over$state), c("s1", "s2", "s3", "s4"))
  expect_true(all(over$n_communities == 2L))
  shares <- tapply(over$share, over$state, sum)
  expect_equal(as.numeric(shares), rep(1, 4), tolerance = 1e-12)
})

test_that("trials are seeded, reproducible, and restore the caller's RNG", {
  hon <- hon(.cm_toy_seqs, max_order = 2L)
  set.seed(42)
  before <- .Random.seed
  a <- hg_communities(hon, trials = 3L, seed = 7L)
  expect_identical(.Random.seed, before)
  b <- hg_communities(hon, trials = 3L, seed = 7L)
  expect_identical(hg_get(a), hg_get(b))
  tr <- hg_get(a, what = "trials")
  expect_named(tr, c("run", "codelength", "n_communities", "ari_to_best",
                   "best"))
  expect_equal(nrow(tr), 3L)
  expect_equal(sum(tr$best), 1L)
})

# ---- accessors and methods --------------------------------------------------

test_that("as.data.frame returns tidy tables for every `what`", {
  hon <- hon(.cm_toy_seqs, max_order = 2L)
  comm <- hg_communities(hon, trials = 2L)
  expect_named(hg_get(comm), c("node", "state", "community", "flow"))
  expect_named(hg_get(comm, what = "physical"),
               c("state", "community", "flow", "share", "n_communities"))
  expect_named(hg_get(comm, what = "modules"),
               c("community", "n_nodes", "n_states", "flow", "exit_flow",
                 "enter_flow"))
  expect_named(hg_get(comm, what = "first_order"),
               c("state", "community", "flow"))
  cl <- hg_get(comm, what = "codelength")
  expect_equal(cl$model, c("memory", "first_order"))
  mods <- hg_get(comm, what = "modules")
  expect_equal(sum(mods$flow), 1, tolerance = 1e-12)
  expect_equal(sum(mods$exit_flow), sum(mods$enter_flow), tolerance = 1e-12)
  one <- hg_get(comm, community = 1L)
  expect_true(all(one$community == 1L))
  expect_equal(nrow(hg_get(comm, what = "physical", overlapping = TRUE)),
               2L)
})

test_that("print and summary are stable", {
  hon <- hon(.cm_toy_seqs, max_order = 2L)
  comm <- hg_communities(hon, trials = 2L)
  expect_snapshot(print(comm))
  out <- summary(comm)
  expect_s3_class(out, "hypernets_summary")
  expect_identical(out$codelength, hg_get(comm, what = "codelength"))
})

# Physical-node centres as drawn: the mean of each node's disc outline, in
# the hypergraph's node order.
.cm_drawn_nodes <- function(p, nodes) {
  disc <- Filter(function(l) inherits(l$geom, "GeomPolygon") &&
                   identical(l$aes_params$fill, "#000000"), p$layers)[[1L]]$data
  data.frame(node = nodes,
             x = as.numeric(tapply(disc$x, disc$id, mean)),
             y = as.numeric(tapply(disc$y, disc$id, mean)))
}

test_that("the physical view is a community hypergraph with flow circles", {
  skip_on_cran()
  hon <- hon(.cm_sim_ring(seed = 1L), max_order = 2L)
  comm <- hg_communities(hon, trials = 3L)
  spec <- hypernets:::.hcm_physical_plot_data(comm)
  phys <- hg_get(comm, what = "physical")
  expect_setequal(spec$hypergraph$nodes, unique(phys$state))
  expect_setequal(colnames(spec$hypergraph$incidence),
                  sprintf("Community %d", 1:4))
  expect_true(is.factor(spec$community))
  expect_identical(levels(spec$community), sprintf("Community %d", 1:4))
  flow <- tapply(phys$flow, phys$state, sum)
  expect_equal(spec$sizes$value, as.numeric(flow[spec$sizes$node]))
  W1 <- comm$physical_weights
  expect_equal(nrow(spec$moves), sum(W1 > 0))
  expect_equal(sum(spec$moves$weight), sum(W1))

  p <- plot(comm)
  expect_s3_class(p, "ggplot")
  built <- ggplot2::ggplot_build(p)
  expect_identical(p$labels$title, "4 communities, 4 shared physical nodes")
  expect_match(p$labels$caption, "Circle area: physical flow", fixed = TRUE)
  expect_match(p$labels$caption, "Triangle: points to the event", fixed = TRUE)
  expect_match(p$labels$caption, "4 zero-flow states not drawn", fixed = TRUE)
  fill_scale <- p$scales$get_scales("fill")
  expect_true(fill_scale$is_discrete())
  expect_identical(fill_scale$name, "Community")
  # the blob look: a title box per community with its flow, legend below,
  # haloed labels, wide margins
  expect_identical(p$theme$legend.position, "bottom")
  expect_equal(as.numeric(p$theme$plot.margin), c(40, 130, 30, 130))
  # the key one sentence to a line, left-aligned under the whole plot
  expect_identical(strsplit(p$labels$caption, "\n", fixed = TRUE)[[1L]][2L],
                   "Triangle: points to the event that most often follows it.")
  expect_identical(p$theme$plot.caption.position, "plot")
  boxes <- Filter(function(l) inherits(l$geom, "GeomLabel"), p$layers)
  expect_length(boxes, 1L)
  module_flow <- tapply(comm$states$flow, comm$states$community, sum)
  # zero-flow modules hold no physical node and are not drawn
  module_flow <- module_flow[module_flow > 0]
  expect_setequal(boxes[[1L]]$data$label,
                  sprintf("Community %d\n%s flow", as.integer(names(module_flow)),
                          hypernets:::.thg_comma(as.numeric(module_flow), count = TRUE)))
  expect_identical(spec$flow,
                   stats::setNames(as.numeric(module_flow),
                                   sprintf("Community %s", names(module_flow))))
  texts <- Filter(function(l) inherits(l$geom, "GeomText"), p$layers)
  expect_length(texts, 9L)

  # every shared node lies inside the pebble of each of its communities; a
  # single-community node's pebble is its community's Okabe-Ito colour
  hulls <- p$layers[[1L]]$data
  hull_fill <- built$data[[1L]]
  drawn <- .cm_drawn_nodes(p, spec$hypergraph$nodes)
  palette <- stats::setNames(hypernets:::.thg_okabe_ito[1:4],
                             sprintf("Community %d", 1:4))
  inside <- function(v, k) {
    ring <- hulls[hulls$hyperedge == sprintf("Community %d", k), ]
    hypernets:::.thg_in_convex(ring$x, ring$y, drawn$x[drawn$node == v],
                               drawn$y[drawn$node == v])
  }
  shared <- phys[phys$n_communities > 1L, ]
  expect_equal(length(unique(shared$state)), 4L)
  expect_true(all(mapply(inside, shared$state, shared$community)))
  single <- phys[phys$n_communities == 1L, ]
  expect_true(all(mapply(inside, single$state, single$community)))
  group_of <- levels(hulls$hyperedge)[hull_fill$group]
  pebble_colour <- vapply(split(hull_fill$fill, group_of), unique,
                          character(1L))
  expect_identical(unname(pebble_colour[sprintf("Community %d", single$community)]),
                   unname(palette[sprintf("Community %d", single$community)]))

  # `...` reaches plot.net_hg()
  outside <- plot(comm, arrow_style = "outside", seed = 2L)
  expect_s3_class(outside, "ggplot")
  expect_error(plot(comm, node_fill = "not a colour"),
               class = "hypernets_bad_input")
})

test_that("the physical view of human_long names its one-node community", {
  skip_on_cran()
  # One trajectory per session, ordered within it (Nestimate's hon()
  # reads sequences, not a long table).
  ord <- order(human_long$session_id, human_long$order_in_session)
  seqs <- split(human_long$code[ord], human_long$session_id[ord])
  human <- hon(seqs, max_order = 2L, min_freq = 10L)
  comm <- hg_communities(human, trials = 3L)
  p <- plot(comm)
  expect_s3_class(p, "ggplot")
  lone <- hg_get(comm, what = "modules")
  lone <- lone$community[lone$n_physical == 1L & lone$flow > 0]
  if (length(lone)) {
    expect_match(p$labels$caption,
                 sprintf("Community %d (", lone[1L]), fixed = TRUE)
    expect_false(sprintf("Community %d", lone[1L]) %in%
                   levels(hypernets:::.hcm_physical_plot_data(comm)$community))
  }
  expect_no_error(ggplot2::ggplot_build(p))
})

test_that("the state view still delegates to cograph::overlay_communities", {
  skip_on_cran()
  hon <- hon(.cm_sim_ring(seed = 1L), max_order = 2L)
  comm <- hg_communities(hon, trials = 3L)
  okabe <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
             "#D55E00", "#CC79A7", "#999999")
  st <- hypernets:::.hcm_state_plot_args(comm, FALSE)
  zero <- hg_get(comm)$flow <= 0
  expect_equal(nrow(st$args$x), sum(!zero))
  expect_length(st$args$communities, 4L)
  expect_false(any(grepl("[", st$args$labels, fixed = TRUE)))
  expect_true(all(st$args$node_fill %in% okabe))
  expect_equal(st$n_hidden, sum(zero))
  all_states <- hypernets:::.hcm_state_plot_args(comm, TRUE)
  expect_equal(nrow(all_states$args$x), nrow(hg_get(comm)))
  expect_length(all_states$args$communities, 8L)
  expect_equal(all_states$n_hidden, 0L)

  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off(), add = TRUE)
  expect_identical(plot(comm, type = "states"), comm)
  expect_identical(plot(comm, type = "states", show_zero_flow = TRUE), comm)
  expect_error(plot(comm, show_zero_flow = NA))
  seen <- NULL
  testthat::local_mocked_bindings(
    overlay_communities = function(x, communities, ...) {
      seen <<- list(...)
      invisible(NULL)
    },
    .package = "cograph")
  plot(comm, type = "states")
  expect_equal(seen$title, "4 communities of state nodes")
})

# ---- error paths ------------------------------------------------------------

test_that("bad input raises classed conditions", {
  hon <- hon(.cm_toy_seqs, max_order = 2L)
  expect_error(hg_communities(list()), class = "hypernets_bad_input")
  expect_error(hg_communities(hon, trials = 0), class = "hypernets_bad_input")
  expect_error(hg_communities(hon, trials = 1.5), class = "hypernets_bad_input")
  expect_error(hg_communities(hon, teleportation = 1),
               class = "hypernets_bad_input")
  expect_error(hg_communities(hon, seed = NA), class = "hypernets_bad_input")
  expect_error(hg_communities(hon, partition = c(a = 1)),
               class = "hypernets_bad_input")
  expect_error(hg_communities(hon, partition = 1:7),
               class = "hypernets_bad_input")
  expect_error(hg_communities(hon, partition = data.frame(x = 1)),
               class = "hypernets_bad_input")
  expect_error(plot(hg_communities(hon, trials = 1L), type = "nope"))
})

test_that("teleportation = 0 on a reducible walk is hypernets_not_ergodic", {
  # two disconnected cycles: no unique stationary flow without teleportation
  s <- c("a", "b", "c", "d")
  W <- matrix(0, 4, 4, dimnames = list(s, s))
  W["a", "b"] <- W["b", "a"] <- W["c", "d"] <- W["d", "c"] <- 1
  expect_error(hg_communities(.cm_fake_hon(W), teleportation = 0),
               class = "hypernets_not_ergodic")
  # with teleportation the flow exists and splits into the two cycles
  comm <- hg_communities(.cm_fake_hon(W), trials = 2L)
  expect_equal(nrow(hg_get(comm, what = "modules")), 2L)
})
