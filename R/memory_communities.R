# ---- Memory-network communities: the map equation on state nodes ----------
#
# Community detection on a higher-order network with the map equation for
# memory networks (Rosvall, Esquivel, Lancichinetti, West & Lambiotte 2014;
# Edler, Bohlin & Rosvall 2017). A HON state "a -> b" is a STATE node of the
# PHYSICAL node b (the last element; same convention as memory_centrality.R).
# State nodes are clustered, but the module codebooks are written over
# physical nodes: state nodes of the same physical node in the same module
# share one codeword (Edler et al. 2017, Eq. 25-26). A physical node whose
# state nodes fall in different modules belongs to all of them, so the
# partition of physical nodes overlaps.
#
# Everything below is written out from the papers; no code is ported.
#
# Flow (Edler et al. 2017, Eq. 8-10, with the teleportation of Lambiotte &
# Rosvall 2012 as implemented in Infomap): the walker follows a link with
# probability 1 - tau and otherwise (and always from a node without
# out-links) restarts at the tail of a link drawn in proportion to its
# weight. Teleportation is UNRECORDED: the map equation encodes only the link
# steps, so link flow is f_uv = pi_u w_uv / w_u renormalised to sum 1, and a
# state's flow is its link in-flow. This is exactly what Infomap 2.15.1
# computes for a directed network (local oracle:
# local_testing_and_equivalence/test-equiv-hon-communities-infomap.R).
#
# Codelength (Edler et al. 2017, Eq. 11, 12, 20-22, 24-26), expanded with
# plogp(x) = x log2 x:
#   L = plogp(sum_m q_m_enter) - sum_m plogp(q_m_enter)
#     + sum_m [ plogp(q_m_exit + sum_{a in m} pi_a) - plogp(q_m_exit)
#               - sum_{i phys} plogp(pi_{i cap m}) ]
#
# Search (Edler et al. 2017, Algorithms 3-5): repeated node aggregation with
# local moves to neighbouring or new modules in random sequential order, then
# fine-tuning (re-running the aggregation from single state nodes started in
# their current modules) until the codelength stops decreasing. Coarse-tuning
# (Algorithm 6) is not implemented; more trials compensate.
#
# The node moves are sequential by construction (each move changes the
# module flows the next move is scored against), so a sweep carries the search
# state through the random node order with Reduce(); the sweeps, levels and
# fine-tuning rounds repeat until the codelength stops falling, each bounded.

.hcm_plogp <- function(x) {
  x <- pmax(x, 0)
  out <- x * log2(x)
  out[x <= 0] <- 0
  out
}

#' Raise a classed bad-input error for the communities verb
#' @noRd
.hcm_bad_input <- function(msg) {
  stop(errorCondition(msg, class = "hypernets_bad_input", call = NULL))
}

#' Flow on a directed state network (Infomap directed flow model)
#'
#' @param W Square non-negative weight matrix (rows = from).
#' @param tau Teleportation probability in [0, 1).
#' @return list(node_flow, from, to, flow): node flow (link in-flow) and the
#'   link flows of every non-zero link, summing to 1.
#' @noRd
.hcm_flow <- function(W, tau) {
  n <- nrow(W)
  w_out <- rowSums(W)
  if (sum(w_out) <= 0) {
    .hcm_bad_input("the network has no links, so there is no flow to cluster")
  }
  dangling <- w_out <= 0
  P <- W
  P[!dangling, ] <- W[!dangling, , drop = FALSE] / w_out[!dangling]
  v <- w_out / sum(w_out)
  # pi = pi G with G = (1 - tau) P + (tau * [non-dangling] + [dangling]) v'
  restart <- ifelse(dangling, 1, tau)
  G <- (1 - tau) * P + outer(restart, v)
  A <- t(G) - diag(n)
  A[n, ] <- 1
  b <- c(rep(0, n - 1L), 1)
  qa <- qr(A)
  if (qa$rank < n) {
    stop(errorCondition(
      paste0("the walk on this network has no unique stationary flow ",
             "(it is reducible); use `teleportation` > 0"),
      class = c("hypernets_not_ergodic", "hypernets_bad_input"), call = NULL))
  }
  pi <- qr.solve(qa, b)
  pi <- pmax(pi, 0)
  pi <- pi / sum(pi)
  ij <- which(W > 0, arr.ind = TRUE)
  f <- pi[ij[, 1L]] * P[ij]
  f <- f / sum(f)
  node_flow <- numeric(n)
  inflow <- rowsum(f, ij[, 2L])
  node_flow[as.integer(rownames(inflow))] <- inflow[, 1L]
  list(node_flow = node_flow, from = unname(ij[, 1L]), to = unname(ij[, 2L]),
       flow = f)
}

#' Two-level map-equation codelength of a state-node partition
#'
#' @param fl Output of .hcm_flow().
#' @param phys Integer physical id of every state node (1..n_phys).
#' @param module Integer module of every state node.
#' @return list(codelength, index, module) in bits.
#' @noRd
.hcm_codelength <- function(fl, phys, module) {
  mods <- sort(unique(module))
  m <- match(module, mods)
  k <- length(mods)
  cross <- m[fl$from] != m[fl$to]
  exit <- vapply(seq_len(k), function(j) sum(fl$flow[cross & m[fl$from] == j]),
                 numeric(1L))
  enter <- vapply(seq_len(k), function(j) sum(fl$flow[cross & m[fl$to] == j]),
                  numeric(1L))
  mod_flow <- vapply(seq_len(k), function(j) sum(fl$node_flow[m == j]),
                     numeric(1L))
  # physical flow inside each module: Edler et al. (2017) Eq. 25
  phys_in_mod <- tapply(fl$node_flow, list(m, phys), sum)
  phys_in_mod[is.na(phys_in_mod)] <- 0
  index <- .hcm_plogp(sum(enter)) - sum(.hcm_plogp(enter))
  modcl <- sum(.hcm_plogp(exit + mod_flow)) - sum(.hcm_plogp(exit)) -
    sum(.hcm_plogp(phys_in_mod))
  list(codelength = index + modcl, index = index, module = modcl)
}

#' Build a search level from state nodes or aggregated modules
#'
#' @param from,to,flow Link flows between level nodes.
#' @param node_flow Flow of each level node.
#' @param node_phys Dense level-node x physical flow matrix.
#' @noRd
.hcm_level <- function(from, to, flow, node_flow, node_phys) {
  n <- length(node_flow)
  keep <- from != to
  from <- from[keep]
  to <- to[keep]
  flow <- flow[keep]
  out_nb <- split(to, factor(from, levels = seq_len(n)))
  out_fl <- split(flow, factor(from, levels = seq_len(n)))
  in_nb <- split(from, factor(to, levels = seq_len(n)))
  in_fl <- split(flow, factor(to, levels = seq_len(n)))
  list(n = n, from = from, to = to, flow = flow, node_flow = node_flow,
       node_phys = node_phys,
       out_nb = lapply(out_nb, as.integer), out_fl = unname(out_fl),
       in_nb = lapply(in_nb, as.integer), in_fl = unname(in_fl),
       out_tot = vapply(out_fl, sum, numeric(1L)),
       in_tot = vapply(in_fl, sum, numeric(1L)))
}

#' Algorithm 4 (optimize) on one level: local moves until no improvement
#'
#' @param lv A level from .hcm_level().
#' @param module Initial integer module per level node (ids in 1..n).
#' @param eps Minimum codelength decrease that counts as a move.
#' @param max_sweeps Bound on sweeps.
#' @return Integer module per level node.
#' @noRd
.hcm_optimize <- function(lv, module, eps = 1e-10, max_sweeps = 100L) {
  n <- lv$n
  cross <- module[lv$from] != module[lv$to]
  exit <- numeric(n)
  enter <- numeric(n)
  ex <- tapply(lv$flow[cross], module[lv$from][cross], sum)
  en <- tapply(lv$flow[cross], module[lv$to][cross], sum)
  exit[as.integer(names(ex))] <- ex
  enter[as.integer(names(en))] <- en
  mflow <- numeric(n)
  mf <- tapply(lv$node_flow, module, sum)
  mflow[as.integer(names(mf))] <- mf
  pm <- matrix(0, nrow = n, ncol = ncol(lv$node_phys))
  pm_used <- rowsum(lv$node_phys, module)
  pm[as.integer(rownames(pm_used)), ] <- pm_used
  state <- list(module = module, exit = exit, enter = enter, mflow = mflow,
                pm = pm, size = tabulate(module, nbins = n),
                total_enter = sum(enter), moved = 0L)
  sweeps <- 0L
  # sweeps until no node moves (bounded by `max_sweeps`); within a sweep the
  # moves are sequential by construction -- each changes the module flows the
  # next node is scored against (Edler et al. 2017, Algorithm 4) -- so the
  # module state is carried through the random order by Reduce()
  repeat {
    sweeps <- sweeps + 1L
    state$moved <- 0L
    state <- Reduce(function(s, x) .hcm_move(lv, s, x, eps), sample.int(n),
                    state)
    if (state$moved == 0L || sweeps >= max_sweeps) break
  }
  state$module
}

#' One local move: node `x` to the neighbouring or new module that most
#' shortens the codelength, if any does by more than `eps`
#'
#' @param lv A level from .hcm_level().
#' @param s Search state: module, exit, enter, mflow, pm (module x physical
#'   flow), size, total_enter, moved.
#' @param x The level node to move.
#' @return The updated state.
#' @noRd
.hcm_move <- function(lv, s, x, eps) {
  a <- s$module[x]
  cand <- unique(s$module[c(lv$out_nb[[x]], lv$in_nb[[x]])])
  cand <- cand[cand != a]
  if (s$size[a] > 1L) {
    empty <- which(s$size == 0L)
    if (length(empty)) cand <- c(cand, empty[1L])
  }
  if (!length(cand)) return(s)
  pl <- .hcm_plogp
  om <- s$module[lv$out_nb[[x]]]
  im <- s$module[lv$in_nb[[x]]]
  o_to <- function(mm) vapply(mm, function(j) sum(lv$out_fl[[x]][om == j]),
                              numeric(1L))
  i_from <- function(mm) vapply(mm, function(j) sum(lv$in_fl[[x]][im == j]),
                                numeric(1L))
  ox <- lv$out_tot[x]
  ix <- lv$in_tot[x]
  fx <- lv$node_flow[x]
  oA <- o_to(a)
  iA <- i_from(a)
  oB <- o_to(cand)
  iB <- i_from(cand)
  exA2 <- s$exit[a] - (ox - oA) + iA
  enA2 <- s$enter[a] - (ix - iA) + oA
  exB2 <- s$exit[cand] + (ox - oB) - iB
  enB2 <- s$enter[cand] + (ix - iB) - oB
  flA2 <- s$mflow[a] - fx
  flB2 <- s$mflow[cand] + fx
  te2 <- s$total_enter - s$enter[a] - s$enter[cand] + enA2 + enB2
  ph <- which(lv$node_phys[x, ] > 0)
  px <- lv$node_phys[x, ph]
  pA <- s$pm[a, ph]
  pB <- s$pm[cand, ph, drop = FALSE]
  d_phys_A <- sum(pl(pA - px) - pl(pA))
  d_phys_B <- rowSums(pl(sweep(pB, 2L, px, `+`)) - pl(pB))
  delta <- (pl(te2) - pl(s$total_enter)) -
    (pl(enA2) + pl(enB2) - pl(s$enter[a]) - pl(s$enter[cand])) +
    (pl(exA2 + flA2) - pl(s$exit[a] + s$mflow[a]) +
       pl(exB2 + flB2) - pl(s$exit[cand] + s$mflow[cand])) -
    (pl(exA2) - pl(s$exit[a]) + pl(exB2) - pl(s$exit[cand])) -
    (d_phys_A + d_phys_B)
  best <- which.min(delta)
  if (delta[best] >= -eps) return(s)
  b <- cand[best]
  s$total_enter <- te2[best]
  s$exit[c(a, b)] <- c(exA2, exB2[best])
  s$enter[c(a, b)] <- c(enA2, enB2[best])
  s$mflow[c(a, b)] <- c(flA2, flB2[best])
  s$pm[a, ph] <- pmax(pA - px, 0)
  s$pm[b, ph] <- pB[best, ] + px
  s$size[c(a, b)] <- s$size[c(a, b)] + c(-1L, 1L)
  s$module[x] <- b
  s$moved <- s$moved + 1L
  s
}

#' Algorithm 3: repeated node aggregation from an initial state partition
#'
#' @param fl Flow list; @param phys physical id per state; @param n_phys.
#' @param init Integer module per state node (ids in 1..n).
#' @return Integer module per state node.
#' @noRd
.hcm_aggregate <- function(fl, phys, n_phys, init, eps = 1e-10,
                           max_levels = 50L) {
  n <- length(fl$node_flow)
  node_phys <- matrix(0, n, n_phys)
  node_phys[cbind(seq_len(n), phys)] <- fl$node_flow
  # state -> current level node
  map <- seq_len(n)
  lv <- .hcm_level(fl$from, fl$to, fl$flow, fl$node_flow, node_phys)
  module <- init
  l_old <- .hcm_codelength(fl, phys, module[map])$codelength
  levels_done <- 0L
  repeat {
    levels_done <- levels_done + 1L
    module <- .hcm_optimize(lv, module, eps = eps)
    state_mod <- module[map]
    l_new <- .hcm_codelength(fl, phys, state_mod)$codelength
    k <- length(unique(module))
    if (l_new - l_old > -eps || k == 1L || k == lv$n ||
        levels_done >= max_levels) {
      break
    }
    l_old <- l_new
    # rebuild the network with modules as nodes
    ids <- sort(unique(module))
    new <- match(module, ids)
    map <- new[map]
    agg <- stats::aggregate(
      lv$flow, by = list(from = new[lv$from], to = new[lv$to]), FUN = sum)
    np <- rowsum(lv$node_phys, new, reorder = TRUE)
    nf <- as.numeric(rowsum(lv$node_flow, new, reorder = TRUE))
    lv <- .hcm_level(agg$from, agg$to, agg$x, nf, matrix(np, ncol = n_phys))
    module <- seq_len(lv$n)
  }
  module[map]
}

#' One trial: aggregation from singletons, then fine-tuning (Algorithms 2-5)
#' @noRd
.hcm_trial <- function(fl, phys, n_phys, eps = 1e-10, max_tune = 10L) {
  n <- length(fl$node_flow)
  mod <- .hcm_aggregate(fl, phys, n_phys, seq_len(n), eps = eps)
  l_old <- .hcm_codelength(fl, phys, mod)$codelength
  tunes <- 0L
  repeat {
    tunes <- tunes + 1L
    init <- match(mod, sort(unique(mod)))
    cand <- .hcm_aggregate(fl, phys, n_phys, init, eps = eps)
    l_new <- .hcm_codelength(fl, phys, cand)$codelength
    if (l_new < l_old) mod <- cand
    if (l_new - l_old > -eps || tunes >= max_tune) break
    l_old <- l_new
  }
  mod
}

#' Relabel modules 1..M by decreasing flow (ties: first state)
#' @noRd
.hcm_relabel <- function(module, node_flow) {
  mf <- tapply(node_flow, module, sum)
  first <- tapply(seq_along(module), module, min)
  ids <- names(mf)[order(-mf, first)]
  match(as.character(module), ids)
}

#' Search: several trials from seeded random orders, keep the best
#' @return list(module, codelength, trials data.frame, partitions)
#' @noRd
.hcm_search <- function(fl, phys, n_phys, trials) {
  parts <- lapply(seq_len(trials), function(t) {
    .hcm_relabel(.hcm_trial(fl, phys, n_phys), fl$node_flow)
  })
  cls <- vapply(parts, function(p) .hcm_codelength(fl, phys, p)$codelength,
                numeric(1L))
  # ties (a zero-flow state placed alone or not codes identically): prefer
  # the partition with fewer modules, then the earlier trial
  n_mod <- vapply(parts, function(p) length(unique(p)), integer(1L))
  tied <- which(cls - min(cls) < 1e-10)
  best <- tied[order(n_mod[tied], tied)][1L]
  list(module = parts[[best]], codelength = cls[best], codelengths = cls,
       partitions = parts, best = best)
}

#' Resolve a user partition onto the state nodes
#' @noRd
.hcm_partition <- function(partition, states) {
  if (is.data.frame(partition)) {
    if (!all(c("state", "module") %in% names(partition))) {
      .hcm_bad_input("a `partition` data.frame needs `state` and `module` columns")
    }
    if (anyDuplicated(partition$state)) {
      .hcm_bad_input("`partition` lists a state more than once")
    }
    partition <- stats::setNames(partition$module,
                                 as.character(partition$state))
  }
  if (!is.atomic(partition) || is.null(names(partition))) {
    .hcm_bad_input(paste0("`partition` must be a named vector (names = state ",
                          "labels) or a data.frame with `state` and `module`"))
  }
  missing_states <- setdiff(states, names(partition))
  if (length(missing_states)) {
    .hcm_bad_input(sprintf("`partition` does not cover state(s): %s",
                           paste(utils::head(missing_states, 5L),
                                 collapse = ", ")))
  }
  lab <- partition[states]
  if (anyNA(lab)) .hcm_bad_input("`partition` has missing module labels")
  match(as.character(lab), unique(as.character(lab)))
}

#' Pairwise adjusted Rand index between trial partitions
#' @noRd
.hcm_trial_ari <- function(parts) {
  k <- length(parts)
  if (k < 2L) return(NA_real_)
  pairs <- utils::combn(k, 2L)
  vapply(seq_len(ncol(pairs)), function(j) {
    .thg_ari(parts[[pairs[1L, j]]], parts[[pairs[2L, j]]])
  }, numeric(1L))
}

#' Memory-network communities with the map equation
#'
#' @description
#' Finds modules in a higher-order network by minimising the map equation for
#' memory networks (Rosvall et al. 2014; Edler, Bohlin & Rosvall 2017). The
#' state nodes of the HON (`"a -> b"` is physical node `b` reached from `a`)
#' are clustered, but each module's codebook is written over **physical**
#' nodes: state nodes of one physical node that land in the same module share
#' a codeword. This is not first-order Infomap run on the state graph, and a
#' physical node whose state nodes fall in different modules belongs to all of
#' them, so modules of physical nodes can overlap.
#'
#' The same search is run on the first-order network obtained by projecting
#' the memory network's link flow onto physical nodes, so the result says
#' whether memory compresses the flow better than a first-order description.
#'
#' @details
#' \strong{Flow.} State node visit rates are the stationary distribution of
#' the walk on the HON transition matrix \eqn{P_{\alpha\beta} =
#' w_{\alpha\beta}/w_\alpha} (Edler et al. 2017, Eq. 8-9). With probability
#' `teleportation` (\eqn{\tau}), and always from a state with no out-links,
#' the walker restarts at the tail of a link drawn in proportion to its weight
#' (teleportation to links; Lambiotte & Rosvall 2012). Teleportation is
#' unrecorded: link flow is \eqn{f_{\alpha\beta} = \pi_\alpha
#' P_{\alpha\beta}} renormalised to sum to one, and a state's flow is its link
#' in-flow. This is Infomap's directed flow model; node flows and
#' codelengths match the Python `infomap` package (2.15.1). `teleportation =
#' 0` uses the plain stationary distribution and needs an irreducible chain.
#'
#' \strong{Codelength} (Edler et al. 2017, Eq. 24-26). With module enter flow
#' \eqn{q^{in}_m}, exit flow \eqn{q^{out}_m}, and physical flow inside a
#' module \eqn{\pi_{i \cap m} = \sum_{\beta_i \in m} \pi_{\beta_i}},
#' \deqn{L(M) = q^{in} H(Q) + \sum_m p_m H(P_m),}
#' where \eqn{q^{in} = \sum_m q^{in}_m}, \eqn{H(Q)} is the entropy of the
#' enter rates \eqn{q^{in}_m / q^{in}}, \eqn{p_m = q^{out}_m + \sum_i
#' \pi_{i \cap m}}, and \eqn{H(P_m)} is the entropy of \eqn{\{q^{out}_m,
#' \pi_{i \cap m}\}} normalised by \eqn{p_m}. The one-module codelength is
#' the entropy of the physical-node flow.
#'
#' \strong{Search} (Edler et al. 2017, Algorithms 2-5): every state node
#' starts in its own module; in random sequential order each node moves to
#' the neighbouring or new module that most decreases \eqn{L}; modules are
#' then aggregated into nodes and the moves repeat, until \eqn{L} stops
#' decreasing (repeated node aggregation, the Louvain machinery). The result
#' is fine-tuned by restarting the aggregation from single state nodes placed
#' in their current modules. Coarse-tuning (Algorithm 6) is not implemented.
#' `trials` independent runs are made from one seeded random stream and the
#' shortest codelength is kept; `as.data.frame(x, what = "trials")` reports
#' every run and its agreement (adjusted Rand index) with the kept one.
#'
#' \strong{First-order comparison.} The link flow of the memory network is
#' summed over state nodes into a physical-node network, \eqn{w_{ij} =
#' \sum_{\alpha_i, \beta_j} f_{\alpha\beta}}, and searched with the same flow
#' model and the same number of trials. Rosvall et al. (2014) compare memory
#' and first-order maps of the same pathway data this way.
#'
#' @param hon A `net_hon` from [build_hon()].
#' @param partition `NULL` (default) to search, or a given state-node
#'   partition to evaluate without searching: a named vector (names = state
#'   labels as in `as.data.frame(hon, what = "nodes")`) or a data.frame with
#'   columns `state` and `module`.
#' @param trials Number of independent search trials (default 10). Ignored
#'   when `partition` is given, except for the first-order search.
#' @param teleportation Teleportation probability \eqn{\tau \in [0, 1)}
#'   (default 0.15, Infomap's default).
#' @param seed Integer seed for the random node orders (default 1). The
#'   caller's random-number stream is restored on exit.
#' @return A `net_hon_communities` object. Read it with [as.data.frame()]:
#'   `what = "states"` (one row per state node), `"physical"` (one row per
#'   physical node x module), `"modules"`, `"trials"`, `"first_order"` and
#'   `"codelength"`; see [as.data.frame.net_hon_communities()].
#' @section Conditions:
#' `hypernets_bad_input` for a non-`net_hon` input, invalid arguments, or a
#' partition that does not cover every state; `hypernets_not_ergodic` when
#' `teleportation = 0` and the walk has no unique stationary flow.
#' @references
#' Rosvall, M., Esquivel, A. V., Lancichinetti, A., West, J. D., & Lambiotte,
#' R. (2014). Memory in network flows and its effects on spreading dynamics
#' and community detection. \emph{Nature Communications}, 5, 4630.
#' \doi{10.1038/ncomms5630}
#'
#' Edler, D., Bohlin, L., & Rosvall, M. (2017). Mapping higher-order network
#' flows in memory and multilayer networks with Infomap. \emph{Algorithms},
#' 10(4), 112. \doi{10.3390/a10040112}
#'
#' Rosvall, M., & Bergstrom, C. T. (2008). Maps of random walks on complex
#' networks reveal community structure. \emph{Proceedings of the National
#' Academy of Sciences}, 105(4), 1118-1123. \doi{10.1073/pnas.0706851105}
#'
#' Lambiotte, R., & Rosvall, M. (2012). Ranking and clustering of nodes in
#' networks with smart teleportation. \emph{Physical Review E}, 85(5),
#' 056107. \doi{10.1103/PhysRevE.85.056107}
#'
#' Hubert, L., & Arabie, P. (1985). Comparing partitions. \emph{Journal of
#' Classification}, 2(1), 193-218. \doi{10.1007/BF01908075}
#' @examples
#' seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
#'              c("c", "h", "d", "c", "h", "d", "c"),
#'              c("a", "h", "b", "a", "h", "b"),
#'              c("c", "h", "d", "c", "h", "d"))
#' hon <- build_hon(seqs, max_order = 2L)
#' comm <- hon_communities(hon, trials = 3L)
#' comm
#' summary(comm)
#' as.data.frame(comm)
#' as.data.frame(comm, what = "physical")
#' @export
hon_communities <- function(hon, partition = NULL, trials = 10L,
                            teleportation = 0.15, seed = 1L) {
  if (!inherits(hon, "net_hon")) {
    .hcm_bad_input("`hon` must be a net_hon from build_hon()")
  }
  whole <- function(x) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1 &&
      abs(x - round(x)) < sqrt(.Machine$double.eps)
  }
  if (!whole(trials)) .hcm_bad_input("`trials` must be one whole number >= 1")
  if (!(is.numeric(teleportation) && length(teleportation) == 1L &&
        is.finite(teleportation) && teleportation >= 0 && teleportation < 1)) {
    .hcm_bad_input("`teleportation` must be a single number in [0, 1)")
  }
  if (!(is.numeric(seed) && length(seed) == 1L && is.finite(seed))) {
    .hcm_bad_input("`seed` must be a single number")
  }
  trials <- as.integer(trials)
  W <- hon$weights
  if (!is.matrix(W) || nrow(W) < 1L) {
    .hcm_bad_input("`hon` has no state network to cluster")
  }
  states <- rownames(W)
  phys_label <- vapply(strsplit(states, " -> ", fixed = TRUE),
                       function(p) p[length(p)], character(1L))
  phys_names <- sort(unique(phys_label))
  phys <- match(phys_label, phys_names)
  n_phys <- length(phys_names)

  # scoped RNG: restore the caller's stream on exit
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE, after = FALSE)
  set.seed(seed)

  fl <- .hcm_flow(unname(W), teleportation)

  if (is.null(partition)) {
    srch <- .hcm_search(fl, phys, n_phys, trials)
    module <- srch$module
    # agreement over states that carry flow: a zero-flow state (reached only
    # by teleportation) can sit anywhere at no codelength cost
    carries <- fl$node_flow > 0
    trial_ari <- vapply(srch$partitions, function(p) {
      .thg_ari(p[carries], module[carries])
    }, numeric(1L))
    trial_tab <- data.frame(
      trial = seq_len(trials),
      codelength = srch$codelengths,
      n_modules = vapply(srch$partitions, function(p) length(unique(p)),
                         integer(1L)),
      ari_to_best = trial_ari,
      best = seq_len(trials) == srch$best,
      stringsAsFactors = FALSE)
    pair_ari <- .hcm_trial_ari(lapply(srch$partitions, function(p) p[carries]))
  } else {
    module <- .hcm_relabel(.hcm_partition(partition, states), fl$node_flow)
    trial_tab <- data.frame(trial = integer(0L), codelength = numeric(0L),
                            n_modules = integer(0L), ari_to_best = numeric(0L),
                            best = logical(0L))
    pair_ari <- NA_real_
  }
  cl <- .hcm_codelength(fl, phys, module)
  one_level <- .hcm_codelength(fl, phys, rep(1L, length(module)))$codelength

  # first-order comparison on the projected physical network
  W1 <- matrix(0, n_phys, n_phys, dimnames = list(phys_names, phys_names))
  proj <- stats::aggregate(fl$flow, by = list(from = phys[fl$from],
                                              to = phys[fl$to]), FUN = sum)
  W1[cbind(proj$from, proj$to)] <- proj$x
  fl1 <- .hcm_flow(W1, teleportation)
  srch1 <- .hcm_search(fl1, seq_len(n_phys), n_phys, trials)
  cl1 <- .hcm_codelength(fl1, seq_len(n_phys), srch1$module)
  one_level1 <- .hcm_codelength(fl1, seq_len(n_phys),
                                rep(1L, n_phys))$codelength

  # tidy tables
  state_tab <- data.frame(state = states, physical = phys_label,
                          module = module, flow = fl$node_flow,
                          stringsAsFactors = FALSE)
  k <- max(module)
  pim <- stats::aggregate(fl$node_flow, by = list(physical = phys_label,
                                                  module = module), FUN = sum)
  names(pim)[names(pim) == "x"] <- "flow"
  # a state with no link in-flow (reached only by teleportation) carries no
  # physical flow, so it does not make its physical node a module member
  pim <- pim[pim$flow > 0, , drop = FALSE]
  tot <- tapply(pim$flow, pim$physical, sum)
  pim$share <- pim$flow / tot[pim$physical]
  pim$share[!is.finite(pim$share)] <- 0
  nmod <- table(pim$physical)
  pim$n_modules <- as.integer(nmod[pim$physical])
  pim <- pim[order(pim$physical, pim$module), , drop = FALSE]
  rownames(pim) <- NULL

  cross <- module[fl$from] != module[fl$to]
  mod_tab <- data.frame(
    module = seq_len(k),
    n_states = tabulate(module, nbins = k),
    n_physical = vapply(seq_len(k), function(j)
      length(unique(phys_label[module == j])), integer(1L)),
    flow = vapply(seq_len(k), function(j) sum(fl$node_flow[module == j]),
                  numeric(1L)),
    exit_flow = vapply(seq_len(k), function(j)
      sum(fl$flow[cross & module[fl$from] == j]), numeric(1L)),
    enter_flow = vapply(seq_len(k), function(j)
      sum(fl$flow[cross & module[fl$to] == j]), numeric(1L)),
    stringsAsFactors = FALSE)

  fo_tab <- data.frame(physical = phys_names, module = srch1$module,
                       flow = fl1$node_flow, stringsAsFactors = FALSE)

  cl_tab <- data.frame(
    model = c("memory", "first_order"),
    codelength = c(cl$codelength, cl1$codelength),
    index_codelength = c(cl$index, cl1$index),
    module_codelength = c(cl$module, cl1$module),
    one_level_codelength = c(one_level, one_level1),
    savings_bits = c(one_level - cl$codelength, one_level1 - cl1$codelength),
    savings_pct = 100 * c(1 - cl$codelength / one_level,
                          1 - cl1$codelength / one_level1),
    n_modules = c(k, length(unique(srch1$module))),
    stringsAsFactors = FALSE)

  structure(list(
    states = state_tab,
    physical = pim,
    modules = mod_tab,
    trials = trial_tab,
    first_order = fo_tab,
    codelengths = cl_tab,
    trial_pairwise_ari = pair_ari,
    weights = W,
    physical_weights = W1,
    teleportation = teleportation,
    n_trials = if (is.null(partition)) trials else 0L,
    searched = is.null(partition),
    seed = seed
  ), class = "net_hon_communities")
}

#' Tidy tables of a memory-network community result
#'
#' @param x A `net_hon_communities` object.
#' @param row.names,optional Ignored (S3 consistency).
#' @param ... Ignored.
#' @param what Which table:
#'   \describe{
#'     \item{`"states"` (default)}{one row per state node: `state`,
#'       `physical`, `module`, `flow`.}
#'     \item{`"physical"`}{one row per physical node x module it appears in:
#'       `physical`, `module`, `flow` (\eqn{\pi_{i \cap m}}), `share` (of that
#'       physical node's flow), `n_modules` (> 1 marks an overlapping node).}
#'     \item{`"modules"`}{one row per module: `module`, `n_states`,
#'       `n_physical`, `flow`, `exit_flow`, `enter_flow`.}
#'     \item{`"trials"`}{one row per search trial: `trial`, `codelength`,
#'       `n_modules`, `ari_to_best` (adjusted Rand index with the kept
#'       partition), `best`.}
#'     \item{`"first_order"`}{one row per physical node: its module in the
#'       best first-order partition and its first-order flow.}
#'     \item{`"codelength"`}{one row per model (`memory`, `first_order`):
#'       `codelength`, `index_codelength`, `module_codelength`,
#'       `one_level_codelength`, `savings_bits`, `savings_pct`, `n_modules`.}
#'   }
#' @param module Integer vector or `NULL`: keep only these modules (tables
#'   `"states"`, `"physical"`, `"modules"`).
#' @param overlapping Logical: for `what = "physical"`, keep only physical
#'   nodes that appear in more than one module. Default `FALSE`.
#' @return A base data.frame as described under `what`.
#' @examples
#' seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
#'              c("c", "h", "d", "c", "h", "d", "c"))
#' comm <- hon_communities(build_hon(seqs, max_order = 2L), trials = 2L)
#' as.data.frame(comm, what = "modules")
#' as.data.frame(comm, what = "physical", overlapping = TRUE)
#' @export
as.data.frame.net_hon_communities <- function(
    x, row.names = NULL, optional = FALSE, ...,
    what = c("states", "physical", "modules", "trials", "first_order",
             "codelength"),
    module = NULL, overlapping = FALSE) {
  what <- match.arg(what)
  stopifnot(
    "`overlapping` must be TRUE or FALSE" =
      is.logical(overlapping) && length(overlapping) == 1L && !is.na(overlapping),
    "`module` must be NULL or whole numbers" =
      is.null(module) || (is.numeric(module) && all(is.finite(module)))
  )
  out <- switch(what,
                states = x$states,
                physical = x$physical,
                modules = x$modules,
                trials = x$trials,
                first_order = x$first_order,
                codelength = x$codelengths)
  if (!is.null(module) && what %in% c("states", "physical", "modules")) {
    out <- out[out$module %in% module, , drop = FALSE]
  }
  if (isTRUE(overlapping) && identical(what, "physical")) {
    out <- out[out$n_modules > 1L, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}

#' Print a memory-network community result
#'
#' @param x A `net_hon_communities` object.
#' @param ... Ignored.
#' @return `x`, invisibly.
#' @examples
#' seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
#'              c("c", "h", "d", "c", "h", "d", "c"))
#' print(hon_communities(build_hon(seqs, max_order = 2L), trials = 2L))
#' @export
print.net_hon_communities <- function(x, ...) {
  cl <- x$codelengths
  mem <- cl[cl$model == "memory", , drop = FALSE]
  fo <- cl[cl$model == "first_order", , drop = FALSE]
  n_over <- length(unique(x$physical$physical[x$physical$n_modules > 1L]))
  cat("Memory-network communities (map equation)\n")
  cat(sprintf("  State nodes: %d | physical nodes: %d | modules: %d\n",
              nrow(x$states), length(unique(x$states$physical)),
              mem$n_modules))
  n_zero <- sum(x$modules$flow <= 0)
  if (n_zero > 0L) {
    cat(sprintf(paste0("  (%d module(s) hold only states with zero flow, ",
                       "reached by teleportation alone)\n"), n_zero))
  }
  cat(sprintf("  Overlapping physical nodes: %d\n", n_over))
  cat(sprintf("  Codelength memory:      %.4f bits (one module %.4f)\n",
              mem$codelength, mem$one_level_codelength))
  cat(sprintf("  Codelength first-order: %.4f bits (one module %.4f; modules: %d)\n",
              fo$codelength, fo$one_level_codelength, fo$n_modules))
  if (isTRUE(x$searched)) {
    cat(sprintf("  Trials: %d (teleportation %.2f, seed %s)\n", x$n_trials,
                x$teleportation, format(x$seed)))
  } else {
    cat("  Partition supplied (no memory-network search)\n")
  }
  cat("  Tables: as.data.frame(x, what = \"states\" | \"physical\" |",
      "\"modules\" | \"trials\" | \"first_order\" | \"codelength\")\n")
  invisible(x)
}

#' Summarise a memory-network community result
#'
#' Prints the codelength comparison (memory vs first-order, bits and percent
#' saved against the one-module code), the module count, the overlap, and
#' trial stability.
#'
#' @param object A `net_hon_communities` object.
#' @param ... Ignored.
#' @return The codelength table (as `as.data.frame(object, what =
#'   "codelength")`), invisibly: one row per model with `codelength`,
#'   `index_codelength`, `module_codelength`, `one_level_codelength`,
#'   `savings_bits`, `savings_pct`, `n_modules`.
#' @examples
#' seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
#'              c("c", "h", "d", "c", "h", "d", "c"))
#' summary(hon_communities(build_hon(seqs, max_order = 2L), trials = 2L))
#' @export
summary.net_hon_communities <- function(object, ...) {
  cl <- as.data.frame(object, what = "codelength")
  mem <- cl[cl$model == "memory", , drop = FALSE]
  fo <- cl[cl$model == "first_order", , drop = FALSE]
  phys <- object$physical
  n_phys <- length(unique(phys$physical))
  n_over <- length(unique(phys$physical[phys$n_modules > 1L]))
  cat("Memory-network communities: summary\n")
  cat(sprintf("  Memory:      %.4f bits, modules: %d, saves %.4f bits (%.1f%%) vs one module\n",
              mem$codelength, mem$n_modules, mem$savings_bits,
              mem$savings_pct))
  cat(sprintf("  First-order: %.4f bits, modules: %d, saves %.4f bits (%.1f%%) vs one module\n",
              fo$codelength, fo$n_modules, fo$savings_bits, fo$savings_pct))
  cat(sprintf("  Memory minus first-order: %+.4f bits\n",
              mem$codelength - fo$codelength))
  cat(sprintf("  Overlap: %d of %d physical nodes in more than one module (mean %.2f modules per node)\n",
              n_over, n_phys,
              mean(as.numeric(table(phys$physical)))))
  tr <- object$trials
  if (nrow(tr)) {
    n_best <- sum(tr$codelength - min(tr$codelength) < 1e-10)
    cat(sprintf("  Trials: %d; %d reached the best codelength; ARI to best: mean %.3f, min %.3f\n",
                nrow(tr), n_best, mean(tr$ari_to_best), min(tr$ari_to_best)))
  }
  invisible(cl)
}

.HCM_PALETTE <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
                  "#D55E00", "#CC79A7", "#999999")

.hcm_col <- function(m) .HCM_PALETTE[(m - 1L) %% length(.HCM_PALETTE) + 1L]

.hcm_node_size <- function(f, n = length(f)) {
  top <- if (n > 40L) 3.5 else 7
  0.35 * top + 0.65 * top * sqrt(f / max(f))
}

#' Arguments handed to cograph::overlay_communities() for the state view
#'
#' Kept separate from the drawing so tests can assert exactly what reaches
#' cograph.
#' @return list(args = <overlay_communities() arguments>, n_hidden =
#'   zero-flow state nodes left out)
#' @noRd
.hcm_state_plot_args <- function(x, show_zero_flow) {
  st <- x$states
  hidden <- if (show_zero_flow) 0L else sum(st$flow <= 0)
  if (!show_zero_flow) st <- st[st$flow > 0, , drop = FALSE]
  st <- st[order(st$module, st$state), , drop = FALSE]
  W <- x$weights[st$state, st$state, drop = FALSE]
  mods <- sort(unique(st$module))
  biggest <- max(tabulate(st$module))
  # cograph's group layout: one ring per community; the ring widens with
  # the community so large ones do not pile their nodes on top of each other
  lay <- cograph::layout_groups(
    cograph::as_cograph(W), st$module,
    inner_radius = min(0.3, max(0.12, 0.006 * biggest)))
  args <- list(
    x = W,
    communities = split(st$state, factor(st$module, levels = mods)),
    blob_colors = .hcm_col(mods), blob_alpha = 0.18,
    layout = as.matrix(lay),
    groups = sprintf("Community %d", st$module),
    legend = TRUE, legend_edge_colors = FALSE,
    node_fill = .hcm_col(st$module),
    node_size = .hcm_node_size(st$flow),
    labels = st$state, label_size = if (nrow(st) > 40L) 0.45 else 0.8,
    title_size = 1.1,
    edge_color = "grey40", edge_alpha = 0.25, threshold = 0.05,
    title = sprintf("%d communit%s of state nodes", length(mods),
                    if (length(mods) == 1L) "y" else "ies"))
  list(args = args, n_hidden = hidden)
}

#' The physical view as a community hypergraph plus its overlays
#'
#' Every community is a hyperedge over the physical nodes it holds, so a
#' shared physical node is a member of each of its communities' pebbles.
#' @return list(hypergraph, community (factor named by hyperedge), sizes
#'   (node, value = physical flow), flow (community flow, named by
#'   hyperedge), moves (from, to, weight = projected link flow), title, notes (caption lines beyond the size/triangle key))
#' @noRd
.hcm_physical_plot_data <- function(x) {
  ph <- x$physical
  mods <- sort(unique(ph$module))
  community_names <- sprintf("Community %d", mods)
  members <- data.frame(physical = ph$physical,
                        community = sprintf("Community %d", ph$module),
                        stringsAsFactors = FALSE)
  hg <- group_hypergraph(members, actor = "physical", group = "community")
  # a community with one physical node has no pebble, so it gets no colour
  # key either (an empty key reads as a missing colour); the caption names it
  counts <- tabulate(match(ph$module, mods), nbins = length(mods))
  lone <- mods[counts == 1L]
  community <- stats::setNames(
    factor(community_names, levels = community_names[counts > 1L]),
    community_names)
  flow <- tapply(ph$flow, ph$physical, sum)
  sizes <- data.frame(node = names(flow), value = as.numeric(flow),
                      stringsAsFactors = FALSE)
  W1 <- x$physical_weights
  ij <- which(W1 > 0, arr.ind = TRUE)
  moves <- data.frame(from = rownames(W1)[ij[, 1L]],
                      to = colnames(W1)[ij[, 2L]], weight = W1[ij],
                      stringsAsFactors = FALSE)
  n_mod <- tapply(ph$module, ph$physical, length)
  n_shared <- sum(n_mod > 1L)
  n_hidden <- sum(x$states$flow <= 0)
  notes <- c(
    if (length(lone)) {
      alone <- ph$physical[match(lone, ph$module)]
      sprintf("%s: a single physical node, drawn without a pebble.",
              paste(sprintf("Community %d (%s)", lone, alone),
                    collapse = ", "))
    },
    if (n_hidden > 0L) {
      sprintf("%d zero-flow state%s not drawn (reached only by teleportation).",
              n_hidden, if (n_hidden == 1L) "" else "s")
    }
  )
  # a community's flow: the visit rate of its state nodes, for its title box
  module_flow <- tapply(x$states$flow, x$states$module, sum)
  community_flow <- stats::setNames(
    as.numeric(module_flow[as.character(mods)]), community_names)
  list(hypergraph = hg, community = community, flow = community_flow,
       sizes = sizes, moves = moves,
       title = sprintf("%d communit%s, %d shared physical node%s",
                       length(mods), if (length(mods) == 1L) "y" else "ies",
                       n_shared, if (n_shared == 1L) "" else "s"),
       notes = notes, n_hidden = n_hidden)
}

#' Plot memory-network communities
#'
#' \describe{
#'   \item{`type = "physical"` (default)}{every community is drawn as a
#'     pebble around the physical nodes it holds, through
#'     [plot.net_hypergraph()] on the community hypergraph (member = physical
#'     node, group = `"Community k"`, from `as.data.frame(x, what =
#'     "physical")`). A physical node shared by several communities lies
#'     inside each of their pebbles. Pebbles are coloured by community
#'     (Okabe-Ito, legend "Community"). Each node is a black circle whose
#'     area follows its physical flow (the memory network's visit rate summed
#'     over its state nodes), and a triangle inside the circle points at the
#'     physical node its link flow most often goes to next. Each community
#'     has a title box with its name and flow (the visit rate of its state
#'     nodes) beside its pebble, drawn by [plot.net_hypergraph()] (haloed
#'     labels, legend below, wide margins); the boxes sit a little further
#'     out than there (`title_gap = 0.14`) to clear the labels of the nodes at
#'     a pebble's rim. A
#'     community holding a single physical node has no pebble of its own;
#'     its title box sits by the node and the caption names it.}
#'   \item{`type = "states"`}{the state nodes of the higher-order network,
#'     drawn by [cograph::overlay_communities()], one ring per community
#'     ([cograph::layout_groups()]), with community blobs, cograph's group
#'     legend, short state labels, and faded edges (transition probability
#'     below 0.05 hidden).}
#' }
#' State nodes with zero flow (reached only by teleportation, see
#' [hon_communities()]) carry no physical flow and are never drawn in the
#' physical view; the state view leaves them and the singleton communities
#' they form out unless `show_zero_flow = TRUE`. A caption says how many were
#' left out. The tables returned by [as.data.frame.net_hon_communities()] are
#' unaffected.
#'
#' @param x A `net_hon_communities` object.
#' @param type `"physical"` (default) or `"states"`.
#' @param show_zero_flow Draw zero-flow state nodes in the state view?
#'   Default `FALSE`.
#' @param ... For `type = "physical"`, passed to [plot.net_hypergraph()]
#'   (e.g. `seed`, `label_size`, `arrow_style = "outside"`, `layout`),
#'   overriding the defaults set here; for `type = "states"`, passed to
#'   [cograph::overlay_communities()] and on to [cograph::splot()].
#' @return For `type = "physical"`, a ggplot object (print it to draw). For
#'   `type = "states"`, `x`, invisibly (cograph draws with base graphics).
#' @section Conditions:
#' `hypernets_bad_input` from [plot.net_hypergraph()] for arguments passed
#' through `...` that it rejects.
#' @examples
#' seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
#'              c("c", "h", "d", "c", "h", "d", "c"))
#' comm <- hon_communities(build_hon(seqs, max_order = 2L), trials = 2L)
#' plot(comm)
#' plot(comm, type = "states")
#' @export
plot.net_hon_communities <- function(x, type = c("physical", "states"),
                                     show_zero_flow = FALSE, ...) {
  type <- match.arg(type)
  stopifnot("`show_zero_flow` must be TRUE or FALSE" =
              is.logical(show_zero_flow) && length(show_zero_flow) == 1L &&
              !is.na(show_zero_flow))
  dots <- list(...)
  if (identical(type, "states")) {
    spec <- .hcm_state_plot_args(x, show_zero_flow)
    args <- spec$args
    args[names(dots)] <- dots
    do.call(cograph::overlay_communities, args)
    if (spec$n_hidden > 0L) {
      # a caption in the bottom margin, clear of the network
      graphics::mtext(sprintf(
        "%d zero-flow state%s not drawn (reached only by teleportation)",
        spec$n_hidden, if (spec$n_hidden == 1L) "" else "s"),
        side = 1L, line = 4, cex = 0.85)
    }
    return(invisible(x))
  }
  spec <- .hcm_physical_plot_data(x)
  # a title box per community with its flow; plot.net_hypergraph() puts the
  # legends below and halos the labels
  args <- list(x = spec$hypergraph,
               color_by = spec$community, legend_title = "Community",
               titles = spec$flow, unit = "flow",
               node_sizes = spec$sizes, direction = spec$moves,
               arrow_style = "inside",
               size_title = "physical flow (visit rate summed over the node's states)",
               title_gap = 0.14)
  args[names(dots)] <- dots
  p <- do.call(plot.net_hypergraph, args)
  # the key one sentence to a line, then the notes, left-aligned under the
  # whole plot: one long line is wider than the figure
  p + ggplot2::labs(
    title = spec$title,
    caption = paste(c(
      "Circle area: physical flow (visit rate summed over the node's states).",
      "Triangle: points to the event that most often follows it.",
      spec$notes), collapse = "\n")
  ) + ggplot2::theme(plot.caption = ggplot2::element_text(hjust = 0),
                     plot.caption.position = "plot")
}
