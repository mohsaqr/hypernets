# Hypergraph-native modularity: the Kaminski et al. (2019) score of a
# partition, and Kumar et al. (2020) iteratively reweighted modularity
# maximisation (IRMM). Both read the hypergraph as a family of node SETS
# (the binary incidence pattern) with one weight per hyperedge. The oracle
# is HyperNetX 2.4.3 (`hypernetx.algorithms`: modularity(), kumar(),
# two_section()); see local_testing_and_equivalence/test-equiv-hg-modularity-*.R.

# ---- shared parts -----------------------------------------------------------

# Binary incidence, hyperedge sizes and weights, with empty hyperedges
# removed (they bind no node, so neither model defines a contribution).
.hg_mod_parts <- function(hg, edge_weights = NULL) {
  incidence <- hg$incidence
  m <- ncol(incidence)
  edge_weights <- .hl_check_edge_weights(edge_weights, m)
  w <- as.numeric(edge_weights %||% hg$window_counts %||% rep(1, m))
  pattern <- if (.thg_is_sparse(hg)) {
    methods::as(incidence > 0, "dMatrix")
  } else {
    (as.matrix(incidence) > 0) * 1
  }
  sizes <- as.numeric(Matrix::colSums(pattern))
  keep <- sizes > 0
  edges <- colnames(incidence) %||% paste0("e", seq_len(m))
  nodes <- rownames(incidence) %||% hg$nodes
  list(
    pattern = pattern[, keep, drop = FALSE],
    sizes = sizes[keep],
    weights = w[keep],
    edges = edges[keep],
    nodes = nodes
  )
}

# Kaminski et al. (2019) / HyperNetX weight w(d, c) of a d-edge whose
# largest part holds c of its nodes. Vectorised over d and c.
.hg_mod_wdc <- function(type) {
  switch(type,
    linear = function(d, c) ifelse(c > d / 2, c / d, 0),
    majority = function(d, c) as.numeric(c > d / 2),
    strict = function(d, c) as.numeric(c == d)
  )
}

# Node x community indicator matrix (dense; communities are few).
.hg_mod_membership <- function(labels) {
  levels <- unique(labels)
  m <- matrix(0, length(labels), length(levels))
  m[cbind(seq_along(labels), match(labels, levels))] <- 1
  colnames(m) <- levels
  m
}

# Per-community edge contribution and degree tax, unnormalised.
# Edge contribution of community A: sum over hyperedges whose strict
# majority lies in A of w_e * wdc(d_e, c_e). Degree tax of A:
# sum_d W_d sum_{c > d/2} wdc(d, c) Binom(c; d, vol(A) / vol(V)), where W_d
# is the total weight of d-edges and vol the weighted degree (Kaminski et
# al. 2019, Eqs. for the weighted general modularity; HyperNetX modularity()).
.hg_mod_fit <- function(parts, labels, type) {
  wdc <- .hg_mod_wdc(type)
  membership <- .hg_mod_membership(labels)
  communities <- colnames(membership)
  counts <- as.matrix(Matrix::crossprod(parts$pattern, membership))
  top <- max.col(counts, ties.method = "first")
  largest <- counts[cbind(seq_len(nrow(counts)), top)]
  contribution <- parts$weights * wdc(parts$sizes, largest)
  edge_contribution <- vapply(seq_along(communities), function(a) {
    sum(contribution[top == a])
  }, numeric(1L))

  degree <- as.numeric(parts$pattern %*% parts$weights)
  volume <- as.numeric(crossprod(membership, degree))
  share <- volume / sum(volume)
  weight_by_size <- tapply(parts$weights, parts$sizes, sum)
  d_values <- as.numeric(names(weight_by_size))
  # every (d, c) cell with c > d/2 and its coefficient W_d * wdc(d, c)
  cells <- do.call(rbind, lapply(seq_along(d_values), function(i) {
    d <- d_values[[i]]
    c_values <- seq.int(floor(d / 2) + 1, d)
    data.frame(d = d, c = c_values,
               coef = weight_by_size[[i]] * wdc(d, c_values))
  }))
  degree_tax <- vapply(share, function(p) {
    sum(cells$coef * stats::dbinom(cells$c, cells$d, p))
  }, numeric(1L))

  data.frame(
    community = communities,
    n_nodes = as.integer(colSums(membership)),
    volume = volume,
    edge_contribution = edge_contribution,
    degree_tax = degree_tax,
    stringsAsFactors = FALSE
  )
}

#' Hypergraph modularity of a partition
#'
#' Scores a partition of the nodes of a hypergraph with the hypergraph
#' modularity of Kaminski et al. (2019): the weighted share of hyperedges
#' that a community "owns" (the edge contribution) minus the share expected
#' under the generalised Chung-Lu null model that preserves weighted degrees
#' and hyperedge sizes (the degree tax). The hypergraph is read as a family
#' of node sets: cell values of the incidence are ignored, only membership
#' and one weight per hyperedge count.
#'
#' A hyperedge of size \eqn{d} whose largest part holds \eqn{c} of its nodes
#' contributes \eqn{w_e\,\omega(d, c)} to the community holding that part
#' when \eqn{c > d/2}. The three `type`s are HyperNetX's `wdc` functions:
#' `"strict"` (\eqn{\omega = 1} only when the whole edge lies in one
#' community, the strict modularity of Kaminski et al. 2019),
#' `"majority"` (\eqn{\omega = 1} when a strict majority does), and
#' `"linear"` (\eqn{\omega = c/d} for a strict majority, the HyperNetX
#' default; Kaminski, Pralat & Theberge 2020). The degree tax of community
#' \eqn{A} is \eqn{\sum_d W_d \sum_{c > d/2} \omega(d, c)\,
#' \mathrm{Binom}(c; d, \mathrm{vol}(A)/\mathrm{vol}(V))}, with \eqn{W_d}
#' the total weight of \eqn{d}-edges and \eqn{\mathrm{vol}} the weighted
#' degree. Both terms are divided by the total hyperedge weight.
#'
#' @param hg A `net_hg`.
#' @param partition The partition to score: a data.frame with `node` and a
#'   label column (`community`, `cluster`, `label` or `predicted`), a named
#'   label vector, an unnamed vector in node order, or an [hg_communities()]
#'   fit (its AMI medoid is scored). Every node must be labelled.
#' @param type Hyperedge weighting \eqn{\omega(d, c)}: `"linear"` (default),
#'   `"majority"` or `"strict"`.
#' @param edge_weights Positive hyperedge weights (one per hyperedge, or one
#'   value recycled). `NULL` uses the window counts of a
#'   [window_hypergraph()], else unit weights.
#' @param what `"score"` (default; one row) or `"communities"` (one row per
#'   community; its `modularity` column sums to the score).
#' @return A base data.frame. For `what = "score"`: one row with `type`,
#'   `modularity`, `edge_contribution`, `degree_tax` (both already divided
#'   by `total_weight`, so `modularity = edge_contribution - degree_tax`),
#'   `n_communities`, `n_edges` (non-empty hyperedges) and `total_weight`.
#'   For `what = "communities"`: one row per community with `community`,
#'   `n_nodes`, `volume` (summed weighted degree), `edge_contribution`,
#'   `degree_tax` and `modularity`, sorted by decreasing `modularity`.
#'   Raises `hypernets_bad_input` for a non-`net_hg` input, a partition that
#'   misses a node or has `NA` labels, or a hypergraph with no non-empty
#'   hyperedge; invalid `edge_weights` fail the shared weight check.
#' @references Kaminski, B., Poulin, V., Pralat, P., Szufel, P., &
#'   Theberge, F. (2019). Clustering via hypergraph modularity. *PLoS ONE*,
#'   14(11), e0224307. \doi{10.1371/journal.pone.0224307}
#'
#'   Kaminski, B., Pralat, P., & Theberge, F. (2020). Community detection
#'   algorithm using hypergraph modularity. In *Complex Networks & Their
#'   Applications IX*, Studies in Computational Intelligence 943, 152-163.
#'   \doi{10.1007/978-3-030-65347-7_13}
#' @seealso [hg_communities()] with `type = "irmm"` to find a partition.
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z",
#'              "c", "x"),
#'   edge = c(rep(paste0("e", 1:4), each = 3), "e5", "e5")
#' )
#' h <- group_hypergraph(dat, "member", "edge")
#' split <- data.frame(node = c("a", "b", "c", "x", "y", "z"),
#'                     community = c(1, 1, 1, 2, 2, 2))
#' hg_modularity(h, split)
#' hg_modularity(h, split, type = "strict")
#' hg_modularity(h, split, what = "communities")
#' @export
hg_modularity <- function(hg, partition,
                          type = c("linear", "majority", "strict"),
                          edge_weights = NULL,
                          what = c("score", "communities")) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  what <- match.arg(what)
  parts <- .hg_mod_parts(hg, edge_weights)
  if (!length(parts$sizes)) {
    .thg_bad_input("hypergraph modularity needs at least one non-empty hyperedge")
  }
  labels <- .thg_partition_vector(partition, parts$nodes)
  if (anyNA(labels)) .thg_bad_input("the partition labels contain NA")
  total <- sum(parts$weights)
  per <- .hg_mod_fit(parts, labels, type)
  per$edge_contribution <- per$edge_contribution / total
  per$degree_tax <- per$degree_tax / total
  per$modularity <- per$edge_contribution - per$degree_tax
  if (identical(what, "communities")) {
    per <- per[order(-per$modularity, per$community), , drop = FALSE]
    rownames(per) <- NULL
    return(per)
  }
  data.frame(
    type = type,
    modularity = sum(per$modularity),
    edge_contribution = sum(per$edge_contribution),
    degree_tax = sum(per$degree_tax),
    n_communities = nrow(per),
    n_edges = length(parts$sizes),
    total_weight = total,
    stringsAsFactors = FALSE
  )
}

# ---- IRMM (Kumar et al. 2020) -----------------------------------------------

# Random-walk clique reduction A = H W (D_e - I)^{-1} H^T with a zero
# diagonal (Kumar et al. 2020, Algorithm 1 lines 4-6; HyperNetX
# two_section()). Hyperedges of size 1 have no pairs and are left out. A
# sparse incidence gives a sparse reduction (igraph reads either).
.hg_irmm_two_section <- function(parts, weights) {
  pairs <- parts$sizes > 1
  scale <- weights[pairs] / (parts$sizes[pairs] - 1)
  h <- parts$pattern[, pairs, drop = FALSE]
  a <- Matrix::tcrossprod(h %*% Matrix::Diagonal(x = scale), h)
  a <- if (methods::is(parts$pattern, "sparseMatrix")) {
    methods::as(methods::as(a, "generalMatrix"), "CsparseMatrix")
  } else {
    as.matrix(a)
  }
  diag(a) <- 0
  dimnames(a) <- list(parts$nodes, parts$nodes)
  a
}

# Kumar et al. (2020) Eq. 6: w'(e) = (1/m) sum_i 1 / (k_i(e) + 1) * (d_e + c),
# the sum running over all c clusters (k_i(e) = |e intersect C_i|, possibly 0)
# and m the number of hyperedges.
.hg_irmm_reweight <- function(parts, labels) {
  membership <- .hg_mod_membership(labels)
  counts <- as.matrix(Matrix::crossprod(parts$pattern, membership))
  n_clusters <- ncol(membership)
  as.numeric(rowSums(1 / (counts + 1)) * (parts$sizes + n_clusters) /
               length(parts$sizes))
}

# Louvain modularity maximisation on the weighted two-section (Blondel et
# al. 2008) through igraph, as HyperNetX does through python-igraph.
# Labels are renumbered by first appearance in node order.
.hg_irmm_louvain <- function(adjacency, iteration) {
  graph <- igraph::graph_from_adjacency_matrix(
    adjacency, mode = "undirected", weighted = TRUE, diag = FALSE
  )
  fit <- igraph::cluster_louvain(graph, weights = igraph::E(graph)$weight)
  labels <- as.integer(igraph::membership(fit))
  match(labels, unique(labels))
}

# One IRMM run from initial weights w0: Louvain, then reweight with the
# moving average w <- (w + w') / 2 and re-cluster until the largest weight
# change is at most `delta` (HyperNetX kumar(); the paper's ||W - W_prev||
# with the max norm) or `max_iter` reweighting passes were made. `cluster`
# is the graph-clustering step, a function(adjacency, iteration); the oracle
# tests replace it with recorded HyperNetX memberships.
.hg_irmm_fit <- function(parts, delta, max_iter, cluster = .hg_irmm_louvain) {
  weights <- parts$weights
  labels <- cluster(.hg_irmm_two_section(parts, weights), 0L)
  history <- list(weights)
  iteration <- 0L
  change <- NA_real_
  converged <- FALSE
  # An iterative fit: each pass needs the previous pass's clustering, so
  # the passes cannot be vectorised; the loop is bounded by `max_iter`.
  while (iteration < max_iter) {
    updated <- 0.5 * weights + 0.5 * .hg_irmm_reweight(parts, labels)
    change <- max(abs(updated - weights))
    weights <- updated
    iteration <- iteration + 1L
    history[[iteration + 1L]] <- weights
    labels <- cluster(.hg_irmm_two_section(parts, weights), iteration)
    if (change <= delta) {
      converged <- TRUE
      break
    }
  }
  list(labels = labels, weights = weights, iterations = iteration,
       converged = converged, max_weight_change = change, history = history)
}

# Pairwise partition similarity over runs, as in hg_communities().
.hg_irmm_similarity <- function(partitions) {
  n_runs <- length(partitions)
  run_names <- paste0("run_", seq_len(n_runs))
  pairs <- if (n_runs > 1L) utils::combn(n_runs, 2L) else matrix(integer(), 2L, 0L)
  metrics <- c(ami = .thg_ami, ari = .thg_ari, nmi = .thg_nmi)
  lapply(metrics, function(metric) {
    m <- diag(1, n_runs, n_runs)
    dimnames(m) <- list(run_names, run_names)
    values <- vapply(seq_len(ncol(pairs)), function(k) {
      metric(partitions[[pairs[1L, k]]], partitions[[pairs[2L, k]]])
    }, numeric(1L))
    m[t(pairs)] <- values
    m[t(pairs[2:1, , drop = FALSE])] <- values
    m
  })
}

#' IRMM communities of a hypergraph (engine of hg_communities(type = "irmm"))
#'
#' Runs Kumar et al. (2020) iteratively reweighted modularity maximisation
#' once per seed and returns the AMI medoid, in the `hg_communities`
#' contract shared with the Infomap path.
#' @return An `hg_communities` object.
#' @noRd
.hg_irmm_communities <- function(hg, n_runs = 50L, seeds = NULL, delta = 0.01,
                                 max_iter = 50L, edge_weights = NULL) {
  .thg_check_hg(hg)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(errorCondition(
      "IRMM needs the suggested package `igraph` for the Louvain step",
      class = "hypernets_missing_dependency", call = NULL
    ))
  }
  whole <- function(x, name) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 1 ||
        abs(x - round(x)) > sqrt(.Machine$double.eps)) {
      .thg_bad_input(sprintf("`%s` must be one positive whole number", name))
    }
    as.integer(x)
  }
  n_runs <- whole(n_runs, "n_runs")
  max_iter <- whole(max_iter, "max_iter")
  if (length(delta) != 1L || !is.numeric(delta) || !is.finite(delta) ||
      delta <= 0) {
    .thg_bad_input("`delta` must be one positive number")
  }
  seeds <- seeds %||% seq_len(n_runs)
  if (length(seeds) != n_runs || any(!is.finite(seeds)) ||
      any(abs(seeds - round(seeds)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`seeds` must contain one whole-number seed per run")
  }
  seeds <- as.integer(seeds)
  parts <- .hg_mod_parts(hg, edge_weights)
  if (!length(parts$sizes)) {
    .thg_bad_input("IRMM needs at least one non-empty hyperedge")
  }

  # Seed each run locally and restore the caller's RNG stream afterwards.
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  on.exit({
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, add = TRUE)
  fits <- lapply(seeds, function(seed) {
    set.seed(seed)
    .hg_irmm_fit(parts, delta = delta, max_iter = max_iter)
  })

  labels <- lapply(fits, function(fit) as.character(fit$labels))
  partitions <- do.call(rbind, lapply(seq_len(n_runs), function(i) {
    data.frame(node = parts$nodes, community = labels[[i]], run = i,
               seed = seeds[[i]], stringsAsFactors = FALSE)
  }))
  runs <- data.frame(
    run = seq_len(n_runs),
    seed = seeds,
    n_communities = vapply(labels, function(l) length(unique(l)), integer(1L)),
    iterations = vapply(fits, `[[`, integer(1L), "iterations"),
    converged = vapply(fits, `[[`, logical(1L), "converged"),
    max_weight_change = vapply(fits, `[[`, numeric(1L), "max_weight_change"),
    modularity = vapply(labels, function(l) {
      per <- .hg_mod_fit(parts, l, "linear")
      sum(per$edge_contribution - per$degree_tax) / sum(parts$weights)
    }, numeric(1L)),
    stringsAsFactors = FALSE
  )
  if (any(!runs$converged)) {
    warning(warningCondition(
      sprintf(paste("IRMM did not converge within %d reweighting passes",
                    "(delta = %g) in %d of %d runs; see",
                    "hg_get(fit, what = \"runs\")"),
              max_iter, delta, sum(!runs$converged), n_runs),
      class = "hypernets_no_converge", call = NULL
    ))
  }

  similarity <- .hg_irmm_similarity(labels)
  medoid_run <- unname(which.max(rowSums(similarity$ami))[1L])
  medoid <- data.frame(node = parts$nodes, community = labels[[medoid_run]],
                       stringsAsFactors = FALSE)
  sizes <- as.data.frame(table(medoid$community), stringsAsFactors = FALSE)
  names(sizes) <- c("community", "n_nodes")
  sizes <- sizes[order(-sizes$n_nodes, sizes$community), , drop = FALSE]
  rownames(sizes) <- NULL
  best <- fits[[medoid_run]]

  structure(list(
    medoid = medoid,
    partitions = partitions,
    similarity = similarity,
    runs = runs,
    sizes = sizes,
    projection = as.matrix(.hg_irmm_two_section(parts, best$weights)),
    weights = data.frame(edge = parts$edges, size = as.integer(parts$sizes),
                         initial_weight = parts$weights,
                         weight = best$weights, stringsAsFactors = FALSE),
    medoid_run = medoid_run,
    params = list(type = "irmm", n_runs = n_runs, seeds = seeds,
                  delta = delta, max_iter = max_iter)
  ), class = "hg_communities")
}
