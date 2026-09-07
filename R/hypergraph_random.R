# Random hypergraph generators. The model definitions follow HyperG's four
# core samplers, but return hypernets' net_hypergraph representation and use
# names that do not mask HyperG when both packages are loaded.

#' Sample Bernoulli-Incidence Random Hypergraphs
#'
#' Generates an Erdos--Renyi-style random hypergraph with `n` nodes and `m`
#' hyperedges by drawing every incidence independently as Bernoulli(`p`). If
#' `m` is omitted it is Poisson with mean `lambda`, or `n * p` when `lambda`
#' is also omitted. Empty and singleton hyperedges are retained because they
#' are valid outcomes of this incidence model.
#'
#' @param n Number of nodes.
#' @param m Number of hyperedges. For `hg_sample_gnp()`, `NULL` draws it from
#'   a Poisson distribution. For the uniform and regular models it is
#'   required.
#' @param p Incidence probability in `[0, 1]`.
#' @param lambda Optional Poisson mean for `m`.
#' @param seed Optional reproducibility seed. The caller's RNG state is
#'   restored on exit.
#' @return A `net_hypergraph` with binary incidence and model parameters in
#'   `$params`.
#' @references
#' Marchette, D. J. (2021). HyperG: Hypergraphs in R. R package version
#' 1.0.0.
#' @examples
#' h <- hg_sample_gnp(n = 20, m = 8, p = 0.2, seed = 1)
#' summary(h)
#' @export
hg_sample_gnp <- function(n, m = NULL, p, lambda = NULL, seed = NULL) {
  n <- .hgr_count(n, "n", minimum = 1L)
  .hgr_probability(p, "p")
  if (!is.null(lambda)) {
    if (length(lambda) != 1L || !is.finite(lambda) || lambda < 0) {
      stop("`lambda` must be one non-negative number.", call. = FALSE)
    }
  }
  .hgr_seed(seed)
  if (is.null(m)) m <- stats::rpois(1L, lambda %||% (n * p))
  m <- .hgr_count(m, "m", minimum = 0L)
  # Draw in edge x node order, then transpose. This is distributionally
  # immaterial but preserves exact seeded parity with HyperG's definition.
  incidence <- if (m == 0L) {
    matrix(integer(0), nrow = n, ncol = 0L)
  } else {
    t(matrix(stats::rbinom(n * m, 1L, p), nrow = m, ncol = n))
  }
  dimnames(incidence) <- list(paste0("V", seq_len(n)),
                              if (m) paste0("h", seq_len(m)) else character())
  .hgr_from_incidence(
    incidence, "hg_sample_gnp",
    list(model = "gnp", n = n, m = m, p = p, lambda = lambda, seed = seed)
  )
}

#' Sample Stochastic-Block-Model Hypergraphs
#'
#' Samples a graph stochastic block model and augments every sampled dyad to
#' a hyperedge. Additional members come from the endpoint blocks; `impurity`
#' members can then be replaced by nodes outside those blocks. This is the
#' hypergraph SBM construction used by HyperG.
#'
#' @param n Total node count; `NULL` uses `sum(block_sizes)`.
#' @param P Symmetric block-to-block dyad probability matrix.
#' @param block_sizes Positive integer community sizes.
#' @param d Hyperedge size, recycled across sampled dyads. With
#'   `variable_size = TRUE`, it is instead the Poisson mean and two is added.
#' @param impurity Non-negative integer number of augmented members replaced
#'   by nodes outside the endpoint blocks when possible.
#' @param variable_size Draw sizes as `2 + Poisson(d)`.
#' @param absolute_purity If true, impurity replacements must come from
#'   outside the endpoint blocks; otherwise any nonmember may be used.
#' @inheritParams hg_sample_gnp
#' @return A `net_hypergraph`; `$blocks` records each node's planted block.
#' @examples
#' P <- matrix(c(.5, .05, .05, .5), 2, 2)
#' h <- hg_sample_sbm(P = P, block_sizes = c(10, 10), d = 3, seed = 1)
#' as.data.frame(h, what = "nodes")
#' @export
hg_sample_sbm <- function(n = NULL, P, block_sizes, d, impurity = 0L,
                          variable_size = FALSE, absolute_purity = TRUE,
                          seed = NULL) {
  if (!is.numeric(block_sizes) || !length(block_sizes) ||
      any(!is.finite(block_sizes)) || any(block_sizes < 1) ||
      any(block_sizes != as.integer(block_sizes))) {
    stop("`block_sizes` must contain positive integers.", call. = FALSE)
  }
  block_sizes <- as.integer(block_sizes)
  n_blocks <- length(block_sizes)
  n <- if (is.null(n)) sum(block_sizes) else .hgr_count(n, "n", 2L)
  if (n != sum(block_sizes)) {
    stop("`n` must equal sum(`block_sizes`).", call. = FALSE)
  }
  if (!is.matrix(P) || !is.numeric(P) ||
      !identical(dim(P), c(n_blocks, n_blocks)) || anyNA(P) ||
      any(P < 0 | P > 1) || !isTRUE(all.equal(P, t(P)))) {
    stop("`P` must be a symmetric probability matrix matching the blocks.",
         call. = FALSE)
  }
  if (!is.numeric(d) || !length(d) || any(!is.finite(d)) || any(d < 0)) {
    stop("`d` must contain non-negative finite values.", call. = FALSE)
  }
  if (!variable_size && any(d < 2 | d != as.integer(d))) {
    stop("fixed `d` values must be integers >= 2.", call. = FALSE)
  }
  impurity <- .hgr_count(impurity, "impurity", 0L)
  if (!is.logical(variable_size) || length(variable_size) != 1L ||
      !is.logical(absolute_purity) || length(absolute_purity) != 1L) {
    stop("`variable_size` and `absolute_purity` must be single logicals.",
         call. = FALSE)
  }
  .hgr_seed(seed)

  blocks <- rep.int(seq_len(n_blocks), block_sizes)
  pairs <- utils::combn(n, 2L)
  probs <- P[cbind(blocks[pairs[1L, ]], blocks[pairs[2L, ]])]
  keep <- stats::rbinom(ncol(pairs), 1L, probs) == 1L
  pairs <- pairs[, keep, drop = FALSE]
  m <- ncol(pairs)
  target <- if (variable_size) {
    2L + stats::rpois(m, rep(d, length.out = m))
  } else {
    rep(as.integer(d), length.out = m)
  }
  if (any(target > n)) stop("sampled hyperedge size exceeds `n`.", call. = FALSE)

  edges <- lapply(seq_len(m), function(j) {
    endpoints <- pairs[, j]
    endpoint_blocks <- unique(blocks[endpoints])
    eligible <- setdiff(which(blocks %in% endpoint_blocks), endpoints)
    need <- target[j] - 2L
    if (need > length(eligible)) {
      stop("not enough nodes in the endpoint blocks for hyperedge size `d`.",
           call. = FALSE)
    }
    added <- if (need) sample(eligible, need) else integer(0)
    edge <- c(endpoints, added)
    replace_n <- min(length(added), impurity)
    # As in HyperG, cross-block dyads in a two-block model cannot receive an
    # outside-block impurity under the absolute-purity rule.
    if (replace_n && !(n_blocks == 2L && length(endpoint_blocks) == 2L)) {
      remove <- sample(added, replace_n)
      candidates <- if (absolute_purity) {
        which(!blocks %in% endpoint_blocks)
      } else {
        setdiff(seq_len(n), edge)
      }
      if (length(candidates) < replace_n) {
        stop("not enough eligible nodes for the requested `impurity`.",
             call. = FALSE)
      }
      edge <- c(setdiff(edge, remove), sample(candidates, replace_n))
    }
    sort(edge)
  })
  incidence <- .hgr_edges_to_incidence(edges, n)
  out <- .hgr_from_incidence(
    incidence, "hg_sample_sbm",
    list(model = "sbm", n = n, P = P, block_sizes = block_sizes, d = d,
         impurity = impurity, variable_size = variable_size,
         absolute_purity = absolute_purity, seed = seed)
  )
  out$blocks <- stats::setNames(blocks, out$nodes)
  out
}

#' Sample Uniform or Regular Random Hypergraphs
#'
#' `hg_sample_uniform()` makes every hyperedge contain exactly `k` distinct
#' nodes. `hg_sample_regular()` makes every node incident to exactly `k`
#' distinct hyperedges. Sampling probabilities can be unequal.
#'
#' @param k Hyperedge size for the uniform model; node degree for the regular
#'   model.
#' @param prob Sampling weights: length `n` for the uniform model and length
#'   `m` for the regular model. `NULL` is uniform.
#' @inheritParams hg_sample_gnp
#' @return A `net_hypergraph`.
#' @examples
#' u <- hg_sample_uniform(20, 8, k = 3, seed = 1)
#' r <- hg_sample_regular(20, 8, k = 2, seed = 1)
#' hg_measures(u, what = "distribution", measure = "size")
#' hg_measures(r, what = "distribution", measure = "hyperdegree")
#' @export
hg_sample_uniform <- function(n, m, k, prob = NULL, seed = NULL) {
  n <- .hgr_count(n, "n", 1L)
  m <- .hgr_count(m, "m", 0L)
  k <- .hgr_count(k, "k", 0L)
  if (k > n) stop("`k` cannot exceed `n` for a uniform hypergraph.",
                  call. = FALSE)
  prob <- .hgr_sampling_prob(prob, n, k, "n")
  .hgr_seed(seed)
  incidence <- matrix(0L, n, m,
                      dimnames = list(paste0("V", seq_len(n)),
                                      if (m) paste0("h", seq_len(m)) else
                                        character()))
  for (j in seq_len(m)) incidence[sample.int(n, k, prob = prob), j] <- 1L
  .hgr_from_incidence(
    incidence, "hg_sample_uniform",
    list(model = "uniform", n = n, m = m, k = k, prob = prob, seed = seed)
  )
}

#' @rdname hg_sample_uniform
#' @export
hg_sample_regular <- function(n, m, k, prob = NULL, seed = NULL) {
  n <- .hgr_count(n, "n", 1L)
  m <- .hgr_count(m, "m", 1L)
  k <- .hgr_count(k, "k", 0L)
  if (k > m) stop("`k` cannot exceed `m` for a regular hypergraph.",
                  call. = FALSE)
  prob <- .hgr_sampling_prob(prob, m, k, "m")
  .hgr_seed(seed)
  incidence <- matrix(0L, n, m,
                      dimnames = list(paste0("V", seq_len(n)),
                                      paste0("h", seq_len(m))))
  for (i in seq_len(n)) incidence[i, sample.int(m, k, prob = prob)] <- 1L
  .hgr_from_incidence(
    incidence, "hg_sample_regular",
    list(model = "regular", n = n, m = m, k = k, prob = prob, seed = seed)
  )
}

# Descriptive aliases coexist with the concise hg_* family without taking
# HyperG's sample_* names.
#' @rdname hg_sample_gnp
#' @export
hypergraph_sample_gnp <- hg_sample_gnp
#' @rdname hg_sample_sbm
#' @export
hypergraph_sample_sbm <- hg_sample_sbm
#' @rdname hg_sample_uniform
#' @export
hypergraph_sample_uniform <- hg_sample_uniform
#' @rdname hg_sample_uniform
#' @export
hypergraph_sample_regular <- hg_sample_regular

#' @noRd
.hgr_count <- function(x, name, minimum) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
      x < minimum || x != as.integer(x)) {
    stop(sprintf("`%s` must be one integer >= %d.", name, minimum),
         call. = FALSE)
  }
  as.integer(x)
}

#' @noRd
.hgr_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x > 1) {
    stop(sprintf("`%s` must be one probability in [0, 1].", name),
         call. = FALSE)
  }
  invisible(x)
}

#' @noRd
.hgr_sampling_prob <- function(prob, size, k, size_name) {
  # Explicit equal weights preserve seeded parity with HyperG's calls to
  # sample(..., prob = rep(1 / size, size)).
  if (is.null(prob)) return(rep(1 / size, size))
  if (!is.numeric(prob) || length(prob) != size || any(!is.finite(prob)) ||
      any(prob < 0) || sum(prob) <= 0 || sum(prob > 0) < k) {
    stop(sprintf(paste0("`prob` must contain %s non-negative finite weights ",
                        "with at least `k` positive values."), size_name),
         call. = FALSE)
  }
  prob
}

# Set a local RNG scope and restore the exact caller state at function exit.
#' @noRd
.hgr_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  seed <- .hgr_count(seed, "seed", 0L)
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = globalenv()) else NULL
  caller <- parent.frame()
  do.call(on.exit, list(substitute({
    if (HAD) {
      assign(".Random.seed", OLD, envir = globalenv())
    } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
      rm(".Random.seed", envir = globalenv())
    }
  }, list(HAD = had_seed, OLD = old_seed)), add = TRUE), envir = caller)
  set.seed(seed)
  invisible(NULL)
}

#' @noRd
.hgr_edges_to_incidence <- function(edges, n) {
  m <- length(edges)
  out <- matrix(0L, n, m,
                dimnames = list(paste0("V", seq_len(n)),
                                if (m) paste0("h", seq_len(m)) else
                                  character()))
  if (m && sum(lengths(edges))) {
    out[cbind(unlist(edges, use.names = FALSE),
              rep.int(seq_len(m), lengths(edges)))] <- 1L
  }
  out
}

#' @noRd
.hgr_from_incidence <- function(incidence, source, params) {
  incidence <- as.matrix(incidence)
  storage.mode(incidence) <- "integer"
  n <- nrow(incidence)
  m <- ncol(incidence)
  rownames(incidence) <- rownames(incidence) %||% paste0("V", seq_len(n))
  colnames(incidence) <- colnames(incidence) %||%
    if (m) paste0("h", seq_len(m)) else character()
  edges <- lapply(seq_len(m), function(j) which(incidence[, j] > 0L))
  sizes <- lengths(edges)
  size_dist <- if (length(sizes)) {
    tab <- table(sizes)
    stats::setNames(as.integer(tab), paste0("size_", names(tab)))
  } else integer(0)
  params$source <- source
  structure(list(
    hyperedges = edges,
    incidence = incidence,
    nodes = rownames(incidence),
    n_nodes = n,
    n_hyperedges = m,
    size_distribution = size_dist,
    params = params
  ), class = "net_hypergraph")
}
