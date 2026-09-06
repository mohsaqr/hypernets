# Degree-preserving null-model test for hypergraph structure. The null
# fixes both margins of the binary membership matrix (every vertex keeps its
# hyperdegree, every hyperedge its size) and randomizes which memberships
# occur, via checkerboard swaps (the bipartite swap null of Gotelli 2000).
# Observed structure beyond that null is structure the degree sequences
# alone cannot explain.

# One MCMC chain of checkerboard swap attempts on a binary matrix. A swap
# picks two rows and two columns showing a [1,0;0,1] or [0,1;1,0] pattern
# and flips it, preserving all row and column sums exactly.
.thg_swap_chain <- function(m, attempts) {
  n_row <- nrow(m)
  n_col <- ncol(m)
  rows_1 <- sample.int(n_row, attempts, replace = TRUE)
  rows_2 <- sample.int(n_row, attempts, replace = TRUE)
  cols_1 <- sample.int(n_col, attempts, replace = TRUE)
  cols_2 <- sample.int(n_col, attempts, replace = TRUE)
  # sequential MCMC: each accepted swap changes the state the next attempt
  # sees, so the loop cannot be vectorized; candidate indices are drawn in
  # one vectorized batch above
  for (i in seq_len(attempts)) {
    r1 <- rows_1[i]; r2 <- rows_2[i]; c1 <- cols_1[i]; c2 <- cols_2[i]
    if (r1 == r2 || c1 == c2) next
    a <- m[r1, c1]; b <- m[r1, c2]; d <- m[r2, c1]; e <- m[r2, c2]
    if (a + e == 2L && b + d == 0L) {
      m[r1, c1] <- 0L; m[r2, c2] <- 0L; m[r1, c2] <- 1L; m[r2, c1] <- 1L
    } else if (a + e == 0L && b + d == 2L) {
      m[r1, c1] <- 1L; m[r2, c2] <- 1L; m[r1, c2] <- 0L; m[r2, c1] <- 0L
    }
  }
  m
}

# One configuration-model draw: vertex stubs (one per membership) are
# matched uniformly at random to hyperedge slots (Chodrow 2020). Margins are
# preserved only up to collapse -- if a vertex's stubs land twice in the same
# hyperedge the repeat collapses, costing that vertex a degree and that
# hyperedge a slot -- which is exactly how this null differs from the swap
# null, where both margins are exact.
.thg_configuration_draw <- function(m) {
  stubs <- rep.int(seq_len(nrow(m)), rowSums(m))
  slots <- rep.int(seq_len(ncol(m)), colSums(m))
  out <- matrix(0L, nrow(m), ncol(m), dimnames = dimnames(m))
  out[cbind(sample(stubs), slots)] <- 1L
  out
}

# The degree-ordered assignment of Coupette et al. (2024, footnote 7): nodes
# in decreasing degree each choose their hyperedges uniformly, without
# replacement, among the hyperedges that still have a free slot, so a node
# never repeats within a hyperedge and both margins hold exactly. A draw
# that strands a node without enough open hyperedges is redrawn.
.thg_assignment_draw <- function(m, max_attempts = 100L) {
  degree <- rowSums(m)
  capacity_all <- colSums(m)
  order_nodes <- order(-degree, seq_along(degree))
  for (attempt in seq_len(max_attempts)) {
    out <- matrix(0L, nrow(m), ncol(m), dimnames = dimnames(m))
    capacity <- capacity_all
    ok <- TRUE
    # sequential by construction: each node's choice depends on the slots
    # the previous nodes left open
    for (i in order_nodes) {
      k <- degree[i]
      if (k == 0) next
      open <- which(capacity > 0)
      if (length(open) < k) {
        ok <- FALSE
        break
      }
      chosen <- if (length(open) == 1L) open else sample(open, k)
      out[i, chosen] <- 1L
      capacity[chosen] <- capacity[chosen] - 1
    }
    if (ok) return(out)
  }
  stop(errorCondition(
    "the degree-ordered assignment stranded a node in every attempt",
    class = "honets_no_converge", call = NULL
  ))
}

# Statistics on the binary membership, computed directly from sparse
# cross-products. The delegated path (rebuild the hypergraph, take
# hg_measures(what = "overlap")) materializes a table quadratic in the
# hyperedge count on EVERY null replicate and dominated the test's
# runtime (measured: 4.4 s per replicate on 165 x 4043; the direct path
# is ~50 ms). Numerical identity with the delegated path is asserted in
# tests/testthat/test-null.R.
.thg_null_statistics <- function(m, statistic) {
  m_sparse <- Matrix::Matrix(m, sparse = TRUE) * 1
  sizes <- as.numeric(Matrix::colSums(m_sparse))
  n_nodes <- nrow(m)
  n_edges <- ncol(m)
  vapply(statistic, \(s) {
    switch(s,
      density = sum(sizes) / (n_nodes * n_edges),
      avg_edge_size = mean(sizes),
      # Hyperedges whose member set already occurred: m minus distinct sets.
      repeated_edges = {
        members <- .thg_edge_members(m_sparse)
        n_edges - length(unique(vapply(members, paste, collapse = "\r",
                                       character(1L))))
      },
      # Member pairs counted with multiplicity minus distinct co-occurring
      # pairs: how many collaborations repeat an earlier one.
      repeated_pairs = {
        co <- methods::as(Matrix::tcrossprod(m_sparse), "TsparseMatrix")
        distinct <- sum(co@i < co@j & co@x != 0)
        sum(choose(sizes, 2)) - distinct
      },
      pairwise_participation = {
        co <- methods::as(Matrix::tcrossprod(m_sparse), "TsparseMatrix")
        sharing <- sum(co@i < co@j & co@x > 0)
        sharing / (n_nodes * (n_nodes - 1) / 2)
      },
      avg_jaccard = {
        co <- methods::as(Matrix::crossprod(m_sparse), "TsparseMatrix")
        keep <- co@i < co@j & co@x > 0
        inter <- co@x[keep]
        s_i <- sizes[co@i[keep] + 1L]
        s_j <- sizes[co@j[keep] + 1L]
        sum(inter / (s_i + s_j - inter)) / (n_edges * (n_edges - 1) / 2)
      },
      stop("unknown statistic: ", s)
    )
  }, numeric(1))
}

#' Degree-preserving null-model test of hypergraph structure
#'
#' Tests whether an observed structural statistic exceeds what the degree
#' sequences alone imply. The null model fixes both margins of the binary
#' membership (every vertex keeps its hyperdegree, every hyperedge its
#' size) and randomizes the memberships by checkerboard swaps (Gotelli
#' 2000) -- a sequential MCMC with burn-in `10 * nnz` swap attempts and
#' thinning `nnz` between samples, `nnz` being the number of memberships.
#' Statistics are evaluated on the binarized hypergraph (weights carry no
#' meaning under this null), through the same delegated measures as
#' [hg_measures()].
#'
#' The permutation p-value is `(1 + extreme) / (n + 1)` (Phipson & Smyth
#' 2010), two-sided by default around the null mean; `null_lo`/`null_hi`
#' are the 2.5% and 97.5% null quantiles -- report them with the observed
#' value, not the p-value alone.
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any honets
#'   `net_hypergraph`.
#' @param statistic Statistics to test; any of `"pairwise_participation"`,
#'   `"density"`, `"avg_edge_size"`, `"avg_jaccard"` (mean pairwise edge
#'   Jaccard), `"repeated_edges"` (hyperedges whose member set already
#'   occurred: the paper's repeated collaboration trios) and
#'   `"repeated_pairs"` (member pairs counted with multiplicity minus the
#'   distinct co-occurring pairs: its repeated collaboration pairs). Several
#'   allowed.
#' @param method Which null. `"swap"` (default) is the checkerboard swap
#'   chain above: both margins are preserved *exactly*, and no vertex may
#'   repeat within a hyperedge. `"configuration"` draws independent
#'   stub-matchings (Chodrow 2020), the standard higher-order configuration
#'   model: vertex stubs are matched uniformly at random to hyperedge slots,
#'   so margins hold only up to collapse (a vertex whose stubs land twice in
#'   one hyperedge loses a degree). `"assignment"` is the degree-ordered
#'   assignment Coupette et al. (2024) use for their repeated-collaboration
#'   test: vertices in decreasing degree each choose their hyperedges without
#'   replacement among those with a free slot, so both margins hold exactly
#'   and no vertex repeats within a hyperedge, at the price of a draw that is
#'   not uniform over configurations. `"swap"` is the stricter null; use
#'   `"configuration"` to compare against the higher-order network
#'   literature, which reports it.
#' @param n Number of null samples (default `199L`).
#' @param seed Seed for the null draws; set it for a reproducible test.
#'   The global RNG state is restored on exit.
#' @param alternative `"two_sided"` (default), `"greater"`, or `"less"`.
#' @return A base `data.frame`, one row per statistic: `statistic`,
#'   `observed`, `null_mean`, `null_lo`, `null_hi`, `z`, `p_value`, `n`,
#'   `method`.
#' @section Conditions: Raises `honets_bad_input` for broken contracts.
#'   `method = "configuration"` signals a `honets_configuration_collapse`
#'   warning when stub matching collapses more than 1% of memberships on
#'   average, which happens whenever hyperedges are large relative to the
#'   vertex set -- the usual case in the document orientation.
#' @references
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence
#' patterns. *Ecology*, 81(9).
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero. *Statistical Applications in Genetics and Molecular Biology*, 9(1).
#'
#' Chodrow, P. S. (2020). Configuration models of random hypergraphs.
#' *Journal of Complex Networks*, 8(3), cnaa018.
#' \doi{10.1093/comnet/cnaa018}
#' @examples
#' hg <- text_hypergraph(c(
#'   a = "salt and soup and onions",
#'   b = "soup and salt",
#'   c = "stars and sky and salt",
#'   d = "stars and sky"
#' ))
#' hg_null_test(hg, statistic = "pairwise_participation", n = 49, seed = 1)
#' @export
hg_null_test <- function(hg,
                         statistic = c("pairwise_participation", "density",
                                       "avg_edge_size", "avg_jaccard",
                                       "repeated_edges", "repeated_pairs"),
                         method = c("swap", "configuration", "assignment"),
                         n = 199L, seed = NULL,
                         alternative = c("two_sided", "greater", "less")) {
  .thg_check_hg(hg)
  if (.thg_is_sparse(hg)) {
    stop(errorCondition(
      "the null test needs the dense representation for now",
      class = "honets_sparse_unsupported", call = NULL
    ))
  }
  statistic <- match.arg(statistic, several.ok = TRUE)
  method <- match.arg(method)
  alternative <- match.arg(alternative)
  stopifnot(
    "`n` must be a single count >= 19" =
      length(n) == 1L && is.finite(n) && n >= 19
  )

  if (!is.null(seed)) {
    old_seed <- if (exists(".Random.seed", envir = globalenv())) {
      get(".Random.seed", envir = globalenv())
    } else {
      NULL
    }
    on.exit({
      if (is.null(old_seed)) {
        rm(".Random.seed", envir = globalenv())
      } else {
        assign(".Random.seed", old_seed, envir = globalenv())
      }
    }, add = TRUE)
    set.seed(seed)
  }

  membership <- (hg$incidence > 0) * 1L
  nnz <- sum(membership)
  observed <- .thg_null_statistics(membership, statistic)

  draws <- if (identical(method, "swap")) {
    state <- .thg_swap_chain(membership, attempts = 10L * nnz)
    vapply(seq_len(n), \(i) {
      state <<- .thg_swap_chain(state, attempts = nnz)
      .thg_null_statistics(state, statistic)
    }, numeric(length(statistic)))
  } else if (identical(method, "assignment")) {
    vapply(seq_len(n), \(i) {
      .thg_null_statistics(.thg_assignment_draw(membership), statistic)
    }, numeric(length(statistic)))
  } else {
    retained <- numeric(n)
    stats <- vapply(seq_len(n), \(i) {
      draw <- .thg_configuration_draw(membership)
      retained[i] <<- sum(draw)
      .thg_null_statistics(draw, statistic)
    }, numeric(length(statistic)))
    # Collapse is not a defect of the draw but it does weaken the null, and
    # it is severe whenever hyperedges are large relative to the vertex set
    # -- the usual case for document-orientation text hypergraphs. Surface it
    # rather than let the margins quietly shrink.
    lost <- 1 - mean(retained) / nnz
    if (lost > 0.01) {
      warning(warningCondition(
        sprintf(paste0("the configuration draws collapse %.1f%% of memberships ",
                       "on average, so vertex degrees and hyperedge sizes are ",
                       "preserved only approximately; method = \"swap\" preserves ",
                       "both margins exactly"), 100 * lost),
        class = "honets_configuration_collapse"
      ))
    }
    stats
  }
  draws <- matrix(draws, nrow = length(statistic))

  rows <- lapply(seq_along(statistic), \(i) {
    null_draws <- draws[i, ]
    null_mean <- mean(null_draws)
    null_sd <- stats::sd(null_draws)
    centered_obs <- abs(observed[i] - null_mean)
    extreme <- switch(alternative,
      two_sided = sum(abs(null_draws - null_mean) >= centered_obs -
                        sqrt(.Machine$double.eps)),
      greater = sum(null_draws >= observed[i] - sqrt(.Machine$double.eps)),
      less = sum(null_draws <= observed[i] + sqrt(.Machine$double.eps))
    )
    data.frame(
      statistic = statistic[i],
      observed = unname(observed[i]),
      null_mean = null_mean,
      null_lo = unname(stats::quantile(null_draws, 0.025)),
      null_hi = unname(stats::quantile(null_draws, 0.975)),
      z = if (null_sd > 0) (observed[i] - null_mean) / null_sd else NA_real_,
      p_value = (1 + extreme) / (n + 1),
      n = as.integer(n),
      method = method
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname hg_null_test
#' @export
hypergraph_null_test <- hg_null_test
