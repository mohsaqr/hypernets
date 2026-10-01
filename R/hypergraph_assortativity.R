# Degree assortativity of a hypergraph (ROADMAP Phase 2, item 3).
#   hg_assortativity()      Chodrow (2020), Sec. 4.2: the correlation of the
#                           degrees of two nodes chosen from each hyperedge by
#                           a choice function (uniform / top-2 / top-bottom)
#   hg_degree_correlation() Lotito et al. (2023): the Pearson correlation of
#                           the order-specific hyperdegree sequences
# Oracles (local_testing_and_equivalence/test-oracle-assortativity.R):
# XGI 0.10.2 degree_assortativity (exact top-2; exact uniform on uniform
# hypergraphs; its Monte Carlo sampler for uniform on non-uniform ones),
# igraph::assortativity_degree on a simple graph, an independent Python
# implementation of Chodrow's top-bottom definition, and hypergraphx 1.8.0
# degree_correlation.

# (node index, edge index) of every membership, dense or sparse.
.hg_membership_pairs <- function(membership) {
  if (methods::is(membership, "sparseMatrix")) {
    triplet <- Matrix::summary(methods::as(membership, "CsparseMatrix"))
    pairs <- data.frame(node = triplet$i, edge = triplet$j)
  } else {
    idx <- which(membership > 0, arr.ind = TRUE)
    pairs <- data.frame(node = idx[, 1L], edge = idx[, 2L])
  }
  pairs[order(pairs$edge, pairs$node), , drop = FALSE]
}

# Correlation of a pair distribution given its moments; NA when a marginal
# has zero variance (the coefficient is undefined, not 0).
.hg_moment_correlation <- function(mean_1, mean_2, second_1, second_2, cross) {
  var_1 <- second_1 - mean_1^2
  var_2 <- second_2 - mean_2^2
  scale <- max(abs(c(second_1, second_2)), 1)
  if (var_1 <= sqrt(.Machine$double.eps) * scale ||
      var_2 <= sqrt(.Machine$double.eps) * scale) {
    return(NA_real_)
  }
  (cross - mean_1 * mean_2) / sqrt(var_1 * var_2)
}

# One assortativity coefficient for one choice function. `x` is the per-node
# score (degree or degree rank); `members` a list of node-index vectors, one
# per hyperedge of size >= 2. Every hyperedge carries equal weight, as in
# Chodrow's averages over edges.
.hg_assort_one <- function(x, members, type) {
  if (identical(type, "uniform")) {
    # Uniform unordered pair per edge, both orders equally likely: per-edge
    # moments in closed form, E[XY] = (S1^2 - S2) / (m (m - 1)).
    m <- lengths(members)
    s1 <- vapply(members, \(v) sum(x[v]), numeric(1))
    s2 <- vapply(members, \(v) sum(x[v]^2), numeric(1))
    first <- mean(s1 / m)
    second <- mean(s2 / m)
    cross <- mean((s1^2 - s2) / (m * (m - 1)))
    return(.hg_moment_correlation(first, first, second, second, cross))
  }
  if (identical(type, "top_2")) {
    # The two largest scores of each edge, as an unordered pair.
    top <- vapply(members, \(v) sort(x[v], decreasing = TRUE)[seq_len(2L)],
                  numeric(2))
    a <- top[1L, ]
    b <- top[2L, ]
    first <- mean((a + b) / 2)
    second <- mean((a^2 + b^2) / 2)
    return(.hg_moment_correlation(first, first, second, second, mean(a * b)))
  }
  # top_bottom: the ordered pair (largest, smallest) score of each edge.
  high <- vapply(members, \(v) max(x[v]), numeric(1))
  low <- vapply(members, \(v) min(x[v]), numeric(1))
  .hg_moment_correlation(mean(high), mean(low), mean(high^2), mean(low^2),
                         mean(high * low))
}

#' Degree assortativity of a hypergraph
#'
#' Measures whether nodes of similar hyperdegree share hyperedges, with the
#' three choice functions of Chodrow (2020, Sec. 4.2). From every hyperedge
#' with at least two members a pair of nodes is chosen, and the coefficient
#' is the correlation of their scores over hyperedges, each hyperedge
#' counting once:
#'
#' * `"uniform"` -- a uniformly random pair of distinct members (computed
#'   exactly, as the expectation over the random choice; the pair is
#'   unordered, so both orders enter with equal weight);
#' * `"top_2"` -- the two members of largest degree (unordered);
#' * `"top_bottom"` -- the ordered pair (largest-degree member,
#'   smallest-degree member), measuring whether the least-connected member
#'   of a group is better connected when its best-connected member is.
#'
#' @details
#' The score of a node is its hyperdegree (the number of hyperedges
#' containing it). With `scale = "rank"` (default) the score is the node's
#' rank by hyperdegree among all nodes (ties averaged), which gives
#' Chodrow's "generalized Spearman" coefficient; `scale = "degree"` uses the
#' raw hyperdegree (a Pearson coefficient, as XGI does).
#'
#' On a graph (every hyperedge of size two) the three choice functions give
#' the same pairs, and `type = "uniform", scale = "degree"` is Newman's
#' degree assortativity (tested against `igraph::assortativity_degree()`).
#'
#' **Convention differences with XGI 0.10.2** (`degree_assortativity(exact =
#' TRUE)`), documented and tested: XGI's exact `"uniform"` weights every
#' hyperedge by its number of member pairs \eqn{m(m-1)} instead of counting
#' each hyperedge once, so the two agree only on uniform hypergraphs (XGI's
#' Monte Carlo sampler, `exact = FALSE`, samples Chodrow's edge-weighted
#' quantity and agrees within sampling error); XGI symmetrizes
#' `"top-bottom"` (both orders), whereas Chodrow's definition keeps the
#' order (largest, smallest); XGI returns 0 where the coefficient is
#' undefined. `"top_2"` agrees with XGI exactly.
#'
#' The coefficient is undefined (`NA`) when either chosen score has zero
#' variance across hyperedges, e.g. on a regular hypergraph.
#'
#' **Zero is not the null value for `"top_2"` and `"top_bottom"`.** Both
#' pair order statistics of a hyperedge (its two largest scores, or its
#' largest and smallest), which are correlated even when members are drawn
#' at random, so random hypergraphs give positive values (Chodrow 2020
#' compares against a null model for this reason). Judge them against
#' coefficients on degree-matched random hypergraphs, not against 0;
#' `"uniform"` is centred near 0 under random mixing.
#'
#' @param hg A hypernets `net_hg`, dense or sparse.
#' @param type One or more of `"uniform"` (default), `"top_2"`,
#'   `"top_bottom"`.
#' @param scale `"rank"` (default; Chodrow's generalized Spearman) or
#'   `"degree"` (Pearson on raw hyperdegrees).
#'
#' @return A base `data.frame`, one row per `type`, with columns `type`,
#'   `scale`, `assortativity` (in `[-1, 1]`, `NA` when undefined) and
#'   `n_edges` (the hyperedges of size at least two that entered).
#'
#' @section Conditions: Raises `hypernets_bad_input` when `hg` is not a
#'   `net_hg` or has no hyperedge with at least two members.
#'
#' @references
#' Chodrow, P. S. (2020). Configuration models of random hypergraphs.
#' \emph{Journal of Complex Networks}, 8(3), cnaa018.
#' \doi{10.1093/comnet/cnaa018}
#'
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
#' Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. \emph{Physics Reports}, 874, 1-92.
#' \doi{10.1016/j.physrep.2020.05.004}
#'
#' @seealso [hg_degree_correlation()], [hg_transitivity()].
#'
#' @examples
#' groups <- data.frame(
#'   member = c("a", "b", "c", "a", "d", "b", "d", "e", "c", "d", "e", "f"),
#'   group  = c("g1", "g1", "g1", "g2", "g2", "g3", "g3", "g3",
#'              "g4", "g4", "g4", "g4")
#' )
#' hg <- group_hypergraph(groups, actor = "member", group = "group")
#' hg_assortativity(hg, type = c("uniform", "top_2", "top_bottom"))
#' @export
hg_assortativity <- function(hg, type = "uniform",
                             scale = c("rank", "degree")) {
  .thg_check_hg(hg)
  type <- match.arg(type, choices = c("uniform", "top_2", "top_bottom"),
                    several.ok = TRUE)
  scale <- match.arg(scale)
  membership <- .hg_membership(hg)
  degree <- as.numeric(Matrix::rowSums(membership))
  pairs <- .hg_membership_pairs(membership)
  members <- split(pairs$node, factor(pairs$edge,
                                      levels = seq_len(ncol(membership))))
  members <- members[lengths(members) >= 2L]
  if (length(members) == 0L) {
    stop(errorCondition(
      paste("degree assortativity needs at least one hyperedge",
            "with two or more members"),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  x <- if (identical(scale, "rank")) {
    rank(degree, ties.method = "average")
  } else {
    degree
  }
  data.frame(
    type = type,
    scale = scale,
    assortativity = vapply(type, \(t) .hg_assort_one(x, members, t),
                           numeric(1)),
    n_edges = length(members),
    row.names = NULL
  )
}

#' Correlation of hyperdegrees across interaction orders
#'
#' For every pair of hyperedge sizes \eqn{s < t} present in the hypergraph,
#' the Pearson correlation across nodes between the size-\eqn{s} hyperdegree
#' (the number of size-\eqn{s} hyperedges containing the node) and the
#' size-\eqn{t} hyperdegree, as reported by Lotito et al. (2023) for the
#' high-school contact hypergraph (\eqn{\rho_{2,3}}, \eqn{\rho_{2,4}},
#' \eqn{\rho_{3,4}}). A high value means that nodes active in pairs are also
#' active in larger groups.
#'
#' @details
#' Every node of `hg` enters every correlation (with degree 0 at the sizes
#' it does not take part in). Sizes below two are ignored, as in
#' hypergraphx. A pair is `NA` when one of the two sequences is constant.
#' hypergraphx 1.8.0 `degree_correlation()` returns the full
#' `(max_size - 1)`-square matrix, including sizes with no hyperedge (all
#' `NaN`); hypernets returns one row per pair of sizes that occur. Tested
#' against hypergraphx at `1e-12`.
#'
#' @param hg A hypernets `net_hg`, dense or sparse.
#'
#' @return A base `data.frame`, one row per pair of present hyperedge sizes
#'   `size_1 < size_2`, with columns `size_1`, `size_2`, `correlation`
#'   (Pearson, `NA` when undefined) and `n_nodes`.
#'
#' @section Conditions: Raises `hypernets_bad_input` when `hg` is not a
#'   `net_hg` or has fewer than two distinct hyperedge sizes of at least
#'   two.
#'
#' @references
#' Lotito, Q. F., Contisciani, M., De Bacco, C., Di Gaetano, L., Gallo, L.,
#' Montresor, A., Musciotto, F., Ruggeri, N., & Battiston, F. (2023).
#' Hypergraphx: a library for higher-order network analysis. \emph{Journal
#' of Complex Networks}, 11(3), cnad019. \doi{10.1093/comnet/cnad019}
#'
#' @seealso [hg_assortativity()].
#'
#' @examples
#' groups <- data.frame(
#'   member = c("a", "b", "a", "c", "b", "c", "a", "d",
#'              "a", "b", "c", "a", "b", "d", "c", "d", "e"),
#'   group  = c("p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4",
#'              "t1", "t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3")
#' )
#' hg <- group_hypergraph(groups, actor = "member", group = "group")
#' hg_degree_correlation(hg)
#' @export
hg_degree_correlation <- function(hg) {
  .thg_check_hg(hg)
  membership <- .hg_membership(hg)
  sizes <- as.numeric(Matrix::colSums(membership))
  present <- sort(unique(sizes[sizes >= 2]))
  if (length(present) < 2L) {
    stop(errorCondition(
      paste("cross-order degree correlation needs at least two distinct",
            "hyperedge sizes >= 2"),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  by_size <- vapply(present, \(s) {
    as.numeric(membership %*% as.numeric(sizes == s))
  }, numeric(nrow(membership)))
  size_pairs <- utils::combn(length(present), 2L)
  correlation <- vapply(seq_len(ncol(size_pairs)), \(p) {
    a <- by_size[, size_pairs[1L, p]]
    b <- by_size[, size_pairs[2L, p]]
    if (length(unique(a)) < 2L || length(unique(b)) < 2L) return(NA_real_)
    stats::cor(a, b)
  }, numeric(1))
  data.frame(
    size_1 = as.integer(present[size_pairs[1L, ]]),
    size_2 = as.integer(present[size_pairs[2L, ]]),
    correlation = correlation,
    n_nodes = nrow(membership),
    row.names = NULL
  )
}
