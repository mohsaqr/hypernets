# Local clustering coefficients of a hypergraph (ROADMAP Phase 2, item 3).
# Three published families, each named after its definition:
#   projection      Watts & Strogatz (1998) on the unweighted clique expansion
#                   (Battiston et al. 2020, Sec. III.C: paths through
#                   hyperedges, triangles closed by any hyperedge)
#   extra_overlap   Zhou & Nakhleh (2011) as stated by Klimm et al. (2021)
#   two_node_*      Latapy, Magnien & Del Vecchio (2008), the two-mode
#                   clustering cc(u, v) with union / min / max denominators,
#                   averaged over the node's hypergraph neighbours; applied to
#                   hypernetworks by Gallagher & Goldberg (2013)
# Oracle: XGI 0.10.2 clustering_coefficient / local_clustering_coefficient /
# two_node_clustering_coefficient (local_testing_and_equivalence/
# test-oracle-transitivity-xgi.R). XGI returns 0 where the coefficient is
# undefined; hypernets returns NA (see @details).
#
# Estrada & Rodriguez-Velazquez's (2006) global C2(H) (hyper-triangles over
# 2-paths) is NOT implemented: it has no independent oracle, and the worked
# example in the paper (its Fig. 1, value 0.25 "of the 18 triples") cannot be
# reproduced from the text alone.

# Binary membership (node x hyperedge) as a numeric matrix, dense or
# dgCMatrix according to the stored incidence.
.hg_membership <- function(hg) {
  incidence <- hg$incidence
  if (.thg_is_sparse(hg)) {
    membership <- methods::as(incidence > 0, "CsparseMatrix") * 1
  } else {
    membership <- (incidence > 0) * 1
  }
  dimnames(membership) <- dimnames(incidence)
  membership
}

# Unweighted clique-expansion adjacency (zero diagonal) from a membership.
.hg_projection_adjacency <- function(membership) {
  co <- Matrix::tcrossprod(membership)
  adjacency <- (co > 0) * 1
  Matrix::diag(adjacency) <- 0
  adjacency
}

# Watts-Strogatz local clustering on the unweighted projection:
# closed 2-paths / C(k, 2); NA for k < 2.
.hg_trans_projection <- function(membership) {
  adjacency <- .hg_projection_adjacency(membership)
  k <- as.numeric(Matrix::rowSums(adjacency))
  two_step <- adjacency %*% adjacency
  closed <- as.numeric(Matrix::rowSums(two_step * adjacency)) / 2
  out <- rep(NA_real_, length(k))
  defined <- k >= 2
  out[defined] <- closed[defined] / choose(k[defined], 2)
  out
}

# Extra overlap EO(e, f) for every pair of intersecting hyperedges
# (Klimm et al. 2021, eq. for EO): with D1 = e \ f and D2 = f \ e,
#   EO = (|N(D1) & D2| + |N(D2) & D1|) / (|D1| + |D2|),  EO = 0 when e = f.
# Returned as a symmetric edge x edge matrix with zero diagonal (dense: the
# node-level sum needs every intersecting pair).
.hg_extra_overlap <- function(membership) {
  adjacency <- .hg_projection_adjacency(membership)
  n_edges <- ncol(membership)
  eo <- matrix(0, n_edges, n_edges)
  overlap <- as.matrix(Matrix::crossprod(membership))
  pairs <- which(overlap > 0 & upper.tri(overlap), arr.ind = TRUE)
  if (nrow(pairs) == 0L) return(eo)
  values <- vapply(seq_len(nrow(pairs)), \(p) {
    in_e <- as.numeric(membership[, pairs[p, 1L]]) > 0
    in_f <- as.numeric(membership[, pairs[p, 2L]]) > 0
    d1 <- in_e & !in_f
    d2 <- in_f & !in_e
    size <- sum(d1) + sum(d2)
    if (size == 0) return(0)
    reach_1 <- as.numeric(adjacency %*% d1) > 0
    reach_2 <- as.numeric(adjacency %*% d2) > 0
    (sum(reach_1 & d2) + sum(reach_2 & d1)) / size
  }, numeric(1))
  eo[pairs] <- values
  eo + t(eo)
}

# Node-level extra-overlap clustering: the mean EO over unordered pairs of
# hyperedges containing the node; 0 for hyperdegree 1, NA for hyperdegree 0.
.hg_trans_extra_overlap <- function(membership) {
  eo <- .hg_extra_overlap(membership)
  dense <- as.matrix(membership)
  k <- rowSums(dense)
  # b_i' EO b_i sums EO over ordered pairs of incident hyperedges
  pair_sum <- rowSums((dense %*% eo) * dense)
  out <- rep(NA_real_, length(k))
  out[k == 1] <- 0
  many <- k >= 2
  out[many] <- pair_sum[many] / (k[many] * (k[many] - 1))
  out
}

# Two-node clustering averaged over hypergraph neighbours:
#   cc(u, v) = |M(u) & M(v)| / denom(|M(u)|, |M(v)|), denom = union/min/max
# NA for a node with no neighbour.
.hg_trans_two_node <- function(membership, kind) {
  shared <- Matrix::tcrossprod(membership)
  neighbour <- (shared > 0) * 1
  Matrix::diag(neighbour) <- 0
  shared <- as.matrix(shared)
  neighbour <- as.matrix(neighbour)
  degree <- diag(shared)
  denom <- switch(kind,
    union = outer(degree, degree, `+`) - shared,
    min = outer(degree, degree, pmin),
    max = outer(degree, degree, pmax)
  )
  ratio <- ifelse(neighbour > 0, shared / denom, 0)
  n_neighbour <- rowSums(neighbour)
  out <- rep(NA_real_, length(degree))
  has <- n_neighbour > 0
  out[has] <- rowSums(ratio)[has] / n_neighbour[has]
  out
}

#' Local clustering coefficients of a hypergraph
#'
#' Computes one or more published local clustering coefficients for every
#' node of a [net_hg][network_hypergraph]. Each variant is named after its
#' definition:
#'
#' * `"projection"` -- the Watts & Strogatz (1998) local clustering of the
#'   unweighted clique expansion: the share of pairs of a node's neighbours
#'   that are themselves neighbours. A path through hyperedges is closed by
#'   *any* hyperedge, so a single hyperedge of size three already closes a
#'   triangle (the path-based generalization reviewed by Battiston et al.
#'   2020, Sec. III.C).
#' * `"extra_overlap"` -- Zhou & Nakhleh (2011), in the statement of Klimm,
#'   Deane & Reinert (2021): for every unordered pair of hyperedges \eqn{e,
#'   f} containing the node, the *extra overlap*
#'   \deqn{EO(e, f) = \frac{|N(D_{e,f}) \cap D_{f,e}| + |N(D_{f,e}) \cap
#'   D_{e,f}|}{|D_{e,f}| + |D_{f,e}|},\qquad D_{e,f} = e \setminus f,}
#'   is the share of the non-shared members that are linked to the other
#'   side by a *third* hyperedge (\eqn{N(S)} is the union of the
#'   hypergraph neighbourhoods of \eqn{S}); \eqn{EO = 0} when \eqn{e = f}.
#'   The node's coefficient is the mean EO over its \eqn{k(k-1)/2} pairs of
#'   hyperedges. It reduces to the ordinary clustering coefficient on a graph
#'   and does not count the triangles *inside* a hyperedge.
#' * `"two_node_union"`, `"two_node_min"`, `"two_node_max"` -- the two-mode
#'   clustering of Latapy, Magnien & Del Vecchio (2008), applied to
#'   hypernetworks by Gallagher & Goldberg (2013): for nodes \eqn{u, v} with
#'   hyperedge sets \eqn{M(u), M(v)},
#'   \eqn{cc(u, v) = |M(u) \cap M(v)| / |M(u) \cup M(v)|} (`union`), or with
#'   \eqn{\min} / \eqn{\max(|M(u)|, |M(v)|)} in the denominator; the node's
#'   coefficient averages \eqn{cc(u, v)} over its hypergraph neighbours
#'   \eqn{v} (Latapy et al.'s \eqn{N(N(u))}).
#'
#' @details
#' **Undefined values are `NA`.** A node with fewer than two projection
#' neighbours has no `"projection"` coefficient, a node in no hyperedge has no
#' `"extra_overlap"` coefficient (hyperdegree 1 gives 0, as Klimm et al.
#' define it), and a node without neighbours has no `"two_node_*"`
#' coefficient. XGI 0.10.2 reports 0 in all of these cases; the values agree
#' everywhere else (the oracle test maps `NA` to 0 before comparing).
#'
#' The global clustering coefficient \eqn{C_2(H)} of Estrada &
#' Rodriguez-Velazquez (2006) is not provided (no independent oracle).
#'
#' Duplicate hyperedges are distinct hyperedges (as in XGI): they raise the
#' hyperdegree and, for `"extra_overlap"`, contribute pairs with EO = 0.
#' `"extra_overlap"` forms a dense edge-by-edge matrix and visits every pair
#' of intersecting hyperedges, so it is meant for hypergraphs with at most a
#' few thousand hyperedges; the other variants use sparse products when the
#' incidence is sparse.
#'
#' @param hg A hypernets `net_hg` ([network_hypergraph()],
#'   [group_hypergraph()], [text_hypergraph()], ...), dense or sparse.
#' @param type One or more of `"projection"` (default), `"extra_overlap"`,
#'   `"two_node_union"`, `"two_node_min"`, `"two_node_max"`.
#' @param what `"nodes"` (default; one row per node) or `"summary"` (one row
#'   per `type` with the mean local coefficient over the nodes where it is
#'   defined, the average local clustering used by Chodrow 2020, and the
#'   number of nodes where it is undefined).
#' @param sort_by For `what = "nodes"`, optional `type` to sort by,
#'   descending (ties broken by node name, `NA` last).
#' @param n Return only the first `n` node rows after sorting.
#'
#' @return A base `data.frame`. For `what = "nodes"`: one row per node with
#'   column `node` and one numeric column per requested `type` (in `[0, 1]`,
#'   `NA` where undefined). For `what = "summary"`: one row per `type` with
#'   columns `type`, `mean_local` and `n_undefined`.
#'
#' @section Conditions: Raises `hypernets_bad_input` when `hg` is not a
#'   `net_hg` or `n` is not a count.
#'
#' @references
#' Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
#' 'small-world' networks. \emph{Nature}, 393(6684), 440-442.
#' \doi{10.1038/30918}
#'
#' Zhou, W., & Nakhleh, L. (2011). Properties of metabolic graphs:
#' biological organization or representation artifacts? \emph{BMC
#' Bioinformatics}, 12, 132. \doi{10.1186/1471-2105-12-132}
#'
#' Klimm, F., Deane, C. M., & Reinert, G. (2021). Hypergraphs for predicting
#' essential genes using multiprotein complex data. \emph{Journal of Complex
#' Networks}, 9(2), cnaa028. \doi{10.1093/comnet/cnaa028}
#'
#' Latapy, M., Magnien, C., & Del Vecchio, N. (2008). Basic notions for the
#' analysis of large two-mode networks. \emph{Social Networks}, 30(1),
#' 31-48. \doi{10.1016/j.socnet.2007.04.006}
#'
#' Gallagher, S. R., & Goldberg, D. S. (2013). Clustering coefficients in
#' protein interaction hypernetworks. \emph{Proceedings of the ACM
#' Conference on Bioinformatics, Computational Biology and Biomedical
#' Informatics (BCB '13)}, 552-560. \doi{10.1145/2506583.2506635}
#'
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
#' Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. \emph{Physics Reports}, 874, 1-92.
#' \doi{10.1016/j.physrep.2020.05.004}
#'
#' Chodrow, P. S. (2020). Configuration models of random hypergraphs.
#' \emph{Journal of Complex Networks}, 8(3), cnaa018.
#' \doi{10.1093/comnet/cnaa018}
#'
#' @seealso [hg_assortativity()], [hg_measures()], [hg_clique_expansion()].
#'
#' @examples
#' groups <- data.frame(
#'   member = c("a", "b", "c", "a", "d", "b", "d", "c", "d", "e"),
#'   group  = c("g1", "g1", "g1", "g2", "g2", "g3", "g3", "g4", "g4", "g4")
#' )
#' hg <- group_hypergraph(groups, actor = "member", group = "group")
#' hg_transitivity(hg, type = c("projection", "extra_overlap"))
#' hg_transitivity(hg, type = "two_node_union", what = "summary")
#' @export
hg_transitivity <- function(hg,
                            type = "projection",
                            what = c("nodes", "summary"),
                            sort_by = NULL, n = Inf) {
  .thg_check_hg(hg)
  type <- match.arg(type, choices = c("projection", "extra_overlap",
                                      "two_node_union", "two_node_min",
                                      "two_node_max"), several.ok = TRUE)
  what <- match.arg(what)
  if (!(length(n) == 1L && is.numeric(n) &&
        (is.infinite(n) || (is.finite(n) && n >= 1)))) {
    stop(errorCondition("`n` must be a single count >= 1",
                        class = "hypernets_bad_input", call = NULL))
  }
  membership <- .hg_membership(hg)
  values <- lapply(type, \(t) {
    switch(t,
      projection = .hg_trans_projection(membership),
      extra_overlap = .hg_trans_extra_overlap(membership),
      two_node_union = .hg_trans_two_node(membership, "union"),
      two_node_min = .hg_trans_two_node(membership, "min"),
      two_node_max = .hg_trans_two_node(membership, "max")
    )
  })
  names(values) <- type

  if (identical(what, "summary")) {
    return(data.frame(
      type = type,
      mean_local = vapply(values, \(v) {
        if (all(is.na(v))) NA_real_ else mean(v[!is.na(v)])
      }, numeric(1)),
      n_undefined = vapply(values, \(v) sum(is.na(v)), integer(1)),
      row.names = NULL
    ))
  }

  out <- data.frame(node = rownames(hg$incidence), values,
                    row.names = NULL, check.names = FALSE,
                    stringsAsFactors = FALSE)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, choices = type)
    out <- out[order(-out[[sort_by]], out$node, na.last = TRUE), ,
               drop = FALSE]
    rownames(out) <- NULL
  }
  if (is.finite(n) && n < nrow(out)) out <- out[seq_len(n), , drop = FALSE]
  out
}
