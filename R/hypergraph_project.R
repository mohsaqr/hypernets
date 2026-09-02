# The projection tier: hypergraph -> weighted graph, and hypergraph ->
# s-line graph. Both are pure incidence algebra. Results expose graph-shaped
# matrices and edge lists at the cograph boundary: honets owns these
# hypergraph transformations; cograph owns downstream graph analysis and
# plotting.

# Membership pattern of an incidence matrix, sparse or dense, weights dropped.
.thg_binary <- function(x) (x != 0) * 1

# Scale columns of an incidence matrix by `coef`, sparse or dense. diag() with
# a length-one `coef` would build an identity of that size, hence `nrow=`.
.thg_scale_cols <- function(x, coef) {
  if (methods::is(x, "sparseMatrix")) {
    x %*% Matrix::Diagonal(x = coef)
  } else {
    x %*% diag(x = coef, nrow = length(coef))
  }
}

# Upper triangle of a symmetric weight matrix as a tidy edge list, sorted by
# (from, to) so the row order never depends on storage or input order.
.thg_tidy_pairs <- function(w, labels) {
  if (methods::is(w, "sparseMatrix")) {
    triplet <- methods::as(w, "TsparseMatrix")
    keep <- triplet@i < triplet@j & triplet@x != 0
    from <- labels[triplet@i[keep] + 1L]
    to <- labels[triplet@j[keep] + 1L]
    weight <- triplet@x[keep]
  } else {
    nz <- which(w != 0 & upper.tri(w), arr.ind = TRUE)
    from <- labels[nz[, "row"]]
    to <- labels[nz[, "col"]]
    weight <- as.numeric(w[nz])
  }
  out <- data.frame(from = from, to = to, weight = as.numeric(weight),
                    row.names = NULL)
  out <- out[order(out$from, out$to), , drop = FALSE]
  row.names(out) <- NULL
  out
}

#' Project a hypergraph onto a weighted graph
#'
#' Collapses every hyperedge into pairwise vertex relations, giving the
#' weighted graph that a graph engine (paths, betweenness, communities) can
#' consume. Two weightings are available. `"clique"` is plain co-occurrence,
#' the clique expansion: a hyperedge of size \eqn{|e|} contributes the same
#' amount to each of its \eqn{|e|(|e|-1)/2} pairs, so a single large hyperedge
#' can dominate every downstream measure. `"association"` divides each
#' hyperedge's contribution by \eqn{|e|-1}, following Coupette et al. (2024):
#'
#' \deqn{w(\{u,v\}) = \sum_{e \supseteq \{u,v\}} \frac{1}{|e| - 1}}
#'
#' so the total weight a hyperedge adds around any one of its members is
#' exactly 1, and each vertex's weighted degree in the projection equals the
#' number of hyperedges of size at least two that contain it. That
#' normalisation is what makes the projection a legitimate random-walk
#' operator rather than an arbitrary co-occurrence count.
#'
#' A hyperedge of size one has no pairs, so it contributes nothing and is
#' skipped rather than dividing by zero. Such hyperedges are common in text:
#' a word used in exactly one document is a singleton hyperedge in the
#' document orientation, which is why the degree identity above counts only
#' hyperedges of size at least two.
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any honets
#'   `net_hypergraph`.
#' @param method Weighting. `"clique"` (default) sums incidence products;
#'   `"association"` applies the \eqn{1/(|e|-1)} normalisation above.
#' @param weighted `method = "clique"` only. `TRUE` (default) uses the
#'   incidence weights, `FALSE` their membership pattern. Setting it together
#'   with `method = "association"` is an error, because the association
#'   weighting is defined on hyperedge cardinality and never on the incidence
#'   weights.
#' @param duplicate_edges For `method = "association"`, `"count"` (default)
#'   lets repeated hyperedges contribute repeatedly (the paper's
#'   multi-hypergraph representation); `"collapse"` lets each distinct member
#'   set contribute once (its binary-hypergraph representation).
#' @param self_association Add the paper's source-to-member association term?
#'   Requires one source vertex per hyperedge through `edge_source` or
#'   `hg$edge_data$source`. For source `u`, every membership occurrence of
#'   target `v` contributes `1 / sum_e |e|` over hyperedges sourced by `u`, so
#'   the added incident weight from all of `u`'s citations sums to one.
#' @param edge_source Optional source identifiers: a vector of length
#'   `n_hyperedges`, a named vector keyed by hyperedge, or a two-column data
#'   frame named `edge` and `source`. Used only with `self_association = TRUE`.
#' @param what `"edges"` (default) for the tidy edge list, or `"matrix"` for
#'   the symmetric weight matrix to hand to a graph engine.
#' @return With `what = "edges"`, a base data.frame with one row per
#'   unordered vertex pair of non-zero weight, sorted by `from` then `to`,
#'   with columns `from`, `to` (vertex names) and `weight` (numeric). Vertices
#'   sharing no hyperedge do not appear, so an edgeless hypergraph yields a
#'   zero-row data.frame with those columns. With `what = "matrix"`, the
#'   symmetric `n_nodes` x `n_nodes` weight matrix with zero diagonal and
#'   vertex names as dimnames, sparse if `hg`'s incidence is sparse.
#' @references
#' Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
#' *Philosophical Transactions of the Royal Society A*, 382(2270), 20230141.
#' \doi{10.1098/rsta.2023.0141}
#'
#' Zhou, D., Huang, J., & Schoelkopf, B. (2006). Learning with hypergraphs:
#' clustering, classification, and embedding. *NeurIPS 19*, 1601-1608.
#' @seealso [hg_line_graph()] for the dual projection, onto hyperedges.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
#'                         c = "stars and salt"))
#' hg_project(hg)
#' hg_project(hg, method = "association")
#' @export
hg_project <- function(hg, method = c("clique", "association"),
                       weighted = TRUE, what = c("edges", "matrix"),
                       duplicate_edges = c("count", "collapse"),
                       self_association = FALSE, edge_source = NULL) {
  .thg_check_hg(hg)
  method <- match.arg(method)
  what <- match.arg(what)
  duplicate_edges <- match.arg(duplicate_edges)
  stopifnot("`weighted` must be TRUE or FALSE" =
              length(weighted) == 1L && is.logical(weighted) &&
              !is.na(weighted),
            "`self_association` must be TRUE or FALSE" =
              length(self_association) == 1L && is.logical(self_association) &&
              !is.na(self_association))
  if (identical(method, "association") && !missing(weighted)) {
    stop(errorCondition(
      "`weighted` applies to `method = \"clique\"` only; the association weighting is defined on hyperedge cardinality, not on incidence weights",
      class = "honets_bad_input", call = NULL
    ))
  }
  if (identical(method, "clique") && !identical(duplicate_edges, "count")) {
    .thg_bad_input("`duplicate_edges` applies to `method = \"association\"` only")
  }
  if (identical(method, "clique") && self_association) {
    .thg_bad_input("`self_association` requires `method = \"association\"`")
  }
  incidence <- hg$incidence
  w <- if (identical(method, "clique")) {
    tcrossprod(if (weighted) incidence else .thg_binary(incidence))
  } else {
    b <- .thg_binary(incidence)
    b_projection <- b
    if (identical(duplicate_edges, "collapse") && ncol(b) > 1L) {
      signatures <- vapply(seq_len(ncol(b)), function(j) {
        paste(rownames(b)[which(b[, j] != 0)], collapse = "\r")
      }, character(1L))
      b_projection <- b[, !duplicated(signatures), drop = FALSE]
    }
    cardinality <- Matrix::colSums(b_projection)
    # A singleton hyperedge has no pairs: contribute 0 rather than divide by 0.
    coef <- ifelse(cardinality > 1, 1 / (cardinality - 1), 0)
    tcrossprod(.thg_scale_cols(b_projection, coef), b_projection)
  }
  diag(w) <- 0
  nodes <- rownames(incidence)
  dimnames(w) <- list(nodes, nodes)

  if (self_association) {
    sources <- .thg_resolve_edge_source(hg, edge_source)
    all_nodes <- union(nodes, unique(sources))
    node_index <- match(nodes, all_nodes)
    source_index <- match(sources, all_nodes)
    b_original <- .thg_binary(incidence)
    nz <- which(b_original != 0, arr.ind = TRUE)
    from <- source_index[nz[, "col"]]
    to <- match(nodes[nz[, "row"]], all_nodes)
    keep <- from != to
    from <- from[keep]
    to <- to[keep]
    source_total <- rowsum(rep(1, nrow(nz)), sources[nz[, "col"]],
                           reorder = FALSE)
    denom <- stats::setNames(source_total[, 1L], rownames(source_total))
    contribution <- 1 / unname(denom[sources[nz[, "col"]][keep]])

    sparse_result <- methods::is(w, "sparseMatrix")
    full <- Matrix::Matrix(0, nrow = length(all_nodes), ncol = length(all_nodes),
                           sparse = TRUE,
                           dimnames = list(all_nodes, all_nodes))
    full[node_index, node_index] <- w
    if (length(from)) {
      self_w <- Matrix::sparseMatrix(
        i = c(from, to), j = c(to, from),
        x = c(contribution, contribution),
        dims = c(length(all_nodes), length(all_nodes)),
        dimnames = list(all_nodes, all_nodes)
      )
      full <- full + self_w
    }
    w <- if (sparse_result) full else as.matrix(full)
    nodes <- all_nodes
  }

  diag(w) <- 0
  dimnames(w) <- list(nodes, nodes)
  if (identical(what, "matrix")) return(w)
  .thg_tidy_pairs(w, nodes)
}

.thg_resolve_edge_source <- function(hg, edge_source) {
  edges <- colnames(hg$incidence) %||% paste0("h", seq_len(hg$n_hyperedges))
  source <- edge_source
  if (is.null(source) && !is.null(hg$edge_data) &&
      all(c("edge", "source") %in% names(hg$edge_data))) {
    source <- stats::setNames(as.character(hg$edge_data$source),
                              as.character(hg$edge_data$edge))
  }
  if (is.data.frame(source)) {
    if (!all(c("edge", "source") %in% names(source))) {
      .thg_bad_input("an `edge_source` data.frame needs `edge` and `source` columns")
    }
    source <- stats::setNames(as.character(source$source),
                              as.character(source$edge))
  }
  if (is.null(source)) {
    .thg_bad_input(paste0("self-association needs one source per hyperedge; ",
                          "supply `edge_source` or construct a temporal ",
                          "hypergraph with `source =`"))
  }
  if (!is.atomic(source)) .thg_bad_input("`edge_source` must be a vector or data.frame")
  if (!is.null(names(source))) {
    missing_edges <- setdiff(edges, names(source))
    if (length(missing_edges)) {
      .thg_bad_input("the named `edge_source` vector does not cover every hyperedge")
    }
    source <- source[edges]
  } else if (length(source) != length(edges)) {
    .thg_bad_input("an unnamed `edge_source` must have one value per hyperedge")
  }
  source <- as.character(source)
  if (anyNA(source) || any(!nzchar(source))) {
    .thg_bad_input("every hyperedge must have a non-missing source")
  }
  unname(source)
}

#' The s-line graph of a hypergraph
#'
#' Turns the hyperedges into vertices: two hyperedges are adjacent when they
#' share at least `s` vertices, and the edge weight is how many they share.
#' `s = 1` is the ordinary line graph (any overlap connects); raising `s`
#' keeps only substantial overlaps, which is how walk-based measures on
#' hyperedges are parametrised (Aksoy et al. 2020).
#'
#' This is the projection of the dual: `hg_line_graph(hg, s = 1)` returns the
#' same graph as `hg_project(dual_hypergraph(hg), weighted = FALSE)`, which
#' the tests assert.
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any honets
#'   `net_hypergraph`.
#' @param s Minimum number of shared vertices for two hyperedges to be
#'   adjacent. A single integer of at least 1; `1` (default) is the ordinary
#'   line graph.
#' @param what `"edges"` (default) for the tidy edge list, or `"matrix"` for
#'   the symmetric overlap matrix to hand to a graph engine.
#' @return With `what = "edges"`, a base data.frame with one row per
#'   unordered hyperedge pair sharing at least `s` vertices, sorted by `from`
#'   then `to`, with columns `from`, `to` (hyperedge names) and `weight` (the
#'   number of shared vertices). A zero-row data.frame with those columns when
#'   no pair reaches `s`. With `what = "matrix"`, the symmetric
#'   `n_hyperedges` x `n_hyperedges` overlap matrix, zero on the diagonal and
#'   wherever the overlap is below `s`, sparse if `hg`'s incidence is sparse.
#' @references
#' Aksoy, S. G., Joslyn, C., Ortiz Marrero, C., Praggastis, B., & Purvine, E.
#' (2020). Hypernetwork science via high-order hypergraph walks.
#' *EPJ Data Science*, 9(1), 16. \doi{10.1140/epjds/s13688-020-00231-0}
#' @seealso [hg_project()] for the projection onto vertices,
#'   [dual_hypergraph()] for the role swap itself.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
#'                         c = "stars and salt"))
#' hg_line_graph(hg)
#' hg_line_graph(hg, s = 2)
#' @export
hg_line_graph <- function(hg, s = 1, what = c("edges", "matrix")) {
  .thg_check_hg(hg)
  what <- match.arg(what)
  stopifnot(
    "`s` must be a single whole number of at least 1" =
      length(s) == 1L && is.numeric(s) && is.finite(s) && s >= 1 &&
      isTRUE(all.equal(s, round(s)))
  )
  s <- as.integer(round(s))
  b <- .thg_binary(hg$incidence)
  overlap <- crossprod(b)
  diag(overlap) <- 0
  # Multiply by the threshold mask rather than index-assigning, which would
  # densify a sparse overlap matrix.
  overlap <- overlap * (overlap >= s)
  if (methods::is(overlap, "sparseMatrix")) overlap <- Matrix::drop0(overlap)
  edges <- colnames(hg$incidence)
  dimnames(overlap) <- list(edges, edges)
  if (identical(what, "matrix")) return(overlap)
  .thg_tidy_pairs(overlap, edges)
}

#' @rdname hg_project
#' @export
hypergraph_project <- hg_project

#' @rdname hg_line_graph
#' @export
hypergraph_line_graph <- hg_line_graph
