# Sub-hypergraphs: keep chosen hyperedges, the hyperedges inside a node set,
# or the hyperedges of one source (the citation blocks of one decision).

# Rebuild the derived fields of a net_hypergraph from a new incidence matrix
# (dense or sparse), keeping edge-level metadata aligned with the columns.
.thg_rebuild <- function(hg, incidence, keep_edges) {
  edge_names <- colnames(incidence)
  hg$incidence <- incidence
  hg$nodes <- rownames(incidence)
  hg$n_nodes <- nrow(incidence)
  hg$n_hyperedges <- ncol(incidence)
  hg$hyperedges <- .thg_edge_members(incidence)
  sizes <- lengths(hg$hyperedges)
  hg$size_distribution <- if (length(sizes)) {
    size_tab <- table(sizes)
    stats::setNames(as.integer(size_tab), paste0("size_", names(size_tab)))
  } else {
    integer(0L)
  }
  if (!is.null(hg$edge_multiplicity)) {
    hg$edge_multiplicity <- hg$edge_multiplicity[keep_edges]
  }
  if (!is.null(hg$edge_data)) {
    hg$edge_data <- hg$edge_data[
      match(edge_names, as.character(hg$edge_data$edge)), , drop = FALSE
    ]
    rownames(hg$edge_data) <- NULL
  }
  hg
}

# Member index vectors of every column, without touching one column at a
# time (which is slow on a sparse matrix with tens of thousands of columns).
.thg_edge_members <- function(incidence) {
  m <- ncol(incidence)
  if (m == 0L) return(list())
  if (methods::is(incidence, "sparseMatrix")) {
    triplet <- methods::as(methods::as(incidence, "generalMatrix"),
                           "TsparseMatrix")
    nz <- triplet@x != 0
    i <- triplet@i[nz] + 1L
    j <- triplet@j[nz] + 1L
  } else {
    nz <- which(incidence != 0, arr.ind = TRUE)
    i <- nz[, 1L]
    j <- nz[, 2L]
  }
  members <- split(i, factor(j, levels = seq_len(m)))
  lapply(unname(members), function(v) sort(unique(v)))
}

#' Sub-hypergraph by hyperedges, nodes or hyperedge attributes
#'
#' Keeps part of a hypergraph and returns it as a `net_hypergraph` with the
#' same incidence weights, edge metadata and multiplicities. Three selectors
#' combine by intersection. `edges` names the hyperedges to keep. `nodes`
#' keeps the hyperedges whose members all lie in the set, the induced
#' sub-hypergraph of HypergraphX. `where` keeps the hyperedges whose
#' attributes take given values, so `where = c(citing = "153-001")` on a
#' citation-block hypergraph is the hypergraph of one citing decision
#' (Coupette et al. 2024, Figure 3).
#'
#' @param hg A `net_hypergraph` (from [group_hypergraph()],
#'   [hypergraph_snapshot()], [text_hypergraph()], ...).
#' @param edges Hyperedge names to keep: a character vector, or a data.frame
#'   with an `edge` column such as the table [hg_edge_centrality()] or
#'   [hg_edges()] returns, so a ranking can be passed straight through.
#' @param nodes Character vector of node names. Only hyperedges contained in
#'   this set are kept, and the node set of the result is exactly `nodes`.
#' @param where Hyperedge attribute values to keep: a named vector or list,
#'   one element per attribute column of the edge metadata (the columns
#'   [temporal_hypergraph()] keeps because they are constant within a
#'   hyperedge), each holding the value or values to keep.
#' @param size Hyperedge sizes to keep, as a vector of member counts (distinct
#'   members): `size = 3` keeps the hyperedges with exactly three members, the
#'   3-uniform hypergraph that [hg_motifs()] needs.
#' @param drop_isolated Drop nodes that belong to no retained hyperedge?
#'   Default `TRUE`. Ignored when `nodes` is given.
#' @return A `net_hypergraph` whose incidence matrix is the selected
#'   sub-matrix of the input, sparse if the input is sparse. Edge metadata
#'   (`edge_data`) and duplicate multiplicities (`edge_multiplicity`) are
#'   subset alongside.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "b", "c", "d", "d", "e"),
#'   event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
#'   kind = c("x", "x", "x", "x", "x", "x", "y", "y")
#' )
#' hg <- group_hypergraph(dat, actor = "member", group = "event")
#' hg_subset(hg, edges = c("e1", "e2"))
#' hg_subset(hg, nodes = c("b", "c", "d", "e"))
#' hg_subset(hg, where = c(kind = "y"))
#' @export
hg_subset <- function(hg, edges = NULL, nodes = NULL, where = NULL,
                      size = NULL, drop_isolated = TRUE) {
  .thg_check_hg(hg)
  stopifnot(
    "`drop_isolated` must be TRUE or FALSE" =
      is.logical(drop_isolated) && length(drop_isolated) == 1L &&
      !is.na(drop_isolated)
  )
  if (is.null(edges) && is.null(nodes) && is.null(where) && is.null(size)) {
    .thg_bad_input("supply at least one of `edges`, `nodes`, `where` or `size`")
  }
  if (!is.null(size) && (!is.numeric(size) || anyNA(size) || any(size < 1))) {
    .thg_bad_input("`size` must be a vector of positive member counts")
  }
  edge_names <- colnames(hg$incidence)
  keep_edge <- rep(TRUE, hg$n_hyperedges)

  if (!is.null(edges)) {
    if (is.data.frame(edges)) {
      if (!"edge" %in% names(edges)) {
        .thg_bad_input("an `edges` data.frame needs an `edge` column")
      }
      edges <- edges$edge
    }
    edges <- unique(as.character(edges))
    unknown <- setdiff(edges, edge_names)
    if (length(unknown)) {
      .thg_bad_input(sprintf("unknown hyperedge(s): %s",
                             paste(utils::head(unknown, 5L), collapse = ", ")))
    }
    keep_edge <- keep_edge & edge_names %in% edges
  }
  if (!is.null(where)) {
    where <- as.list(where)
    if (is.null(names(where)) || any(!nzchar(names(where)))) {
      .thg_bad_input("`where` must be named by hyperedge attribute")
    }
    unknown <- setdiff(names(where), setdiff(names(hg$edge_data), "edge"))
    if (length(unknown)) {
      .thg_bad_input(sprintf("no hyperedge attribute named %s",
                             paste(unknown, collapse = ", ")))
    }
    row <- match(edge_names, as.character(hg$edge_data$edge))
    for (attribute in names(where)) {
      values <- hg$edge_data[[attribute]][row]
      keep_edge <- keep_edge & !is.na(values) & values %in% where[[attribute]]
    }
  }
  if (!is.null(size)) {
    members <- as.numeric(Matrix::colSums(.thg_binary(hg$incidence)))
    keep_edge <- keep_edge & members %in% size
  }
  if (!is.null(nodes)) {
    nodes <- as.character(nodes)
    unknown <- setdiff(nodes, hg$nodes)
    if (length(unknown)) {
      .thg_bad_input(sprintf("unknown node(s): %s",
                             paste(utils::head(unknown, 5L), collapse = ", ")))
    }
    inside <- hg$nodes %in% nodes
    b <- .thg_binary(hg$incidence)
    outside_count <- as.numeric(Matrix::colSums(b[!inside, , drop = FALSE]))
    keep_edge <- keep_edge & outside_count == 0
  }

  incidence <- hg$incidence[, keep_edge, drop = FALSE]
  keep_node <- if (!is.null(nodes)) {
    hg$nodes %in% nodes
  } else if (drop_isolated) {
    as.numeric(Matrix::rowSums(.thg_binary(incidence))) > 0
  } else {
    rep(TRUE, hg$n_nodes)
  }
  incidence <- incidence[keep_node, , drop = FALSE]
  out <- .thg_rebuild(hg, incidence, keep_edge)
  out$params$subset <- list(edges = edges, nodes = nodes, where = where,
                            size = size, drop_isolated = drop_isolated)
  out
}

#' @rdname hg_subset
#' @export
hypergraph_subset <- hg_subset
