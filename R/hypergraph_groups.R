# ---- Bipartite group hypergraph (EG-6) -----------------------------------
# Direct constructor: long-format event data with member + group columns
# becomes a net_hypergraph where each group is a hyperedge spanning all
# members that appeared in it.

#' Hypergraph from bipartite group / event data
#'
#' Constructs a [net_hypergraph][build_hypergraph] from long-format event
#' data in which each row records a `member` participating in a `group`
#' (a session, team, project, transaction, or any group context). Each
#' unique group becomes one hyperedge spanning the members that appeared in
#' it. Optional `weight` column produces a weighted incidence matrix.
#'
#' @param data Data frame in long format. Must contain `member` and `group`
#'   columns; optionally a `weight` column.
#' @param member Character. Name of the column whose values become the
#'   hypergraph's nodes (members, participants, actors).
#' @param group Character. Name of the column whose values become the
#'   hypergraph's hyperedges (groups, sessions, teams).
#' @param weight Character or `NULL`. If supplied, the column is summed per
#'   `(member, group)` pair to produce a weighted incidence matrix. Default
#'   `NULL` produces a 0/1 binary incidence matrix.
#' @param nodes Optional vector giving the complete node universe. This keeps
#'   nodes with no observed group memberships as zero-incidence rows, which is
#'   needed for representations such as citation hypergraphs where every
#'   decision is a node but some decisions are never cited.
#' @param sparse Logical. Store incidence as a sparse `Matrix`? Use this for
#'   large, sparse event data such as the full GFCC citation-block corpus.
#'
#' @return A `net_hypergraph` object with the same structure produced by
#'   [build_hypergraph()] (`hyperedges`, `incidence`, `nodes`, `n_nodes`,
#'   `n_hyperedges`, `size_distribution`, `params`). The `params` list
#'   records `source = "group_hypergraph"` and the original column names.
#'
#' @details
#' The bipartite representation preserves the full group structure without
#' projecting to a pairwise network. A group of three members A, B, C
#' produces a single 3-hyperedge containing all three, not three pairwise
#' edges AB, AC, BC. This avoids information loss when group interactions are
#' the primary unit of analysis (Perc et al. 2013).
#'
#' Unlike [build_hypergraph()] (which derives hyperedges from a network's
#' clique structure), `group_hypergraph()` takes group memberships
#' directly. The two functions are complementary:
#' \itemize{
#'   \item `group_hypergraph()` - when group membership is observed
#'     (sessions, transactions, co-authorships).
#'   \item `build_hypergraph()` - when only pairwise interactions are
#'     observed and triadic structure must be inferred from triangles.
#' }
#'
#' Rows with `NA` in either the `member` or `group` column (or, when
#' supplied, the `weight` column) are dropped silently.
#'
#' @seealso [build_hypergraph()] for the clique-based constructor.
#'
#' @examples
#' df <- data.frame(
#'   person = c("Alice", "Bob", "Carol", "Alice", "Bob",
#'              "Dave", "Carol", "Dave", "Eve"),
#'   session = c("S1", "S1", "S1", "S2", "S2",
#'               "S3", "S3", "S3", "S3")
#' )
#' hg <- group_hypergraph(df, member = "person", group = "session")
#' print(hg)
#' summary(hg)
#'
#' @references
#' Perc, M., Gomez-Gardenes, J., Szolnoki, A., Floria, L. M., & Moreno, Y.
#' (2013). Evolutionary dynamics of group interactions on structured
#' populations: a review. \emph{Journal of the Royal Society Interface}
#' 10(80), 20120997. \doi{10.1098/rsif.2012.0997}
#'
#' @note Dense and sparse paths are tested for exact equality. The sparse
#'   path is additionally exercised by the full GFCC reproduction from the
#'   Legal Hypergraphs Zenodo archive (3,618 nodes, 46,165 hyperedges and
#'   77,187 nonzero incidences).
#'
#' @export
group_hypergraph <- function(data, member, group, weight = NULL, nodes = NULL,
                             sparse = FALSE) {
  stopifnot(
    is.data.frame(data),
    is.character(member), length(member) == 1L,
    is.character(group),  length(group)  == 1L,
    member %in% names(data),
    group  %in% names(data),
    is.null(weight) ||
      (is.character(weight) && length(weight) == 1L && weight %in% names(data)),
    is.null(nodes) || is.atomic(nodes),
    is.logical(sparse), length(sparse) == 1L, !is.na(sparse)
  )

  cols <- c(member, group, weight)
  d <- data[, cols, drop = FALSE]
  d <- d[stats::complete.cases(d), , drop = FALSE]
  if (nrow(d) == 0L) {
    stop("No complete observations after dropping NAs.", call. = FALSE)
  }

  d[[member]] <- as.character(d[[member]])
  d[[group]]  <- as.character(d[[group]])

  if (!is.null(nodes)) {
    nodes <- as.character(nodes)
    if (anyNA(nodes) || any(!nzchar(nodes)) || anyDuplicated(nodes)) {
      stop("`nodes` must contain unique, non-missing node names.", call. = FALSE)
    }
    missing_members <- setdiff(unique(d[[member]]), nodes)
    if (length(missing_members)) {
      stop("Every observed member must occur in `nodes`.", call. = FALSE)
    }
  }
  member_levels <- sort(if (is.null(nodes)) unique(d[[member]]) else nodes)
  group_levels  <- sort(unique(d[[group]]))
  n_members <- length(member_levels)
  n_groups  <- length(group_levels)

  # Map values to row/col indices, then accumulate into the flat cell index.
  # A (member, group) pair may repeat across rows, so duplicated cells must be
  # summed BEFORE assignment -- assigning by an index vector keeps the last
  # write, not the total.
  mi <- match(d[[member]], member_levels)
  gj <- match(d[[group]],  group_levels)
  if (sparse) {
    incidence <- Matrix::sparseMatrix(
      i = mi, j = gj,
      x = if (is.null(weight)) rep.int(1, length(mi)) else as.numeric(d[[weight]]),
      dims = c(n_members, n_groups),
      dimnames = list(member_levels, group_levels)
    )
    if (is.null(weight) && length(incidence@x)) incidence@x[] <- 1
    incidence <- Matrix::drop0(incidence)
  } else {
    cell <- (gj - 1L) * n_members + mi
    if (is.null(weight)) {
      counts <- tabulate(cell, nbins = n_members * n_groups)
      incidence <- matrix(as.integer(counts > 0L), n_members, n_groups,
                          dimnames = list(member_levels, group_levels))
    } else {
      incidence <- matrix(0, n_members, n_groups,
                          dimnames = list(member_levels, group_levels))
      acc <- rowsum(as.numeric(d[[weight]]), cell, reorder = FALSE)
      incidence[as.integer(rownames(acc))] <- acc[, 1L]
    }
  }

  # Drop hyperedges that ended up empty (e.g. all-zero weight)
  he_sizes_pre <- if (sparse) Matrix::colSums(incidence > 0) else
    colSums(incidence > 0)
  keep <- he_sizes_pre > 0
  incidence <- incidence[, keep, drop = FALSE]
  group_levels <- group_levels[keep]
  n_groups <- length(group_levels)

  hyperedges <- lapply(seq_len(n_groups), function(j) {
    sort(which(incidence[, j] > 0))
  })

  he_sizes <- vapply(hyperedges, length, integer(1L))
  size_dist <- if (length(he_sizes)) {
    tab <- table(he_sizes)
    out <- as.integer(tab)
    names(out) <- paste0("size_", names(tab))
    out
  } else {
    integer(0L)
  }

  structure(
    list(
      hyperedges        = hyperedges,
      incidence         = incidence,
      nodes             = member_levels,
      n_nodes           = n_members,
      n_hyperedges      = n_groups,
      size_distribution = size_dist,
      params = list(
        source         = "group_hypergraph",
        member         = member,
        group          = group,
        weight         = weight,
        nodes          = nodes,
        sparse         = sparse,
        n_observations = nrow(d)
      )
    ),
    class = "net_hypergraph"
  )
}
