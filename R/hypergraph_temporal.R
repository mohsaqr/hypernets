# Temporal hypergraphs as edge-membership spells. This is deliberately an
# honets object, not a Dynet adapter: hyperedges remain first-class throughout.

.thg_bad_input <- function(message) {
  stop(errorCondition(message, class = "honets_bad_input", call = NULL))
}

.thg_selector <- function(data, selector, arg, required = FALSE) {
  if (is.null(selector)) {
    if (required) .thg_bad_input(sprintf("`%s` must name a column", arg))
    return(NULL)
  }
  if (!is.character(selector) || length(selector) != 1L ||
      !selector %in% names(data)) {
    .thg_bad_input(sprintf("`%s` must name one column in `data`", arg))
  }
  selector
}

#' Temporal hypergraph from an edge list or co-occurrence data
#'
#' Builds a temporal hypergraph the way a network is defined from data: one
#' row per relation, with the columns that name its ends and its time. Two
#' input shapes are accepted. An **edge list** names `from` and `to`, and
#' every row is a hyperedge of size two, so an ordinary temporal network is
#' the same object. **Co-occurrence data** name an `actor` and the column it
#' co-occurs `by`: every actor sharing one value of `cooccur_by` (a case, a
#' citation block, a session) belongs to one hyperedge. The clock is either
#' a single `time`, after which a hyperedge stays present (a growing
#' hypergraph, the point-aggregation model of Coupette et al. 2024), or a
#' `start` and an `end`, between which it is active (the interval-event
#' model). A missing `end` means the hyperedge stays active through the end
#' of observation.
#'
#' Every other column that is constant within a hyperedge is kept as a
#' hyperedge attribute in the edge metadata, where `plot()` can colour by it
#' and where a `source` column names the decision a citation block belongs
#' to, as [hg_project()] needs for self-association. Columns that vary
#' within a hyperedge, such as the seat an arbitrator held, are not
#' attributes of the hyperedge and are left out.
#'
#' @param data A data frame with one row per relation.
#' @param from,to Column names of a pairwise edge list.
#' @param actor,cooccur_by Column names of co-occurrence data: the node, and
#'   the grouping whose shared values bind nodes into one hyperedge.
#' @param time Column with the time a hyperedge appears (growing evolution).
#' @param start,end Columns with the interval on which a hyperedge is active
#'   (interval evolution). `start` alone is the same as `time`.
#' @param weight Optional membership-weight column.
#' @param nodes Optional node universe: a character vector of node names, or
#'   a data frame whose first column holds the names and whose `start`,
#'   `time` or `date` column, if present, gives the time each node enters. Nodes that
#'   never appear in a hyperedge are then kept as zero-degree nodes of every
#'   snapshot, and [hg_growth()] counts a node from its own start. Every
#'   observed node must be in the universe.
#' @param sparse Store every snapshot's incidence as a sparse `Matrix`?
#'   Default `FALSE`.
#' @return A `net_temporal_hypergraph` holding the membership table (`node`,
#'   `edge`, `start`, `end`, `weight`), the edge metadata (`edge`, `start`,
#'   `end` and the hyperedge attributes), the node universe with entry
#'   times, the sorted event times, and `evolution` (`"growing"` or
#'   `"interval"`).
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' seats <- data.frame(
#'   case = c("A", "A", "A", "B", "B", "B"),
#'   arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
#'   constituted = c(1, 1, 1, 2, 2, 2), concluded = c(4, 4, 4, 5, 5, 5),
#'   sector = c("oil", "oil", "oil", "gas", "gas", "gas")
#' )
#' thg <- temporal_hypergraph(seats, actor = "arbitrator", cooccur_by = "case",
#'                            start = "constituted", end = "concluded")
#' hypergraph_snapshot(thg, at = 3)
#'
#' contacts <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
#'                        time = 1:3)
#' temporal_hypergraph(contacts, from = "from", to = "to", time = "time")
#' @export
temporal_hypergraph <- function(data, from = NULL, to = NULL, actor = NULL,
                                cooccur_by = NULL, time = NULL, start = NULL,
                                end = NULL, weight = NULL, nodes = NULL,
                                sparse = FALSE) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    .thg_bad_input("`data` must be a non-empty data.frame")
  }
  if (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse)) {
    .thg_bad_input("`sparse` must be TRUE or FALSE")
  }
  edge_list <- !is.null(from) || !is.null(to)
  if (edge_list) {
    if (!is.null(actor) || !is.null(cooccur_by)) {
      .thg_bad_input("name either `from` and `to` or `actor` and `cooccur_by`, not both")
    }
    from <- .thg_selector(data, from, "from", required = TRUE)
    to <- .thg_selector(data, to, "to", required = TRUE)
  } else {
    actor <- .thg_selector(data, actor, "actor", required = TRUE)
    cooccur_by <- .thg_selector(data, cooccur_by, "cooccur_by", required = TRUE)
  }
  if (!is.null(start) && !is.null(time)) {
    .thg_bad_input("supply only one of `start` and `time`")
  }
  clock <- .thg_selector(data, start %||% time, "time", required = TRUE)
  end <- .thg_selector(data, end, "end")
  weight <- .thg_selector(data, weight, "weight")
  evolution <- if (is.null(end)) "growing" else "interval"

  # One row per membership. An edge list is unpivoted to its two ends; the
  # remaining columns ride along as candidate hyperedge attributes.
  used <- c(from, to, actor, cooccur_by, clock, end, weight)
  candidates <- setdiff(names(data), used)
  if (edge_list) {
    edge_id <- paste0("e", seq_len(nrow(data)))
    memberships <- data.frame(
      member = c(as.character(data[[from]]), as.character(data[[to]])),
      edge = c(edge_id, edge_id),
      stringsAsFactors = FALSE
    )
    rows <- c(seq_len(nrow(data)), seq_len(nrow(data)))
  } else {
    memberships <- data.frame(
      member = as.character(data[[actor]]),
      edge = as.character(data[[cooccur_by]]),
      stringsAsFactors = FALSE
    )
    rows <- seq_len(nrow(data))
  }
  memberships$start <- data[[clock]][rows]
  memberships$end <- if (is.null(end)) rep(NA, length(rows)) else data[[end]][rows]
  memberships$weight <- if (is.null(weight)) rep(1, length(rows)) else
    as.numeric(data[[weight]][rows])
  for (column in candidates) memberships[[column]] <- data[[column]][rows]

  keep <- !is.na(memberships$member) & nzchar(memberships$member) &
    !is.na(memberships$edge) & nzchar(memberships$edge) &
    !is.na(memberships$start) & !is.na(memberships$weight)
  memberships <- memberships[keep, , drop = FALSE]
  rownames(memberships) <- NULL
  if (nrow(memberships) == 0L) {
    .thg_bad_input("no complete memberships remain after dropping missing rows")
  }
  if (any(!is.finite(memberships$weight)) || any(memberships$weight < 0)) {
    .thg_bad_input("membership weights must be finite and non-negative")
  }

  # Times must not depend on which membership row carried them; attribute
  # columns that vary within a hyperedge are not hyperedge attributes.
  constant_within_edge <- function(column, allow_na = FALSE) {
    by_edge <- split(memberships[[column]], memberships$edge)
    !any(vapply(by_edge, function(x) {
      x <- if (allow_na) x[!is.na(x)] else x
      length(unique(x)) > 1L
    }, logical(1L)))
  }
  if (!constant_within_edge("start")) .thg_bad_input("`start`/`time` must be constant within each hyperedge")
  if (!is.null(end) && !constant_within_edge("end", allow_na = TRUE)) {
    .thg_bad_input("`end` must be constant within each hyperedge")
  }
  attributes <- Filter(function(column) constant_within_edge(column, allow_na = TRUE),
                       candidates)

  edge_names <- sort(unique(memberships$edge))
  first <- match(edge_names, memberships$edge)
  edge_data <- memberships[first, c("edge", "start", "end", attributes), drop = FALSE]
  rownames(edge_data) <- NULL
  if (identical(evolution, "interval")) {
    bad_interval <- !is.na(edge_data$end) & edge_data$end < edge_data$start
    if (any(bad_interval)) .thg_bad_input("every interval must satisfy `end >= start`")
  }
  memberships <- memberships[, c("member", "edge", "start", "end", "weight"), drop = FALSE]

  observed_nodes <- sort(unique(memberships$member))
  node_data <- .thg_node_universe(nodes, observed_nodes)

  times <- sort(unique(c(edge_data$start, edge_data$end[!is.na(edge_data$end)])))
  structure(
    list(
      memberships = memberships,
      edge_data = edge_data,
      node_data = node_data,
      nodes = node_data$node,
      edges = edge_names,
      times = times,
      evolution = evolution,
      params = list(from = from, to = to, actor = actor, cooccur_by = cooccur_by,
                    time = clock, end = end, weight = weight,
                    attributes = attributes, sparse = sparse,
                    nodes_given = !is.null(nodes))
    ),
    class = "net_temporal_hypergraph"
  )
}

# The node universe as a sorted `node` / `start` table. `start` is NA when the
# caller gave no entry times.
.thg_node_universe <- function(nodes, observed_nodes) {
  if (is.null(nodes)) {
    return(data.frame(node = observed_nodes, start = rep(NA, length(observed_nodes)),
                      stringsAsFactors = FALSE))
  }
  if (is.data.frame(nodes)) {
    if (ncol(nodes) == 0L) .thg_bad_input("a `nodes` data.frame needs a column of node names")
    node <- as.character(nodes[[1L]])
    clock <- intersect(c("start", "time", "date"), names(nodes))
    start <- if (length(clock)) nodes[[clock[[1L]]]] else rep(NA, length(node))
  } else if (is.atomic(nodes)) {
    node <- as.character(nodes)
    start <- rep(NA, length(node))
  } else {
    .thg_bad_input("`nodes` must be a character vector or a data.frame")
  }
  if (anyNA(node) || any(!nzchar(node)) || anyDuplicated(node)) {
    .thg_bad_input("`nodes` must contain unique, non-missing node names")
  }
  missing_members <- setdiff(observed_nodes, node)
  if (length(missing_members)) {
    .thg_bad_input(sprintf("%d observed node(s) are not in `nodes`",
                           length(missing_members)))
  }
  ord <- order(node)
  data.frame(node = node[ord], start = start[ord], stringsAsFactors = FALSE)
}

.thg_empty_hypergraph <- function(nodes = character(), sparse = FALSE) {
  n <- length(nodes)
  incidence <- if (sparse) {
    Matrix::Matrix(0, n, 0L, sparse = TRUE, dimnames = list(nodes, NULL))
  } else {
    matrix(0, n, 0L, dimnames = list(nodes, NULL))
  }
  structure(list(
    hyperedges = list(), incidence = incidence, nodes = nodes,
    n_nodes = n, n_hyperedges = 0L, size_distribution = integer(),
    params = list(source = "group_hypergraph", member = "member",
                  group = "edge", weight = NULL)
  ), class = "net_hypergraph")
}

.thg_collapse_duplicate_edges <- function(hg) {
  if (hg$n_hyperedges < 2L) {
    hg$edge_multiplicity <- rep.int(1L, hg$n_hyperedges)
    return(hg)
  }
  members <- .thg_edge_members(hg$incidence)
  signatures <- vapply(members, paste, collapse = "\r", character(1L))
  unique_sig <- unique(signatures)
  first <- match(unique_sig, signatures)
  multiplicity <- tabulate(match(signatures, unique_sig), length(unique_sig))
  incidence <- .thg_binary(hg$incidence[, first, drop = FALSE])
  keep <- seq_len(hg$n_hyperedges) %in% first
  hg <- .thg_rebuild(hg, incidence, keep)
  hg$edge_multiplicity <- as.integer(multiplicity)
  hg
}

#' Extract one snapshot from a temporal hypergraph
#'
#' @param x A [temporal_hypergraph()].
#' @param at One time value. Not used for `mode = "all"`.
#' @param mode `"active"` applies the object's growing or interval contract;
#'   `"cumulative"` retains every edge begun by `at`; `"all"` ignores time
#'   and returns the static aggregate.
#' @param multiedges Keep distinct edge identities with identical member sets?
#'   `TRUE` matches a multi-hypergraph; `FALSE` collapses them to one edge and
#'   records their counts in `$edge_multiplicity`.
#' @return A static `net_hypergraph` usable by every honets hypergraph verb.
#' @export
hypergraph_snapshot <- function(x, at = NULL,
                                mode = c("active", "cumulative", "all"),
                                multiedges = TRUE) {
  if (!inherits(x, "net_temporal_hypergraph")) {
    .thg_bad_input("`x` must come from temporal_hypergraph()")
  }
  mode <- match.arg(mode)
  if (!is.logical(multiedges) || length(multiedges) != 1L || is.na(multiedges)) {
    .thg_bad_input("`multiedges` must be TRUE or FALSE")
  }
  if (!identical(mode, "all") && (is.null(at) || length(at) != 1L || is.na(at))) {
    .thg_bad_input("`at` must be one time value unless `mode = \"all\"`")
  }
  ed <- x$edge_data
  keep <- switch(mode,
    all = rep(TRUE, nrow(ed)),
    cumulative = ed$start <= at,
    active = if (identical(x$evolution, "growing")) {
      ed$start <= at
    } else {
      ed$start <= at & (is.na(ed$end) | ed$end >= at)
    }
  )
  active_edges <- ed$edge[keep]
  d <- x$memberships[x$memberships$edge %in% active_edges, , drop = FALSE]
  sparse <- isTRUE(x$params$sparse)
  universe <- .thg_snapshot_nodes(x, at, mode, d$member)
  if (nrow(d) == 0L) return(.thg_empty_hypergraph(universe, sparse))
  hg <- group_hypergraph(d, member = "member", group = "edge", weight = "weight",
                         nodes = universe, sparse = sparse)
  hg$edge_data <- ed[match(colnames(hg$incidence), ed$edge), , drop = FALSE]
  rownames(hg$edge_data) <- NULL
  hg$params$temporal_mode <- mode
  hg$params$at <- at
  hg$params$evolution <- x$evolution
  if (!multiedges) hg <- .thg_collapse_duplicate_edges(hg)
  hg
}

# Node universe of a snapshot. Without a universe, the nodes are the members
# of the active hyperedges (as before). With one, every node is kept; when
# nodes carry entry times and the snapshot is time-bound, only nodes entered
# by `at` plus any active member are kept.
.thg_snapshot_nodes <- function(x, at, mode, active_members) {
  if (!isTRUE(x$params$nodes_given)) return(NULL)
  node_data <- x$node_data
  if (identical(mode, "all") || all(is.na(node_data$start))) {
    return(node_data$node)
  }
  entered <- !is.na(node_data$start) & node_data$start <= at
  sort(union(node_data$node[entered], unique(active_members)))
}

#' Extract a sequence of temporal-hypergraph snapshots
#'
#' @inheritParams hypergraph_snapshot
#' @param at Time values. `NULL` uses every event time stored by the temporal
#'   hypergraph.
#' @return A named list of `net_hypergraph` objects with class
#'   `net_hypergraph_snapshots`.
#' @export
hypergraph_snapshots <- function(x, at = NULL,
                                 mode = c("active", "cumulative", "all"),
                                 multiedges = TRUE) {
  if (!inherits(x, "net_temporal_hypergraph")) {
    .thg_bad_input("`x` must come from temporal_hypergraph()")
  }
  mode <- match.arg(mode)
  if (identical(mode, "all")) {
    out <- list(all = hypergraph_snapshot(x, mode = "all",
                                          multiedges = multiedges))
  } else {
    at <- at %||% x$times
    if (length(at) == 0L || anyNA(at)) .thg_bad_input("`at` contains no valid times")
    out <- lapply(seq_along(at), function(i) {
      hypergraph_snapshot(x, at = at[i], mode = mode, multiedges = multiedges)
    })
    names(out) <- make.unique(as.character(at))
  }
  structure(out, class = c("net_hypergraph_snapshots", "list"))
}

#' @export
print.net_temporal_hypergraph <- function(x, ...) {
  cat(sprintf("Temporal hypergraph: %d nodes, %d hyperedges, %d event times\n",
              length(x$nodes), length(x$edges), length(x$times)))
  cat(sprintf("Evolution: %s\n", x$evolution))
  invisible(x)
}

#' @export
summary.net_temporal_hypergraph <- function(object, ...) {
  memberships_per_edge <- table(object$memberships$edge)
  duration <- if (all(is.na(object$edge_data$end))) NA_real_ else
    as.numeric(object$edge_data$end - object$edge_data$start)
  data.frame(
    n_nodes = length(object$nodes),
    n_hyperedges = length(object$edges),
    n_event_times = length(object$times),
    first_time = min(object$times),
    last_time = max(object$times),
    n_memberships = nrow(object$memberships),
    mean_edge_size = mean(as.numeric(memberships_per_edge)),
    median_edge_size = stats::median(as.numeric(memberships_per_edge)),
    mean_duration = if (all(is.na(duration))) NA_real_ else
      mean(as.numeric(duration), na.rm = TRUE),
    evolution = object$evolution,
    row.names = NULL
  )
}

#' @export
as.data.frame.net_temporal_hypergraph <- function(x, row.names = NULL,
                                                   optional = FALSE, ...) {
  x$memberships
}

#' Plot a temporal-hypergraph snapshot through cograph
#'
#' @param x A [temporal_hypergraph()].
#' @param at Snapshot time; defaults to the latest event time.
#' @param mode Snapshot mode passed to [hypergraph_snapshot()].
#' @param method Projection weighting passed to [hypergraph_project()].
#' @param ... Additional arguments passed to [cograph::splot()].
#' @return The cograph plot object, invisibly when rendered interactively.
#' @export
plot.net_temporal_hypergraph <- function(x, at = NULL,
                                         mode = c("active", "cumulative", "all"),
                                         method = c("association", "clique"), ...) {
  mode <- match.arg(mode)
  method <- match.arg(method)
  if (is.null(at) && !identical(mode, "all")) at <- utils::tail(x$times, 1L)
  hg <- hypergraph_snapshot(x, at = at, mode = mode)
  if (hg$n_nodes == 0L) .thg_bad_input("the selected snapshot has no active nodes")
  projection <- if (identical(method, "association")) {
    hg_project(hg, method = "association", what = "matrix")
  } else {
    hg_project(hg, method = "clique", what = "matrix")
  }
  cograph::splot(as.matrix(projection), directed = FALSE, ...)
}

# Compact aliases retain identical bodies and formals.
#' @rdname hypergraph_snapshot
#' @export
hg_snapshot <- hypergraph_snapshot

#' @rdname hypergraph_snapshots
#' @export
hg_snapshots <- hypergraph_snapshots
