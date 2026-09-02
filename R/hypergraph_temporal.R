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

#' Construct a temporal hypergraph from membership spells
#'
#' Preserves group interactions as temporal hyperedges. Each input row gives
#' one member of one hyperedge; `member` may instead name several columns for
#' wide, fixed-cardinality data such as three-person tribunals. Two evolution
#' contracts cover the paper's examples: `"growing"` edges enter at `start`
#' and remain present, while `"interval"` edges are active on their closed
#' `[start, end]` interval.
#'
#' @param data A data frame containing memberships and times.
#' @param member One or more column names containing hyperedge members.
#' @param edge Column identifying the hyperedge or event.
#' @param start,time Column containing the edge start or point time. `time` is
#'   an alias for `start`; supply exactly one of them.
#' @param end End-time column for `evolution = "interval"`. Missing end values
#'   mean that the edge remains active through the end of observation.
#' @param source Optional column identifying the source/colour of each
#'   hyperedge, such as the decision containing a citation block. It is used
#'   by [hypergraph_project()] when `self_association = TRUE`.
#' @param weight Optional membership-weight column.
#' @param evolution `"growing"` for point-aggregation or `"interval"` for
#'   interval events.
#' @return A `net_temporal_hypergraph` retaining a tidy membership-spell table,
#'   edge metadata, the node universe and the event-time sequence.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' tribunals <- data.frame(
#'   case = c("A", "B"), president = c("p1", "p2"),
#'   arbitrator_1 = c("a1", "a1"), arbitrator_2 = c("a2", "a3"),
#'   constituted = c(1, 2), concluded = c(4, 5)
#' )
#' thg <- temporal_hypergraph(
#'   tribunals, member = c("president", "arbitrator_1", "arbitrator_2"),
#'   edge = "case", start = "constituted", end = "concluded",
#'   evolution = "interval"
#' )
#' hypergraph_snapshot(thg, at = 3)
#' @export
temporal_hypergraph <- function(data, member, edge, start = NULL, end = NULL,
                                time = NULL, source = NULL, weight = NULL,
                                evolution = c("growing", "interval")) {
  if (!is.data.frame(data) || nrow(data) == 0L) {
    .thg_bad_input("`data` must be a non-empty data.frame")
  }
  evolution <- match.arg(evolution)
  if (!is.character(member) || length(member) < 1L ||
      any(!member %in% names(data))) {
    .thg_bad_input("`member` must contain one or more column names in `data`")
  }
  edge <- .thg_selector(data, edge, "edge", required = TRUE)
  if (!is.null(start) && !is.null(time)) {
    .thg_bad_input("supply only one of `start` and `time`")
  }
  start <- .thg_selector(data, start %||% time, "start", required = TRUE)
  end <- .thg_selector(data, end, "end")
  source <- .thg_selector(data, source, "source")
  weight <- .thg_selector(data, weight, "weight")
  if (identical(evolution, "interval") && is.null(end)) {
    .thg_bad_input("`end` is required for `evolution = \"interval\"`")
  }

  # Pivot wide member columns with base R. Repeating the edge metadata here
  # is intentional: the normalized representation is one membership per row.
  normalized <- lapply(member, function(member_col) {
    data.frame(
      member = as.character(data[[member_col]]),
      edge = as.character(data[[edge]]),
      start = data[[start]],
      end = if (is.null(end)) rep(NA, nrow(data)) else data[[end]],
      source = if (is.null(source)) rep(NA_character_, nrow(data)) else
        as.character(data[[source]]),
      weight = if (is.null(weight)) rep(1, nrow(data)) else
        as.numeric(data[[weight]]),
      stringsAsFactors = FALSE
    )
  })
  memberships <- do.call(rbind, normalized)
  memberships <- memberships[
    !is.na(memberships$member) & nzchar(memberships$member) &
      !is.na(memberships$edge) & nzchar(memberships$edge) &
      !is.na(memberships$start) & !is.na(memberships$weight), , drop = FALSE
  ]
  rownames(memberships) <- NULL
  if (nrow(memberships) == 0L) {
    .thg_bad_input("no complete temporal memberships remain after dropping missing rows")
  }
  if (any(!is.finite(memberships$weight)) || any(memberships$weight < 0)) {
    .thg_bad_input("membership weights must be finite and non-negative")
  }

  # Edge-level time and source metadata must not depend on which membership
  # row happened to carry it.
  consistency <- function(column, allow_na = FALSE) {
    by_edge <- split(memberships[[column]], memberships$edge)
    bad <- vapply(by_edge, function(x) {
      x <- if (allow_na) x[!is.na(x)] else x
      length(unique(x)) > 1L
    }, logical(1L))
    if (any(bad)) {
      .thg_bad_input(sprintf("`%s` must be constant within each hyperedge", column))
    }
  }
  consistency("start")
  if (!is.null(end)) consistency("end", allow_na = TRUE)
  if (!is.null(source)) consistency("source", allow_na = TRUE)

  edge_names <- sort(unique(memberships$edge))
  first <- match(edge_names, memberships$edge)
  edge_data <- memberships[first, c("edge", "start", "end", "source"),
                           drop = FALSE]
  rownames(edge_data) <- NULL
  if (identical(evolution, "interval")) {
    bad_interval <- !is.na(edge_data$end) & edge_data$end < edge_data$start
    if (any(bad_interval)) .thg_bad_input("every interval must satisfy `end >= start`")
  }

  times <- sort(unique(c(edge_data$start, edge_data$end[!is.na(edge_data$end)])))
  structure(
    list(
      memberships = memberships,
      edge_data = edge_data,
      nodes = sort(unique(memberships$member)),
      edges = edge_names,
      times = times,
      evolution = evolution,
      params = list(member = member, edge = edge, start = start, end = end,
                    source = source, weight = weight)
    ),
    class = "net_temporal_hypergraph"
  )
}

.thg_empty_hypergraph <- function() {
  structure(list(
    hyperedges = list(), incidence = matrix(0, 0, 0), nodes = character(),
    n_nodes = 0L, n_hyperedges = 0L, size_distribution = integer(),
    params = list(source = "group_hypergraph", member = "member",
                  group = "edge", weight = NULL)
  ), class = "net_hypergraph")
}

.thg_collapse_duplicate_edges <- function(hg) {
  if (hg$n_hyperedges < 2L) {
    hg$edge_multiplicity <- rep.int(1L, hg$n_hyperedges)
    return(hg)
  }
  b <- hg$incidence != 0
  signatures <- vapply(seq_len(ncol(b)), function(j) {
    paste(rownames(b)[which(b[, j])], collapse = "\r")
  }, character(1L))
  unique_sig <- unique(signatures)
  first <- match(unique_sig, signatures)
  multiplicity <- tabulate(match(signatures, unique_sig), length(unique_sig))
  incidence <- b[, first, drop = FALSE] * 1
  colnames(incidence) <- colnames(hg$incidence)[first]
  hg$incidence <- incidence
  hg$hyperedges <- lapply(seq_len(ncol(incidence)), function(j) which(incidence[, j] != 0))
  hg$n_hyperedges <- ncol(incidence)
  sizes <- lengths(hg$hyperedges)
  size_tab <- table(sizes)
  hg$size_distribution <- stats::setNames(as.integer(size_tab),
                                           paste0("size_", names(size_tab)))
  hg$edge_multiplicity <- as.integer(multiplicity)
  if (!is.null(hg$edge_data)) hg$edge_data <- hg$edge_data[first, , drop = FALSE]
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
  if (nrow(d) == 0L) return(.thg_empty_hypergraph())
  hg <- group_hypergraph(d, member = "member", group = "edge", weight = "weight")
  hg$edge_data <- ed[match(colnames(hg$incidence), ed$edge), , drop = FALSE]
  hg$params$temporal_mode <- mode
  hg$params$at <- at
  hg$params$evolution <- x$evolution
  if (!multiedges) hg <- .thg_collapse_duplicate_edges(hg)
  hg
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
  duration <- object$edge_data$end - object$edge_data$start
  data.frame(
    n_nodes = length(object$nodes),
    n_hyperedges = length(object$edges),
    n_event_times = length(object$times),
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
