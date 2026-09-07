# Temporal hypergraphs as hyperedge spells. This is deliberately an hypernets
# object, not a Dynet adapter: hyperedges remain first-class throughout. The
# vocabulary is Dynet's (`../temporal`): column aliases, time parsing,
# observation bounds and the step/window measurement grid follow its
# definitions so a co-presence log reads the same way in both packages.

.thg_bad_input <- function(message, class = character()) {
  stop(errorCondition(message, class = c(class, "hypernets_bad_input"), call = NULL))
}

.thg_deprecated <- function(old, new, fn) {
  warning(warningCondition(
    sprintf("`%s` is deprecated in %s(); use `%s` instead", old, fn, new),
    class = c("hypernets_deprecated", "deprecatedWarning"), call = NULL
  ))
}

# ---------------------------------------------------------------------------
# Column detection: Dynet's alias table (R/detect.R), checked case-
# insensitively after stripping non-alphanumerics; order is priority order.
# ---------------------------------------------------------------------------

.thg_aliases <- list(
  from     = c("from", "source", "sender", "tail", "ego", "actor1", "node1",
               "sourceid", "fromid", "i"),
  to       = c("to", "target", "receiver", "head", "alter", "actor2", "node2",
               "targetid", "toid", "j"),
  start    = c("start", "onset", "begin", "starttime", "startdate", "tstart",
               "fromtime"),
  end      = c("end", "terminus", "finish", "endtime", "enddate", "tend",
               "totime", "stop"),
  duration = c("duration", "dur", "elapsed", "length", "lasted"),
  time     = c("time", "timestamp", "datetime", "date", "when", "t", "occurred",
               "createdat", "posttime"),
  thread   = c("thread", "threadid", "discussion", "discussiontitle", "topic",
               "conversation", "parent", "postid", "root"),
  session  = c("session", "period", "wave", "phase", "cohort", "course"),
  group    = c("group", "event", "context", "room", "class", "meeting",
               "venue", "team", "channel"),
  actor    = c("actor", "person", "member", "student", "participant", "user",
               "id", "name"),
  weight   = c("weight", "weights", "strength")
)

.thg_norm_name <- function(x) gsub("[^a-z0-9]", "", tolower(x))

.thg_match_column <- function(data, role, exclude = character()) {
  aliases <- .thg_aliases[[role]]
  cand <- setdiff(names(data), exclude)
  hit <- match(aliases, .thg_norm_name(cand))
  hit <- hit[!is.na(hit)]
  if (length(hit) == 0L) NULL else cand[hit[1L]]
}

# One column: an explicit name must exist as given; otherwise the first alias
# match, or NULL.
.thg_resolve_column <- function(data, given, role, exclude = character(),
                                arg = role) {
  if (!is.null(given)) {
    if (!is.character(given) || length(given) != 1L || is.na(given)) {
      .thg_bad_input(sprintf("`%s` must be a single column name", arg))
    }
    if (!given %in% names(data)) {
      .thg_bad_input(
        sprintf("column `%s` (given as `%s`) is not in the data; available: %s",
                given, arg, paste(names(data), collapse = ", ")),
        class = "hypernets_missing_column"
      )
    }
    return(given)
  }
  .thg_match_column(data, role, exclude)
}

# ---------------------------------------------------------------------------
# Time parsing: Dynet's rule. Numeric times pass through with unit "step";
# Date, POSIXct and character date-times become elapsed time since an origin
# in a unit chosen for the span (or requested).
# ---------------------------------------------------------------------------

.thg_time_formats <- c(
  "%Y-%m-%d %H:%M:%S", "%Y-%m-%dT%H:%M:%S", "%Y-%m-%d %H:%M", "%Y-%m-%d",
  "%Y/%m/%d %H:%M:%S", "%Y/%m/%d %H:%M", "%Y/%m/%d",
  "%d/%m/%Y %H:%M:%S", "%d/%m/%Y %H:%M", "%d/%m/%Y",
  "%m/%d/%Y %H:%M:%S", "%m/%d/%Y %H:%M", "%m/%d/%Y",
  "%d-%m-%Y %H:%M:%S", "%d-%m-%Y %H:%M", "%d-%m-%Y",
  "%d %b %Y %H:%M", "%d %B %Y %H:%M", "%d %b %Y", "%d %B %Y"
)

.thg_unit_seconds <- function(time_unit) {
  switch(time_unit, seconds = 1, minutes = 60, hours = 3600, days = 86400,
         weeks = 604800)
}

.thg_auto_unit <- function(span_seconds) {
  if (!is.finite(span_seconds) || span_seconds <= 0) return("seconds")
  if (span_seconds < 2 * 60) "seconds"
  else if (span_seconds < 3 * 3600) "minutes"
  else if (span_seconds < 3 * 86400) "hours"
  else "days"
}

.thg_parse_datetime_strings <- function(x) {
  present <- !is.na(x)
  trial <- lapply(.thg_time_formats, function(fmt) {
    as.POSIXct(x, format = fmt, tz = "UTC")
  })
  n_ok <- vapply(trial, function(p) sum(!is.na(p[present])), integer(1L))
  best <- which.max(n_ok)
  if (n_ok[best] == 0L) {
    .thg_bad_input(
      sprintf("could not parse time strings such as '%s'; supply numeric, Date or POSIXct times",
              x[present][1L]),
      class = "hypernets_unparsed_time"
    )
  }
  trial[[best]]
}

# Parse one or more time vectors on a shared clock. `columns` is a named
# list; NA entries are allowed (an open end, a node without an entry time)
# and stay NA. Returns the parsed columns plus `unit` and `origin`.
.thg_parse_clock <- function(columns, time_unit = "auto") {
  if (!is.character(time_unit) || length(time_unit) != 1L ||
      !time_unit %in% c("auto", "seconds", "minutes", "hours", "days", "weeks")) {
    .thg_bad_input("`time_unit` must be \"auto\", \"seconds\", \"minutes\", \"hours\", \"days\" or \"weeks\"")
  }
  columns <- lapply(columns, function(x) if (is.factor(x)) as.character(x) else x)
  kinds <- vapply(columns, function(x) {
    if (is.numeric(x)) "numeric"
    else if (inherits(x, "POSIXct") || inherits(x, "Date")) "calendar"
    else if (is.character(x)) "character"
    else "unsupported"
  }, character(1L))
  if (any(kinds == "unsupported")) {
    bad <- names(columns)[kinds == "unsupported"][1L]
    .thg_bad_input(sprintf(
      "time column `%s` has unsupported type <%s>; supply numeric, Date, POSIXct or character times",
      bad, paste(class(columns[[bad]]), collapse = "/")))
  }
  if (all(kinds == "numeric")) {
    return(c(lapply(columns, as.numeric), list(unit = "step", origin = 0)))
  }
  if (any(kinds == "numeric")) {
    .thg_bad_input("time columns must all be dates or all be numbers; they cannot be mixed")
  }
  parsed <- lapply(columns, function(x) {
    if (is.character(x)) .thg_parse_datetime_strings(x) else as.POSIXct(x, tz = "UTC")
  })
  all_times <- do.call(c, lapply(parsed, function(p) p[!is.na(p)]))
  if (length(all_times) == 0L) .thg_bad_input("no time value could be read")
  span <- as.numeric(difftime(max(all_times), min(all_times), units = "secs"))
  unit <- if (identical(time_unit, "auto")) .thg_auto_unit(span) else time_unit
  origin <- min(all_times)
  values <- lapply(parsed, function(p) {
    as.numeric(difftime(p, origin, units = "secs")) / .thg_unit_seconds(unit)
  })
  c(values, list(unit = unit, origin = origin))
}

# A user-supplied time (`at`, `start`, `observation_end`, ...) on the stored
# clock: numbers pass through, dates are converted for a calendar hypergraph
# and refused for a numeric one. Vectorised.
.thg_as_time <- function(v, x, arg) {
  if (is.null(v)) return(NULL)
  if (is.factor(v)) v <- as.character(v)
  calendar <- inherits(x$origin, "POSIXt")
  if (inherits(v, "Date") || inherits(v, "POSIXt") || is.character(v)) {
    if (!calendar) {
      .thg_bad_input(sprintf(
        "`%s` was given as a date, but this hypergraph's times are plain numbers",
        arg))
    }
    parsed <- if (is.character(v)) .thg_parse_datetime_strings(v) else
      as.POSIXct(v, tz = "UTC")
    if (anyNA(parsed)) .thg_bad_input(sprintf("`%s` contains missing times", arg))
    return(as.numeric(difftime(parsed, x$origin, units = "secs")) /
             .thg_unit_seconds(x$time_unit))
  }
  if (!is.numeric(v) || anyNA(v) || any(!is.finite(v))) {
    .thg_bad_input(sprintf("`%s` must be finite numbers or dates", arg))
  }
  as.numeric(v)
}

# The calendar value of stored times, for axes and labels; numeric clocks
# return the numbers unchanged.
.thg_calendar <- function(t, origin, time_unit) {
  if (!inherits(origin, "POSIXt")) return(t)
  out <- origin + t * .thg_unit_seconds(time_unit)
  if (identical(time_unit, "days") || identical(time_unit, "weeks")) as.Date(out) else out
}

.thg_clock_label <- function(origin, time_unit) {
  if (!inherits(origin, "POSIXt")) return("time (numeric steps)")
  sprintf("%s since %s", time_unit, format(origin, if (identical(time_unit, "days") ||
                                                          identical(time_unit, "weeks"))
    "%Y-%m-%d" else "%Y-%m-%d %H:%M:%S"))
}

#' Temporal hypergraph from an edge list or co-occurrence data
#'
#' Builds a temporal hypergraph the way a temporal network is defined from
#' data in this ecosystem (Dynet's `dynet()`): one row per relation, with the
#' columns that name its ends and its clock. Two input shapes are accepted.
#' An **edge list** names `from` and `to`, and every row is a hyperedge of
#' size two, so an ordinary temporal network is the same object.
#' **Co-presence data** name an `actor` and a `group`: every actor sharing one
#' value of `group` (a case, a citation block, a seminar) belongs to one
#' hyperedge. Two clocks are understood, and the one you name selects the
#' format:
#'
#' \describe{
#'   \item{interval}{`start` and `end`: the hyperedge is active on the closed
#'     interval between them. A missing `end` means it stays active through
#'     the end of observation.}
#'   \item{contact}{`time`: the hyperedge is an instantaneous event, as in a
#'     contact log -- a citation, a message, a meeting on one day. It is
#'     present at that instant only. The cumulative view in which every
#'     hyperedge stays once it has appeared (the point-aggregation model of
#'     Coupette et al. 2024) is a snapshot `mode`, not a property of the
#'     data: ask for it with `mode = "cumulative"` in [hypergraph_snapshot()],
#'     [hg_growth()] and [hg_edges()].}
#' }
#'
#' A third shape is a **sequence table**: one row per session and one column
#' per position holding the state at that step (the wide format of tna and
#' TraMineR), a list of character vectors, or a `tna` / `netobject` model
#' built from one. It is recognised when no relational column is named or
#' detected and every column is categorical. Each session is one hyperedge
#' whose members are the states it contains, and a state's membership is a
#' contact at its position, `1` to the session's length, on a `"step"`
#' clock: the simple co-occurrence reading in which time is order.
#'
#' In any shape a membership may carry its own time, as when a log has one
#' row per attendance rather than one time per group. The hyperedge then
#' spans from its first to its last membership, each membership is present
#' on its own spell only, and a snapshot keeps the memberships present in
#' its window: `mode = "cumulative"` at step `t` is what each session had
#' shown by `t`, and `window = 3` at `t` is what it showed on steps `t` to
#' `t + 2`. When every membership carries its hyperedge's time, as a
#' tribunal or a citation block does, nothing changes.
#'
#' Column names are resolved case-insensitively from the same alias table
#' Dynet uses, so `Sender`/`Receiver`, `source`/`target`, `onset`/`terminus`
#' and `timestamp` are understood without being spelled out; a name you give
#' explicitly must exist as written. Times may be numeric, `Date`, `POSIXct`
#' or character date-time strings. Numeric times are kept as they are and
#' reported in `"step"` units; calendar input is converted to elapsed time
#' since the earliest time in the data, in a unit chosen for the span
#' (`"days"` beyond three days) or given as `time_unit`, and the unit and
#' origin are stored in the object and shown by `print()`. Every time you
#' pass later -- `at`, `start`, `end`, the observation bounds -- may be a
#' date for a calendar hypergraph or a number on that clock, and every
#' `time` column in a result is on that clock.
#'
#' `observation_start` and `observation_end` declare the study window when it
#' is known independently of the log, with Dynet's meaning: they bound the
#' snapshot times and the measurement grid, and an open-ended hyperedge is
#' active through `observation_end`, but the stored memberships are never
#' rewritten and `as.data.frame()` returns the original spells. Without them
#' the window is the span of the data.
#'
#' Every other column that is constant within a hyperedge is kept as a
#' hyperedge attribute in the edge metadata, where `plot()` can colour by it
#' and where [hg_project()] finds the source a citation block belongs to.
#' Columns that vary within a hyperedge, such as the seat an arbitrator held,
#' are not attributes of the hyperedge and are left out.
#'
#' @param data A data frame with one row per relation, or a sequence table
#'   (a wide data frame of states, a list of character vectors, or a `tna` /
#'   `netobject` model); see Details.
#' @param from,to Column names of a pairwise edge list. Detected from the
#'   alias table when neither is given and `actor`/`group` are not named.
#' @param actor,group Column names of co-presence data: the node, and the
#'   grouping whose shared values bind nodes into one hyperedge. Naming
#'   either selects the co-presence format; the other is then detected by
#'   alias if not given.
#' @param time Column with the instant of a contact hyperedge.
#' @param start,end Columns with the interval on which a hyperedge is active.
#'   `start` without an `end` column is read as a contact clock.
#' @param weight Optional membership-weight column.
#' @param nodes Optional node universe: a character vector of node names, or
#'   a data frame whose first column holds the names and whose `start`,
#'   `time` or `date` column, if present, gives the time each node enters (on
#'   the same clock as `data`). Nodes that never appear in a hyperedge are
#'   then kept as zero-degree nodes of every snapshot, and [hg_growth()]
#'   counts a node from its own start. Every observed node must be in the
#'   universe.
#' @param time_unit Unit for calendar times: `"auto"` (default),
#'   `"seconds"`, `"minutes"`, `"hours"`, `"days"` or `"weeks"`. Numeric
#'   times are left alone and reported as `"step"`.
#' @param observation_start,observation_end Optional bounds of the
#'   observation window, as numbers on the stored clock or as dates for a
#'   calendar hypergraph. Either may be omitted; the corresponding limit of
#'   the data is then used.
#' @param sparse Store every snapshot's incidence as a sparse `Matrix`?
#'   Default `FALSE`.
#' @param cooccur_by Deprecated name of `group`; using it warns with a
#'   `hypernets_deprecated` condition.
#' @return A `net_temporal_hypergraph` holding the membership table (`node`,
#'   `edge`, `start`, `end`, `weight`), the edge metadata (`edge`, `start`,
#'   `end` and the hyperedge attributes), the node universe with entry
#'   times, the sorted event times, `format` (`"interval"` or `"contact"`),
#'   `time_unit`, `origin` and the `observation` bounds. `as.data.frame(x,
#'   what = "memberships" | "edges" | "nodes")` returns the three tables;
#'   `summary()` the one-row description.
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
#' thg <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
#'                            start = "constituted", end = "concluded")
#' thg
#' hypergraph_snapshot(thg, at = 3)
#'
#' # a contact log on a calendar: instants at each date, cumulative on request
#' contacts <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
#'                        date = as.Date(c("2024-01-01", "2024-01-05", "2024-01-09")))
#' calls <- temporal_hypergraph(contacts, from = "from", to = "to", time = "date")
#' calls
#' hypergraph_snapshot(calls, at = as.Date("2024-01-05"))
#' hypergraph_snapshot(calls, at = as.Date("2024-01-05"), mode = "cumulative")
#' @export
temporal_hypergraph <- function(data, from = NULL, to = NULL, actor = NULL,
                                group = NULL, time = NULL, start = NULL,
                                end = NULL, weight = NULL, nodes = NULL,
                                time_unit = "auto",
                                observation_start = NULL, observation_end = NULL,
                                sparse = FALSE, cooccur_by = NULL) {
  if (.thg_is_sequence_input(data, from, to, actor, group, time, start, end)) {
    data <- .thg_sequence_memberships(data)
    actor <- "state"
    group <- "sequence"
    time <- "position"
  }
  if (!is.data.frame(data) || nrow(data) == 0L) {
    .thg_bad_input("`data` must be a non-empty data.frame")
  }
  if (!is.logical(sparse) || length(sparse) != 1L || is.na(sparse)) {
    .thg_bad_input("`sparse` must be TRUE or FALSE")
  }
  if (!is.null(cooccur_by)) {
    .thg_deprecated("cooccur_by", "group", "temporal_hypergraph")
    if (is.null(group)) group <- cooccur_by
  }
  roles <- .thg_resolve_roles(data, from, to, actor, group, time, start, end,
                              weight)
  from <- roles$from; to <- roles$to; actor <- roles$actor; group <- roles$group
  clock <- roles$clock; end <- roles$end; weight <- roles$weight
  edge_list <- !is.null(from)
  format <- roles$format

  # One row per membership. An edge list is unpivoted to its two ends; the
  # remaining columns ride along as candidate hyperedge attributes.
  used <- c(from, to, actor, group, clock, end, weight)
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
      edge = as.character(data[[group]]),
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

  # One clock for every time the object holds: hyperedge starts and ends and
  # the node entry times share the origin and the unit.
  observed_nodes <- sort(unique(memberships$member))
  node_data <- .thg_node_universe(nodes, observed_nodes)
  clock_columns <- list(start = memberships$start)
  if (!is.null(end)) clock_columns$end <- memberships$end
  if (!all(is.na(node_data$start))) clock_columns$node_start <- node_data$start
  parsed <- .thg_parse_clock(clock_columns, time_unit)
  memberships$start <- parsed$start
  memberships$end <- if (is.null(end)) rep(NA_real_, nrow(memberships)) else parsed$end
  node_data$start <- if (is.null(parsed$node_start)) rep(NA_real_, nrow(node_data)) else
    parsed$node_start

  # Times must not depend on which membership row carried them; attribute
  # columns that vary within a hyperedge are not hyperedge attributes.
  constant_within_edge <- function(column, allow_na = FALSE) {
    by_edge <- split(memberships[[column]], memberships$edge)
    !any(vapply(by_edge, function(x) {
      x <- if (allow_na) x[!is.na(x)] else x
      length(unique(x)) > 1L
    }, logical(1L)))
  }
  # A hyperedge's time is the hull of its memberships. When every membership
  # carries the hyperedge's time (a tribunal, a citation block) the hull is
  # that time and nothing below changes. When memberships carry their own
  # times (a session read as a sequence, a log with one row per attendance)
  # the hyperedge spans its first to its last membership and each membership
  # is present on its own spell only; a contact is present at its instant.
  membership_times <- !constant_within_edge("start") ||
    (!is.null(end) && !constant_within_edge("end", allow_na = TRUE))
  if (membership_times && identical(format, "contact")) {
    memberships$end <- memberships$start
    format <- "interval"
  }
  attributes <- Filter(function(column) constant_within_edge(column, allow_na = TRUE),
                       candidates)

  edge_names <- sort(unique(memberships$edge))
  first <- match(edge_names, memberships$edge)
  edge_data <- memberships[first, c("edge", "start", "end", attributes), drop = FALSE]
  rownames(edge_data) <- NULL
  if (membership_times) {
    by_edge <- factor(memberships$edge, levels = edge_names)
    edge_data$start <- as.numeric(tapply(memberships$start, by_edge, min))
    edge_data$end <- as.numeric(tapply(memberships$end, by_edge, function(e) {
      if (anyNA(e)) NA_real_ else max(e)
    }))
  }
  if (identical(format, "interval")) {
    bad_interval <- !is.na(edge_data$end) & edge_data$end < edge_data$start
    if (any(bad_interval)) .thg_bad_input("every interval must satisfy `end >= start`")
  }
  memberships <- memberships[, c("member", "edge", "start", "end", "weight"), drop = FALSE]

  times <- sort(unique(c(edge_data$start, edge_data$end[!is.na(edge_data$end)],
                         memberships$start, memberships$end[!is.na(memberships$end)])))
  span <- range(c(times, node_data$start[!is.na(node_data$start)]))

  x <- structure(
    list(
      memberships = memberships,
      edge_data = edge_data,
      node_data = node_data,
      nodes = node_data$node,
      edges = edge_names,
      times = times,
      format = format,
      time_unit = parsed$unit,
      origin = parsed$origin,
      observation = c(start = span[[1L]], end = span[[2L]]),
      params = list(from = from, to = to, actor = actor, group = group,
                    time = clock, end = end, weight = weight,
                    attributes = attributes, sparse = sparse,
                    membership_times = membership_times,
                    nodes_given = !is.null(nodes),
                    observation_explicit = !is.null(observation_start) ||
                      !is.null(observation_end))
    ),
    class = "net_temporal_hypergraph"
  )
  lower <- .thg_as_time(observation_start, x, "observation_start")
  upper <- .thg_as_time(observation_end, x, "observation_end")
  if (!is.null(lower) && length(lower) != 1L) .thg_bad_input("`observation_start` must be one time")
  if (!is.null(upper) && length(upper) != 1L) .thg_bad_input("`observation_end` must be one time")
  x$observation <- c(start = lower %||% span[[1L]], end = upper %||% span[[2L]])
  if (x$observation[["end"]] < x$observation[["start"]]) {
    .thg_bad_input("`observation_end` is earlier than `observation_start`")
  }
  x
}

# A sequence table is states only: no relational column is named or detected,
# no clock column is detected, and every column is categorical. Anything
# else takes the relational route, where a table without ends still fails
# with the message that names them.
.thg_is_sequence_input <- function(data, from, to, actor, group, time, start, end) {
  named <- !is.null(from) || !is.null(to) || !is.null(actor) || !is.null(group) ||
    !is.null(time) || !is.null(start) || !is.null(end)
  if (named) return(FALSE)
  if (inherits(data, c("tna", "netobject", "cograph_network"))) return(TRUE)
  if (is.matrix(data)) return(!is.numeric(data))
  if (!is.data.frame(data)) return(is.list(data))
  if (ncol(data) == 0L) return(FALSE)
  categorical <- vapply(data, function(column) {
    is.character(column) || is.factor(column) || all(is.na(column))
  }, logical(1L))
  if (!all(categorical)) return(FALSE)
  detected <- function(role) !is.null(.thg_match_column(data, role))
  relational <- (detected("from") && detected("to")) ||
    (detected("actor") && detected("group"))
  !relational && !detected("time") && !detected("start")
}

# One row per state occurrence: the session is the hyperedge, the state the
# member, and the position 1..n its contact time. Sessions are named as
# window_hypergraph() names them, so the two readings of one table agree.
.thg_sequence_memberships <- function(data) {
  if (inherits(data, c("tna", "netobject", "cograph_network"))) {
    data <- .coerce_sequence_input(data)
  }
  trajectories <- .wh_parse_input(data, action = NULL, actor = NULL, time = NULL)
  trajectories <- trajectories[lengths(trajectories) > 0L]
  if (!length(trajectories)) .thg_bad_input("`data` holds no non-empty sequence")
  n <- lengths(trajectories)
  data.frame(
    state = unlist(trajectories, use.names = FALSE),
    sequence = rep(names(trajectories), n),
    position = unlist(lapply(n, seq_len), use.names = FALSE),
    stringsAsFactors = FALSE
  )
}

# Which columns play which role, and therefore the input shape and the clock.
.thg_resolve_roles <- function(data, from, to, actor, group, time, start, end,
                               weight) {
  copresence <- !is.null(actor) || !is.null(group)
  if (copresence && (!is.null(from) || !is.null(to))) {
    .thg_bad_input("name either `from` and `to` or `actor` and `group`, not both")
  }
  claimed <- character()
  if (copresence) {
    actor <- .thg_resolve_column(data, actor, "actor", arg = "actor")
    group <- .thg_resolve_column(data, group, "group", exclude = actor, arg = "group")
    if (is.null(actor) || is.null(group)) {
      .thg_bad_input("co-presence data need both `actor` and `group`; name the one that could not be detected")
    }
    claimed <- c(actor, group)
  } else {
    from <- .thg_resolve_column(data, from, "from", arg = "from")
    to <- .thg_resolve_column(data, to, "to", exclude = from, arg = "to")
    if (is.null(from) || is.null(to)) {
      .thg_bad_input("name `from` and `to` for an edge list, or `actor` and `group` for co-presence data")
    }
    claimed <- c(from, to)
  }
  if (!is.null(start) && !is.null(time)) {
    .thg_bad_input("supply only one of `start` and `time`")
  }
  weight <- .thg_resolve_column(data, weight, "weight", exclude = claimed, arg = "weight")
  claimed <- c(claimed, weight)
  if (!is.null(time)) {
    clock <- .thg_resolve_column(data, time, "time", arg = "time")
    if (!is.null(end)) .thg_bad_input("`end` belongs to the interval format; name `start` with it, not `time`")
    format <- "contact"
    end <- NULL
  } else if (!is.null(start)) {
    clock <- .thg_resolve_column(data, start, "start", arg = "start")
    end <- .thg_resolve_column(data, end, "end", exclude = c(claimed, clock), arg = "end")
    format <- if (is.null(end)) "contact" else "interval"
  } else {
    detected_start <- .thg_match_column(data, "start", claimed)
    detected_end <- .thg_match_column(data, "end", c(claimed, detected_start))
    if (!is.null(end)) {
      end <- .thg_resolve_column(data, end, "end", arg = "end")
      if (is.null(detected_start)) .thg_bad_input("`end` needs a `start` column; name it")
      clock <- detected_start
      format <- "interval"
    } else if (!is.null(detected_start) && !is.null(detected_end)) {
      clock <- detected_start
      end <- detected_end
      format <- "interval"
    } else {
      clock <- .thg_match_column(data, "time", claimed) %||% detected_start
      if (is.null(clock)) {
        .thg_bad_input("no clock found: name `time` for a contact log or `start` and `end` for intervals")
      }
      end <- NULL
      format <- "contact"
    }
  }
  list(from = from, to = to, actor = actor, group = group, clock = clock,
       end = end, weight = weight, format = format)
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

# ---------------------------------------------------------------------------
# The measurement grid (Dynet's `.window_spec`): `start`/`end` bound the
# measured period, `step` is how often to look, `window` how much time each
# look covers (0 = a point; "all" = one window over the whole period), and
# `at` names instants directly. Without `step` and `at`, the grid is the
# event times of the hypergraph.
# ---------------------------------------------------------------------------

.thg_check_temporal <- function(x) {
  if (!inherits(x, "net_temporal_hypergraph")) {
    .thg_bad_input("`x` must come from temporal_hypergraph()")
  }
}

.thg_grid <- function(x, start = NULL, end = NULL, step = NULL, window = NULL,
                      at = NULL) {
  whole <- is.character(window) && length(window) == 1L && !is.na(window) &&
    identical(tolower(window), "all")
  if (!whole && !is.null(window) &&
      (!is.numeric(window) || length(window) != 1L || !is.finite(window) || window < 0)) {
    .thg_bad_input("`window` must be a single non-negative number, or the string \"all\"")
  }
  if (!is.null(step) && (!is.numeric(step) || length(step) != 1L ||
                         !is.finite(step) || step <= 0)) {
    .thg_bad_input("`step` must be a single positive number")
  }
  if (whole && !is.null(step)) {
    .thg_bad_input("`step` has no meaning with `window = \"all\"`, which measures one window")
  }
  if (!is.null(at) && !is.null(step)) {
    .thg_bad_input("supply either `at` (instants) or `step` (a grid), not both")
  }
  obs <- x$observation
  explicit <- isTRUE(x$params$observation_explicit)
  start <- .thg_as_time(start, x, "start")
  end <- .thg_as_time(end, x, "end")
  if (!is.null(start) && length(start) != 1L) .thg_bad_input("`start` must be one time")
  if (!is.null(end) && length(end) != 1L) .thg_bad_input("`end` must be one time")
  if (!is.null(start) && !is.null(end) && end < start) {
    .thg_bad_input(sprintf("`end` (%s) is earlier than `start` (%s)", end, start))
  }
  if (explicit) {
    query_start <- start %||% obs[["start"]]
    query_end <- end %||% obs[["end"]]
    if (query_end < obs[["start"]] || query_start > obs[["end"]]) {
      .thg_bad_input("the requested measurement range does not intersect the observation window",
                     class = "hypernets_outside_observation")
    }
    if (!is.null(start)) start <- max(start, obs[["start"]])
    if (!is.null(end)) end <- min(end, obs[["end"]])
  }
  start <- start %||% obs[["start"]]
  end <- end %||% obs[["end"]]

  if (!is.null(at)) {
    times <- .thg_as_time(at, x, "at")
    if (length(times) == 0L) .thg_bad_input("`at` contains no valid times")
    if (explicit && any(times < obs[["start"]] | times > obs[["end"]])) {
      .thg_bad_input("`at` lies outside the observation window",
                     class = "hypernets_outside_observation")
    }
    times <- sort(unique(times))
    window <- if (whole) end - start else as.numeric(window %||% 0)
  } else if (whole) {
    times <- start
    window <- end - start
  } else if (!is.null(step)) {
    times <- seq(start, end, by = step)
    window <- as.numeric(window %||% step)
  } else {
    times <- x$times[x$times >= start & x$times <= end]
    if (length(times) == 0L) times <- end
    window <- as.numeric(window %||% 0)
  }
  list(times = times, window = window, closed = whole)
}

# Logical over the edge metadata: which hyperedges a window from `t`
# measures. A point (window 0) is closed: an interval ending exactly at `t`
# is still active, a contact at `t` is present. A positive window covers
# [t, t + window), closed on the right only for the whole-period window.
# Which spells [start, end] (end NA = open) are in the window [t, t + window)
# -- closed on the right for a point or a closed grid -- or, cumulatively,
# have begun by its end.
.thg_spells_in_window <- function(start, end, t, window, mode, closed = FALSE) {
  upper <- t + window
  begun <- if (window > 0 && !closed) start < upper else start <= upper
  if (identical(mode, "cumulative")) return(begun)
  begun & (is.na(end) | end >= t)
}

.thg_edges_in_window <- function(x, t, window, mode, closed = FALSE) {
  ed <- x$edge_data
  if (identical(x$format, "contact") && !identical(mode, "cumulative")) {
    upper <- t + window
    begun <- if (window > 0 && !closed) ed$start < upper else ed$start <= upper
    return(begun & ed$start >= t)
  }
  .thg_spells_in_window(ed$start, ed$end, t, window, mode, closed)
}

# The memberships of the active hyperedges; when memberships carry their own
# times, only those present in the window.
.thg_memberships_in_window <- function(x, t, window, mode, closed = FALSE) {
  keep <- .thg_edges_in_window(x, t, window, mode, closed)
  mem <- x$memberships
  present <- mem$edge %in% x$edge_data$edge[keep]
  if (isTRUE(x$params$membership_times)) {
    present <- present & .thg_spells_in_window(mem$start, mem$end, t, window, mode, closed)
  }
  present
}

# Node universe of a snapshot. Without a universe, the nodes are the members
# of the active hyperedges. With one, every node is kept; when nodes carry
# entry times, only nodes entered by the end of the window plus any active
# member are kept.
.thg_snapshot_nodes <- function(x, t_end, active_members) {
  if (!isTRUE(x$params$nodes_given)) return(NULL)
  node_data <- x$node_data
  if (all(is.na(node_data$start))) return(node_data$node)
  entered <- !is.na(node_data$start) & node_data$start <= t_end
  sort(union(node_data$node[entered], unique(active_members)))
}

.thg_snapshot_at <- function(x, t, window, mode, multiedges, closed = FALSE) {
  ed <- x$edge_data
  d <- x$memberships[.thg_memberships_in_window(x, t, window, mode, closed), , drop = FALSE]
  sparse <- isTRUE(x$params$sparse)
  universe <- .thg_snapshot_nodes(x, t + window, d$member)
  if (nrow(d) == 0L) {
    hg <- .thg_empty_hypergraph(universe, sparse)
  } else {
    hg <- group_hypergraph(d, actor = "member", group = "edge", weight = "weight",
                           nodes = universe, sparse = sparse)
    hg$edge_data <- ed[match(colnames(hg$incidence), ed$edge), , drop = FALSE]
    rownames(hg$edge_data) <- NULL
  }
  hg$params$temporal_mode <- mode
  hg$params$at <- t
  hg$params$window <- window
  hg$params$format <- x$format
  hg$params$time_unit <- x$time_unit
  hg$params$origin <- x$origin
  if (!multiedges) hg <- .thg_collapse_duplicate_edges(hg)
  hg
}

.thg_check_mode <- function(mode, fn, arg = "mode") {
  if (identical(mode, "all")) {
    .thg_deprecated(sprintf("%s = \"all\"", arg),
                    sprintf("%s = \"cumulative\" (with `at` unset for the whole period)", arg),
                    fn)
    return("cumulative")
  }
  match.arg(mode, c("active", "cumulative"))
}

.thg_check_multiedges <- function(multiedges) {
  if (!is.logical(multiedges) || length(multiedges) != 1L || is.na(multiedges)) {
    .thg_bad_input("`multiedges` must be TRUE or FALSE")
  }
}

#' Extract one snapshot from a temporal hypergraph
#'
#' A snapshot is a static `net_hypergraph` holding the hyperedges a window
#' measures. `mode = "active"` follows the format: an interval hyperedge is
#' active on its closed interval, a contact hyperedge at its instant.
#' `mode = "cumulative"` keeps every hyperedge begun by the end of the
#' window, the growing view of a contact log; with `at` unset it is the
#' static aggregate of everything observed.
#'
#' @param x A [temporal_hypergraph()].
#' @param at One time value, as a number on the hypergraph's clock or a date
#'   for a calendar hypergraph. `NULL` (default) is the end of observation.
#' @param window How much time the snapshot covers, in the hypergraph's time
#'   unit: `0` (default) evaluates the instant `at`; a positive width covers
#'   `[at, at + window)`; `"all"` covers the whole observation period.
#' @param mode `"active"` (default) or `"cumulative"`, see Description. The
#'   former `"all"` is deprecated: it is `"cumulative"` with `at` unset.
#' @param multiedges Keep distinct edge identities with identical member sets?
#'   `TRUE` matches a multi-hypergraph; `FALSE` collapses them to one edge and
#'   records their counts in `edge_multiplicity`.
#' @return A static `net_hypergraph` usable by every hypernets hypergraph verb;
#'   its `params` record `at`, `window`, `temporal_mode` and the clock.
#' @export
hypergraph_snapshot <- function(x, at = NULL, window = 0,
                                mode = c("active", "cumulative"),
                                multiedges = TRUE) {
  .thg_check_temporal(x)
  mode <- .thg_check_mode(mode, "hypergraph_snapshot")
  .thg_check_multiedges(multiedges)
  if (!is.null(at) && length(at) != 1L) {
    .thg_bad_input("`at` must be one time value; use hypergraph_snapshots() for several")
  }
  whole <- is.character(window) && length(window) == 1L && !is.na(window) &&
    identical(tolower(window), "all")
  if (whole && !is.null(at)) {
    .thg_bad_input("`at` has no meaning with `window = \"all\"`, which covers the whole period")
  }
  grid <- .thg_grid(x, window = window,
                    at = if (whole) NULL else at %||% x$observation[["end"]])
  .thg_snapshot_at(x, grid$times, grid$window, mode, multiedges, grid$closed)
}

#' Extract a sequence of temporal-hypergraph snapshots
#'
#' Snapshots on a measurement grid with Dynet's meaning: `start` and `end`
#' bound the measured period (default: the observation window), `step` is
#' how often to look, `window` how much time each look covers, and `at`
#' names instants directly. Without `step` and `at`, the grid is the event
#' times of the hypergraph, each looked at as a point.
#'
#' @inheritParams hypergraph_snapshot
#' @param start,end Bounds of the measured period, as numbers on the
#'   hypergraph's clock or dates for a calendar hypergraph. Default: the
#'   observation window.
#' @param step How often to measure, in the hypergraph's time unit. `NULL`
#'   (default) measures at every event time.
#' @param window How much time each measurement covers. Defaults to `step`
#'   when a step is given (a partition of the period) and to `0` otherwise
#'   (a point sample); larger than `step` gives a rolling window; `"all"`
#'   measures the whole period as one window.
#' @param at Instants to measure instead of a grid; may be a vector.
#' @return A named list of `net_hypergraph` objects with class
#'   `net_hypergraph_snapshots`, named by the window starts on the
#'   hypergraph's clock.
#' @examples
#' seats <- data.frame(
#'   case = rep(c("A", "B", "C"), each = 3),
#'   arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
#'   constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3)
#' )
#' thg <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
#'                            start = "constituted", end = "concluded")
#' yearly <- hypergraph_snapshots(thg, step = 2)
#' names(yearly)
#' @export
hypergraph_snapshots <- function(x, start = NULL, end = NULL, step = NULL,
                                 window = NULL, mode = c("active", "cumulative"),
                                 at = NULL, multiedges = TRUE) {
  .thg_check_temporal(x)
  mode <- .thg_check_mode(mode, "hypergraph_snapshots")
  .thg_check_multiedges(multiedges)
  grid <- .thg_grid(x, start, end, step, window, at)
  out <- lapply(grid$times, function(t) {
    .thg_snapshot_at(x, t, grid$window, mode, multiedges, grid$closed)
  })
  names(out) <- make.unique(as.character(grid$times))
  structure(out, class = c("net_hypergraph_snapshots", "list"),
            time_unit = x$time_unit, origin = x$origin, window = grid$window)
}

#' @export
print.net_temporal_hypergraph <- function(x, ...) {
  cat(sprintf("Temporal hypergraph: %d nodes, %d hyperedges, %d event times\n",
              length(x$nodes), length(x$edges), length(x$times)))
  cat(if (isTRUE(x$params$membership_times)) {
    paste0("Format: interval; memberships carry their own times ",
           "(a hyperedge spans its first to its last membership)\n")
  } else if (identical(x$format, "interval")) {
    "Format: interval (a hyperedge is active from its start to its end)\n"
  } else {
    paste0("Format: contact (a hyperedge is an instantaneous event; ",
           "mode = \"cumulative\" keeps every hyperedge once it appears)\n")
  })
  cat(sprintf("Time: %s; observed from %s to %s%s\n",
              .thg_clock_label(x$origin, x$time_unit),
              format(x$observation[["start"]]), format(x$observation[["end"]]),
              if (isTRUE(x$params$observation_explicit)) " (declared)" else ""))
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
    format = object$format,
    time_unit = object$time_unit,
    row.names = NULL
  )
}

#' Tidy tables of a temporal hypergraph
#'
#' @param x A [temporal_hypergraph()].
#' @param row.names,optional Ignored; present for compatibility.
#' @param what `"memberships"` (default) for one row per node-in-hyperedge
#'   spell (`member`, `edge`, `start`, `end`, `weight`), `"edges"` for one
#'   row per hyperedge with its clock and attributes, `"nodes"` for the node
#'   universe with entry times. Times are on the hypergraph's clock.
#' @param ... Ignored.
#' @return A base `data.frame`.
#' @export
as.data.frame.net_temporal_hypergraph <- function(x, row.names = NULL,
                                                   optional = FALSE,
                                                   what = c("memberships", "edges",
                                                            "nodes"), ...) {
  what <- match.arg(what)
  switch(what, memberships = x$memberships, edges = x$edge_data, nodes = x$node_data)
}

#' Plot a temporal-hypergraph snapshot through cograph
#'
#' @param x A [temporal_hypergraph()].
#' @param at Snapshot time; defaults to the end of observation.
#' @param mode Snapshot mode passed to [hypergraph_snapshot()].
#' @param method Projection weighting passed to [hg_project()].
#' @param ... Additional arguments passed to [cograph::splot()].
#' @return The cograph plot object, invisibly when rendered interactively.
#' @export
plot.net_temporal_hypergraph <- function(x, at = NULL,
                                         mode = c("active", "cumulative"),
                                         method = c("association", "clique"), ...) {
  mode <- .thg_check_mode(mode, "plot")
  method <- match.arg(method)
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
