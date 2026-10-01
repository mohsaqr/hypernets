# ---- Idiographic higher-order networks: one BuildHON model per unit -------
#
# build_hon(data, by = "<unit>") fits an INDEPENDENT BuildHON model to each
# unit's own trajectories (design (a) of ROADMAP-analytics.md Tier C). Each
# per-unit model is literally build_hon() on that unit's rows, so it is
# identical() to fitting the unit alone -- the test suite asserts this for
# every unit. The alternative (b), a pooled higher-order node set with
# per-unit probabilities on it, is not built: it imposes the pooled
# dependency structure on every person, which is the nomothetic assumption
# idiographic analysis exists to test (Molenaar 2004), and no published
# method defines it for BuildHON.
#
# Support guard: a unit is fitted only when its first-order transitions reach
# the number of free parameters of a first-order transition matrix on the
# pooled state space, S * (S - 1) (Anderson & Goodman 1957: the ML estimate
# of each row is n_ij / n_i, so a unit with fewer transitions than free
# parameters cannot estimate its own chain, let alone a higher-order one).

#' Idiographic higher-order networks: one model per unit
#'
#' @name net_hon_group
#' @description
#' `build_hon(data, by = "<unit>")` fits one Higher-Order Network per unit
#' (an actor: a student, a team, a project) with the BuildHON algorithm (Xu,
#' Wickramarathne & Chawla 2016), each on that unit's own trajectories alone,
#' and returns them together as a `net_hon_group`. Every per-unit model is
#' exactly the `net_hon` that [build_hon()] returns for that unit's data, so
#' the container adds bookkeeping, never a different estimate.
#'
#' @section Why one independent model per unit:
#' Pooling trajectories across people and reading the pooled model as a
#' description of each person assumes the process is the same in everyone.
#' Molenaar (2004) shows that inferences from inter-individual variation
#' transfer to the individual only under ergodicity, which learning and
#' behavioural processes rarely satisfy; in learning analytics, Saqr (2024)
#' found that group-level engagement networks poorly reflect individual
#' students' processes. Each unit therefore gets its own higher-order
#' structure: the dependency orders BuildHON detects are part of what is
#' allowed to differ between people, as in person-specific network modelling
#' (Beltz et al. 2016).
#'
#' A consequence is that two units' networks need not share nodes: a memory
#' node `"A -> B"` appears for a unit only where that unit's data support it.
#' [as.data.frame()] therefore reports each unit's rules in one long table,
#' and a rule absent for a unit is absent, not zero.
#'
#' @section Support guard:
#' Most units in learning data have short sequences. A unit's support is its
#' number of first-order transitions (after `collapse_repeats`). By default a
#' unit is fitted only when that number reaches
#' \eqn{S(S-1)}, the number of free parameters of a first-order transition
#' matrix over the \eqn{S} states observed in the pooled data: the maximum
#' likelihood estimate of each row of a Markov chain is \eqn{n_{ij}/n_i}
#' (Anderson & Goodman 1957), and a unit with fewer transitions than free
#' parameters cannot estimate its own first-order chain, let alone a
#' higher-order one. Pass `min_transitions` to state a different floor.
#'
#' Units below the floor raise one warning of class
#' `hypernets_insufficient_support` naming them (the condition carries the
#' full list in its `actors` field) and are left out of the container. With
#' `keep_sparse = TRUE` they are fitted anyway and flagged
#' `sufficient = FALSE` in `as.data.frame(x, what = "actors")`; the warning
#' still fires. A unit with no transition at all can never be fitted and is
#' recorded with `fitted = FALSE`. If no unit can be fitted, the same class
#' is raised as an error.
#'
#' @section Not included:
#' Distances, clustering and outcome models over the per-unit networks are
#' not provided: the per-unit networks have different node sets, and this
#' package ships no published distance between BuildHON networks of
#' different structure.
#'
#' @section Input:
#' `by` names a column of a wide `data.frame` (removed before the states are
#' read) or of a long `data.frame` given with `action` (and optionally
#' `actor`, `time`), or, for a list of trajectories, is a vector with one
#' entry per trajectory. `tna` / `netobject` input carries no unit column and
#' is refused.
#'
#' @section Conditions:
#' * `hypernets_bad_input` -- `by` missing from the data, of the wrong
#'   length, containing `NA`, given with model-object input; `min_transitions`
#'   not a single whole number \eqn{\ge 1}; `min_transitions` or
#'   `keep_sparse` given without `by`.
#' * `hypernets_insufficient_support` -- warning when some units fall below
#'   the floor; error when no unit can be fitted.
#'
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances}, 2(5),
#' e1600028. \doi{10.1126/sciadv.1600028}
#'
#' Saebi, M., Xu, J., Kaplan, L. M., Ribeiro, B., & Chawla, N. V. (2020).
#' Efficient modeling of higher-order dependencies in networks: from
#' algorithm to application for anomaly detection. \emph{EPJ Data Science},
#' 9(1), 15.
#'
#' Molenaar, P. C. M. (2004). A manifesto on psychology as idiographic
#' science: Bringing the person back into scientific psychology, this time
#' forever. \emph{Measurement: Interdisciplinary Research and Perspectives},
#' 2(4), 201--218. \doi{10.1207/s15366359mea0204_1}
#'
#' Beltz, A. M., Wright, A. G. C., Sprague, B. N., & Molenaar, P. C. M.
#' (2016). Bridging the nomothetic and idiographic approaches to the analysis
#' of clinical data. \emph{Assessment}, 23(4), 447--458.
#' \doi{10.1177/1073191116648209}
#'
#' Saqr, M. (2024). Group-level analysis of engagement poorly reflects
#' individual students' processes: Why we need idiographic learning
#' analytics. \emph{Computers in Human Behavior}, 150, 107991.
#' \doi{10.1016/j.chb.2023.107991}
#'
#' Anderson, T. W., & Goodman, L. A. (1957). Statistical inference about
#' Markov chains. \emph{The Annals of Mathematical Statistics}, 28(1),
#' 89--110. \doi{10.1214/aoms/1177707039}
#'
#' @seealso [build_hon()], [as.data.frame.net_hon_group()],
#'   [plot.net_hon_group()]
#' @examples
#' fits <- build_hon(human_long, action = "code", actor = "session_id",
#'                   time = "timestamp", by = "project", max_order = 2L)
#' fits
#' as.data.frame(fits, what = "actors")
#' as.data.frame(fits, actor = "Project_7", order_min = 2L)
NULL

# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

#' Build a net_hon_group (engine behind build_hon(by =))
#' @noRd
.hon_group_build <- function(data, by, max_order, min_freq, collapse_repeats,
                             method, action, actor, time, min_transitions,
                             keep_sparse) {
  .hon_group_check_args(min_transitions, keep_sparse)
  units <- .hon_group_units(data, by)
  ids <- sort(unique(units))
  subsets <- lapply(ids, function(id) {
    .hon_group_subset(data, by, units, id, drop_by = is.null(action))
  })
  names(subsets) <- ids

  # Support is counted on exactly the trajectories build_hon() would fit.
  trajectories <- lapply(subsets, .hon_group_trajectories, action = action,
                         actor = actor, time = time,
                         collapse_repeats = collapse_repeats)
  pooled_states <- sort(unique(unlist(trajectories, use.names = FALSE)))
  n_pooled <- length(pooled_states)
  floor_used <- as.integer(min_transitions %||%
                             max(n_pooled * (n_pooled - 1L), 1L))

  n_transitions <- vapply(trajectories,
                          \(tr) as.integer(sum(lengths(tr) - 1L)), integer(1L))
  support <- data.frame(
    actor = ids,
    n_sequences = vapply(trajectories, length, integer(1L)),
    n_transitions = n_transitions,
    n_states = vapply(trajectories,
                      \(tr) length(unique(unlist(tr, use.names = FALSE))),
                      integer(1L)),
    min_transitions = rep.int(floor_used, length(ids)),
    sufficient = n_transitions >= floor_used,
    stringsAsFactors = FALSE
  )
  support$fitted <- n_transitions >= 1L &
    (support$sufficient | isTRUE(keep_sparse))
  rownames(support) <- NULL

  sparse <- support$actor[!support$sufficient]
  if (!any(support$fitted)) {
    stop(errorCondition(
      sprintf(paste0("No unit of `%s` has enough transitions to fit a model ",
                     "(floor %d; largest unit has %d). Lower ",
                     "`min_transitions` or pool the data."),
              .hon_group_by_label(by), floor_used, max(n_transitions)),
      actors = sparse,
      class = c("hypernets_insufficient_support", "hypernets_bad_input"),
      call = NULL))
  }
  if (length(sparse)) .hon_group_warn_sparse(sparse, floor_used, keep_sparse)

  fit_ids <- support$actor[support$fitted]
  models <- lapply(fit_ids, function(id) {
    build_hon(subsets[[id]], max_order = max_order, min_freq = min_freq,
              collapse_repeats = collapse_repeats, method = method,
              action = action, actor = actor, time = time)
  })
  names(models) <- fit_ids

  structure(
    list(
      models = models,
      support = support,
      by = .hon_group_by_label(by),
      pooled_states = pooled_states,
      min_transitions = floor_used,
      floor_rule = if (is.null(min_transitions)) "S(S-1)" else "user",
      keep_sparse = isTRUE(keep_sparse),
      max_order = as.integer(max_order),
      min_freq = as.integer(min_freq),
      method = method
    ),
    class = "net_hon_group"
  )
}

#' Validate the support-guard arguments
#' @noRd
.hon_group_check_args <- function(min_transitions, keep_sparse) {
  ok_floor <- is.null(min_transitions) ||
    (is.numeric(min_transitions) && length(min_transitions) == 1L &&
       is.finite(min_transitions) && min_transitions >= 1 &&
       min_transitions == round(min_transitions))
  if (!ok_floor) {
    stop(errorCondition(
      "`min_transitions` must be NULL or a single whole number >= 1.",
      class = "hypernets_bad_input", call = NULL))
  }
  if (!(is.logical(keep_sparse) && length(keep_sparse) == 1L &&
        !is.na(keep_sparse))) {
    stop(errorCondition("`keep_sparse` must be TRUE or FALSE.",
                        class = "hypernets_bad_input", call = NULL))
  }
  invisible(NULL)
}

#' The unit label of every row (data.frame) or trajectory (list)
#' @noRd
.hon_group_units <- function(data, by) {
  if (inherits(data, c("tna", "netobject", "cograph_network"))) {
    stop(errorCondition(
      paste0("`by` needs the unit of every trajectory, which a model object ",
             "does not carry; pass the sequence data.frame instead."),
      class = "hypernets_bad_input", call = NULL))
  }
  units <- if (is.data.frame(data)) {
    if (!(is.character(by) && length(by) == 1L && by %in% names(data))) {
      stop(errorCondition(
        "`by` must name one column of `data`.",
        class = "hypernets_bad_input", call = NULL))
    }
    data[[by]]
  } else if (is.list(data)) {
    if (length(by) != length(data)) {
      stop(errorCondition(
        sprintf(paste0("For list input `by` must have one entry per ",
                       "trajectory (%d), not %d."),
                length(data), length(by)),
        class = "hypernets_bad_input", call = NULL))
    }
    by
  } else {
    stop(errorCondition(
      "`data` must be a data.frame or a list of trajectories when `by` is given.",
      class = "hypernets_bad_input", call = NULL))
  }
  if (anyNA(units)) {
    stop(errorCondition(
      sprintf("`by` has %d missing unit label(s); drop or label them first.",
              sum(is.na(units))),
      class = "hypernets_bad_input", call = NULL))
  }
  as.character(units)
}

#' One unit's data, in the shape build_hon() takes
#' @noRd
.hon_group_subset <- function(data, by, units, id, drop_by) {
  keep <- units == id
  if (!is.data.frame(data)) return(data[keep])
  # A wide frame reads every column as a time step, so the unit column goes.
  cols <- if (drop_by) setdiff(names(data), by) else names(data)
  sub <- data[keep, cols, drop = FALSE]
  rownames(sub) <- NULL
  sub
}

#' One unit's trajectories, parsed exactly as build_hon() parses them
#' @noRd
.hon_group_trajectories <- function(sub, action, actor, time,
                                    collapse_repeats) {
  if (!is.null(action)) {
    return(.hi_parse(sub, action = action, actor = actor, time = time,
                     collapse_repeats = collapse_repeats, allow_empty = TRUE))
  }
  .hon_parse_input(sub, collapse_repeats = collapse_repeats,
                   verb = "build_hon")
}

#' Printable name of the unit variable
#' @noRd
.hon_group_by_label <- function(by) {
  if (is.character(by) && length(by) == 1L) by else "unit"
}

#' Raise the classed sparse-support warning
#' @noRd
.hon_group_warn_sparse <- function(sparse, floor_used, keep_sparse) {
  shown <- utils::head(sparse, 10L)
  more <- if (length(sparse) > 10L) {
    sprintf(" and %d more", length(sparse) - 10L)
  } else {
    ""
  }
  fate <- if (isTRUE(keep_sparse)) {
    "fitted anyway (keep_sparse = TRUE) and flagged `sufficient = FALSE`"
  } else {
    "left out of the result"
  }
  warning(warningCondition(
    sprintf("%d unit(s) have fewer than %d transitions and were %s: %s%s.",
            length(sparse), floor_used, fate,
            paste(shown, collapse = ", "), more),
    actors = sparse,
    class = "hypernets_insufficient_support", call = NULL))
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Tidy tables of a net_hon_group
#'
#' @description
#' The supported way into a [net_hon_group]: every unit's higher-order rules
#' in one long table, or one row per unit with its support.
#'
#' @param x A `net_hon_group`, from `build_hon(..., by = )`.
#' @param row.names,optional,... Unused; S3 consistency.
#' @param what `"edges"` (default): one row per unit per rule edge. `"actors"`:
#'   one row per unit, including units that were not fitted.
#' @param actor `NULL` (every unit) or a character vector of unit labels to
#'   keep.
#' @param order_min Integer or `NULL`. For `what = "edges"`, keep only rules
#'   whose context is at least this order (`2` for the memory rules).
#' @param sort_by `NULL` (unit, then construction order), `"count"` or
#'   `"probability"`: largest first, ties broken by `actor`, `from`, `to`.
#'   For `what = "edges"` only.
#' @param top Integer or `NULL`: keep the first `top` rows, applied last.
#' @return A base `data.frame`.
#'   For `what = "edges"`, one row per fitted unit per rule edge, columns
#'   `actor`, `path`, `from`, `to`, `count`, `probability`, `from_order`,
#'   `to_order` -- the columns of [as.data.frame.net_hon()] with the unit
#'   first. Within a unit, `probability` sums to 1 over the rules leaving
#'   each `from` context.
#'   For `what = "actors"`, one row per unit with `actor`, `n_sequences`,
#'   `n_transitions`, `n_states`, `min_transitions` (the floor applied),
#'   `sufficient`, `fitted`, and for fitted units `n_nodes`, `n_edges`,
#'   `max_order_observed` (`NA` otherwise).
#' @section Conditions:
#' Raises `hypernets_bad_input` for an `actor` label that is not a unit of
#' `x`.
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances}, 2(5),
#' e1600028.
#' @examples
#' fits <- build_hon(human_long, action = "code", actor = "session_id",
#'                   time = "timestamp", by = "project", max_order = 2L)
#' as.data.frame(fits, what = "actors")
#' as.data.frame(fits, order_min = 2L, sort_by = "count", top = 10L)
#' @export
as.data.frame.net_hon_group <- function(x, row.names = NULL, optional = FALSE,
                                        ..., what = c("edges", "actors"),
                                        actor = NULL, order_min = NULL,
                                        sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  .hon_group_check_actor(x, actor)
  if (what == "actors") {
    out <- .hon_group_actor_table(x)
    if (!is.null(actor)) out <- out[out$actor %in% actor, , drop = FALSE]
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  ids <- names(x$models)
  if (!is.null(actor)) ids <- ids[ids %in% actor]
  parts <- lapply(ids, function(id) {
    rules <- as.data.frame(x$models[[id]], order_min = order_min)
    cbind(actor = rep.int(id, nrow(rules)), rules, stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, parts)
  if (is.null(out) || !nrow(out)) {
    out <- data.frame(actor = character(0), path = character(0),
                      from = character(0), to = character(0),
                      count = integer(0), probability = numeric(0),
                      from_order = integer(0), to_order = integer(0),
                      stringsAsFactors = FALSE)
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("count", "probability"))
    out <- out[order(-out[[sort_by]], out$actor, out$from, out$to), ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Refuse unknown unit labels
#' @noRd
.hon_group_check_actor <- function(x, actor) {
  if (is.null(actor)) return(invisible(NULL))
  unknown <- setdiff(as.character(actor), x$support$actor)
  if (!is.character(actor) || length(unknown)) {
    stop(errorCondition(
      sprintf("`actor` names no unit of this model: %s.",
              paste(unknown, collapse = ", ")),
      class = "hypernets_bad_input", call = NULL))
  }
  invisible(NULL)
}

#' Support table plus per-model size
#' @noRd
.hon_group_actor_table <- function(x) {
  out <- x$support
  model_stat <- function(field) {
    vapply(out$actor, function(id) {
      m <- x$models[[id]]
      if (is.null(m)) NA_integer_ else as.integer(m[[field]])
    }, integer(1L), USE.NAMES = FALSE)
  }
  out$n_nodes <- model_stat("n_nodes")
  out$n_edges <- model_stat("n_edges")
  out$max_order_observed <- model_stat("max_order_observed")
  out
}

#' Print a net_hon_group
#'
#' @param x A `net_hon_group`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' fits <- build_hon(human_long, action = "code", actor = "session_id",
#'                   time = "timestamp", by = "project", max_order = 2L)
#' print(fits)
#' @export
print.net_hon_group <- function(x, ...) {
  s <- x$support
  cat("Idiographic Higher-Order Networks (one BuildHON model per unit)\n")
  cat(sprintf("  Units (%s):   %d fitted of %d\n", x$by, sum(s$fitted),
              nrow(s)))
  cat(sprintf("  Support floor: %d transitions (%s, %d pooled states)\n",
              x$min_transitions,
              if (identical(x$floor_rule, "S(S-1)")) "S(S-1)" else "stated",
              length(x$pooled_states)))
  below <- s$actor[!s$sufficient]
  if (length(below)) {
    cat(sprintf("  Below floor:   %d (%s)\n", length(below),
                if (x$keep_sparse) "kept, flagged" else "excluded"))
  }
  cat(sprintf("  Max order:     %d requested | Min freq: %d | Method: %s\n",
              x$max_order, x$min_freq, x$method))
  cat("  Tables: as.data.frame(x) for rules, what = \"actors\" for support.\n")
  invisible(x)
}

#' Summarise a net_hon_group
#'
#' @param object A `net_hon_group`.
#' @param ... Unused.
#' @return A base `data.frame`, one row per unit: the table of
#'   `as.data.frame(object, what = "actors")` (support, fit flags and
#'   per-model size).
#' @examples
#' fits <- build_hon(human_long, action = "code", actor = "session_id",
#'                   time = "timestamp", by = "project", max_order = 2L)
#' summary(fits)
#' @export
summary.net_hon_group <- function(object, ...) {
  as.data.frame(object, what = "actors")
}

#' Plot a net_hon_group
#'
#' @description
#' Three views, each delegated to cograph:
#' \describe{
#'   \item{`what = "edges"` (default)}{[cograph::plot_heatmap()] of units
#'     (rows) by rule paths (columns), cells the per-unit transition
#'     probability. A blank cell means the rule is absent from that unit's
#'     network, not zero. Columns are the `top` paths with the largest count
#'     summed over units (default 30); `order_min = 2` shows memory rules
#'     only.}
#'   \item{`what = "support"`}{[cograph::plot_centrality()] of each unit's
#'     `n_transitions`, units ordered by it; the units at or above the
#'     support floor are highlighted.}
#'   \item{`what = "network"`}{[cograph::plot_tna()] of one unit's
#'     higher-order network (`actor =` names it).}
#' }
#' @param x A `net_hon_group`.
#' @param what `"edges"`, `"support"` or `"network"`.
#' @param actor For `what = "network"`, the single unit to draw (required);
#'   for `what = "edges"`, optional units to keep as rows.
#' @param order_min For `what = "edges"`: keep rules of at least this order.
#' @param top For `what = "edges"`: number of paths shown (default 30).
#' @param ... Passed to the cograph function, overriding the defaults set
#'   here.
#' @return For `"edges"` and `"support"`, the `ggplot` from cograph; for
#'   `"network"`, `x` invisibly (cograph draws with base graphics).
#' @section Conditions:
#' Raises `hypernets_bad_input` for an unknown `what`, for `what =
#' "network"` without exactly one fitted `actor`, and when the edge view has
#' no rule left to draw.
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances}, 2(5),
#' e1600028.
#' @examples
#' fits <- build_hon(human_long, action = "code", actor = "session_id",
#'                   time = "timestamp", by = "project", max_order = 2L)
#' plot(fits, top = 20L)
#' plot(fits, what = "support")
#' plot(fits, what = "network", actor = "Project_7")
#' @export
plot.net_hon_group <- function(x, what = c("edges", "support", "network"),
                               actor = NULL, order_min = NULL, top = 30L,
                               ...) {
  views <- c("edges", "support", "network")
  if (missing(what)) what <- "edges"
  if (!is.character(what) || length(what) != 1L || !(what %in% views)) {
    stop(errorCondition(
      sprintf("`what` must be one of %s.",
              paste0("\"", views, "\"", collapse = ", ")),
      class = "hypernets_bad_input", call = NULL))
  }
  switch(what,
         edges = .hon_group_plot_edges(x, actor, order_min, top, ...),
         support = .hon_group_plot_support(x, ...),
         network = .hon_group_plot_network(x, actor, ...))
}

#' Units x paths probability heatmap through cograph::plot_heatmap()
#' @noRd
.hon_group_plot_edges <- function(x, actor, order_min, top, ...) {
  rules <- as.data.frame(x, actor = actor, order_min = order_min)
  if (!nrow(rules)) {
    stop(errorCondition("No rule left to draw for these units / `order_min`.",
                        class = "hypernets_bad_input", call = NULL))
  }
  total <- tapply(rules$count, rules$path, sum)
  paths <- names(total)[order(-total, names(total))]
  paths <- utils::head(paths, as.integer(top %||% length(paths)))
  units <- unique(rules$actor)
  M <- matrix(NA_real_, length(units), length(paths),
              dimnames = list(units, paths))
  shown <- rules[rules$path %in% paths, , drop = FALSE]
  M[cbind(shown$actor, shown$path)] <- shown$probability
  defaults <- list(
    x = M, colors = c("#FFFFFF", "#4A6FE3"), limits = c(0, 1),
    na_color = "grey92", legend_title = "Probability",
    show_values = FALSE, show_diagonal = TRUE, aspect_ratio = NULL,
    xlab = "rule (context -> next state)", ylab = x$by,
    title = sprintf("Per-%s transition probabilities", x$by),
    subtitle = sprintf(paste0("%d of %d rules by pooled count; blank = rule ",
                              "absent from that %s's network"),
                       length(paths), length(total), x$by))
  do.call(cograph::plot_heatmap, utils::modifyList(defaults, list(...)))
}

#' Per-unit support through cograph::plot_centrality()
#' @noRd
.hon_group_plot_support <- function(x, ...) {
  s <- x$support
  # Only the floor measure: plot_centrality() highlights the top values of
  # every panel, so a second panel would highlight units unrelated to the
  # floor.
  table <- data.frame(node = s$actor, n_transitions = s$n_transitions)
  defaults <- list(
    x = table, style = "dot", order_by = "n_transitions",
    highlight = sum(s$sufficient),
    title = sprintf("Support per %s", x$by),
    subtitle = sprintf("Highlighted: the %d unit(s) at or above the floor of %d transitions",
                       sum(s$sufficient), x$min_transitions))
  p <- do.call(cograph::plot_centrality,
               utils::modifyList(defaults, list(...)))
  p + ggplot2::labs(x = "first-order transitions")
}

#' One unit's HON through cograph::plot_tna()
#' @noRd
.hon_group_plot_network <- function(x, actor, ...) {
  if (!(is.character(actor) && length(actor) == 1L &&
        actor %in% names(x$models))) {
    stop(errorCondition(
      "`what = \"network\"` needs `actor =` naming one fitted unit.",
      class = "hypernets_bad_input", call = NULL))
  }
  m <- x$models[[actor]]
  defaults <- list(x = m$weights, minimum = 0.05, label_size = 0.7,
                   title = sprintf("%s: %d nodes, max order %d", actor,
                                   m$n_nodes, m$max_order_observed))
  do.call(cograph::plot_tna, utils::modifyList(defaults, list(...)))
  invisible(x)
}
