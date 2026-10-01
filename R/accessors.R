# ---------------------------------------------------------------------------
# Shared accessor plumbing
#
# Every hypernets accessor that can return many rows takes `top =`. It is
# applied LAST -- after `what`, after every filter (`order_min`, `min_count`,
# `dim`, `k`, `dimension`, `significant`, ...) and after `sort_by` -- so
# `sort_by` and `top` compose: `top = n` means "the first n rows of the table
# as it would otherwise have been returned". Semantics match
# the path counters' `top =`, which shipped first.
#
# This exists so that a caller never has to write head() on the public
# surface. A vignette line that subsets, sorts or filters a returned table is
# a defect; the accessor owns the argument instead.
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# One vocabulary for column names
#
# Every table hg_get() returns names the same quantity the same way:
# `count` (observations of the row's unit), `expected`, `z`, `p_value`,
# `p_adj`, `significant` (a test), `ci_lower` / `ci_upper`,
# `log_likelihood`, `df`, `aic`, `bic`, `community`, `run`, `dimension`,
# `node` (one node), `members` (a set of nodes, comma separated), `path`
# (an ordered sequence of states). Tables that come from an imported
# estimator are renamed on the way out with .ho_rename().
# ---------------------------------------------------------------------------

#' Rename columns of a data frame by a map of old to new names
#'
#' @param d A data.frame.
#' @param map Named character vector, `c(old = "new")`; names absent from
#'   `d` are ignored.
#' @return `d` with the mapped columns renamed.
#' @noRd
.ho_rename <- function(d, map) {
  hit <- names(d) %in% names(map)
  names(d)[hit] <- map[names(d)[hit]]
  d
}

#' Print the first rows of a result's default table
#'
#' Every print method ends here: the header lines are the method's own, the
#' rows are those of `hg_get(x)`.
#'
#' @param x A result object with an `hg_get()` method.
#' @param n Number of rows to show.
#' @param ... Passed to `hg_get()`.
#' @return `NULL`, invisibly.
#' @noRd
.ho_print_table <- function(x, n = 10L, ...) {
  stopifnot("`n` must be a single whole number >= 1" =
              is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1)
  table <- hg_get(x, ...)
  if (nrow(table) == 0L) {
    cat("(no rows)\n")
    return(invisible(NULL))
  }
  shown <- utils::head(table, as.integer(n))
  print(shown, row.names = FALSE)
  rest <- nrow(table) - nrow(shown)
  if (rest > 0L) {
    cat(sprintf("... %d more row%s\n", rest, if (rest == 1L) "" else "s"))
  }
  invisible(NULL)
}

#' Truncate a tidy accessor result to its first `top` rows
#'
#' @param x A data.frame, already filtered and ordered.
#' @param top Integer or NULL. NULL returns `x` unchanged.
#' @return `x`, or its first `top` rows with row names reset.
#' @noRd
.ho_top <- function(x, top) {
  if (is.null(top)) return(x)
  stopifnot(
    "`top` must be a single whole number >= 1" =
      is.numeric(top) && length(top) == 1L && is.finite(top) &&
      top == round(top) && top >= 1
  )
  out <- utils::head(x, as.integer(top))
  rownames(out) <- NULL
  out
}

# ---------------------------------------------------------------------------
# hg_get(): the one reader
#
# Every hypernets object -- a hypergraph, a memory network, a simplicial
# complex, an inference or community result -- is read through hg_get().
# Each class registers an hg_get.<class> method; `what =` picks a secondary
# table and every filter, ordering and `top` is a named argument of that
# method. hypernets registers no as.data.frame() method.
# ---------------------------------------------------------------------------

#' Read a table from any hypernets object
#'
#' `hg_get()` is the single reader of the package. Every object a
#' constructor or a verb returns -- a hypergraph, a memory network, a
#' simplicial complex, an inference or a community result -- hands over its
#' tables through `hg_get()`, never through `$`. `what =` selects the table;
#' `what = NULL` (the default) returns the object's primary table. Filters,
#' orderings and `top =` are named arguments of the method for that class.
#'
#' @section Tables by class:
#' \describe{
#'   \item{`net_hg` (hypergraphs)}{`"edges"` (default), `"nodes"`,
#'     `"memberships"`, `"sets"`, `"state_counts"`, `"node_data"`,
#'     `"edge_data"`, `"incidence_data"`; see [hg_get.net_hg()].}
#'   \item{`text_hypergraph`}{`"weights"` (default), `"documents"`,
#'     `"vocabulary"`, `"sentences"`; see [hg_get.text_hypergraph()].}
#'   \item{`net_temporal_hypergraph`}{`"memberships"` (default), `"edges"`,
#'     `"nodes"`.}
#'   \item{`net_hon` ([hon()])}{`"rules"` (default), `"nodes"`,
#'     `"pathways"`.}
#'   \item{`net_honem` ([honem()])}{`"embeddings"` (default), `"variance"`.}
#'   \item{`net_hypa` ([hypa()] on sequences)}{`"scores"` (default),
#'     `"over"`, `"under"`, `"pathways"`.}
#'   \item{`net_mogen` ([mogen()])}{`"orders"` (default), `"transitions"`,
#'     `"paths"`, `"pathways"`.}
#'   \item{`net_markov_order` ([markov_order()])}{`"orders"` (default),
#'     `"null"`.}
#'   \item{`net_path_dependence` ([memory()])}{`"contexts"`.}
#'   \item{`net_markov_stability` ([hg_markov_stability()])}{`"states"`
#'     (default), `"passage_time"`, `"stationary"`.}
#'   \item{`net_hon_boot` ([hg_bootstrap()])}{`"edges"`.}
#'   \item{`net_hon_compare` ([hg_compare()])}{`"edges"`.}
#'   \item{`net_hon_communities` ([hg_communities()] on a memory
#'     network)}{`"states"` (default), `"physical"`, `"modules"`,
#'     `"codelength"`, `"trials"`.}
#'   \item{`simplicial_complex` ([simplicial()])}{`"simplices"` (default),
#'     `"f_vector"`, `"degree"`.}
#'   \item{`q_analysis` ([hg_qanalysis()])}{`"q_levels"` (default),
#'     `"nodes"`.}
#'   \item{`persistent_homology` ([hg_homology()])}{`"persistence"`
#'     (default), `"betti"`.}
#'   \item{`persistence_landscape` ([hg_landscape()])}{`"landscape"`.}
#'   \item{`hg_communities` ([hg_communities()] on a hypergraph)}{see
#'     [hg_get.hg_communities()].}
#'   \item{`hypernets_community_comparison`, `hypernets_motifs`,
#'     `net_hg_mmsbm`}{see their methods.}
#' }
#'
#' @param x A hypernets object.
#' @param what `NULL` (default) for the object's primary table, or the name
#'   of one of its tables (see *Tables by class*).
#' @param ... Arguments of the method for `class(x)`: filters (`order_min`,
#'   `min_count`, `dim`, `k`, `dimension`, `significant`, ...), `sort_by`,
#'   and `top` (the first `top` rows, applied last).
#' @return A base `data.frame`, one row per observation of the selected
#'   table. An object of a class without a method raises
#'   `hypernets_bad_input`, naming the class.
#' @examples
#' hg <- window_hypergraph(list(s1 = c("a", "b", "a", "c")), window = 2L)
#' hg_get(hg)
#' hg_get(hg, what = "nodes")
#'
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' hg_get(hon(seqs, max_order = 2), what = "nodes")
#' @export
hg_get <- function(x, what = NULL, ...) {
  UseMethod("hg_get")
}

#' @rdname hg_get
#' @export
hg_get.default <- function(x, what = NULL, ...) {
  stop(errorCondition(
    sprintf(paste0("hg_get() has no table for an object of class %s. ",
                   "It reads the objects hypernets' constructors and verbs ",
                   "return."),
            paste(sprintf("`%s`", class(x)), collapse = "/")),
    class = "hypernets_bad_input", call = NULL))
}

#' Drop a result table's subclasses, keeping a plain data.frame
#'
#' Several verbs return a data.frame with a subclass (for its print or plot
#' method). Internally, when the plain table is needed, the subclasses in
#' front of `"data.frame"` are dropped -- what base `hg_get()` does
#' for a data.frame, spelled out so no package object is ever read through
#' `hg_get()`.
#'
#' @param x A data.frame, possibly subclassed.
#' @return `x` with class `"data.frame"`.
#' @noRd
.ho_plain <- function(x) {
  cl <- oldClass(x)
  i <- match("data.frame", cl)
  if (!is.na(i) && i > 1L) class(x) <- cl[-seq_len(i - 1L)]
  x
}

# ---------------------------------------------------------------------------
# Results of the imported estimators
#
# hon(), honem(), mogen(), markov_order(), memory(), hg_markov_stability(),
# hg_homology(), hg_landscape() and hg_qanalysis() return the estimator's
# object with one class in front, `hypernets_result`, so that it prints as
# every hypernets result does: a header line and its default table. Every
# other method (summary, plot, the estimator's own verbs) is the
# estimator's, reached through the next class.
# ---------------------------------------------------------------------------

#' Mark an estimator's result as a hypernets result
#' @param x The estimator's object.
#' @return `x` with class `hypernets_result` in front.
#' @noRd
.ho_result <- function(x) {
  if (inherits(x, "hypernets_result")) return(x)
  structure(x, class = c("hypernets_result", class(x)))
}

#' The estimator's object without the hypernets class
#' @noRd
.ho_unresult <- function(x) {
  class(x) <- setdiff(class(x), "hypernets_result")
  x
}

#' One line naming a result and its main settings
#' @noRd
.ho_result_header <- function(x) {
  cls <- setdiff(class(x), "hypernets_result")[1L]
  switch(cls,
    net_hon = sprintf(paste0(
      "Higher-order network: %d states, %d nodes, %d rules ",
      "(highest order %d, min_freq %d, %d sequences)"),
      length(x$first_order_states), x$n_nodes, nrow(x$ho_edges),
      x$max_order_observed, x$min_freq, x$n_trajectories),
    net_honem = sprintf("Higher-order network embedding: %d nodes in %d dimensions",
                        x$n_nodes, x$dim),
    net_mogen = sprintf(paste0(
      "Multi-order generative model: order %d selected by %s ",
      "(orders %d to %d, %d paths)"),
      x$optimal_order, toupper(x$criterion), min(x$orders), max(x$orders),
      x$n_paths),
    net_markov_order = sprintf(paste0(
      "Markov order test: order %d selected (%d sequences, ",
      "%d permutations)"),
      x$optimal_order, x$n_sequences, x$n_perm),
    net_markov_order_group = sprintf("Markov order tests by group: %d groups",
                                     length(x)),
    net_path_dependence = sprintf(paste0(
      "Memory of order %d: %d contexts, %d change the most likely next ",
      "state"),
      x$order, nrow(x$contexts), sum(x$contexts$flips)),
    net_markov_stability = sprintf("Markov stability of %d states",
                                   nrow(x$stability)),
    persistent_homology = sprintf(
      "Persistent homology: %d features in dimension%s %s",
      nrow(x$persistence),
      if (length(unique(x$persistence$dimension)) == 1L) "" else "s",
      paste(sort(unique(x$persistence$dimension)), collapse = ", ")),
    persistence_landscape = sprintf(
      "Persistence landscape of dimension %d: %d landscape%s on %d points",
      x$dimension, x$k_max, if (x$k_max == 1L) "" else "s",
      length(x$t_grid)),
    q_analysis = sprintf("Q-analysis: highest q = %d", x$max_q),
    sprintf("A %s result", cls)
  )
}

#' @rdname print.hypernets_result
#' @export
hg_get.hypernets_result <- function(x, ...) NextMethod()

#' Print a result of an imported estimator
#'
#' Prints one line naming the result and its main settings, then the first
#' rows of its default table, the table [hg_get()] returns.
#'
#' @param x A result of [hon()], [honem()], [mogen()], [markov_order()],
#'   [memory()], [hg_markov_stability()], [hg_homology()], [hg_landscape()]
#'   or [hg_qanalysis()].
#' @param n Number of rows of the table to show. Default `10`.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' print(hon(seqs, max_order = 2), n = 3)
#' @export
print.hypernets_result <- function(x, n = 10L, ...) {
  cat(.ho_result_header(x), "\n", sep = "")
  .ho_print_table(x, n)
  invisible(x)
}

# ---------------------------------------------------------------------------
# summary(): every table of a result as a list of data frames
#
# summary(x) returns a `hypernets_summary`, a list whose elements are the
# data frames hg_get(x, what = ) returns, named after `what`, plus the
# overall figures of a result as further named data frames. The result
# itself is left untouched, so every field the estimator's own methods read
# keeps its shape; the summary is where `$` gives a clear table.
# ---------------------------------------------------------------------------

#' The `what` values of the hg_get() method that reads `x`
#' @noRd
.ho_what_values <- function(x) {
  probe <- if (inherits(x, c("net_hon_group", "net_hon_boot_group",
                             "net_markov_order_group")) && length(x)) {
    x[[1L]]
  } else {
    x
  }
  for (cls in class(probe)) {
    method <- utils::getS3method("hg_get", cls, optional = TRUE)
    if (!is.null(method) && "what" %in% names(formals(method))) {
      what <- formals(method)$what
      return(unique(if (is.call(what)) eval(what) else what))
    }
  }
  NULL
}

#' Build a hypernets summary from a result's tables
#'
#' @param x A result with an hg_get() method.
#' @param extra Named list of further data frames (overall figures).
#' @return A `hypernets_summary`.
#' @noRd
.ho_summary <- function(x, extra = list()) {
  whats <- .ho_what_values(x)
  tables <- if (is.null(whats)) {
    list(table = hg_get(x))
  } else {
    # a table the object does not hold (e.g. sentences of a bag-of-words
    # corpus) is refused with hypernets_bad_input and left out
    stats::setNames(lapply(whats, \(w) {
      tryCatch(hg_get(x, what = w),
               hypernets_bad_input = function(e) NULL)
    }), whats)
  }
  tables <- c(Filter(Negate(is.null), tables), extra)
  structure(tables, class = "hypernets_summary",
            header = utils::capture.output(print(x, n = 1L))[1L])
}

#' @rdname result-summary
#' @param x A `hypernets_summary`.
#' @param what For [hg_get()] on a summary: the name of one of its tables.
#'   Default: the first.
#' @export
hg_get.hypernets_summary <- function(x, what = names(x)[1L], ...) {
  what <- match.arg(what, names(x))
  x[[what]]
}

#' @param n For `print()`: number of rows of each table to show. Default
#'   `5`.
#' @param ... Unused.
#' @rdname result-summary
#' @export
print.hypernets_summary <- function(x, n = 5L, ...) {
  cat(attr(x, "header"), "\n", sep = "")
  cat(sprintf("Tables: %s\n", paste(names(x), collapse = ", ")))
  lapply(names(x), \(name) {
    table <- x[[name]]
    cat(sprintf("\n%s (%d row%s)\n", name, nrow(table),
                if (nrow(table) == 1L) "" else "s"))
    if (nrow(table) > 0L) {
      print(utils::head(table, as.integer(n)), row.names = FALSE)
    }
  })
  invisible(x)
}

#' Summaries of results
#'
#' `summary()` on a result returns every table of the result as a list of
#' data frames, named after the `what` values of [hg_get()], so
#' `summary(x)$validation` is the same data frame as
#' `hg_get(x, what = "validation")`. Results with overall figures add them
#' as further one-row or per-group tables. Tables that a particular result
#' does not hold are left out. The result itself is not changed.
#'
#' @name result-summary
#' @param object A result: from [hon()], [honem()], [mogen()],
#'   [markov_order()], [memory()], [hg_markov_stability()], [hypa()] on
#'   sequences, [hg_bootstrap()], [hg_compare()], [hg_communities()],
#'   [hg_compare_communities()], [hg_mmsbm()], [hg_motifs()], [simplicial()],
#'   [hg_homology()], [hg_landscape()] or [hg_qanalysis()].
#' @return A `hypernets_summary`: a named list of data frames. Its print
#'   method shows the first rows of each.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' fit_summary <- summary(memory(seqs, order = 2, min_count = 1))
#' fit_summary
#' fit_summary$contexts
NULL

#' @rdname result-summary
#' @export
summary.hypernets_result <- function(object, ...) {
  extra <- if (inherits(object, "net_path_dependence")) {
    chain <- object$chain
    list(overall = data.frame(
      n_contexts = as.integer(chain$n_contexts),
      n_flips = as.integer(chain$n_flips),
      kl_weighted = as.numeric(chain$KL_weighted),
      entropy_drop_weighted = as.numeric(chain$H_drop_weighted)))
  } else if (inherits(object, "net_markov_order")) {
    list(overall = data.frame(
      optimal_order = as.integer(object$optimal_order),
      aic_order = as.integer(object$aic_order),
      bic_order = as.integer(object$bic_order),
      n_sequences = as.integer(object$n_sequences),
      n_permutations = as.integer(object$n_perm)))
  } else {
    list()
  }
  .ho_summary(object, extra)
}

#' @rdname print.hypernets_result
#' @param y Unused.
#' @param what For a Markov stability result: `"states"` (default), one
#'   panel of bars per measure, or `"passage_time"`, the mean first-passage
#'   times as a heatmap (rows are the starting state).
#' @return `plot()` returns a ggplot for a persistence landscape (the
#'   landscapes that are not zero everywhere, one colour and line type
#'   each) and for a Markov stability result; every other result is plotted
#'   by its estimator's method.
#' @export
plot.hypernets_result <- function(x, y, what = c("states", "passage_time"),
                                  ...) {
  if (inherits(x, "net_markov_stability")) {
    return(.hms_plot(x, match.arg(what), ...))
  }
  if (!inherits(x, "persistence_landscape")) return(NextMethod())
  curves <- hg_get(x)
  drawn <- unique(curves$k[curves$value > 0])
  curves <- curves[curves$k %in% drawn, , drop = FALSE]
  curves$landscape <- factor(curves$k, levels = drawn)
  ggplot2::ggplot(curves, ggplot2::aes(x = .data$t, y = .data$value,
                                       colour = .data$landscape,
                                       linetype = .data$landscape)) +
    ggplot2::geom_line(linewidth = 0.9) +
    ggplot2::scale_colour_manual(
      values = rep_len(.thg_okabe_ito, length(drawn))) +
    ggplot2::labs(x = "filtration value", y = "landscape value",
                  colour = "landscape", linetype = "landscape",
                  title = sprintf("Persistence landscapes, dimension %d",
                                  x$dimension)) +
    ggplot2::theme_minimal(base_size = 12)
}

# Markov stability: one panel of bars per measure, states in the order of
# their stationary share (the axis names them, so the bars share one
# colour); or the first-passage matrix as a heatmap through cograph.
.hms_plot <- function(x, what, ...) {
  if (identical(what, "passage_time")) {
    passage <- hg_get(x, what = "passage_time")
    states <- hg_get(x)$state
    steps <- matrix(NA_real_, length(states), length(states),
                    dimnames = list(states, states))
    steps[cbind(match(passage$from, states), match(passage$to, states))] <-
      passage$steps
    return(cograph::plot_heatmap(steps, colors = c("#FFFFFF", "#0072B2"),
                                 show_values = TRUE, value_digits = 1,
                                 legend_title = "steps", xlab = "to",
                                 ylab = "from", axis_text_angle = 45, ...))
  }
  table <- hg_get(x)
  measures <- c(persistence = "persistence", stationary = "stationary share",
                return_time = "return time", sojourn_time = "sojourn time",
                avg_time_to_others = "steps to the others",
                avg_time_from_others = "steps from the others")
  long <- do.call(rbind, lapply(names(measures), \(m) {
    data.frame(state = table$state, measure = measures[[m]],
               value = table[[m]], stringsAsFactors = FALSE)
  }))
  long$state <- factor(long$state,
                       levels = table$state[order(table$stationary)])
  long$measure <- factor(long$measure, levels = unname(measures))
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$state)) +
    ggplot2::geom_col(fill = "#0072B2", width = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(.data$measure), scales = "free_x",
                        ncol = 3) +
    ggplot2::labs(x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12)
}
