# ---- Group models of memory networks ----
#
# hon(group =) splits the sequences by a grouping and estimates one memory
# network per group. The group model keeps each group's data and the
# estimation settings, so hg_compare() and hg_bootstrap() work from it
# without the caller splitting the data.

# Split sequence data by a grouping: a column name of a data frame, or one
# label per sequence (per row of a data frame, per element of a list).
.hon_split_groups <- function(data, group) {
  if (is.data.frame(data) && is.character(group) && length(group) == 1L &&
      group %in% names(data)) {
    labels <- data[[group]]
    parts <- split(data, as.character(labels))
    # a wide data frame would read the group column as a state
    parts <- lapply(parts, \(d) d[setdiff(names(d), group)])
  } else {
    n_seq <- if (is.data.frame(data)) nrow(data) else length(data)
    if (length(group) != n_seq) {
      .thg_bad_input(paste0(
        "`group` must name a column of `data` or give one label per ",
        "sequence (", n_seq, " here)"))
    }
    labels <- group
    parts <- split(data, as.character(labels))
  }
  if (anyNA(labels)) .thg_bad_input("`group` has missing labels")
  if (length(parts) < 2L) .thg_bad_input("`group` defines fewer than 2 groups")
  parts
}

.hon_group <- function(data, group, max_order, min_freq, collapse_repeats,
                       method, action, actor, time, session, time_threshold,
                       timezone) {
  parts <- .hon_split_groups(data, group)
  args <- list(max_order = max_order, min_freq = min_freq,
               collapse_repeats = collapse_repeats, method = method,
               action = action, actor = actor, time = time,
               session = session, time_threshold = time_threshold,
               timezone = timezone)
  fits <- lapply(parts, \(d) do.call(hon, c(list(d), args)))
  structure(fits, class = c("net_hon_group", "list"), data = parts,
            args = args,
            group = if (is.character(group) && length(group) == 1L) group)
}

#' @rdname hg_get.net_hon_group
#' @param x A `net_hon_group` from [hon()] with `group`.
#' @param n Number of rows of the default table to print. Default `10`.
#' @param ... Unused.
#' @return `print()` returns `x` invisibly.
#' @export
print.net_hon_group <- function(x, n = 10L, ...) {
  cat(sprintf("Memory networks by %s: %d groups (%s)\n",
              attr(x, "group") %||% "group", length(x),
              paste(sprintf("%s: %d sequences", names(x),
                            vapply(x, \(fit) as.integer(fit$n_trajectories),
                                   integer(1L))),
                    collapse = ", ")))
  .ho_print_table(x, n)
  invisible(x)
}

#' Tables of a group model of memory networks
#'
#' Returns the requested table of every group's memory network, stacked,
#' with a `group` column first.
#'
#' @param x A `net_hon_group` from [hon()] with `group`.
#' @param ... Passed to [hg_get()] on each group's network (`what`,
#'   `order_min`, `sort_by`, `top`, ...).
#' @return A data.frame: `group`, then the columns of the requested table.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b"), c("a", "b", "c", "a", "b"),
#'              c("x", "b", "d", "x", "b"), c("x", "b", "d", "x", "b"))
#' by_group <- hon(seqs, group = c("g1", "g1", "g2", "g2"), max_order = 2)
#' by_group
#' hg_get(by_group, top = 3)
#' @export
hg_get.net_hon_group <- function(x, ...) {
  tables <- lapply(names(x), \(g) {
    t <- hg_get(x[[g]], ...)
    data.frame(group = rep(g, nrow(t)), t, stringsAsFactors = FALSE,
               check.names = FALSE)
  })
  out <- do.call(rbind, tables)
  rownames(out) <- NULL
  out
}

#' Permutation comparison of the rules of two groups
#'
#' Tests whether two groups of sequences differ in their higher-order rule
#' probabilities. The rule set is extracted from the pooled data of the two
#' groups; for every pooled rule the statistic is the absolute difference of
#' the two groups' conditional probabilities, and its null distribution
#' comes from permuting the group labels over sequences. Per-rule p-values
#' are Benjamini-Hochberg adjusted; a global test aggregates the rule
#' differences weighted by pooled counts.
#'
#' Per-sequence counts are precomputed once and every permutation is a
#' weighted aggregation; permutations are drawn before any parallel work, so
#' `parallel = TRUE` reproduces the serial result under the same `seed`.
#' Rules whose context is unobserved in a group under some permutation
#' contribute only their valid permutations (`n_perm_used`).
#'
#' @param x A group model from [hon()] with `group`. The estimation settings
#'   (`max_order`, `min_freq`, `collapse_repeats`) and the data of each group
#'   are taken from it.
#' @param groups The two groups to compare, as a character vector of length
#'   2. May be omitted when the model has exactly two groups.
#' @param n_perm Integer >= 2. Label permutations. Default `1000`.
#' @param alpha Significance level for the `significant` flag on the
#'   adjusted p-values. Default `0.05`.
#' @param parallel,n_cores,seed As in [hg_bootstrap()].
#' @return An object of class `net_hon_compare`: a list with `edges` (one
#'   row per pooled rule: `from`, `to`, `order`, `count`, the two groups'
#'   counts and probabilities (columns named after the groups), `diff`
#'   (probability difference, first minus second), `p_value`, `p_adj` (BH),
#'   `significant`, `n_perm_used`), `global` (`statistic`, the
#'   pooled-count-weighted mean absolute difference, and `p_value`),
#'   `names`, `n_perm`, `alpha`, `max_order`, `min_freq`, `n_trajectories`
#'   (per group) and `seed`. Has `print`, `summary` and `plot` methods;
#'   [hg_get()] returns the rule table (`significant = TRUE` restricts it).
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances} 2(5),
#' e1600028. \doi{10.1126/sciadv.1600028}
#'
#' Good, P. (2005). \emph{Permutation, Parametric and Bootstrap Tests of
#' Hypotheses} (3rd ed.). Springer.
#'
#' Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
#' rate. \emph{Journal of the Royal Statistical Society B}, 57(1), 289-300.
#' \doi{10.1111/j.2517-6161.1995.tb02031.x}
#' @examples
#' set.seed(1)
#' first_order <- replicate(6, sample(c("a", "b", "c"), 12, replace = TRUE),
#'                          simplify = FALSE)
#' second_order <- replicate(6, rep(c("a", "b", "c", "b"), 3),
#'                           simplify = FALSE)
#' by_group <- hon(c(first_order, second_order),
#'                 group = rep(c("first", "second"), each = 6), max_order = 2)
#' comparison <- hg_compare(by_group, n_perm = 99, seed = 1)
#' comparison
#' hg_get(comparison, sort_by = "p_adj", top = 6)
#' @seealso [hon()], [hg_bootstrap()], [markov_order()]
#' @export
hg_compare <- function(x, groups = NULL, n_perm = 1000L, alpha = 0.05,
                       parallel = FALSE, n_cores = 2L, seed = NULL) {
  if (!inherits(x, "net_hon_group")) {
    .thg_bad_input(paste0("`x` must be a group model: ",
                          "hon(data, ..., group = <column or labels>)"))
  }
  available <- names(x)
  if (is.null(groups)) {
    if (length(available) != 2L) {
      .thg_bad_input(sprintf(
        "the model has %d groups (%s); choose two with `groups =`",
        length(available), paste(available, collapse = ", ")))
    }
    groups <- available
  }
  groups <- as.character(groups)
  if (length(groups) != 2L || anyDuplicated(groups) ||
      !all(groups %in% available)) {
    .thg_bad_input(sprintf("`groups` must be two of: %s",
                           paste(available, collapse = ", ")))
  }
  parts <- attr(x, "data")
  args <- attr(x, "args")
  .hg_compare_pair(parts[[groups[1L]]], parts[[groups[2L]]],
                   n_perm = n_perm, alpha = alpha,
                   max_order = args$max_order, min_freq = args$min_freq,
                   collapse_repeats = args$collapse_repeats,
                   action = args$action, actor = args$actor,
                   time = args$time, session = args$session,
                   time_threshold = args$time_threshold,
                   timezone = args$timezone,
                   names = groups, parallel = parallel, n_cores = n_cores,
                   seed = seed)
}

#' @rdname hg_get.net_hon_group
#' @export
print.net_hon_boot_group <- function(x, n = 10L, ...) {
  cat(sprintf("HON bootstrap by group: %d groups, %d replicates each\n",
              length(x), x[[1L]]$n_boot))
  .ho_print_table(x, n)
  invisible(x)
}

#' @rdname result-summary
#' @export
summary.net_hon_boot_group <- function(object, ...) .ho_summary(object)

#' @rdname result-summary
#' @export
summary.net_hon_group <- function(object, ...) .ho_summary(object)

#' @rdname hg_get.net_hon_group
#' @export
hg_get.net_hon_boot_group <- function(x, ...) hg_get.net_hon_group(x, ...)
