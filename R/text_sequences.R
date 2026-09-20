# ---------------------------------------------------------------------------
# hg_sequences(): the text family -> memory family bridge.
#
# A clustered corpus is already a sequence: each document carries the actor
# who produced it and the moment it was produced, and clustering gives it a
# state. This verb performs the join, the ordering and the reshaping that the
# caller would otherwise write by hand, and returns the canonical long
# sequence table (`actor` / `time` / `action`) that the memory family reads.
# ---------------------------------------------------------------------------

#' Per-actor state sequences from a clustered corpus
#'
#' Turns document-level topic assignments into the canonical long sequence
#' table the memory family consumes. Each document of a [text_hypergraph()]
#' carries whatever metadata came with the corpus (an author, a turn, a
#' timestamp); [hg_cluster()] gives each document a topic. `hg_sequences()`
#' joins the two, orders the documents within each actor and returns one row
#' per document as `actor` / `time` / `action` -- the long shape
#' [build_hon()], [bootstrap_hon()], [compare_hon()] and
#' [window_hypergraph()] read through their own `action` / `actor` / `time`
#' arguments, so no coercion is written at the call site:
#'
#' ```r
#' seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
#' build_hon(seqs, action = "action", actor = "actor", time = "time")
#' ```
#'
#' Those three column names are fixed, which is what makes the second call
#' the same every time. The memory-family verbs that do not (yet) take the
#' long form -- [build_mogen()], [build_hypa()], [markov_order_test()] --
#' read a list of trajectories instead; pass the result of `hg_sequences()`
#' to them only through one of the verbs above, since a long table handed to
#' them bare is read as a *wide* one (one trajectory per row) and silently
#' gives the wrong model.
#'
#' This closes the one missing edge in the cross-family design: `pathways()`
#' bridges the memory family to the others and [window_hypergraph()] bridges
#' sequences to hypergraphs, but until now the text family dead-ended at
#' clustering and the caller had to assemble the sequence table themselves.
#'
#' @param hg A [text_hypergraph()] -- more generally, any `net_hypergraph`
#'   carrying a documents table, the one reached with
#'   `as.data.frame(hg, what = "documents")`.
#' @param clusters A partition of the documents, as returned by
#'   [hg_cluster()]: a data.frame with one row per node, a `node` column
#'   holding document identifiers and a state column (`cluster` by default,
#'   or the column named by `state`). `NULL` (the default) takes the state
#'   from the documents table instead, in which case `state` is required.
#' @param actor Name of the documents-table column identifying who produced
#'   each document. One sequence is returned per distinct value.
#' @param order_by Name of the documents-table column ordering the documents
#'   within an actor -- a turn number, an index, a date or a timestamp.
#'   Anything `order()` can sort is accepted, and the values are carried
#'   through to the `time` column unchanged.
#' @param state Name of the column holding each document's state. With
#'   `clusters` supplied it names a column of `clusters` and defaults to
#'   `"cluster"`; with `clusters = NULL` it names a column of the documents
#'   table and has no default.
#'
#' @return A base `data.frame` with one row per document of `hg`, ordered by
#'   `actor` and then `time`, with columns
#'   \describe{
#'     \item{actor}{character. The value of the `actor` column; one sequence
#'       per distinct value.}
#'     \item{time}{The value of the `order_by` column, type unchanged, so a
#'       numeric turn stays numeric and a `Date` stays a `Date`.}
#'     \item{action}{character. The document's state -- its cluster label
#'       under `clusters`, or the `state` column of the documents table.}
#'   }
#'
#' @section Conditions:
#' A `hypernets_bad_input` error is raised when `hg` is not a hypergraph,
#' when it carries no documents table, when `actor`, `order_by` or `state`
#' does not name a column of the table it is looked up in, when `clusters`
#' lacks a `node` column or assigns a document more than one state, when
#' some document of `hg` has no row in `clusters`, or when the actor, time
#' or state values contain `NA` (which would silently shorten a sequence).
#'
#' @references
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances}, 2(5),
#' e1600028.
#'
#' Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
#' clustering, classification, and embedding. \emph{Advances in Neural
#' Information Processing Systems}, 19, 1601-1608.
#'
#' @seealso [hg_cluster()] for the partition, [build_hon()] and
#'   [bootstrap_hon()] for what to do with the sequences,
#'   [window_hypergraph()] to read them back into a hypergraph.
#'
#' @examples
#' posts <- data.frame(
#'   student = rep(c("ana", "bo"), each = 4),
#'   turn = rep(1:4, times = 2),
#'   phase = rep(c("early", "early", "late", "late"), times = 2),
#'   text = c(
#'     "simmer the soup with onions and carrots",
#'     "this soup recipe needs salt on a cold night",
#'     "the telescope revealed a distant galaxy and stars",
#'     "astronomers aimed the telescope at the stars all night",
#'     "a galaxy of stars seen through the telescope",
#'     "salt and onions make the soup",
#'     "the soup simmered all night with carrots",
#'     "astronomers watched the distant galaxy"
#'   ),
#'   stringsAsFactors = FALSE
#' )
#' hg <- text_hypergraph(posts, column = "text",
#'                       stop_words = c("the", "with", "and", "a", "this",
#'                                      "at", "on", "all", "of", "through"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_sequences(hg, topics, actor = "student", order_by = "turn")
#'
#' # A state that is already a column of the corpus needs no clustering.
#' hg_sequences(hg, actor = "student", order_by = "turn", state = "phase")
#'
#' # Straight into the memory family, no further coercion.
#' seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
#' build_hon(seqs, action = "action", actor = "actor", time = "time",
#'           max_order = 2L)
#' @export
hg_sequences <- function(hg, clusters = NULL, actor, order_by, state = NULL) {
  .thg_check_hg(hg)
  stopifnot(
    "`actor` must be a single column name" =
      is.character(actor) && length(actor) == 1L && !is.na(actor),
    "`order_by` must be a single column name" =
      is.character(order_by) && length(order_by) == 1L && !is.na(order_by),
    "`state` must be NULL or a single column name" =
      is.null(state) ||
      (is.character(state) && length(state) == 1L && !is.na(state)),
    "`clusters` must be NULL or a data.frame (the result of `hg_cluster()`)" =
      is.null(clusters) || is.data.frame(clusters)
  )

  docs <- .thg_documents(hg)
  .thg_require_columns(docs, c(actor, order_by), "the documents table")

  action <- if (is.null(clusters)) {
    if (is.null(state)) {
      .thg_bad_input(paste0(
        "`state` must name a column of the documents table when `clusters` ",
        "is NULL; supply `clusters` (from `hg_cluster()`) or `state`"
      ))
    }
    .thg_require_columns(docs, state, "the documents table")
    as.character(docs[[state]])
  } else {
    .thg_cluster_states(docs, clusters, state = state %||% "cluster")
  }

  actor_values <- as.character(docs[[actor]])
  time_values <- docs[[order_by]]
  if (!is.atomic(time_values)) {
    .thg_bad_input(sprintf(
      "`order_by` must name a sortable column; `%s` is a %s",
      order_by, class(time_values)[1L]
    ))
  }
  .thg_no_missing(actor_values, sprintf("`actor` column \"%s\"", actor))
  .thg_no_missing(time_values, sprintf("`order_by` column \"%s\"", order_by))
  .thg_no_missing(action, "the state column")

  # The document identifier is the deterministic tie-break, so two documents
  # by the same actor at the same time never depend on input order.
  tie <- as.character(docs[["doc"]] %||% seq_len(nrow(docs)))
  ord <- order(actor_values, time_values, tie)
  out <- data.frame(
    actor = actor_values[ord],
    time = time_values[ord],
    action = action[ord],
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  stopifnot(
    "internal: `hg_sequences()` must return one row per document" =
      nrow(out) == nrow(docs)
  )
  out
}

# The documents table of a hypergraph, or a classed error saying it has none.
.thg_documents <- function(hg) {
  docs <- hg[["text"]][["documents"]]
  if (!is.data.frame(docs) || nrow(docs) == 0L) {
    .thg_bad_input(paste0(
      "`hg` carries no documents table; `hg_sequences()` needs a ",
      "text_hypergraph() built from a data.frame whose columns identify the ",
      "actor and the order of each document"
    ))
  }
  docs
}

# Every name in `cols` must be a column of `df`; otherwise a classed error
# naming the missing column and what is on offer.
.thg_require_columns <- function(df, cols, what) {
  missing_cols <- setdiff(cols, names(df))
  if (length(missing_cols) > 0L) {
    .thg_bad_input(sprintf(
      "%s has no column %s; it has: %s",
      what, paste0("\"", missing_cols, "\"", collapse = ", "),
      paste(names(df), collapse = ", ")
    ))
  }
  invisible(df)
}

# `NA` in an actor, a time or a state would silently shorten or reorder a
# sequence, so it is refused rather than dropped.
.thg_no_missing <- function(x, what) {
  if (anyNA(x)) {
    .thg_bad_input(sprintf(
      "%s contains %d missing value(s); drop or impute those documents first",
      what, sum(is.na(x))
    ))
  }
  invisible(x)
}

# Align an `hg_cluster()` partition onto the documents table, by document id.
.thg_cluster_states <- function(docs, clusters, state) {
  .thg_require_columns(clusters, c("node", state), "`clusters`")
  .thg_require_columns(docs, "doc", "the documents table")
  nodes <- as.character(clusters[["node"]])
  if (anyDuplicated(nodes) > 0L) {
    .thg_bad_input(
      "`clusters` must have one row per node; `node` has duplicates"
    )
  }
  idx <- match(as.character(docs[["doc"]]), nodes)
  if (anyNA(idx)) {
    .thg_bad_input(sprintf(
      paste0("`clusters` covers %d of the %d documents; it must assign a ",
             "state to every document of `hg` (is it a partition of the ",
             "words rather than the documents?)"),
      sum(!is.na(idx)), nrow(docs)
    ))
  }
  as.character(clusters[[state]][idx])
}
