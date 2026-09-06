# ---- Windowed sequence hyperedges ----
#
# Sequence-to-hypergraph construction after Ding et al. (2020): a fixed-size
# window slides (or strides) over each sequence, and the distinct states
# inside one window form one hyperedge. Windows with identical state sets
# collapse into a single hyperedge whose weight is the window count; the
# incidence cells hold within-window occurrence totals, i.e. edge-dependent
# vertex weights (Chitra & Raphael 2019) that hypergraph_cluster(type =
# "random_walk") consumes directly.

# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------

#' Parse window_hypergraph input into a list of character trajectories
#'
#' Long format (action given): one trajectory per actor, ordered by time
#' (row order within actor when time is NULL). Otherwise wide data.frame /
#' character matrix (one trajectory per row, trailing NAs stripped) or a
#' list of vectors.
#'
#' @param data data.frame, character matrix, or list.
#' @param action,actor,time Long-format column names or NULL.
#' @return Named list of character vectors.
#' @noRd
.wh_parse_input <- function(data, action, actor, time) {
  if (!is.null(action)) {
    stopifnot(
      "`data` must be a data.frame when `action` is given" = is.data.frame(data),
      "`action` must name a column of `data`" =
        is.character(action) && length(action) == 1L && action %in% names(data),
      "`actor` must be NULL or name a column of `data`" =
        is.null(actor) ||
        (is.character(actor) && length(actor) == 1L && actor %in% names(data)),
      "`time` must be NULL or name a column of `data`" =
        is.null(time) ||
        (is.character(time) && length(time) == 1L && time %in% names(data))
    )
    a <- as.character(data[[action]])
    g <- if (is.null(actor)) rep("sequence_1", length(a)) else
      as.character(data[[actor]])
    keep <- !is.na(g)
    a <- a[keep]
    g <- g[keep]
    # order() is stable, so with time = NULL the row order within each
    # actor is preserved.
    o <- if (is.null(time)) order(g) else order(g, data[[time]][keep])
    return(split(a[o], g[o]))
  }
  stopifnot(
    "`actor` requires `action`" = is.null(actor),
    "`time` requires `action`"  = is.null(time)
  )
  if (is.matrix(data) && !is.numeric(data)) {
    data <- as.data.frame(data, stringsAsFactors = FALSE)
  }
  if (is.data.frame(data)) {
    stopifnot(
      "`data` must have at least one row"    = nrow(data) >= 1L,
      "`data` must have at least one column" = ncol(data) >= 1L
    )
    trajectories <- lapply(seq_len(nrow(data)), function(i) {
      row_vals <- as.character(unlist(data[i, ], use.names = FALSE))
      non_na <- which(!is.na(row_vals))
      if (length(non_na) == 0L) return(character(0L))
      row_vals[seq_len(max(non_na))]  # strip trailing NAs, keep internal
    })
    names(trajectories) <- paste0("sequence_", seq_along(trajectories))
    return(trajectories)
  }
  if (is.list(data)) {
    trajectories <- lapply(data, as.character)
    if (is.null(names(trajectories))) {
      names(trajectories) <- paste0("sequence_", seq_along(trajectories))
    }
    return(trajectories)
  }
  stop("`data` must be a data.frame, character matrix, or list of vectors.",
       call. = FALSE)
}

#' Enumerate the full windows of one trajectory
#'
#' @param traj Character vector.
#' @param window,step Integers.
#' @return NULL when the trajectory is shorter than `window`, else a
#'   data.frame with one row per window position: `win` (window id within
#'   this trajectory) and `state`.
#' @noRd
.wh_windows_one <- function(traj, window, step) {
  n <- length(traj)
  if (n < window) return(NULL)
  starts <- seq.int(1L, n - window + 1L, by = step)
  idx <- outer(starts, 0L:(window - 1L), "+")  # n_windows x window positions
  data.frame(
    win   = as.vector(row(idx)),
    state = as.character(traj[idx]),
    stringsAsFactors = FALSE
  )
}

# ---------------------------------------------------------------------------
# Constructor
# ---------------------------------------------------------------------------

#' Windowed Sequence Hyperedges
#'
#' Builds a hypergraph from categorical sequence data: a fixed-size window
#' moves over each sequence and the distinct states observed inside one
#' window form one hyperedge (Ding et al. 2020). Windows whose distinct
#' state sets coincide collapse into a single hyperedge whose weight is the
#' number of such windows (`window_counts`); the incidence cells hold the
#' total within-window occurrences of each state, so the incidence matrix
#' carries edge-dependent vertex weights (Chitra & Raphael 2019) that
#' [hypergraph_cluster()] with `type = "random_walk"` uses directly. The
#' window counts are the default hyperedge weights of the whole Laplacian
#' family ([hypergraph_laplacian()], [hypergraph_cluster()],
#' [hypergraph_transduction()]).
#'
#' `step = 1` (default) gives sliding windows; `step = window` gives
#' tumbling (non-overlapping) windows. Only full windows are formed: a
#' sequence shorter than `window` contributes none and is counted in
#' `params$n_short_sequences`. `NA` states are excluded from a window's
#' set; a window containing only `NA`s is skipped and counted in
#' `params$n_empty_windows`.
#'
#' Whole-sequence hyperedges (each complete sequence as one hyperedge) are
#' the special case already covered by [group_hypergraph()] on long-format
#' data with `member = action` and `group = actor`; use this verb when the
#' hyperedges should be local in time.
#'
#' @param data Sequence data: a wide data.frame or character matrix (one
#'   sequence per row, trailing `NA`s stripped), a list of character
#'   vectors, or a long data.frame together with `action` (and optionally
#'   `actor`, `time`).
#' @param window Integer >= 2. Window size in sequence positions.
#' @param step Integer >= 1. Offset between consecutive window starts:
#'   `1` slides, `window` tumbles.
#' @param action Character or NULL. Long format only: column holding the
#'   categorical state of each event.
#' @param actor Character or NULL. Long format only: column grouping events
#'   into sequences (one sequence per actor). `NULL` treats all rows as one
#'   sequence.
#' @param time Character or NULL. Long format only: column ordering events
#'   within each actor. `NULL` keeps row order.
#' @param min_size Integer >= 1. Drop hyperedges with fewer distinct states
#'   after collapsing. The default `1` keeps everything.
#' @param min_weight Integer >= 1. Drop hyperedges observed in fewer than
#'   `min_weight` windows, keeping only recurrent state combinations (the
#'   role `min_freq` plays in rule extraction). The default `1` keeps
#'   everything. The total dropped by `min_size` and `min_weight` together
#'   is recorded in `params$n_dropped`; with both at their defaults,
#'   `sum(window_counts)` equals the number of non-empty full windows.
#'
#' @return A `net_hypergraph` object (as from [build_hypergraph()] and
#'   [group_hypergraph()]): a list with `hyperedges` (list of sorted node
#'   index vectors), `incidence` (numeric node x hyperedge matrix of
#'   within-window occurrence totals), `nodes`, `n_nodes`, `n_hyperedges`,
#'   `window_counts` (integer, one weight per hyperedge: the number of
#'   windows collapsed into it), `size_distribution`, and `params`
#'   (`source = "window_hypergraph"`, `window`, `step`, `min_size`,
#'   `n_sequences`, `n_short_sequences`, `n_windows`, `n_empty_windows`,
#'   `n_dropped`). Use [as.data.frame()] for the tidy one-row-per-hyperedge
#'   table.
#'
#' @references
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
#' Hypergraph attention networks for inductive text classification.
#' \emph{Proceedings of EMNLP 2020}, 4927-4936.
#' \doi{10.18653/v1/2020.emnlp-main.399}
#'
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. \emph{Proceedings of the 36th
#' International Conference on Machine Learning}, PMLR 97, 1172-1181.
#'
#' @examples
#' hg <- window_hypergraph(human_long, action = "code",
#'                         actor = "session_id", time = "timestamp",
#'                         window = 3L)
#' hg
#' edges <- as.data.frame(hg)
#' head(edges)
#'
#' # Tumbling windows over wide-format sequences
#' wide <- data.frame(
#'   t1 = c("plan", "code"), t2 = c("code", "test"),
#'   t3 = c("test", "code"), t4 = c("plan", "debug")
#' )
#' window_hypergraph(wide, window = 2L, step = 2L)
#'
#' @seealso [build_hypergraph()], [group_hypergraph()],
#'   [hypergraph_measures()], [hypergraph_cluster()], [clique_expansion()]
#'
#' @export
window_hypergraph <- function(data, window = 3L, step = 1L,
                              action = NULL, actor = NULL, time = NULL,
                              min_size = 1L, min_weight = 1L) {
  stopifnot(
    "`window` must be a single integer >= 2" =
      is.numeric(window) && length(window) == 1L && is.finite(window) &&
      window >= 2 && window == as.integer(window),
    "`step` must be a single integer >= 1" =
      is.numeric(step) && length(step) == 1L && is.finite(step) &&
      step >= 1 && step == as.integer(step),
    "`min_size` must be a single integer >= 1" =
      is.numeric(min_size) && length(min_size) == 1L && is.finite(min_size) &&
      min_size >= 1 && min_size == as.integer(min_size),
    "`min_weight` must be a single integer >= 1" =
      is.numeric(min_weight) && length(min_weight) == 1L &&
      is.finite(min_weight) && min_weight >= 1 &&
      min_weight == as.integer(min_weight)
  )
  window <- as.integer(window)
  step <- as.integer(step)
  min_size <- as.integer(min_size)
  min_weight <- as.integer(min_weight)

  trajectories <- .wh_parse_input(data, action, actor, time)
  n_sequences <- length(trajectories)

  per <- lapply(trajectories, .wh_windows_one, window = window, step = step)
  short <- vapply(per, is.null, logical(1L))
  n_short <- sum(short)
  per <- per[!short]
  if (length(per) == 0L) {
    stop("No windows: every sequence is shorter than `window` (",
         window, ").", call. = FALSE)
  }

  # Stack per-trajectory window tables with globally unique window ids
  n_win_per <- vapply(per, function(d) max(d$win), integer(1L))
  offsets <- cumsum(c(0L, n_win_per[-length(n_win_per)]))
  long <- do.call(rbind, Map(function(d, o) {
    d$win <- d$win + o
    d
  }, per, offsets))
  n_windows <- sum(n_win_per)

  long <- long[!is.na(long$state), , drop = FALSE]
  if (nrow(long) == 0L) {
    stop("No windows: all window positions are NA.", call. = FALSE)
  }

  # One set key per window; windows that were all-NA drop out of split()
  key_by_win <- vapply(
    split(long$state, factor(long$win)),
    function(s) paste(sort(unique(s)), collapse = "\x01"),
    character(1L)
  )
  n_empty <- n_windows - length(key_by_win)

  states_levels <- sort(unique(long$state))
  keys_levels <- sort(unique(key_by_win))
  n_states <- length(states_levels)
  n_keys <- length(keys_levels)

  # Occurrence-weighted incidence: cell [s, h] = total occurrences of state
  # s across all windows collapsed into hyperedge h
  occ_key <- key_by_win[as.character(long$win)]
  si <- match(long$state, states_levels)
  ki <- match(occ_key, keys_levels)
  counts <- tabulate((ki - 1L) * n_states + si, nbins = n_states * n_keys)
  incidence <- matrix(as.numeric(counts), n_states, n_keys,
                      dimnames = list(states_levels,
                                      paste0("h", seq_len(n_keys))))
  window_counts <- as.integer(table(factor(key_by_win, levels = keys_levels)))
  hyperedges <- lapply(strsplit(keys_levels, "\x01", fixed = TRUE),
                       function(s) sort(match(s, states_levels)))

  # min_size / min_weight filters (recorded, never silent)
  he_sizes <- lengths(hyperedges)
  keep <- he_sizes >= min_size & window_counts >= min_weight
  n_dropped <- sum(!keep)
  if (!any(keep)) {
    stop("No hyperedges left: all ", n_keys, " hyperedges fall below ",
         "`min_size` (", min_size, ") or `min_weight` (", min_weight, ").",
         call. = FALSE)
  }
  incidence <- incidence[, keep, drop = FALSE]
  colnames(incidence) <- paste0("h", seq_len(ncol(incidence)))
  hyperedges <- hyperedges[keep]
  window_counts <- window_counts[keep]

  he_sizes <- lengths(hyperedges)
  tab <- table(he_sizes)
  size_dist <- as.integer(tab)
  names(size_dist) <- paste0("size_", names(tab))

  structure(
    list(
      hyperedges        = hyperedges,
      incidence         = incidence,
      nodes             = states_levels,
      n_nodes           = n_states,
      n_hyperedges      = length(hyperedges),
      window_counts     = window_counts,
      size_distribution = size_dist,
      params = list(
        source            = "window_hypergraph",
        window            = window,
        step              = step,
        min_size          = min_size,
        min_weight        = min_weight,
        n_sequences       = n_sequences,
        n_short_sequences = n_short,
        n_windows         = n_windows,
        n_empty_windows   = n_empty,
        n_dropped         = n_dropped
      )
    ),
    class = "net_hypergraph"
  )
}
