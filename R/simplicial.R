# ---- hg_get() for simplicial_complex and q_analysis ----
#
# simplicial() and hg_qanalysis() wrap the Nestimate estimators and return
# their objects unchanged; their tables are read here.

#' Tables of a simplicial complex
#'
#' @param x A `simplicial_complex` object from
#'   [simplicial()].
#' @param ... Additional arguments (ignored).
#' @param what `"simplices"` (default) for one row per simplex,
#'   `"f_vector"` for the face counts by dimension, `"betti"` for the Betti
#'   numbers by dimension, `"degree"` for the
#'   simplicial degree of every node (the table of [hg_degree()]), or
#'   `"validation"` for the tests of a window complex built with
#'   `validate = TRUE`.
#' @param dimension Integer or `NULL`. For `what = "simplices"`: keep only
#'   simplices of this dimension (`0` vertices, `1` edges, `2` triangles).
#'   For `what = "validation"`: keep only the sets of `dimension + 1` actions.
#' @param normalized Logical. Only for `what = "degree"`: divide each count
#'   by its largest possible value. Default `FALSE`.
#' @return A data.frame. For `what = "simplices"`, one row per simplex:
#'   `simplex` (its number), `dimension`, `size`, `members` (the node names,
#'   comma separated). For `what = "f_vector"`, one row per dimension:
#'   `dimension`, `count`. For `what = "betti"`, one row per dimension:
#'   `dimension`, `betti` (the number of connected parts, loops, voids, ...).
#'   For `what = "degree"`, one row per node, most
#'   connected first: `node`, `d0`, `d1`, ... and `total`. For
#'   `what = "validation"`, one row per tested set of actions, by size, then
#'   p-value, then `z`: `members` (the actions, comma separated), `size`,
#'   `count` (windows with exactly these actions), `n_windows` (windows with
#'   this many distinct actions), `expected` (the mean count under the null
#'   model), `z` (the count minus `expected`, in null standard deviations),
#'   `p_value`, `p_adj` (Benjamini-Hochberg over `n_tests` possible sets),
#'   `n_tests`, `significant`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' sc <- simplicial(adj, type = "clique")
#' hg_get(sc)
#' hg_get(sc, what = "f_vector")
#' @export
hg_get.simplicial_complex <- function(x, what = c("simplices", "f_vector",
                                                  "betti", "degree",
                                                  "validation"),
                                      ..., dimension = NULL,
                                      normalized = FALSE,
                                      top = NULL) {
  what <- match.arg(what)
  if (what == "betti") {
    betti <- as.integer(hg_betti(x))
    return(.ho_top(data.frame(dimension = seq_along(betti) - 1L,
                              betti = betti), top))
  }
  if (what == "validation") {
    if (is.null(x$validation)) {
      .thg_bad_input(paste0("this complex has no validation; build it with ",
                            "simplicial(type = \"window\", validate = TRUE)"))
    }
    out <- x$validation
    if (!is.null(dimension)) {
      stopifnot("`dimension` must be a single integer >= 1" =
                  is.numeric(dimension) && length(dimension) == 1L &&
                  dimension >= 1)
      out <- out[out$size == as.integer(dimension) + 1L, , drop = FALSE]
      rownames(out) <- NULL
    }
    return(.ho_top(out, top))
  }
  if (what == "degree") {
    return(.ho_top(.sc_degree_table(x, normalized = normalized), top))
  }
  if (what == "f_vector") {
    fv <- x$f_vector
    return(.ho_top(data.frame(dimension = seq_along(fv) - 1L,
                              count = as.integer(fv),
                              stringsAsFactors = FALSE), top))
  }
  simplices <- x$simplices
  sizes <- vapply(simplices, length, integer(1L))
  out <- data.frame(
    simplex   = seq_along(simplices),
    dimension = sizes - 1L,
    size      = sizes,
    members = vapply(simplices,
                     function(s) paste(x$nodes[s], collapse = ", "),
                     character(1L)),
    stringsAsFactors = FALSE
  )
  if (!is.null(dimension)) {
    stopifnot("`dimension` must be a single integer >= 0" =
                is.numeric(dimension) && length(dimension) == 1L &&
                dimension >= 0)
    out <- out[out$dimension == as.integer(dimension), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' Tables of a q-analysis
#'
#' @param x A `q_analysis` object from [hg_qanalysis()].
#' @param ... Additional arguments (ignored).
#' @param what `"q_levels"` (default) for the number of q-connected
#'   components at each level, or `"nodes"` for each node's highest
#'   q-level.
#' @return A data.frame. For `what = "q_levels"`, one row per level: `q`,
#'   `components`. For `what = "nodes"`, one row per node: `node`, `max_q`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' qa <- hg_qanalysis(simplicial(adj, type = "clique"))
#' hg_get(qa)
#' hg_get(qa, what = "nodes")
#' @export
hg_get.q_analysis <- function(x, what = c("q_levels", "nodes"), ...,
                              top = NULL) {
  what <- match.arg(what)
  if (what == "nodes") {
    sv <- x$structure_vector
    return(.ho_top(data.frame(node = names(sv), max_q = as.integer(sv),
                              stringsAsFactors = FALSE), top))
  }
  qv <- x$q_vector
  # Levels are named "q_<k>"; unnamed or otherwise named vectors fall back
  # to the descending order hg_qanalysis() builds them in.
  labelled <- !is.null(names(qv)) && all(grepl("^q_[0-9]+$", names(qv)))
  q <- if (labelled) {
    as.integer(sub("^q_", "", names(qv)))
  } else {
    rev(seq_along(qv)) - 1L
  }
  out <- data.frame(q = q, components = as.integer(qv),
                    stringsAsFactors = FALSE)
  rownames(out) <- NULL
  .ho_top(out, top)
}

