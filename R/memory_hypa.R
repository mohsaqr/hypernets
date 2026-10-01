# ---- hg_get() for net_hypa ----
#
# hypa() on sequences returns the path-anomaly object of the HYPA builder
# unchanged; its tables are read here.

#' Tables of a path-anomaly (HYPA) fit
#'
#' @param x A `net_hypa` object from [hypa()] on sequences.
#' @param what `"anomalies"` (default) for the anomalous paths selected by
#'   `type`, `"scores"` for every scored path, `"over"` or `"under"` for the
#'   significantly over- or under-represented paths with all their columns,
#'   or `"pathways"` for the anomalous paths written as pathways
#'   (`"a b -> c"`: the context states, then the next state), the form
#'   [simplicial()] turns into simplices.
#' @param ... Unused.
#' @param type For `"anomalies"` and `"pathways"`: `"all"`, `"over"` or
#'   `"under"`. Default: the `type` given to [hypa()], else `"all"`.
#' @param order_by `"sig"` (smallest p-value of the path's own direction
#'   first), `"ratio"` (largest deviation first: highest ratio for
#'   over-represented paths, lowest for under-represented ones, largest
#'   absolute log ratio when both directions are listed), `"freq"` (most
#'   observed first) or `"path"` (alphabetical). Default: the `order_by`
#'   given to [hypa()], else `"sig"`. For `"scores"`, `"sig"` sorts by
#'   `p_value`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `order_by`. Default `NULL`
#'   returns every row.
#' @return A data.frame. For `"anomalies"`, one row per anomalous path:
#'   `path`, `count`, `expected`, `ratio`, `p_adj` (the adjusted p-value
#'   of the path's own direction), `direction`, and `order` when several
#'   orders were scored. For `"scores"`, `"over"` and `"under"`, one row per
#'   path with every column of the fit: `path`, `from`, `to`, `count`,
#'   `expected`, `ratio`, the p-values (`p_value`, `p_under`, `p_over`,
#'   `p_adj_under`, `p_adj_over`), the `direction` (`"over"`, `"under"` or
#'   `"normal"`) and the path `order`. For `"pathways"`, one row per anomalous path with the
#'   column `pathway`.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' fit <- hypa(seqs, order = 2, min_count = 1)
#' hg_get(fit, type = "over", order_by = "ratio", top = 5)
#' @export
hg_get.net_hypa <- function(x, what = c("anomalies", "scores", "over",
                                        "under", "pathways"), ...,
                            type = NULL, order_by = NULL, top = NULL) {
  what <- match.arg(what)
  view <- attr(x, "hypa_view")
  type <- match.arg(type %||% view$type %||% "all",
                    c("all", "over", "under"))
  order_by <- match.arg(order_by %||% view$order_by %||% "sig",
                        c("sig", "ratio", "freq", "path"))
  if (what == "pathways") {
    return(.ho_top(.ho_pathway_table(Nestimate::pathways(x, type = type)),
                   top))
  }
  if (what %in% c("over", "under")) {
    out <- .hypa_rename(.hypa_sort(x[[what]], what, order_by))
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  if (what == "scores") {
    out <- x$scores
    key <- switch(order_by,
                  sig = order(out$p_value, out$path),
                  ratio = order(-abs(log(out$ratio)), out$path),
                  freq = order(-out$observed, out$path),
                  path = order(out$path))
    out <- .hypa_rename(out[key, , drop = FALSE])
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  directions <- switch(type, all = c("over", "under"), over = "over",
                       under = "under")
  parts <- lapply(directions, \(dir) {
    d <- .hypa_sort(x[[dir]], dir, order_by)
    if (is.null(d) || !nrow(d)) return(NULL)
    p_adj <- if (dir == "over") d$p_adjusted_over else d$p_adjusted_under
    cols <- intersect(c("order", "path", "observed", "expected", "ratio"),
                      names(d))
    cbind(.hypa_rename(d[cols]), p_adj = p_adj, direction = dir,
          stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, Filter(Negate(is.null), parts))
  if (is.null(out)) {
    out <- data.frame(path = character(0L), count = numeric(0L),
                      expected = numeric(0L), ratio = numeric(0L),
                      p_adj = numeric(0L), direction = character(0L),
                      stringsAsFactors = FALSE)
  } else if (type == "all") {
    key <- switch(order_by,
                  sig = order(out$p_adj, out$path),
                  ratio = order(-abs(log(out$ratio)), out$path),
                  freq = order(-out$count, out$path),
                  path = order(out$path))
    out <- out[key, , drop = FALSE]
  }
  if (length(unique(out$order)) <= 1L) out$order <- NULL
  rownames(out) <- NULL
  .ho_top(out, top)
}

# The estimator's column names in the package vocabulary.
.hypa_rename <- function(d) {
  if (is.null(d)) return(d)
  .ho_rename(d, c(observed = "count", anomaly = "direction",
                  p_adjusted_under = "p_adj_under",
                  p_adjusted_over = "p_adj_over"))
}

# Sort one direction's table as the HYPA summary does: significance and
# ratio are read in that direction's own tail.
.hypa_sort <- function(d, dir, order_by) {
  if (is.null(d) || !nrow(d)) return(d)
  key <- switch(order_by,
                sig = if (dir == "over") order(d$p_over, d$path) else
                  order(d$p_under, d$path),
                ratio = if (dir == "over") order(-d$ratio, d$path) else
                  order(d$ratio, d$path),
                freq = order(-d$observed, d$path),
                path = order(d$path))
  d[key, , drop = FALSE]
}

#' Wrap a character vector of pathways as a one-column table
#'
#' @param pathways Character vector, as the pathway extractors return it.
#' @return A data.frame with one row per pathway and the column `pathway`.
#' @noRd
.ho_pathway_table <- function(pathways) {
  data.frame(pathway = as.character(pathways), stringsAsFactors = FALSE)
}
