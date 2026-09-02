# Walk-based hyperedge centralities from the paper. honets constructs the
# s-line graph; cograph supplies the ordinary graph centrality kernel.

#' s-betweenness and s-closeness centrality of hyperedges
#'
#' Hyperedges become vertices in the unweighted s-line graph: two are
#' adjacent when they share at least `s` nodes. Betweenness and closeness are
#' then ordinary shortest-path vertex centralities on that graph. At `s = 1`
#' these are the paper's hyperedge 1-betweenness and 1-closeness.
#'
#' @param hg A static `net_hypergraph` or a [temporal_hypergraph()].
#' @param s One or more positive integer intersection thresholds.
#' @param measure Any of `"betweenness"` and `"closeness"`.
#' @param normalized Normalize the graph centralities where the cograph
#'   measure supports normalization.
#' @param top Optional number of highest-scoring hyperedges to retain per
#'   `(time, s, measure)` group.
#' @param at,snapshot_mode,multiedges Temporal snapshot arguments passed to
#'   [hypergraph_snapshots()] when `hg` is temporal.
#' @return A tidy data frame with `edge`, `s`, `measure`, and `value`; temporal
#'   input adds a leading `time` column.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#' hypergraphs. *Philosophical Transactions of the Royal Society A*,
#' 382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' h <- group_hypergraph(
#'   data.frame(member = c("a", "b", "c", "b", "c", "d"),
#'              event = rep(c("e1", "e2"), each = 3)),
#'   "member", "event"
#' )
#' hg_edge_centrality(h, s = 1)
#' @export
hg_edge_centrality <- function(hg, s = 1L,
                               measure = c("betweenness", "closeness"),
                               normalized = TRUE, top = NULL, at = NULL,
                               snapshot_mode = c("active", "cumulative", "all"),
                               multiedges = TRUE) {
  measure <- match.arg(measure, several.ok = TRUE)
  snapshot_mode <- match.arg(snapshot_mode)
  if (!is.numeric(s) || length(s) < 1L || any(!is.finite(s)) ||
      any(s < 1) || any(abs(s - round(s)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`s` must contain positive whole numbers")
  }
  if (!is.logical(normalized) || length(normalized) != 1L || is.na(normalized)) {
    .thg_bad_input("`normalized` must be TRUE or FALSE")
  }
  if (!is.null(top) && (length(top) != 1L || !is.finite(top) || top < 1 ||
                       abs(top - round(top)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`top` must be NULL or one positive whole number")
  }

  if (inherits(hg, "net_temporal_hypergraph")) {
    snaps <- hypergraph_snapshots(hg, at = at, mode = snapshot_mode,
                                  multiedges = multiedges)
    rows <- lapply(seq_along(snaps), function(i) {
      ans <- hg_edge_centrality(snaps[[i]], s = s, measure = measure,
                                normalized = normalized, top = top)
      data.frame(time = names(snaps)[i], ans, row.names = NULL,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }

  .thg_check_hg(hg)
  if (hg$n_hyperedges == 0L) {
    return(data.frame(edge = character(), s = integer(), measure = character(),
                      value = numeric(), stringsAsFactors = FALSE))
  }
  rows <- list()
  cursor <- 1L
  for (ss in as.integer(s)) {
    line <- hg_line_graph(hg, s = ss, what = "matrix")
    line <- (line != 0) * 1
    for (metric in measure) {
      value <- switch(metric,
        betweenness = cograph::centrality_betweenness(
          line, weighted = FALSE, directed = FALSE, normalized = normalized
        ),
        closeness = cograph::centrality_closeness(
          line, weighted = FALSE, directed = FALSE, normalized = normalized
        )
      )
      value[!is.finite(value)] <- 0
      tab <- data.frame(
        edge = colnames(line), s = ss, measure = metric,
        value = as.numeric(value[colnames(line)]), stringsAsFactors = FALSE
      )
      tab <- tab[order(-tab$value, tab$edge), , drop = FALSE]
      if (!is.null(top)) tab <- utils::head(tab, as.integer(top))
      rows[[cursor]] <- tab
      cursor <- cursor + 1L
    }
  }
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

#' @rdname hg_edge_centrality
#' @export
hypergraph_edge_centrality <- hg_edge_centrality
