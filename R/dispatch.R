# ---- One verb per idea: generics that dispatch on the input class ----
#
# hg_centrality(), hg_communities() and hypa() each name one idea and
# pick the estimator from the class of `x`: a memory network (net_hon) or a
# hypergraph (net_hg); hypa() also reads raw sequences. Each method keeps
# its own arguments and defaults. An argument that only the other method
# takes arrives in `...` and is refused by .ho_no_dots(), so a misplaced
# argument is an error, never silently ignored.

#' Refuse arguments that the chosen method does not take
#'
#' @param ... The method's `...`.
#' @param .for Human-readable description of the input the method serves.
#' @return `NULL`, invisibly; raises `hypernets_bad_input` when `...` is not
#'   empty.
#' @noRd
.ho_no_dots <- function(..., .for) {
  n <- ...length()
  if (n == 0L) return(invisible(NULL))
  nm <- ...names() %||% rep("", n)
  nm[!nzchar(nm)] <- "<unnamed>"
  .ho_bad_input(sprintf(
    "%s: not an argument for %s",
    paste0("`", nm, "`", collapse = ", "), .for))
}

#' Refuse an input class a dispatching verb has no method for
#'
#' @param x The input.
#' @param verb The verb name.
#' @param accepted What the verb accepts, for the message.
#' @noRd
.ho_no_method <- function(x, verb, accepted) {
  .ho_bad_input(sprintf(
    "%s() takes %s, not an object of class %s", verb, accepted,
    paste(sprintf("`%s`", class(x)), collapse = "/")))
}

#' Centralities of a hypergraph or a memory network
#'
#' One verb, two estimators, chosen by the class of `x`:
#' \describe{
#'   \item{a hypergraph (`net_hg`, e.g. from [text_hypergraph()],
#'     [window_hypergraph()])}{eigenvector-style hypergraph centralities --
#'     clique-motif, Z- and H-eigenvector, EDVW PageRank, subhypergraph and
#'     Katz; see [hg_centrality.net_hg()].}
#'   \item{a memory network (`net_hon`, from [hon()])}{PageRank,
#'     betweenness and closeness of the higher-order topology, projected
#'     onto the first-order states; see [hg_centrality.net_hon()].}
#' }
#' Each method keeps its own arguments; passing an argument that only the
#' other method takes raises `hypernets_bad_input`.
#'
#' @param x A `net_hg` or a `net_hon`.
#' @param ... Arguments of the method for `class(x)`.
#' @return A base `data.frame`, one row per node (or per state), with one
#'   column per requested centrality. Any other input raises
#'   `hypernets_bad_input`.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
#' hg_centrality(hg, type = "clique")
#'
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' hg_centrality(hon(seqs, max_order = 2), type = "pagerank")
#' @export
hg_centrality <- function(x, ...) {
  UseMethod("hg_centrality")
}

#' @rdname hg_centrality
#' @export
hg_centrality.default <- function(x, ...) {
  .ho_no_method(x, "hg_centrality",
                "a hypergraph (net_hg) or a memory network (net_hon)")
}

#' Communities of a hypergraph or a memory network
#'
#' One verb, two estimators, chosen by the class of `x`:
#' \describe{
#'   \item{a hypergraph (`net_hg`)}{an ensemble of Infomap runs on a
#'     projection (or IRMM with `type = "irmm"`), with the adjusted-mutual-
#'     information medoid; see [hg_communities.net_hg()].}
#'   \item{a memory network (`net_hon`, from [hon()])}{the map equation for
#'     memory networks: overlapping modules of the physical states, compared
#'     with the first-order map; see [hg_communities.net_hon()].}
#' }
#' Each method keeps its own arguments; passing an argument that only the
#' other method takes raises `hypernets_bad_input`.
#'
#' @param x A `net_hg` or a `net_hon`.
#' @param ... Arguments of the method for `class(x)`.
#' @return An `hg_communities` object (hypergraph) or a
#'   `net_hon_communities` object (memory network); read either with
#'   [hg_get()]. Any other input raises `hypernets_bad_input`.
#' @examples
#' seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
#'              c("c", "h", "d", "c", "h", "d", "c"))
#' comm <- hg_communities(hon(seqs, max_order = 2L), trials = 2L)
#' hg_get(comm)
#' @export
hg_communities <- function(x, ...) {
  UseMethod("hg_communities")
}

#' @rdname hg_communities
#' @export
hg_communities.default <- function(x, ...) {
  .ho_no_method(x, "hg_communities",
                "a hypergraph (net_hg) or a memory network (net_hon)")
}

#' Hypergeometric anomalies of paths or of co-occurring node pairs
#'
#' One verb, two estimators, chosen by the class of `x`:
#' \describe{
#'   \item{a hypergraph (`net_hg`)}{co-occurring node pairs scored against a
#'     hypergeometric null built from the hyperdegrees; returns a table. See
#'     [hypa.net_hg()].}
#'   \item{sequences (a long, wide or list input, or a model object carrying
#'     sequences; see [sequence-input])}{paths of a De Bruijn graph scored
#'     against the same null (HYPA); returns a `net_hypa` object. See
#'     [hypa.default()].}
#' }
#' Each method keeps its own arguments; passing an argument that only the
#' other method takes raises `hypernets_bad_input`.
#'
#' @param x A `net_hg`, or sequences.
#' @param ... Arguments of the method for `class(x)`.
#' @return See the methods.
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies in
#' time series data on networks. *Proceedings of the 2020 SIAM International
#' Conference on Data Mining*, 460-468. \doi{10.1137/1.9781611976236.52}
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' hg_get(hypa(seqs, order = 2, min_count = 1))
#' @export
hypa <- function(x, ...) {
  UseMethod("hypa")
}

#' Path anomalies (HYPA) in sequences
#'
#' The sequence method of [hypa()]. Builds the De Bruijn graph of order
#' `order` from the sequences and scores every path of `order + 1` states
#' against a hypergeometric null whose propensities are fitted to the
#' observed path counts: paths observed far more often than the null
#' predicts are over-represented, far less often under-represented.
#'
#' @param x Sequences in any form described in [sequence-input]: a long
#'   event table (with `action`), a wide data.frame, a list of vectors, or a
#'   model object carrying its sequences.
#' @param order Integer (or integer vector) of De Bruijn orders. Default
#'   `2`.
#' @param alpha Significance level of the anomaly labels, in (0, 0.5).
#'   Default `0.05`.
#' @param min_count Integer. Paths observed fewer times are always labelled
#'   `"normal"`. Default `5`.
#' @param p_adjust Multiplicity correction: any [stats::p.adjust] method or
#'   `"none"`. Default `"BH"`; the two tails are adjusted separately.
#' @param type Which anomalous paths the printed result lists: `"all"`
#'   (default), `"over"` (paths observed more often than expected) or
#'   `"under"` (less often).
#' @param order_by How the printed paths are ordered: `"sig"` (default,
#'   smallest p-value of the path's own direction first), `"ratio"` (largest
#'   deviation first: highest observed/expected ratio for over-represented
#'   paths, lowest for under-represented ones), `"freq"` (most observed
#'   first) or `"path"` (alphabetical).
#' @param n Number of paths printed per direction. Default `10`.
#' @param action,actor,time,session,time_threshold,timezone Long-format
#'   arguments (see [sequence-input]); leave the column names `NULL` for
#'   wide, list or model input.
#' @param ... Must be empty: an argument that only the hypergraph method
#'   takes (`top`) raises `hypernets_bad_input`.
#' @return A `net_hypa` object (also a `cograph_network`) that prints the
#'   selected anomalous paths. `type`, `order_by` and `n` change only what is
#'   printed and the default of [hg_get()]; every path is scored. Read the
#'   tables with [hg_get()]: the anomalous paths (`what = "anomalies"`, the
#'   default), every scored path (`"scores"`), the over- or under-represented
#'   ones (`"over"`, `"under"`), or the anomalous paths as pathways
#'   (`"pathways"`).
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies in
#' time series data on networks. *Proceedings of the 2020 SIAM International
#' Conference on Data Mining*, 460-468. \doi{10.1137/1.9781611976236.52}
#' @seealso [hypa.net_hg()] for the same null on hypergraph co-occurrence.
#' @examples
#' seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
#'              c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
#' fit <- hypa(seqs, order = 2, min_count = 1, type = "over")
#' fit
#' hg_get(fit)
#' @export
hypa.default <- function(x, order = 2L, alpha = 0.05, min_count = 5L,
                            p_adjust = "BH",
                            type = c("all", "over", "under"),
                            order_by = c("sig", "ratio", "freq", "path"),
                            n = 10L, action = NULL, actor = NULL,
                            time = NULL, session = NULL,
                            time_threshold = 900, timezone = "UTC", ...) {
  .ho_no_dots(..., .for = "sequences")
  type <- match.arg(type)
  order_by <- match.arg(order_by)
  stopifnot("`n` must be a single count >= 1" =
              is.numeric(n) && length(n) == 1L && is.finite(n) && n >= 1)
  data <- .ho_sequence_input(x, action = action, actor = actor, time = time,
                             session = session,
                             time_threshold = time_threshold,
                             timezone = timezone)
  fit <- Nestimate::build_hypa(data, order = order, alpha = alpha,
                               min_count = min_count, p_adjust = p_adjust)
  structure(fit, class = c("hypernets_hypa", class(fit)),
            hypa_view = list(type = type, order_by = order_by,
                             n = as.integer(n)))
}

#' @rdname result-summary
#' @export
summary.hypernets_hypa <- function(object, ...) .ho_summary(object)

#' @rdname hypa.default
#' @export
print.hypernets_hypa <- function(x, n = NULL, ...) {
  view <- attr(x, "hypa_view")
  cat(sprintf(paste0("Path anomalies (HYPA, order %d): %d of %d paths ",
                     "anomalous at alpha = %s (%d over, %d under, %s)\n"),
              x$order, x$n_anomalous, x$n_edges, format(x$alpha), x$n_over,
              x$n_under, x$p_adjust))
  .ho_print_table(x, n %||% view$n %||% 10L)
  invisible(x)
}

#' @rdname hg_get.net_hypa
#' @export
hg_get.hypernets_hypa <- function(x, ...) NextMethod()
