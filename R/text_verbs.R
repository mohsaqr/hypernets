# Tidy text-facing analysis verbs. Each delegates to the matching engine in
# the hypergraph family (hypergraph_measures, hypergraph_centrality,
# hypergraph_cluster, hypergraph_transduction) and returns a base data.frame.

.thg_check_hg <- function(hg) {
  if (!inherits(hg, "net_hypergraph")) {
    stop(errorCondition(
      "`hg` must be a net_hypergraph (text_hypergraph, knn_hypergraph, group_hypergraph, ...)",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  invisible(hg)
}

#' Structural measures of a hypergraph, as tidy tables
#'
#' Delegates to [hypergraph_measures()] and returns the requested
#' slice as a tidy data.frame.
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hypergraph`).
#' @param what Which table: `"nodes"` (default; one row per node with
#'   `hyperdegree`, `strength`, `max_edge_size` and `n_neighbors`, the
#'   distinct nodes it shares a hyperedge with), `"edges"` (one row per
#'   hyperedge with its `size`), `"overlap"` (one row per hyperedge pair with
#'   `overlap`, `overlap_coefficient`, `jaccard`), `"summary"` (one row per
#'   scalar measure), `"distribution"` (the empirical distribution of
#'   `measure`, as a `hypernets_distribution` table whose `plot()` draws the
#'   CCDF), or `"components"` (one row per connected component through shared
#'   hyperedges, with its `n_nodes`, `n_edges`, `share` of nodes and
#'   `diameter`).
#' @param measure For `what = "distribution"`: `"hyperdegree"` (default),
#'   `"strength"`, `"n_neighbors"` (node measures) or `"size"` (hyperedge
#'   cardinality).
#' @return A base `data.frame`, one row per node, edge, edge pair, measure,
#'   distinct value, or component according to `what`.
#' @examples
#' hg <- text_hypergraph(c(
#'   a = "salt and soup and onions",
#'   b = "soup and salt",
#'   c = "stars and sky"
#' ))
#' hg_measures(hg)
#' hg_measures(hg, what = "summary")
#' hg_measures(hg, what = "components")
#' @export
hg_measures <- function(hg, what = c("nodes", "edges", "overlap", "summary",
                                     "distribution", "components"),
                        measure = c("hyperdegree", "strength", "n_neighbors",
                                    "size")) {
  .thg_check_hg(hg)
  what <- match.arg(what)
  measure <- match.arg(measure)
  if (identical(what, "distribution")) {
    values <- if (identical(measure, "size")) {
      hg_measures(hg, what = "edges")$size
    } else {
      hg_measures(hg, what = "nodes")[[measure]]
    }
    out <- .thg_distribution(values)
    class(out) <- c("hypernets_distribution", "data.frame")
    attr(out, "measure") <- measure
    return(out)
  }
  if (identical(what, "components")) return(.thg_component_table(hg))
  if (.thg_is_sparse(hg)) {
    return(.thg_sparse_measures(hg, what))
  }
  m <- hypergraph_measures(hg)
  edges <- colnames(hg$incidence)
  switch(what,
    nodes = data.frame(
      node = names(m$hyperdegree),
      hyperdegree = as.integer(m$hyperdegree),
      strength = as.numeric(m$node_strength),
      max_edge_size = as.integer(m$max_edge_size),
      n_neighbors = as.integer(rowSums(m$co_degree > 0)),
      row.names = NULL
    ),
    edges = data.frame(
      edge = edges,
      size = as.integer(m$edge_sizes),
      row.names = NULL
    ),
    overlap = {
      pair <- which(upper.tri(m$edge_pairwise_overlap), arr.ind = TRUE)
      data.frame(
        edge_1 = edges[pair[, "row"]],
        edge_2 = edges[pair[, "col"]],
        overlap = as.numeric(m$edge_pairwise_overlap[pair]),
        overlap_coefficient = as.numeric(m$overlap_coefficient[pair]),
        jaccard = as.numeric(m$jaccard[pair]),
        row.names = NULL
      )
    },
    summary = data.frame(
      measure = c("n_nodes", "n_hyperedges", "density", "avg_edge_size",
                  "pairwise_participation"),
      value = c(m$n_nodes, m$n_hyperedges, m$density, m$avg_edge_size,
                m$pairwise_participation),
      row.names = NULL
    )
  )
}

#' Hypergraph node centralities, as a tidy table
#'
#' Delegates to [hypergraph_centrality()]: clique-expansion
#' eigenvector centrality and the tensor Z- and H-eigenvector centralities.
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hypergraph`).
#' @param type Centralities to compute; any of `"clique"`, `"Z"`, `"H"`
#'   (default: all three).
#' @param sort_by Optional centrality name to sort by, descending (ties broken
#'   by node name); default keeps node order.
#' @param n Return only the first `n` rows after sorting (default all) --
#'   e.g. `sort_by = "clique", n = 10` for the ten most central nodes.
#' @param max_iter,tol,normalize Passed to
#'   [hypergraph_centrality()].
#' @return A base `data.frame`, one row per node (or the `n` requested rows),
#'   with one column per requested centrality.
#' @examples
#' hg <- text_hypergraph(c(
#'   a = "salt and soup and onions",
#'   b = "soup and salt",
#'   c = "stars and salt"
#' ))
#' hg_centrality(hg, type = "clique")
#' @export
hg_centrality <- function(hg, type = c("clique", "Z", "H"),
                          sort_by = NULL, n = Inf,
                          max_iter = 1000L, tol = 1e-8, normalize = TRUE) {
  .thg_check_hg(hg)
  type <- match.arg(type, several.ok = TRUE)
  if (.thg_is_sparse(hg)) {
    stop(errorCondition(
      "tensor centralities need the dense representation; use hg_pagerank() at scale",
      class = "hypernets_sparse_unsupported", call = NULL
    ))
  }
  stopifnot(
    "`n` must be a single count >= 1" =
      length(n) == 1L && (is.infinite(n) || (is.finite(n) && n >= 1))
  )
  out <- hypergraph_centrality(
    hg, type = type, max_iter = max_iter, tol = tol, normalize = normalize
  )
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, choices = type)
    out <- out[order(-out[[sort_by]], out$node), , drop = FALSE]
    rownames(out) <- NULL
  }
  if (is.finite(n) && n < nrow(out)) {
    out <- out[seq_len(n), , drop = FALSE]
  }
  out
}

#' Spectral clustering of a hypergraph, as a tidy table
#'
#' Calls the in-package [hypergraph_cluster()] engine (Zhou et al. 2006 normalized
#' Laplacian, or the Hayashi et al. 2020 random-walk Laplacian with
#' edge-dependent vertex weights -- the natural choice for tf-idf-weighted
#' text hypergraphs).
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hypergraph`).
#' @param k Number of clusters (explicit by design; there is no correct
#'   default).
#' @param type `"zhou"` or `"random_walk"`, as in
#'   [hypergraph_cluster()].
#' @param algorithm `"spectral"` (RDC-Spec) or `"symnmf"` (RDC-Sym), as in
#'   [hypergraph_cluster()]. SymNMF currently requires a dense incidence
#'   matrix because its paper objective factorizes a dense node similarity.
#' @param edge_weights Hyperedge weights for the Laplacian: `NULL`
#'   (default, every hyperedge counts once), a positive numeric vector with
#'   one entry per hyperedge, or `"idf"` to weight each word hyperedge by
#'   its smoothed inverse document frequency (needs `weight = "tfidf"`), so
#'   a word used in most documents binds them less than a rare shared word.
#'   Note that the default `type = "zhou"` reads only which documents a
#'   hyperedge binds; the incidence cells (`weight =`) enter the
#'   `"random_walk"` type only.
#' @param seed Random seed passed to the selected solver; set it for a
#'   reproducible partition.
#' @param nstart Number of solver starts (default `25L`).
#' @param max_iter,tol SymNMF convergence controls passed to
#'   [hypergraph_cluster()].
#' @param what What to return: `"clusters"` (default) for the partition,
#'   `"embedding"` for the partition plus the row-normalized spectral
#'   embedding used by k-means (`dim1..dimk` -- plot these to map the
#'   corpus) and the stationary weight `pi`, or `"eigenvalues"` for the
#'   Laplacian spectrum (dense engines return all `n` values; the sparse
#'   engine returns the `k + 1` it computed -- raise `k` for an eigengap
#'   scan).
#' @return A base `data.frame`. For `what = "clusters"`: one row per node,
#'   columns `node` and `cluster`. For `what = "embedding"`: `node`,
#'   `cluster`, `pi`, `dim1..dimk`. For `what = "eigenvalues"`: one row
#'   per eigenvalue, columns `index`, `value` (ascending) and `gap` (the
#'   distance to the next eigenvalue -- large gaps indicate supported
#'   cluster counts; `NA` on the last row).
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' hg_cluster(hg, k = 2, type = "random_walk", seed = 1)
#' hg_cluster(hg, k = 2, seed = 1, what = "embedding")
#' hg_cluster(hg, k = 2, seed = 1, what = "eigenvalues")
#' @export
hg_cluster <- function(hg, k, type = c("zhou", "random_walk"),
                       edge_weights = NULL, seed = NULL, nstart = 25L,
                       what = c("clusters", "embedding", "eigenvalues"),
                       algorithm = c("spectral", "symnmf"),
                       max_iter = 500L, tol = 1e-6) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  algorithm <- match.arg(algorithm)
  what <- match.arg(what)
  edge_weights <- .thg_edge_weights(hg, edge_weights)
  if (.thg_is_sparse(hg) && algorithm == "symnmf") {
    stop(errorCondition(
      paste0("`algorithm = \"symnmf\"` requires a dense incidence matrix; ",
             "RDC-Sym factorizes a dense n_nodes x n_nodes similarity."),
      class = "hypernets_dense_required", call = NULL
    ))
  }
  fit <- if (.thg_is_sparse(hg)) {
    .thg_sparse_cluster(hg, k = k, type = type, edge_weights = edge_weights,
                        nstart = nstart, seed = seed)
  } else {
    hypergraph_cluster(hg, k = k, type = type, edge_weights = edge_weights,
                       algorithm = algorithm, seed = seed, nstart = nstart,
                       max_iter = max_iter, tol = tol)
  }
  if (identical(what, "eigenvalues")) {
    values <- as.numeric(fit$eigenvalues)
    return(data.frame(index = seq_along(values), value = values,
                      gap = c(diff(values), NA_real_)))
  }
  out <- if (identical(what, "embedding")) {
    as.data.frame(fit)
  } else {
    fit$clusters
  }
  rownames(out) <- NULL
  out
}

# Resolve `edge_weights`: NULL, a numeric vector (one per hyperedge), or
# "idf" -- the tf-idf vocabulary's idf, aligned to the incidence columns.
.thg_edge_weights <- function(hg, edge_weights) {
  if (is.null(edge_weights)) {
    return(NULL)
  }
  if (identical(edge_weights, "idf")) {
    vocab <- hg$text$vocabulary
    ok <- !is.null(vocab) && "idf" %in% names(vocab) &&
      identical(hg$text$nodes, "doc")
    if (!ok) {
      stop(errorCondition(
        paste0("`edge_weights = \"idf\"` needs a text_hypergraph(weight = ",
               "\"tfidf\", nodes = \"doc\"), whose hyperedges are words"),
        class = "hypernets_bad_input", call = NULL
      ))
    }
    idf <- stats::setNames(vocab$idf, vocab$word)[colnames(hg$incidence)]
    stopifnot("internal: idf must cover every hyperedge" = !anyNA(idf))
    return(as.numeric(idf))
  }
  if (!is.numeric(edge_weights) || anyNA(edge_weights) ||
      any(edge_weights <= 0) ||
      !length(edge_weights) %in% c(1L, hg$n_hyperedges)) {
    stop(errorCondition(
      "`edge_weights` must be NULL, \"idf\", or positive numbers, one per hyperedge",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  rep_len(as.numeric(edge_weights), hg$n_hyperedges)
}

#' Characteristic words (keywords) per cluster
#'
#' For a clustered hypergraph, ranks each cluster's hyperedges -- on a
#' `nodes = "doc"` [text_hypergraph()], its words -- by one of four scores
#' computed from the hypergraph, or by a score table supplied from another
#' fit. This is the topic-description step that turns a partition from
#' [hg_cluster()] into named topics.
#'
#' * `"mass"` (default): the summed incidence weight a cluster's documents
#'   place on the word (tf-idf mass under `weight = "tfidf"`, counts under
#'   `weight = "n"`). Works on any hypergraph.
#' * `"frequency"`: the raw token count of the word in the cluster's
#'   documents, whatever the incidence weighting.
#' * `"ctfidf"`: class-based tf-idf (Grootendorst 2022), the topic
#'   representation BERTopic uses: the cluster's L1-normalised term
#'   frequency times `log(A / f_w + 1)`, where `f_w` is the word's corpus
#'   count and `A` the average token count per cluster (truncated to an
#'   integer, as in the reference implementation).
#' * `"centrality"`: the word's centrality within the cluster's own word
#'   hypergraph (words as nodes, the cluster's documents as hyperedges,
#'   incidence = the stored weights), from [hypergraph_centrality()] with
#'   the measure named in `centrality`.
#'
#' `"frequency"`, `"ctfidf"` and `"centrality"` need the token-level layer
#' of a bag-of-words document hypergraph (`construction = "bag"`,
#' `nodes = "doc"`); other hypergraphs raise `hypernets_bad_input`.
#'
#' **Sentence scope.** When `hg` is a `text_hypergraph(construction =
#' "sentence")` of the documents, `clusters` still names documents, and the
#' scores are computed over the cluster's *sentences*: two words are then
#' related by sharing a sentence rather than a document. `"mass"` and
#' `"frequency"` sum the in-sentence counts; `"centrality"` is the word's
#' centrality among the cluster's sentence hyperedges.
#'
#' **External scores.** `scores` takes a data.frame with columns `node`,
#' `word` and one numeric column (for instance the `attention` table from
#' `hg_hypergat(what = "attention")`). Its values are summed per cluster
#' and reported as their own block, labelled by the numeric column's name.
#'
#' @param hg The hypergraph the clustering was computed on, or the
#'   sentence hypergraph of the same documents (see Sentence scope).
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels.
#' @param n Keywords per cluster (default `10`); `Inf` returns all.
#' @param type Which score ranks the words: any of `"mass"`,
#'   `"frequency"`, `"ctfidf"`, `"centrality"`. Several at once give one
#'   block of rows per score. Default `"mass"`, or none when only
#'   `scores` is wanted.
#' @param sort_by `"score"` (default) ranks a cluster's words by the
#'   selected score, the heavy vocabulary; `"share"` ranks them by the
#'   fraction of the word's total score that falls in the cluster, the
#'   distinctive vocabulary. Ties are broken by the other column, then
#'   alphabetically.
#' @param min_docs Keep only words that occur in at least this many of the
#'   cluster's documents (default `1`). The support floor that stops
#'   `sort_by = "share"` from surfacing words a cluster owns because they
#'   occur once.
#' @param centrality For `type = "centrality"`, the [hypergraph_centrality()]
#'   measure: `"pagerank"` (default; it reads the incidence weights and does
#'   not tie on words present in every document, as the clique measure
#'   does), `"clique"`, `"Z"` or `"H"`.
#' @param scores An external score table (columns `node`, `word` and one
#'   numeric column), summed per cluster and added as its own block.
#' @param collapse If `TRUE`, return one row per type and cluster with the
#'   words joined into a single comma-separated string -- a display table,
#'   not tidy data. Default `FALSE`.
#' @param x For `plot()` and `print()`: a table returned by `hg_keywords()`.
#' @param value For `plot()`: which column the bars show, `"score"`
#'   (default) or `"share"`; match it to `sort_by`.
#' @param label For `plot()`: print each bar's value at its end
#'   (default `TRUE`).
#' @param ncol For `plot()`: panels per row; default one column per score
#'   type, so clusters run down and types across. With a single type and
#'   many clusters, pass e.g. `ncol = 4`.
#' @param ... Unused; for S3 consistency.
#' @return A base `data.frame` of class `hypernets_keywords`, one row per
#'   type-cluster-keyword triple, columns `type`, `cluster`, `size` (the
#'   cluster's documents), `rank`, `word`, `score` (the selected score),
#'   `share` (`score` divided by the word's summed score over all
#'   clusters) and `n_docs` (the cluster's documents containing the word),
#'   ranked by `sort_by` within type and cluster; only words with a
#'   positive score and at least `min_docs` documents appear. With
#'   `collapse = TRUE`: one row per type and cluster, columns `type`,
#'   `cluster`, `size` and `words`. The print method shows the collapsed
#'   view, truncated to the console width; `as.data.frame()` is the long
#'   form. Raises `hypernets_bad_input` for unknown node names, a `type` that
#'   needs the token layer on a hypergraph without one, or a malformed
#'   `scores` table.
#'
#'   `plot()` returns a ggplot: one panel per cluster (rows) and score type
#'   (columns), each with its own word axis, horizontal bars of `value` per
#'   word, Okabe-Ito fill by type. It needs the long form and raises
#'   `hypernets_bad_input` on a collapsed table.
#' @references
#' Grootendorst, M. (2022). BERTopic: Neural topic modeling with a
#' class-based TF-IDF procedure. arXiv:2203.05794.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_keywords(hg, topics, n = 3)
#' hg_keywords(hg, topics, n = 3, type = "ctfidf")
#' words <- hg_keywords(hg, topics, n = 4,
#'                      type = c("frequency", "ctfidf", "centrality"))
#' plot(words)
#' # distinctive rather than heavy vocabulary: rank by share, with support
#' hg_keywords(hg, topics, n = 3, sort_by = "share", min_docs = 1)
#' @export
hg_keywords <- function(hg, clusters, n = 10L, type = NULL,
                        sort_by = c("score", "share"), min_docs = 1L,
                        centrality = c("pagerank", "clique", "Z", "H"),
                        scores = NULL, collapse = FALSE) {
  .thg_check_hg(hg)
  sort_by <- match.arg(sort_by)
  choices <- c("mass", "frequency", "ctfidf", "centrality")
  type <- type %||% (if (is.null(scores)) "mass" else character(0))
  if (!is.character(type) || !all(type %in% choices) ||
      anyDuplicated(type) > 0L) {
    stop(errorCondition(
      paste0("`type` must be distinct values from: ",
             paste(choices, collapse = ", ")),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  centrality <- match.arg(centrality)
  stopifnot(
    "`n` must be a single positive number" =
      length(n) == 1L && is.numeric(n) && n >= 1,
    "`min_docs` must be a single number >= 1" =
      length(min_docs) == 1L && is.numeric(min_docs) && min_docs >= 1,
    "`collapse` must be TRUE or FALSE" =
      isTRUE(collapse) || isFALSE(collapse)
  )
  assignment <- .thg_labels_input(clusters)
  stopifnot(
    "`clusters` must be a data.frame or a named vector" =
      !is.null(names(assignment))
  )
  scope <- .thg_kw_scope(hg)
  unknown <- setdiff(names(assignment), scope$docs)
  if (length(unknown) > 0L) {
    stop(errorCondition(
      paste0("Unknown node names in `clusters`: ",
             paste(unknown, collapse = ", ")),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  groups <- factor(as.character(assignment)[match(scope$docs,
                                                  names(assignment))])
  groups <- factor(groups, levels = .thg_kw_natural(levels(groups)))
  sizes <- table(groups)

  blocks <- lapply(type, \(one) {
    list(label = one,
         scores = switch(
           one,
           mass = .thg_kw_aggregate(groups, scope$mass),
           frequency = .thg_kw_aggregate(groups, .thg_kw_counts(hg, scope)),
           ctfidf = .thg_kw_ctfidf(.thg_kw_aggregate(groups,
                                                     .thg_kw_counts(hg, scope))),
           centrality = .thg_kw_centrality(hg, scope, groups, centrality)
         ),
         support = .thg_kw_aggregate(groups, (scope$mass != 0) * 1))
  })
  if (!is.null(scores)) {
    external <- .thg_kw_external(scope, groups, scores)
    blocks <- c(blocks, list(external))
  }
  if (length(blocks) == 0L) {
    stop(errorCondition(
      "nothing to rank: give at least one `type`, or `scores`",
      class = "hypernets_bad_input", call = NULL
    ))
  }

  per_block <- lapply(blocks, \(block) {
    total <- colSums(block$scores)
    per_cluster <- lapply(rownames(block$scores), \(cl) {
      row <- block$scores[cl, ]
      share <- ifelse(total > 0, row / total, 0)
      support <- block$support[cl, names(row)]
      ord <- if (identical(sort_by, "share")) {
        order(-share, -row, names(row))
      } else {
        order(-row, -share, names(row))
      }
      eligible <- row[ord] > 0 & support[ord] >= min_docs
      keep <- utils::head(ord[eligible], n)
      data.frame(
        type = block$label, cluster = cl, size = as.integer(sizes[[cl]]),
        rank = seq_along(keep), word = names(row)[keep],
        score = as.numeric(row[keep]), share = as.numeric(share[keep]),
        n_docs = as.numeric(support[keep]),
        stringsAsFactors = FALSE
      )
    })
    do.call(rbind, per_cluster)
  })
  out <- do.call(rbind, per_block)
  rownames(out) <- NULL
  if (isTRUE(collapse)) {
    out <- .thg_kw_collapse(out)
  }
  class(out) <- c("hypernets_keywords", "data.frame")
  out
}

# One row per type and cluster with the words joined, in the long form's
# order; the display table behind collapse = TRUE and the print method.
.thg_kw_collapse <- function(long) {
  key <- paste(long$type, long$cluster, sep = "\r")
  first <- !duplicated(key)
  words <- vapply(split(long$word, factor(key, levels = unique(key))),
                  paste, character(1), collapse = ", ")
  out <- data.frame(
    type = long$type[first], cluster = long$cluster[first],
    size = long$size[first], words = as.character(words),
    stringsAsFactors = FALSE
  )
  rownames(out) <- NULL
  out
}

#' @rdname hg_keywords
#' @export
print.hypernets_keywords <- function(x, ...) {
  shown <- if ("words" %in% names(x)) as.data.frame(x) else
    .thg_kw_collapse(as.data.frame(x))
  # fit the words column to the console: the other columns plus separators
  # take a fixed width, the rest goes to the words, cut with a marker
  fixed <- max(nchar(shown$type), 4L) + max(nchar(shown$cluster), 7L) +
    max(nchar(shown$size), 4L) + 6L
  room <- max(20L, getOption("width", 80L) - fixed)
  long_words <- nchar(shown$words) > room
  shown$words[long_words] <- paste0(substr(shown$words[long_words], 1L,
                                           room - 3L), "...")
  print(shown, right = FALSE, row.names = FALSE)
  if (!"words" %in% names(x)) {
    cat(sprintf("%d rows in the long form (rank, score, share, n_docs): as.data.frame()\n",
                nrow(x)))
  }
  invisible(x)
}

#' @rdname hg_keywords
#' @export
plot.hypernets_keywords <- function(x, value = c("score", "share"),
                                 label = TRUE, ncol = NULL, ...) {
  value <- match.arg(value)
  stopifnot("`ncol` must be NULL or a single positive number" =
              is.null(ncol) || (length(ncol) == 1L && is.numeric(ncol) &&
                                  ncol >= 1))
  if (!all(c("type", "cluster", "rank", "word", "score", "share") %in%
           names(x))) {
    stop(errorCondition(
      "plot() needs the long form of hg_keywords() (collapse = FALSE)",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  stopifnot("`label` must be TRUE or FALSE" = isTRUE(label) || isFALSE(label))
  d <- as.data.frame(x)
  types <- unique(d$type)
  d$type <- factor(d$type, levels = types)
  d$cluster <- factor(d$cluster, levels = unique(d$cluster))
  # one bar per (type, cluster, word); the invisible suffix keeps the same
  # word orderable independently in every panel
  key <- paste(d$word, d$type, d$cluster, sep = "\r")
  ord <- order(d$type, d$cluster, -d$rank)
  d$key <- factor(key, levels = unique(key[ord]))
  d$value <- d[[value]]
  d$label <- format(signif(d$value, 3), scientific = FALSE, trim = TRUE,
                    drop0trailing = TRUE)
  # Okabe-Ito, yellow last: it is too faint on a white panel to lead with
  okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                 "#CC79A7", "#999999", "#000000", "#F0E442")
  fill <- stats::setNames(okabe_ito[seq_along(types)], types)
  p <- ggplot2::ggplot(
    d, ggplot2::aes(x = .data$value, y = .data$key, fill = .data$type)
  ) +
    ggplot2::geom_col(width = 0.75) +
    # facet_wrap, not facet_grid: each panel needs its own word axis
    ggplot2::facet_wrap(ggplot2::vars(.data$cluster, .data$type),
                        scales = "free", ncol = ncol %||% length(types),
                        labeller = ggplot2::label_value) +
    ggplot2::scale_y_discrete(labels = \(k) sub("\r.*$", "", k)) +
    ggplot2::scale_fill_manual(values = fill, guide = "none") +
    # room for the value label at the end of the longest bar
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, .4))) +
    ggplot2::labs(x = value, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank(),
                   panel.spacing.x = ggplot2::unit(1.2, "lines"))
  if (isTRUE(label)) {
    p <- p + ggplot2::geom_text(ggplot2::aes(label = .data$label),
                                hjust = -0.15, size = 3)
  }
  p
}

# Cluster labels in natural order: by their number when every label carries
# one ("Cluster 2" before "Cluster 10"), alphabetically otherwise.
.thg_kw_natural <- function(labels) {
  numbers <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", labels)))
  if (length(labels) > 0L && !anyNA(numbers) && anyDuplicated(numbers) == 0L) {
    labels[order(numbers)]
  } else {
    sort(labels)
  }
}

# The document universe and the document x word mass matrix a keyword
# scope is computed on. A document hypergraph gives its own nodes and
# incidence; a sentence hypergraph gives its documents (from the sentences
# table) and the in-sentence counts summed per document. `edges` carries the
# per-hyperedge weights (edge, word, weight, doc) for centrality.
.thg_kw_scope <- function(hg) {
  layer <- hg$text
  if (!is.null(layer) && identical(layer$construction, "sentence")) {
    sentences <- layer$sentences
    weights <- layer$weights
    weights$doc <- sentences$doc[match(weights$edge, sentences$edge)]
    docs <- unique(sentences$doc)
    words <- sort(unique(weights$word))
    mass <- Matrix::sparseMatrix(
      i = match(weights$doc, docs), j = match(weights$word, words),
      x = as.numeric(weights$weight), dims = c(length(docs), length(words)),
      dimnames = list(docs, words)
    )
    return(list(docs = docs, mass = mass, token = weights,
                edges = weights[, c("edge", "word", "weight", "doc")],
                sentence = TRUE))
  }
  token <- NULL
  if (!is.null(layer) && identical(layer$construction, "bag") &&
      identical(layer$nodes, "doc") && is.data.frame(layer$weights) &&
      all(c("doc", "word", "n", "weight") %in% names(layer$weights))) {
    token <- layer$weights
  }
  list(docs = hg$nodes, mass = hg$incidence, token = token, edges = NULL,
       sentence = FALSE)
}

# cluster x hyperedge indicator-weighted sums of a node x hyperedge matrix;
# unassigned nodes (NA group) contribute nothing.
.thg_kw_aggregate <- function(groups, m) {
  assigned <- which(!is.na(groups))
  indicator <- Matrix::sparseMatrix(
    i = as.integer(groups[assigned]), j = assigned, x = 1,
    dims = c(nlevels(groups), nrow(m))
  )
  out <- as.matrix(indicator %*% m)
  dimnames(out) <- list(levels(groups), colnames(m))
  out
}

# document x word raw token counts (the token layer), or the in-sentence
# counts for a sentence scope
.thg_kw_counts <- function(hg, scope) {
  if (isTRUE(scope$sentence)) {
    return(scope$mass)
  }
  if (is.null(scope$token)) {
    stop(errorCondition(
      paste0("this `type` needs the token layer of a bag-of-words document ",
             "hypergraph (text_hypergraph(construction = \"bag\", ",
             "nodes = \"doc\")) or a sentence hypergraph; use type = \"mass\" otherwise"),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  weights <- scope$token
  Matrix::sparseMatrix(
    i = match(weights$doc, scope$docs),
    j = match(weights$word, colnames(scope$mass)),
    x = as.numeric(weights$n),
    dims = dim(scope$mass), dimnames = dimnames(scope$mass)
  )
}

# Class-based tf-idf exactly as BERTopic's ClassTfidfTransformer computes it
# (Grootendorst 2022): rows L1-normalised, idf = log(A / f_w + 1) with A the
# integer-truncated mean class token count and f_w the word's total count.
.thg_kw_ctfidf <- function(counts) {
  class_total <- rowSums(counts)
  word_total <- colSums(counts)
  avg_tokens <- trunc(mean(class_total))
  tf <- counts / ifelse(class_total > 0, class_total, 1)
  idf <- ifelse(word_total > 0, log(avg_tokens / word_total + 1), 0)
  sweep(tf, 2L, idf, `*`)
}

# One word hypergraph per cluster (words = nodes; the cluster's documents,
# or its sentences under sentence scope, = hyperedges; stored weights =
# incidence); rank by hypergraph_centrality().
.thg_kw_centrality <- function(hg, scope, groups, centrality) {
  edges <- if (isTRUE(scope$sentence)) {
    scope$edges
  } else {
    if (is.null(scope$token)) {
      stop(errorCondition(
        paste0("`type = \"centrality\"` needs the token layer of a ",
               "bag-of-words document hypergraph or a sentence hypergraph"),
        class = "hypernets_bad_input", call = NULL
      ))
    }
    data.frame(edge = scope$token$doc, word = scope$token$word,
               weight = scope$token$weight, doc = scope$token$doc,
               stringsAsFactors = FALSE)
  }
  words <- colnames(scope$mass)
  rows <- lapply(levels(groups), \(cl) {
    docs <- scope$docs[!is.na(groups) & groups == cl]
    sub <- edges[edges$doc %in% docs, , drop = FALSE]
    row <- stats::setNames(numeric(length(words)), words)
    if (nrow(sub) == 0L) {
      return(row)
    }
    word_hg <- group_hypergraph(sub, actor = "word", group = "edge",
                                weight = "weight")
    values <- hypergraph_centrality(word_hg, type = centrality)
    row[values$node] <- values[[centrality]]
    row
  })
  out <- do.call(rbind, rows)
  dimnames(out) <- list(levels(groups), words)
  out
}

# An external (node, word, value) table summed per cluster, as its own
# block; support counts the cluster's documents with a row for the word.
.thg_kw_external <- function(scope, groups, scores) {
  ok <- is.data.frame(scores) && all(c("node", "word") %in% names(scores)) &&
    ncol(scores) >= 3L
  value_col <- if (ok) {
    setdiff(names(scores)[vapply(scores, is.numeric, logical(1))],
            c("node", "word"))
  } else {
    character(0)
  }
  if (!ok || length(value_col) == 0L) {
    stop(errorCondition(
      paste0("`scores` must be a data.frame with columns `node`, `word` and ",
             "one numeric score column (for instance the table from ",
             "hg_hypergat(what = \"attention\"))"),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  value_col <- value_col[[1L]]
  scores <- scores[scores$node %in% scope$docs, , drop = FALSE]
  if (nrow(scores) == 0L) {
    stop(errorCondition(
      "no rows of `scores` refer to documents of `hg`",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  words <- sort(unique(as.character(scores$word)))
  cell <- function(values) {
    Matrix::sparseMatrix(
      i = match(scores$node, scope$docs), j = match(scores$word, words),
      x = values, dims = c(length(scope$docs), length(words)),
      dimnames = list(scope$docs, words)
    )
  }
  list(label = value_col,
       scores = .thg_kw_aggregate(groups, cell(as.numeric(scores[[value_col]]))),
       support = .thg_kw_aggregate(groups, cell(rep(1, nrow(scores)))))
}

#' Relations between topics: a weighted network of shared vocabulary
#'
#' For a clustered document hypergraph, builds the topic-by-topic
#' co-occurrence network of the bibliometric kind: the strength between
#' two topics is the sum, over all words, of the product of the number of
#' documents in each topic that contain the word (full counting, the
#' aggregation bibnets uses for keyword co-occurrence). The raw sum can be
#' normalised by the same similarity measures as bibnets' `normalize()`,
#' with the diagonal of the co-occurrence matrix as each topic's total.
#' The result is an edge list with `source`, `target` and `weight`, the
#' three columns Gephi reads, or a `cograph_network` for plotting.
#'
#' With `D_i` the diagonal entry of topic i and `A_ij` the raw sum:
#' `"association"` is `A_ij / (D_i D_j)`, `"cosine"` is
#' `A_ij / sqrt(D_i D_j)`, `"jaccard"` is `A_ij / (D_i + D_j - A_ij)`,
#' `"inclusion"` is `A_ij / min(D_i, D_j)` and `"equivalence"` is
#' `A_ij^2 / (D_i D_j)`.
#'
#' @param hg The document hypergraph the clustering was computed on.
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels.
#' @param similarity `"none"` (default, the raw sum), `"association"`,
#'   `"cosine"`, `"jaccard"`, `"inclusion"` or `"equivalence"`.
#' @param what `"edges"` (default) for the edge list, `"network"` for a
#'   `cograph_network` built from it with [cograph::as_cograph()], whose
#'   node table carries each topic's `size`.
#' @return For `what = "edges"`: a base `data.frame`, one row per pair of
#'   topics with a positive weight, columns `source`, `target`, `weight`,
#'   pairs in the topics' natural order. For `what = "network"`: a
#'   `cograph_network` with one node per topic (`label`, `name`, `size`)
#'   and one undirected weighted edge per pair. Raises `hypernets_bad_input`
#'   for unknown node names.
#' @references
#' van Eck, N. J., & Waltman, L. (2009). How to normalize cooccurrence
#' data? An analysis of some well-known similarity measures. *Journal of
#' the American Society for Information Science and Technology*, 60(8),
#' 1635--1651.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_relations(hg, topics)
#' hg_relations(hg, topics, similarity = "cosine")
#' @export
hg_relations <- function(hg, clusters,
                         similarity = c("none", "association", "cosine",
                                        "jaccard", "inclusion",
                                        "equivalence"),
                         what = c("edges", "network")) {
  .thg_check_hg(hg)
  similarity <- match.arg(similarity)
  what <- match.arg(what)
  assignment <- .thg_labels_input(clusters)
  stopifnot(
    "`clusters` must be a data.frame or a named vector" =
      !is.null(names(assignment))
  )
  unknown <- setdiff(names(assignment), hg$nodes)
  if (length(unknown) > 0L) {
    stop(errorCondition(
      paste0("Unknown node names in `clusters`: ",
             paste(unknown, collapse = ", ")),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  groups <- factor(as.character(assignment)[match(hg$nodes,
                                                  names(assignment))])
  groups <- factor(groups, levels = .thg_kw_natural(levels(groups)))
  # topic x word document counts (full counting), then the co-occurrence
  # sum over words
  counts <- .thg_kw_aggregate(groups, (hg$incidence != 0) * 1)
  co <- tcrossprod(counts)
  total <- diag(co)
  total[total == 0] <- 1
  weight <- switch(
    similarity,
    none = co,
    association = co / outer(total, total),
    cosine = co / sqrt(outer(total, total)),
    jaccard = co / (outer(total, total, `+`) - co),
    inclusion = co / outer(total, total, pmin),
    equivalence = co^2 / outer(total, total)
  )
  weight[!is.finite(weight)] <- 0
  labels <- levels(groups)
  pairs <- which(upper.tri(weight), arr.ind = TRUE)
  pairs <- pairs[order(pairs[, 1L], pairs[, 2L]), , drop = FALSE]
  edges <- data.frame(
    source = labels[pairs[, 1L]],
    target = labels[pairs[, 2L]],
    weight = as.numeric(weight[pairs]),
    stringsAsFactors = FALSE
  )
  edges <- edges[edges$weight > 0, , drop = FALSE]
  rownames(edges) <- NULL
  if (identical(what, "edges")) {
    return(edges)
  }
  net <- cograph::as_cograph(edges, directed = FALSE)
  sizes <- table(groups)
  net$nodes$size <- as.integer(sizes[net$nodes$name])
  net
}

#' Transductive label spreading on a hypergraph, as a tidy table
#'
#' Calls the in-package [hypergraph_transduction()] engine (Zhou et al. 2006):
#' labels known for a few nodes spread over the hypergraph structure to
#' classify every node.
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hypergraph`).
#' @param labels The known labels: a named character vector (names are
#'   node identifiers -- documents under `nodes = "doc"` -- values their
#'   class labels), or a tidy data.frame with a `node` column and a
#'   `label`, `cluster` or `predicted` column.
#' @param xi,type Passed to [hypergraph_transduction()].
#' @param normalization Decision rule for turning spread scores into
#'   predictions: `"none"` (default, the raw Zhou 2006 argmax) or
#'   `"class_mass"` (class-mass normalization, Zhu et al. 2003). Use
#'   `"class_mass"` when the labeled seeds are class-imbalanced -- the raw
#'   rule can collapse every prediction onto the majority class.
#' @return A base `data.frame`, one row per node, with columns `node`,
#'   `label` (the given label or `NA`), `predicted`, `score`, and `margin`.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' hg_classify(hg, labels = c(cooking_1 = "cooking", space_1 = "space"))
#' @export
hg_classify <- function(hg, labels, xi = 0.99,
                        type = c("zhou", "random_walk"),
                        normalization = c("none", "class_mass")) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  normalization <- match.arg(normalization)
  labels <- .thg_labels_input(labels)
  fit <- if (.thg_is_sparse(hg)) {
    .thg_sparse_transduction(hg, labels = labels, xi = xi, type = type,
                             edge_weights = NULL,
                             normalization = normalization)
  } else {
    hypergraph_transduction(hg, labels = labels, xi = xi, type = type,
                            normalization = normalization)
  }
  out <- fit$predictions
  rownames(out) <- NULL
  out
}

# Long-form aliases for text-facing verbs whose hypergraph names are not
# already occupied by lower-level engines.
#' @rdname hg_keywords
#' @export
hypergraph_keywords <- hg_keywords

#' @rdname hg_relations
#' @export
hypergraph_relations <- hg_relations

#' @rdname hg_classify
#' @export
hypergraph_classify <- hg_classify
