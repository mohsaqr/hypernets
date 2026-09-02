# Tidy text-facing analysis verbs. Each delegates to the matching engine in
# the hypergraph family (hypergraph_measures, hypergraph_centrality,
# hypergraph_cluster, hypergraph_transduction) and returns a base data.frame.

.thg_check_hg <- function(hg) {
  if (!inherits(hg, "net_hypergraph")) {
    stop(errorCondition(
      "`hg` must be a net_hypergraph (text_hypergraph, knn_hypergraph, group_hypergraph, ...)",
      class = "honets_bad_input", call = NULL
    ))
  }
  invisible(hg)
}

#' Structural measures of a hypergraph, as tidy tables
#'
#' Delegates to [hypergraph_measures()] and returns the requested
#' slice as a tidy data.frame.
#'
#' @param hg A [text_hypergraph()] (or any honets `net_hypergraph`).
#' @param what Which table: `"nodes"` (default; one row per node with
#'   `hyperdegree`, `strength`, `max_edge_size`), `"edges"` (one row per
#'   hyperedge with its `size`), `"overlap"` (one row per hyperedge pair with
#'   `overlap`, `overlap_coefficient`, `jaccard`), or `"summary"` (one row per
#'   scalar measure).
#' @return A base `data.frame`, one row per node, edge, edge pair, or measure
#'   according to `what`.
#' @examples
#' hg <- text_hypergraph(c(
#'   a = "salt and soup and onions",
#'   b = "soup and salt",
#'   c = "stars and sky"
#' ))
#' hg_measures(hg)
#' hg_measures(hg, what = "summary")
#' @export
hg_measures <- function(hg, what = c("nodes", "edges", "overlap", "summary")) {
  .thg_check_hg(hg)
  what <- match.arg(what)
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
#' @param hg A [text_hypergraph()] (or any honets `net_hypergraph`).
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
      class = "honets_sparse_unsupported", call = NULL
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
#' @param hg A [text_hypergraph()] (or any honets `net_hypergraph`).
#' @param k Number of clusters (explicit by design; there is no correct
#'   default).
#' @param type `"zhou"` or `"random_walk"`, as in
#'   [hypergraph_cluster()].
#' @param seed Random seed passed to the engine's k-means step; set it for a
#'   reproducible partition.
#' @param nstart Number of k-means starts (default `25L`).
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
hg_cluster <- function(hg, k, type = c("zhou", "random_walk"), seed = NULL,
                       nstart = 25L,
                       what = c("clusters", "embedding", "eigenvalues")) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  what <- match.arg(what)
  fit <- if (.thg_is_sparse(hg)) {
    .thg_sparse_cluster(hg, k = k, type = type, edge_weights = NULL,
                        nstart = nstart, seed = seed)
  } else {
    hypergraph_cluster(hg, k = k, type = type, seed = seed,
                       nstart = nstart)
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

#' Characteristic hyperedges (keywords) per cluster
#'
#' For a clustered hypergraph, ranks each cluster's hyperedges by the
#' incidence weight mass its member nodes place on them. On a
#' `nodes = "doc"` [text_hypergraph()] this is per-topic keyword
#' extraction: a word's score in a cluster is the summed (tf-idf) weight
#' it receives from that cluster's documents, and `share` is the fraction
#' of the word's total corpus mass concentrated in the cluster --
#' `score` finds the cluster's heavy vocabulary, `share` its distinctive
#' vocabulary.
#'
#' @param hg The hypergraph the clustering was computed on.
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels.
#' @param n Keywords per cluster (default `10`); `Inf` returns all.
#' @param collapse If `TRUE`, return one row per cluster with its
#'   keywords joined into a single comma-separated string -- the compact
#'   display for comparing many clusters at once. Default `FALSE`.
#' @return A base `data.frame`. Default: one row per cluster-keyword
#'   pair, columns `cluster`, `rank`, `word` (the hyperedge name),
#'   `score` (in-cluster weight mass) and `share` (`score` divided by
#'   the hyperedge's total mass), ranked by `score` within cluster with
#'   ties broken alphabetically. With `collapse = TRUE`: one row per
#'   cluster, columns `cluster` and `words`.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_keywords(hg, topics, n = 3)
#' @export
hg_keywords <- function(hg, clusters, n = 10L, collapse = FALSE) {
  .thg_check_hg(hg)
  stopifnot(
    "`n` must be a single positive number" =
      length(n) == 1L && is.numeric(n) && n >= 1,
    "`collapse` must be TRUE or FALSE" =
      isTRUE(collapse) || isFALSE(collapse)
  )
  assignment <- if (is.data.frame(clusters)) {
    stopifnot(
      "`clusters` needs `node` and `cluster` columns (from hg_cluster())" =
        all(c("node", "cluster") %in% names(clusters))
    )
    stats::setNames(clusters$cluster, clusters$node)
  } else {
    stopifnot(
      "`clusters` must be a data.frame or a named vector" =
        !is.null(names(clusters))
    )
    clusters
  }
  unknown <- setdiff(names(assignment), hg$nodes)
  if (length(unknown) > 0L) {
    stop(errorCondition(
      paste0("Unknown node names in `clusters`: ",
             paste(unknown, collapse = ", ")),
      class = "honets_bad_input", call = NULL
    ))
  }
  groups <- factor(as.character(assignment)[match(hg$nodes,
                                                  names(assignment))])
  assigned <- which(!is.na(groups))
  indicator <- Matrix::sparseMatrix(
    i = as.integer(groups[assigned]), j = assigned, x = 1,
    dims = c(nlevels(groups), length(hg$nodes))
  )
  mass <- as.matrix(indicator %*% hg$incidence)
  dimnames(mass) <- list(levels(groups), colnames(hg$incidence))
  total <- as.numeric(Matrix::colSums(hg$incidence))
  per_cluster <- lapply(levels(groups), \(cl) {
    scores <- mass[cl, ]
    ord <- order(-scores, names(scores))
    keep <- utils::head(ord[scores[ord] > 0], n)
    data.frame(
      cluster = cl, rank = seq_along(keep), word = names(scores)[keep],
      score = as.numeric(scores[keep]),
      share = as.numeric(scores[keep] / total[keep]),
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, per_cluster)
  rownames(out) <- NULL
  if (isTRUE(collapse)) {
    out <- stats::aggregate(word ~ cluster, data = out, paste,
                            collapse = ", ")
    names(out) <- c("cluster", "words")
  }
  out
}

#' Transductive label spreading on a hypergraph, as a tidy table
#'
#' Calls the in-package [hypergraph_transduction()] engine (Zhou et al. 2006):
#' labels known for a few nodes spread over the hypergraph structure to
#' classify every node.
#'
#' @param hg A [text_hypergraph()] (or any honets `net_hypergraph`).
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

#' @rdname hg_classify
#' @export
hypergraph_classify <- hg_classify
