# Hypergraph PageRank with edge-dependent vertex weights (Chitra & Raphael
# 2019), as a tidy verb with personalization and sparse support. The dense
# EDVW transition matrix is built here (verified at machine precision against
# the stationary distribution of the random-walk Laplacian in
# hypergraph_laplacian.R, itself HyperNetX-parity-tested) and against
# hypergraph_centrality(type = "pagerank"); PageRank adds damping and
# personalization on top.

# EDVW transition matrix: from vertex v, pick an incident hyperedge e with
# probability proportional to omega(e), then land on w in e with
# probability proportional to the incidence weight gamma_e(w). Row-
# stochastic by construction.
.thg_transition <- function(hg, edge_weights = NULL) {
  incidence <- hg$incidence
  membership <- (incidence > 0) * 1
  if (is.null(edge_weights)) {
    # Hayashi et al. (2020) heuristic, the HyperNetX default:
    # population SD of the edge's non-zero vertex weights, plus one.
    # Reduces to unit weights on a binary incidence.
    edge_weights <- apply(incidence, 2, \(col) {
      x <- col[col > 0]
      sqrt(mean((x - mean(x))^2)) + 1
    })
  } else {
    stopifnot(
      "`edge_weights` must be positive and one per hyperedge" =
        is.numeric(edge_weights) &&
        length(edge_weights) == ncol(incidence) &&
        all(is.finite(edge_weights)) && all(edge_weights > 0)
    )
  }
  delta <- colSums(incidence)
  vertex_degree <- as.numeric(membership %*% edge_weights)
  if (any(vertex_degree <= 0)) {
    stop(errorCondition(
      "every vertex must belong to at least one hyperedge",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  transition <- (membership %*% (t(incidence) * (edge_weights / delta))) /
    vertex_degree
  list(transition = transition, edge_weights = edge_weights)
}

#' Hypergraph PageRank with edge-dependent vertex weights
#'
#' Ranks the vertices of a weighted hypergraph by the stationary importance
#' of the Chitra & Raphael (2019) random walk: from a vertex, pick an
#' incident hyperedge with probability proportional to its edge weight,
#' then land on a member vertex with probability proportional to that
#' vertex's *edge-dependent* incidence weight (for a text hypergraph, e.g.,
#' the tf-idf of the word in the document). Damping and an optional
#' personalization vector turn the stationary distribution into PageRank
#' (Page et al. 1999): with probability `damping` the walker follows the
#' hypergraph walk, otherwise it teleports.
#'
#' Chitra & Raphael prove that with edge-*independent* vertex weights the
#' walk collapses to a random walk on a weighted clique-expansion graph —
#' so the hypergraph brings genuinely new information exactly when the
#' incidence weights differ across the hyperedges a vertex belongs to
#' (tested against the closed-form graph stationary distribution). With
#' `damping = 1` and default `edge_weights`, the result equals the
#' stationary distribution of the Hayashi et al. (2020) EDVW walk used by
#' [hg_cluster()] (tested at `1e-12` against [hypergraph_laplacian()]), and
#' with uniform teleportation it equals
#' `hypergraph_centrality(type = "pagerank")` (tested).
#'
#' @param hg A [text_hypergraph()], [knn_hypergraph()], or any hypernets
#'   `net_hypergraph` (connected when `damping = 1`).
#' @param damping Probability of following the hypergraph walk (default
#'   `0.85`); `1 - damping` is the teleport probability. Must be in
#'   `(0, 1]`; `damping = 1` gives the pure stationary distribution and
#'   requires a connected hypergraph to converge.
#' @param personalized Optional restart preference: a character vector of
#'   vertex names (uniform teleport over exactly those vertices), or a
#'   named non-negative vector of teleport weights over (a subset of) the
#'   vertex names; unnamed vertices get teleport probability 0. `NULL`
#'   (default) teleports uniformly.
#' @param edge_weights Positive hyperedge weights (one per hyperedge), or
#'   `NULL` (default) for the Hayashi et al. heuristic used by the
#'   Laplacian engines: the population standard deviation of each edge's
#'   non-zero vertex weights plus one, which reduces to unit weights on a
#'   binary incidence.
#' @param sort_by `NULL` (default, vertex order) or `"pagerank"` to sort
#'   descending (ties broken by vertex name).
#' @param n Return only the first `n` rows after sorting (default all).
#' @param max_iter,tol Power-iteration cap and L1 convergence tolerance.
#'
#' @return A base `data.frame`, one row per vertex, with columns `node` and
#'   `pagerank` (non-negative, summing to 1).
#'
#' @section Conditions: Raises `hypernets_bad_input` for broken contracts and
#'   warns with `hypernets_no_converge` (returning the last iterate) when
#'   `max_iter` is reached before `tol`.
#'
#' @references
#' Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
#' edge-dependent vertex weights. *ICML 2019*.
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. *CIKM 2020*.
#' \doi{10.1145/3340531.3412034}
#'
#' Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank
#' citation ranking: Bringing order to the web. Stanford InfoLab.
#'
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"),
#'    weight = "tfidf")
#' hg_pagerank(hg, sort_by = "pagerank")
#' hg_pagerank(hg, personalized = c(cooking_1 = 1), sort_by = "pagerank")
#' @export
hg_pagerank <- function(hg, damping = 0.85, personalized = NULL,
                        edge_weights = NULL, sort_by = NULL, n = Inf,
                        max_iter = 1000L, tol = 1e-12) {
  .thg_check_hg(hg)
  stopifnot(
    "`damping` must be a single value in (0, 1]" =
      length(damping) == 1L && is.finite(damping) &&
      damping > 0 && damping <= 1,
    "`n` must be a single count >= 1" =
      length(n) == 1L && (is.infinite(n) || (is.finite(n) && n >= 1)),
    "`max_iter` must be a single count >= 1" =
      length(max_iter) == 1L && is.finite(max_iter) && max_iter >= 1,
    "`tol` must be a single positive value" =
      length(tol) == 1L && is.finite(tol) && tol > 0
  )

  if (.thg_is_sparse(hg)) {
    ops <- .thg_walk_operators(hg, edge_weights = edge_weights)
    step <- ops$left
    ids <- rownames(hg$incidence)
  } else {
    walk <- .thg_transition(hg, edge_weights = edge_weights)
    step <- function(v) as.numeric(v %*% walk$transition)
    ids <- rownames(walk$transition)
  }

  if (is.character(personalized)) {
    personalized <- stats::setNames(rep(1, length(personalized)),
                                    personalized)
  }
  if (is.null(personalized)) {
    teleport <- rep(1 / length(ids), length(ids))
  } else {
    if (is.null(names(personalized)) || anyNA(personalized) ||
        !is.numeric(personalized) || any(personalized < 0) ||
        sum(personalized) <= 0 ||
        !all(names(personalized) %in% ids)) {
      stop(errorCondition(
        "`personalized` must be a named non-negative vector with a positive sum, over vertex names of `hg`",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    teleport <- rep(0, length(ids))
    teleport[match(names(personalized), ids)] <- personalized
    teleport <- teleport / sum(teleport)
  }

  rank <- teleport
  converged <- FALSE
  iter <- 0L
  # power iteration: each step depends on the previous iterate, so the loop
  # cannot be vectorized away; bounded by max_iter with convergence surfaced
  while (iter < max_iter) {
    iter <- iter + 1L
    updated <- damping * step(rank) + (1 - damping) * teleport
    if (sum(abs(updated - rank)) < tol) {
      rank <- updated
      converged <- TRUE
      break
    }
    rank <- updated
  }
  if (!converged) {
    warning(warningCondition(
      sprintf("PageRank did not converge in %d iterations (returning the last iterate)",
              as.integer(max_iter)),
      class = "hypernets_no_converge"
    ))
  }

  out <- data.frame(node = ids, pagerank = rank / sum(rank),
                    row.names = NULL)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, choices = "pagerank")
    out <- out[order(-out$pagerank, out$node), , drop = FALSE]
    rownames(out) <- NULL
  }
  if (is.finite(n) && n < nrow(out)) {
    out <- out[seq_len(n), , drop = FALSE]
  }
  out
}

#' @rdname hg_pagerank
#' @export
hypergraph_pagerank <- hg_pagerank
