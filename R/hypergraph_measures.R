# ---- Hypergraph structural measures (HON-9) ------------------------------
# Node-level, hyperedge-level, and global statistics for a net_hypergraph.
# All measures are fast matrix operations on the incidence matrix B.

#' Structural measures for a hypergraph
#'
#' Computes a comprehensive structural-statistics suite for a
#' [net_hypergraph][build_hypergraph]: node-level, hyperedge-level, and
#' global measures. All measures are derived in a few BLAS calls on the
#' incidence matrix.
#'
#' @param hg A `net_hypergraph` (from [build_hypergraph()] or
#'   [group_hypergraph()]).
#'
#' @return An object of class `hypergraph_measures` (a named list) with
#'   components:
#'
#' \describe{
#'   \item{**Node-level (length `n_nodes`)**}{}
#'   \item{`hyperdegree`}{Number of hyperedges containing each node.}
#'   \item{`node_strength`}{Total participation: for node \eqn{i},
#'     \eqn{\sum_{e \ni i} |e|}. A node in many large hyperedges has
#'     high strength.}
#'   \item{`max_edge_size`}{Size of the largest hyperedge containing
#'     each node.}
#'   \item{`co_degree`}{`n_nodes` x `n_nodes` matrix:
#'     `co_degree[i, j] = |{e : i, j in e}|` (number of hyperedges
#'     co-containing nodes \eqn{i} and \eqn{j}). Diagonal is zero.}
#'   \item{**Hyperedge-level (length `n_hyperedges` or m x m)**}{}
#'   \item{`edge_sizes`}{Hyperedge sizes \eqn{|e|}.}
#'   \item{`edge_pairwise_overlap`}{`m` x `m` matrix:
#'     `|e_i intersect e_j|`. Diagonal is zero.}
#'   \item{`overlap_coefficient`}{`m` x `m`:
#'     `|e_i and e_j| / min(|e_i|, |e_j|)`. Measures how much the
#'     smaller hyperedge is contained in the larger.}
#'   \item{`jaccard`}{`m` x `m`: symmetric overlap index
#'     `|e_i and e_j| / |e_i union e_j|`.}
#'   \item{**Global (scalars)**}{}
#'   \item{`density`}{For `k`-uniform hypergraphs: `m / choose(n, k)`.
#'     For mixed sizes: `sum(|e|) / (n * m)` (mean fraction of nodes
#'     per hyperedge).}
#'   \item{`avg_edge_size`}{Mean of `edge_sizes`.}
#'   \item{`size_distribution`}{Tabulation of hyperedge sizes (passed
#'     through from `hg`).}
#'   \item{`intersection_profile`}{Distribution of pairwise hyperedge
#'     intersection sizes - useful for spotting whether hyperedges
#'     overlap mostly trivially or share substantial cores
#'     (Do et al. 2020).}
#'   \item{`pairwise_participation`}{Fraction of node pairs co-appearing
#'     in at least one hyperedge.}
#'   \item{`n_nodes`, `n_hyperedges`}{Convenience scalars.}
#' }
#'
#' @details
#' All measures are computed via standard matrix operations on the binary
#' incidence \eqn{B = (b_{ij})} where \eqn{b_{ij} = 1} iff node \eqn{i}
#' is in hyperedge \eqn{j}:
#' \itemize{
#'   \item `hyperdegree = rowSums(B)`,
#'         `edge_sizes = colSums(B)`
#'   \item `co_degree = tcrossprod(B)` (with zero diagonal)
#'   \item `edge_pairwise_overlap = crossprod(B)` (with zero diagonal)
#'   \item `overlap_coefficient[i, j] = overlap[i, j] /
#'         min(edge_sizes[i], edge_sizes[j])`
#'   \item `jaccard[i, j] = overlap[i, j] /
#'         (edge_sizes[i] + edge_sizes[j] - overlap[i, j])`
#' }
#'
#' Empty hypergraph (`n_hyperedges == 0`) returns trivial zeros and
#' empty matrices.
#'
#' @seealso [build_hypergraph()], [group_hypergraph()],
#'   [clique_expansion()].
#'
#' @examples
#' df <- data.frame(
#'   person  = c("A", "B", "C", "A", "B", "D", "C", "D"),
#'   session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3")
#' )
#' hg <- group_hypergraph(df, "person", "session")
#' m  <- hypergraph_measures(hg)
#' m
#'
#' # tidy accessors: per node, per hyperedge, and the whole hypergraph
#' as.data.frame(m, sort_by = "hyperdegree")
#' as.data.frame(m, what = "edges")
#' as.data.frame(m, what = "global")
#'
#' @references
#' Lee, G., Bu, F., Eliassi-Rad, T., & Shin, K. (2025). A survey on
#' hypergraph mining: patterns, tools, and generators. \emph{ACM Computing
#' Surveys}, 57(8), 203. \doi{10.1145/3719002}
#'
#' Do, M. T., Yoon, S., Hooi, B., & Shin, K. (2020). Structural patterns
#' and generative models of real-world hypergraphs. arXiv:2006.07060.
#'
#' @export
hypergraph_measures <- function(hg) {
  stopifnot(inherits(hg, "net_hypergraph"))

  B <- hg$incidence
  n <- hg$n_nodes
  m <- hg$n_hyperedges

  if (m == 0L) {
    zeros_n <- stats::setNames(numeric(n), hg$nodes)
    empty_mm <- matrix(0, 0L, 0L)
    return(structure(
      list(
        hyperdegree            = stats::setNames(integer(n), hg$nodes),
        node_strength          = zeros_n,
        max_edge_size          = stats::setNames(integer(n), hg$nodes),
        co_degree              = matrix(0, n, n,
                                        dimnames = list(hg$nodes, hg$nodes)),
        edge_sizes             = integer(0),
        edge_pairwise_overlap  = empty_mm,
        overlap_coefficient    = empty_mm,
        jaccard                = empty_mm,
        density                = 0,
        avg_edge_size          = NA_real_,
        size_distribution      = hg$size_distribution,
        intersection_profile   = integer(0),
        pairwise_participation = 0,
        n_nodes                = n,
        n_hyperedges           = 0L
      ),
      class = "net_hypergraph_measures"
    ))
  }

  B_bin <- (B > 0) * 1.0
  edge_sizes <- as.integer(colSums(B_bin))
  edge_names <- colnames(B) %||% paste0("h", seq_len(m))

  # ---- Node-level ----
  hyperdegree <- stats::setNames(as.integer(rowSums(B_bin)), hg$nodes)
  node_strength <- stats::setNames(as.numeric(B_bin %*% edge_sizes), hg$nodes)
  max_edge_size <- stats::setNames(
    vapply(seq_len(n), function(i) {
      mem <- B_bin[i, ] > 0
      if (any(mem)) max(edge_sizes[mem]) else 0L
    }, integer(1L)),
    hg$nodes
  )

  # ---- Pairwise co-degree (n x n) ----
  co_degree <- tcrossprod(B_bin)
  diag(co_degree) <- 0
  dimnames(co_degree) <- list(hg$nodes, hg$nodes)

  # ---- Edge-level pairwise overlap (m x m) ----
  edge_pairwise_overlap <- crossprod(B_bin)
  diag(edge_pairwise_overlap) <- 0
  dimnames(edge_pairwise_overlap) <- list(edge_names, edge_names)

  # ---- Overlap coefficient & Jaccard (m x m) ----
  size_min   <- outer(edge_sizes, edge_sizes, pmin)
  size_union <- outer(edge_sizes, edge_sizes, `+`) - edge_pairwise_overlap
  overlap_coefficient <- ifelse(size_min > 0,
                                edge_pairwise_overlap / size_min, 0)
  jaccard <- ifelse(size_union > 0,
                    edge_pairwise_overlap / size_union, 0)
  diag(overlap_coefficient) <- 0
  diag(jaccard) <- 0
  dimnames(overlap_coefficient) <- list(edge_names, edge_names)
  dimnames(jaccard) <- list(edge_names, edge_names)

  # ---- Global ----
  is_uniform <- length(unique(edge_sizes)) == 1L
  density <- if (is_uniform) {
    k <- edge_sizes[1L]
    if (k > n) 0 else m / choose(n, k)
  } else {
    sum(edge_sizes) / (n * m)
  }
  avg_edge_size <- mean(edge_sizes)

  intersection_profile <- if (m >= 2L) {
    upper <- edge_pairwise_overlap[upper.tri(edge_pairwise_overlap)]
    tab <- table(upper)
    out <- as.integer(tab)
    names(out) <- paste0("overlap_", names(tab))
    out
  } else {
    integer(0L)
  }

  pairwise_participation <- if (n >= 2L) {
    sum(co_degree[upper.tri(co_degree)] > 0) / choose(n, 2L)
  } else {
    0
  }

  structure(
    list(
      hyperdegree            = hyperdegree,
      node_strength          = node_strength,
      max_edge_size          = max_edge_size,
      co_degree              = co_degree,
      edge_sizes             = edge_sizes,
      edge_pairwise_overlap  = edge_pairwise_overlap,
      overlap_coefficient    = overlap_coefficient,
      jaccard                = jaccard,
      density                = density,
      avg_edge_size          = avg_edge_size,
      size_distribution      = hg$size_distribution,
      intersection_profile   = intersection_profile,
      pairwise_participation = pairwise_participation,
      n_nodes                = n,
      n_hyperedges           = m
    ),
    class = "net_hypergraph_measures"
  )
}

#' @param x A `hypergraph_measures` object.
#' @param ... Additional arguments (ignored).
#' @return The input `x` invisibly.
#' @rdname hypergraph_measures
#' @export
print.net_hypergraph_measures <- function(x, ...) {
  cat(sprintf("Hypergraph measures: %d nodes, %d hyperedges\n",
              x$n_nodes, x$n_hyperedges))
  cat(sprintf("  Density:                %.4f\n", x$density))
  cat(sprintf("  Mean edge size:         %.2f\n",
              if (length(x$avg_edge_size)) x$avg_edge_size else NA))
  cat(sprintf("  Pairwise participation: %.4f\n",
              x$pairwise_participation))
  if (length(x$hyperdegree)) {
    cat(sprintf("\n  Hyperdegree:    min=%d  med=%.1f  max=%d\n",
                min(x$hyperdegree),
                stats::median(x$hyperdegree),
                max(x$hyperdegree)))
    cat(sprintf("  Node strength:  min=%.1f  med=%.1f  max=%.1f\n",
                min(x$node_strength),
                stats::median(x$node_strength),
                max(x$node_strength)))
  }
  if (length(x$edge_sizes)) {
    cat(sprintf("  Edge sizes:     min=%d  med=%.1f  max=%d\n",
                min(x$edge_sizes),
                stats::median(x$edge_sizes),
                max(x$edge_sizes)))
  }
  if (length(x$intersection_profile)) {
    cat("\nIntersection profile (pair counts by overlap size):\n")
    for (nm in names(x$intersection_profile)) {
      cat(sprintf("  %-12s: %d\n", nm, x$intersection_profile[[nm]]))
    }
  }
  invisible(x)
}

#' Coerce a net_hypergraph_measures to a tidy table
#'
#' @param x A `net_hypergraph_measures` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"nodes"` (default) for the per-node measures, `"edges"` for
#'   the per-hyperedge measures, or `"global"` for the whole-hypergraph
#'   summary.
#' @param sort_by `NULL` (node/edge order, default), or a column name to
#'   sort by, largest first. Ignored when `what = "global"`.
#' @return A data.frame. For `what = "nodes"`, one row per node: `node`,
#'   `hyperdegree`, `node_strength`, `max_edge_size`. For `what = "edges"`,
#'   one row per hyperedge: `hyperedge`, `size`. For `what = "global"`, one
#'   row per statistic: `measure`, `value`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_hypergraph_measures <- function(
    x, row.names = NULL, optional = FALSE, ...,
    what = c("nodes", "edges", "global"), sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  out <- switch(
    what,
    nodes = data.frame(
      node          = names(x$hyperdegree) %||%
                        as.character(seq_along(x$hyperdegree)),
      hyperdegree   = as.integer(x$hyperdegree),
      node_strength = as.numeric(x$node_strength),
      max_edge_size = as.integer(x$max_edge_size),
      stringsAsFactors = FALSE
    ),
    edges = data.frame(
      hyperedge = names(x$edge_sizes) %||%
                    paste0("h", seq_along(x$edge_sizes)),
      size      = as.integer(x$edge_sizes),
      stringsAsFactors = FALSE
    ),
    global = data.frame(
      measure = c("n_nodes", "n_hyperedges", "density", "avg_edge_size",
                  "pairwise_participation"),
      value   = c(x$n_nodes, x$n_hyperedges, x$density, x$avg_edge_size,
                  x$pairwise_participation),
      stringsAsFactors = FALSE
    )
  )
  if (!is.null(sort_by) && what != "global") {
    sort_by <- match.arg(sort_by, setdiff(names(out), c("node", "hyperedge")))
    key <- if (what == "nodes") out$node else out$hyperedge
    out <- out[order(-out[[sort_by]], key), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
