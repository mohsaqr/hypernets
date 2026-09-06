# ---- Hypergraph construction (EG-4) -------------------------------------
# Build a hypergraph from a network's clique structure, optionally retaining
# pairwise edges. Foundation for higher-order analyses (centrality, walks,
# contagion, public-goods games).

#' Higher-order hypergraph from a network's clique structure
#'
#' Takes a network and produces a hypergraph by promoting k-cliques (k >= 3)
#' to k-hyperedges. Each k-clique is independently included as a k-hyperedge
#' with probability `p`. Optionally retains the underlying pairwise edges as
#' 2-hyperedges. Foundation for higher-order analyses.
#'
#' @param net A `netobject`, `cograph_network`, `net_simplicial`, or
#'   numeric adjacency / weight matrix. Directed inputs are symmetrised by
#'   the underlying clique enumerator.
#' @param p Probability in `[0, 1]` that each k-clique with k >= 3 becomes a
#'   k-hyperedge. Default `1` (deterministic - every found clique is
#'   included).
#' @param type Hyperedge enumeration. `"clique"` (default) promotes
#'   k-cliques in the binarised adjacency to k-hyperedges. A metric
#'   Vietoris-Rips construction is \strong{not implemented}; `"vr"` /
#'   `"rips"` are accepted by `match.arg` but raise an error rather than
#'   silently aliasing `"clique"`.
#' @param include_pairwise Logical. Include 2-edges from the input network as
#'   2-hyperedges. Default `TRUE`. Set `FALSE` for a "fully higher-order"
#'   hypergraph containing only k-hyperedges with k >= 3.
#' @param max_size Integer >= 2. Maximum hyperedge size to extract. Default
#'   `3L` (triangles only). `4L` also includes 4-cliques as 4-hyperedges, etc.
#' @param threshold Numeric. Edge weight cutoff used to binarise the
#'   adjacency for clique enumeration. Default `0` (any non-zero weight is
#'   an edge).
#' @param seed Optional integer for reproducible Bernoulli sampling when
#'   `0 < p < 1`.
#'
#' @return A `net_hypergraph` object: a list with components
#' \describe{
#'   \item{`hyperedges`}{List of integer vectors. Each entry is a hyperedge
#'     given as the sorted node indices it spans.}
#'   \item{`incidence`}{Numeric matrix of size `n_nodes` x `n_hyperedges`.
#'     `incidence[i, j] = 1` iff node i belongs to hyperedge j. Row names
#'     are node names; column names are `h1`, `h2`, ...}
#'   \item{`nodes`}{Character vector of node names.}
#'   \item{`n_nodes`, `n_hyperedges`}{Scalar counts.}
#'   \item{`size_distribution`}{Named integer vector: number of hyperedges
#'     of each size, named `size_2`, `size_3`, ...}
#'   \item{`params`}{Recorded call parameters: `source`, `type`, `p`,
#'     `include_pairwise`, `max_size`, `threshold`, `seed`.}
#' }
#'
#' @details
#' The construction follows Burgio, Matamalas, Gomez & Arenas (2020) on
#' simplicial / hypergraph contagion. For each k-clique with k >= 3 found in
#' the underlying graph (via [build_simplicial()]), an independent
#' Bernoulli(`p`) trial decides whether that clique becomes a k-hyperedge.
#' Underlying pairwise edges are always retained when
#' `include_pairwise = TRUE`, so the resulting hypergraph contains both the
#' original 2-edge structure and the sampled higher-order interactions.
#'
#' At the limits:
#' \itemize{
#'   \item `p = 0` with `include_pairwise = TRUE` reproduces the input
#'     pairwise network as a hypergraph of size-2 edges.
#'   \item `p = 1` with `include_pairwise = FALSE` returns a fully
#'     higher-order hypergraph containing only the k-hyperedges (k >= 3)
#'     found in the network's clique complex.
#' }
#'
#' @seealso [build_simplicial()] (underlying clique enumeration),
#'   \code{Nestimate::build_network()}.
#'
#' @examples
#' set.seed(1)
#' n <- 8
#' adj <- matrix(stats::rbinom(n * n, 1, 0.5), n, n)
#' diag(adj) <- 0
#' adj <- (adj + t(adj)) > 0
#' rownames(adj) <- colnames(adj) <- LETTERS[seq_len(n)]
#' hg <- build_hypergraph(adj, p = 1, max_size = 3L)
#' print(hg)
#' summary(hg)
#'
#' @references
#' Burgio, G., Matamalas, J. T., Gomez, S., & Arenas, A. (2020). Evolution
#' of cooperation in the presence of higher-order interactions: from
#' networks to hypergraphs. \emph{Entropy} 22(7), 744.
#' \doi{10.3390/e22070744}
#'
#' @export
build_hypergraph <- function(net,
                              p = 1,
                              type = c("clique", "vr", "rips"),
                              include_pairwise = TRUE,
                              max_size = 3L,
                              threshold = 0,
                              seed = NULL) {
  type <- match.arg(type)
  if (type == "vr" || type == "rips") {
    stop("type = \"vr\" (Vietoris-Rips) is not implemented. A genuine ",
         "metric Vietoris-Rips construction requires a distance-ball ",
         "graph, which this package does not build. Use ",
         "type = \"clique\" (cliques in the binarised adjacency).",
         call. = FALSE)
  }
  stopifnot(
    is.numeric(p), length(p) == 1L, p >= 0, p <= 1,
    is.logical(include_pairwise), length(include_pairwise) == 1L,
    is.numeric(max_size), length(max_size) == 1L, max_size >= 2L,
    is.numeric(threshold), length(threshold) == 1L
  )
  max_size <- as.integer(max_size)
  if (!is.null(seed)) set.seed(as.integer(seed))

  # ---- Get adjacency + node names from any supported input ------------
  parsed <- .hg_extract_adj(net)
  adj    <- parsed$adj
  nodes  <- parsed$nodes
  n      <- length(nodes)

  # ---- Find all simplices (cliques) up to max_size --------------------
  sc <- build_simplicial(adj, type = type, threshold = threshold,
                         max_dim = max_size - 1L)
  simplices <- sc$simplices

  sizes   <- vapply(simplices, length, integer(1L))
  edges_2 <- simplices[sizes == 2L]

  # k>=3 simplices: sampled with prob p
  hi_simp <- simplices[sizes >= 3L & sizes <= max_size]
  sampled_hi <- if (length(hi_simp) == 0L || p <= 0) {
    list()
  } else if (p >= 1) {
    hi_simp
  } else {
    keep <- as.logical(stats::rbinom(length(hi_simp), 1L, p))
    hi_simp[keep]
  }

  hyperedges <- if (include_pairwise) c(edges_2, sampled_hi) else sampled_hi
  m <- length(hyperedges)

  # ---- Build incidence matrix [n_nodes x n_hyperedges] ----------------
  if (m == 0L) {
    incidence <- matrix(0L, nrow = n, ncol = 0L,
                        dimnames = list(nodes, NULL))
  } else {
    incidence <- matrix(0L, nrow = n, ncol = m,
                        dimnames = list(nodes, paste0("h", seq_len(m))))
    # One (node, hyperedge) cell per membership. Hyperedges hold sorted
    # distinct node indices, so no cell is written twice.
    incidence[cbind(unlist(hyperedges, use.names = FALSE),
                    rep.int(seq_len(m), lengths(hyperedges)))] <- 1L
  }

  # ---- Size distribution ----------------------------------------------
  he_sizes <- vapply(hyperedges, length, integer(1L))
  size_dist <- if (length(he_sizes)) {
    tab <- table(he_sizes)
    out <- as.integer(tab)
    names(out) <- paste0("size_", names(tab))
    out
  } else {
    integer(0L)
  }

  structure(
    list(
      hyperedges        = hyperedges,
      incidence         = incidence,
      nodes             = nodes,
      n_nodes           = n,
      n_hyperedges      = m,
      size_distribution = size_dist,
      params = list(
        source           = "build_hypergraph",
        type             = type,
        p                = p,
        include_pairwise = include_pairwise,
        max_size         = max_size,
        threshold        = threshold,
        seed             = seed
      )
    ),
    class = "net_hypergraph"
  )
}

# ---- Internal helpers --------------------------------------------------

#' Extract symmetric numeric adjacency + node names from any supported input
#' @noRd
.hg_extract_adj <- function(net) {
  if (inherits(net, "net_simplicial")) {
    nodes <- net$nodes
    n     <- length(nodes)
    adj   <- matrix(0, n, n, dimnames = list(nodes, nodes))
    edge_simp <- Filter(function(s) length(s) == 2L, net$simplices)
    for (e in edge_simp) {
      adj[e[1L], e[2L]] <- 1
      adj[e[2L], e[1L]] <- 1
    }
    return(list(adj = adj, nodes = nodes))
  }

  if (inherits(net, "netobject") || inherits(net, "cograph_network")) {
    w <- net$weights
    nm <- if (!is.null(net$nodes$name)) {
      net$nodes$name
    } else if (!is.null(net$nodes$label)) {
      net$nodes$label
    } else if (!is.null(net$nodes$id)) {
      as.character(net$nodes$id)
    } else if (!is.null(rownames(w))) {
      rownames(w)
    } else {
      paste0("V", seq_len(nrow(w)))
    }
    rownames(w) <- colnames(w) <- nm
    return(list(adj = w, nodes = nm))
  }

  if (is.matrix(net) && (is.numeric(net) || is.logical(net))) {
    nm <- rownames(net) %||% paste0("V", seq_len(nrow(net)))
    storage.mode(net) <- "double"
    rownames(net) <- colnames(net) <- nm
    return(list(adj = net, nodes = nm))
  }

  stop("`net` must be a netobject, cograph_network, simplicial_complex, ",
       "or numeric matrix.", call. = FALSE)
}

# ---- S3 methods ---------------------------------------------------------

#' @param x A `net_hypergraph` object (for `print`).
#' @param object A `net_hypergraph` object (for `summary`).
#' @param ... Additional arguments (ignored).
#' @return The input `x` invisibly.
#' @rdname build_hypergraph
#' @export
print.net_hypergraph <- function(x, ...) {
  cat(sprintf("Hypergraph: %d nodes, %d hyperedges\n",
              x$n_nodes, x$n_hyperedges))
  if (length(x$size_distribution)) {
    cat("Size distribution:\n")
    for (nm in names(x$size_distribution)) {
      cat(sprintf("  %-8s : %d\n", nm, x$size_distribution[[nm]]))
    }
  }
  if (identical(x$params$source, "window_hypergraph")) {
    cat(sprintf(
      "Source: windowed sequences, window = %d, step = %d (%d windows from %d sequences)\n",
      x$params$window, x$params$step,
      x$params$n_windows, x$params$n_sequences))
  } else if (identical(x$params$source, "group_hypergraph")) {
    cat(sprintf("Source: group membership (member = %s, group = %s)\n",
                x$params$member, x$params$group))
  } else if (startsWith(x$params$source %||% "", "hg_sample_")) {
    cat(sprintf("Source: random %s model\n", x$params$model))
  } else {
    cat(sprintf(
      "Source: network %ss (p = %.2f, include_pairwise = %s, max_size = %d)\n",
      x$params$type, x$params$p,
      x$params$include_pairwise, x$params$max_size))
  }
  invisible(x)
}

#' Coerce a net_hypergraph to a tidy hyperedge table
#'
#' One row per hyperedge: its identifier, size (number of distinct member
#' states), the member states, and its weight (for hypergraphs built by
#' [window_hypergraph()], the number of windows collapsed into the
#' hyperedge; `NA` for the other constructors, whose hyperedges are
#' unweighted).
#'
#' @param x A `net_hypergraph` object.
#' @param row.names Ignored (present for S3 consistency with the generic).
#' @param optional Ignored (present for S3 consistency with the generic).
#' @param ... Additional arguments (ignored).
#' @param sort_by `NULL` (construction order, default), `"weight"`, or
#'   `"size"` - sort the table by that column, largest first (ties broken
#'   by hyperedge id so the order is deterministic).
#' @param what `"edges"` (default) for one row per hyperedge, `"nodes"`
#'   for one row per node, or `"memberships"` for one row per node-in-
#'   hyperedge cell of the incidence matrix.
#' @return A data.frame. For `what = "edges"`, one row per hyperedge with
#'   columns `hyperedge` (character id), `size` (integer), `states`
#'   (comma-separated member states), and `weight` (numeric window count,
#'   or `NA`). For `what = "nodes"`, one row per node with columns `node`
#'   and `degree` (the number of hyperedges it belongs to), plus `block`
#'   for a hypergraph with planted blocks ([hg_sample_sbm()]); `sort_by =
#'   "degree"` orders it. For `what = "memberships"`, one row per non-zero
#'   incidence cell with columns `node`, `hyperedge` and `weight` (the
#'   incidence value: 1 for a binary hypergraph, the summed weight for a
#'   weighted one), in hyperedge order; `sort_by = "weight"` orders it.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' hg <- window_hypergraph(list(s1 = c("a", "b", "a", "c")), window = 2L)
#' as.data.frame(hg)
#' as.data.frame(hg, sort_by = "weight")
#' as.data.frame(hg, sort_by = "weight", top = 3)
#' as.data.frame(hg, what = "nodes", sort_by = "degree")
#' as.data.frame(hg, what = "memberships")
#' @export
as.data.frame.net_hypergraph <- function(x, row.names = NULL,
                                         optional = FALSE, ...,
                                         what = c("edges", "nodes",
                                                  "memberships"),
                                         sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  if (identical(what, "memberships")) {
    nodes <- x$nodes %||% rownames(x$incidence) %||%
      paste0("n", seq_len(x$n_nodes))
    hyperedges <- colnames(x$incidence) %||%
      paste0("h", seq_len(x$n_hyperedges))
    cells <- Matrix::which(x$incidence != 0, arr.ind = TRUE)
    out <- data.frame(
      node = as.character(nodes[cells[, 1L]]),
      hyperedge = as.character(hyperedges[cells[, 2L]]),
      weight = as.numeric(x$incidence[cells]),
      stringsAsFactors = FALSE
    )
    out <- out[order(cells[, 2L], cells[, 1L]), , drop = FALSE]
    if (!is.null(sort_by)) {
      sort_by <- match.arg(sort_by, "weight")
      out <- out[order(-out$weight, out$hyperedge, out$node), , drop = FALSE]
    }
    rownames(out) <- NULL
    return(.ho_top(out, top))
  }
  if (identical(what, "nodes")) {
    nodes <- x$nodes %||% rownames(x$incidence) %||%
      paste0("n", seq_len(x$n_nodes))
    degree <- if (x$n_hyperedges > 0L) {
      as.integer(Matrix::rowSums(x$incidence > 0))
    } else {
      rep(0L, x$n_nodes)
    }
    out <- data.frame(node = as.character(nodes), degree = degree,
                      stringsAsFactors = FALSE)
    if (!is.null(x$blocks)) out$block <- as.integer(x$blocks)
    if (!is.null(sort_by)) {
      sort_by <- match.arg(sort_by, "degree")
      out <- out[order(-out$degree, out$node), , drop = FALSE]
      rownames(out) <- NULL
    }
    return(.ho_top(out, top))
  }
  out <- data.frame(
    hyperedge = colnames(x$incidence) %||%
      paste0("h", seq_len(x$n_hyperedges)),
    size = lengths(x$hyperedges),
    states = vapply(x$hyperedges,
                    function(idx) paste(x$nodes[idx], collapse = ", "),
                    character(1L)),
    weight = as.numeric(x$window_counts %||% rep(NA_real_, x$n_hyperedges)),
    stringsAsFactors = FALSE
  )
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("weight", "size"))
    out <- out[order(-out[[sort_by]], out$hyperedge), , drop = FALSE]
    rownames(out) <- NULL
  }
  .ho_top(out, top)
}

#' @return A data.frame, one row per node: `node`, `degree` (the number of
#'   hyperedges the node belongs to). Returned **invisibly**: `summary(x)`
#'   prints the summary and nothing else; assign the result to keep the table.
#' @rdname build_hypergraph
#' @export
summary.net_hypergraph <- function(object, ...) {
  he_sizes <- vapply(object$hyperedges, length, integer(1L))
  cat("Hypergraph summary\n")
  cat(sprintf("  Nodes:         %d\n", object$n_nodes))
  cat(sprintf("  Hyperedges:    %d\n", object$n_hyperedges))
  if (length(he_sizes)) {
    cat(sprintf("  Mean size:     %.2f\n", mean(he_sizes)))
    cat(sprintf("  Max size:      %d\n", max(he_sizes)))
  }

  nodes <- if (!is.null(object$nodes)) object$nodes else
    rownames(object$incidence)
  if (is.null(nodes)) nodes <- paste0("n", seq_len(object$n_nodes))
  if (object$n_hyperedges > 0L) {
    deg <- as.integer(rowSums(object$incidence > 0))
  } else {
    deg <- rep(0L, object$n_nodes)
  }
  invisible(data.frame(node = as.character(nodes), degree = deg,
                       stringsAsFactors = FALSE, row.names = NULL))
}
