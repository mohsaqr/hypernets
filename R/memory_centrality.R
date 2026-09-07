# ---- Higher-order centralities with first-order projection ---------------
#
# Centralities computed on the higher-order topology of a net_hon and
# projected back to first-order states (Scholtes, Wider & Garas 2016).
#
# Semantics follow pathpy 2.2.0 (Scholtes' own reference implementation,
# pathpy/algorithms/centralities.py), generalized from its fixed-order
# networks to the variable-order networks build_hon() produces:
#
#   * a higher-order node "a -> b" stands for the first-order path
#     (a, b); its ORDER is the number of states in that path.
#   * a higher-order path (v1, ..., vl) maps to the first-order path
#     "full path of v1, then the last state of each subsequent node"
#     (pathpy: higher_order_path_to_first_order).
#   * first-order distance between states = higher-order hop count plus
#     (order of the source node - 1). pathpy adds the constant `k - 1`
#     because all its nodes share one order; with mixed orders the term
#     is per-node, which reduces to pathpy's formula when the orders
#     coincide.
#
# Graph traversal (BFS, shortest-path enumeration) is inherently
# sequential, so the loops in this file are the justified exception to
# the vectorize-everything rule.

# ---------------------------------------------------------------------------
# Topology helpers
# ---------------------------------------------------------------------------

#' Split higher-order node names into their first-order state paths
#'
#' @param nodes Character vector of node names ("a", "a -> b", ...).
#' @return List of character vectors, one per node.
#' @noRd
.hoc_node_paths <- function(nodes) {
  strsplit(nodes, " -> ", fixed = TRUE)
}

#' Successor adjacency list of a higher-order network
#'
#' @param mat Square weight matrix of the HON.
#' @return List of integer vectors, one per node.
#' @noRd
.hoc_succ <- function(mat) {
  lapply(seq_len(nrow(mat)), function(i) which(mat[i, ] != 0))
}

#' BFS distances from one source over an adjacency list
#'
#' Expands one whole frontier per step, so the loop runs once per BFS
#' level rather than once per node.
#'
#' @param succ Successor list.
#' @param s Integer source index.
#' @param n Integer node count.
#' @return Numeric vector of hop distances (`Inf` when unreachable).
#' @noRd
.hoc_bfs_dist <- function(succ, s, n) {
  dist <- rep(Inf, n)
  dist[s] <- 0
  frontier <- s
  d <- 0
  while (length(frontier) > 0L) {
    d <- d + 1
    nxt <- unique(unlist(succ[frontier], use.names = FALSE))
    nxt <- nxt[is.infinite(dist[nxt])]
    dist[nxt] <- d
    frontier <- nxt
  }
  dist
}

#' All shortest paths from one source, as lists of node-index vectors
#'
#' BFS layering gives the shortest-path predecessor DAG; paths are then
#' expanded from it. The expansion is bounded by `max_paths`: exceeding
#' it raises a classed error rather than silently truncating.
#'
#' @param succ Successor list.
#' @param s Integer source index.
#' @param n Integer node count.
#' @param max_paths Integer cap on enumerated paths.
#' @return List indexed by target node; each element a list of integer
#'   vectors (the shortest paths), or NULL when unreachable.
#' @noRd
.hoc_shortest_paths_from <- function(succ, s, n, max_paths) {
  dist <- .hoc_bfs_dist(succ, s, n)
  reach <- which(is.finite(dist))
  # predecessors on shortest paths: u -> v with dist[v] == dist[u] + 1
  pred <- vector("list", n)
  lapply(reach, function(u) {
    vs <- succ[[u]]
    vs <- vs[is.finite(dist[vs]) & dist[vs] == dist[u] + 1]
    lapply(vs, function(v) pred[[v]] <<- c(pred[[v]], u))
    NULL
  })

  paths <- vector("list", n)
  paths[[s]] <- list(s)
  # targets in increasing distance: every predecessor is already expanded
  ord <- reach[order(dist[reach])]
  n_paths <- 0L
  for (v in ord) {
    # sequential by construction: paths to v are built from paths to its
    # shortest-path predecessors
    if (v == s) next
    ps <- unlist(lapply(pred[[v]], function(u) {
      lapply(paths[[u]], function(p) c(p, v))
    }), recursive = FALSE)
    n_paths <- n_paths + length(ps)
    if (n_paths > max_paths) {
      stop(errorCondition(
        paste0("Shortest-path enumeration exceeded `max_paths` (",
               max_paths, "). The higher-order topology has too many ",
               "equally short paths; raise `max_paths` or reduce the ",
               "network (higher `min_freq`, lower `max_order`)."),
        class = "hypernets_too_many_paths", call = NULL))
    }
    paths[[v]] <- ps
  }
  paths
}

# ---------------------------------------------------------------------------
# PageRank on the higher-order topology
# ---------------------------------------------------------------------------

#' PageRank of a directed (optionally weighted) adjacency matrix
#'
#' @param mat Square weight matrix.
#' @param damping Damping factor.
#' @param weighted Logical. Use edge weights instead of the binary pattern.
#' @param max_iter,tol Power-iteration controls.
#' @return Named numeric vector summing to 1.
#' @noRd
.hoc_pagerank <- function(mat, damping, weighted, max_iter, tol) {
  n <- nrow(mat)
  a <- if (weighted) mat * 1.0 else (mat != 0) * 1.0
  out <- rowSums(a)
  dangling <- out == 0
  # dangling nodes redistribute uniformly (standard PageRank convention)
  a[dangling, ] <- 1
  out[dangling] <- n
  p <- a / out
  x <- rep(1 / n, n)
  converged <- FALSE
  for (iter in seq_len(max_iter)) {
    # power iteration: a sequential fixed-point loop
    y <- damping * as.numeric(crossprod(p, x)) + (1 - damping) / n
    if (sum(abs(y - x)) < tol) {
      x <- y
      converged <- TRUE
      break
    }
    x <- y
  }
  if (!converged) {
    warning(warningCondition(
      sprintf("PageRank did not converge in %d iterations (L1 change > %g).",
              as.integer(max_iter), tol),
      class = "hypernets_no_converge"))
  }
  stats::setNames(x, rownames(mat))
}

# ---------------------------------------------------------------------------
# First-order projection
# ---------------------------------------------------------------------------

#' Project higher-order node values onto first-order states
#'
#' @param values Numeric vector, one per higher-order node.
#' @param node_paths List of first-order state paths, one per node.
#' @param states Character vector of first-order states (output order).
#' @param projection "scaled", "last", "first", or "all" (pathpy's
#'   vocabulary).
#' @return Named numeric vector, one entry per state.
#' @noRd
.hoc_project <- function(values, node_paths, states, projection) {
  contrib <- switch(
    projection,
    scaled = lapply(seq_along(values), function(i) {
      p <- node_paths[[i]]
      list(state = p, value = rep(values[i] / length(p), length(p)))
    }),
    last = lapply(seq_along(values), function(i) {
      p <- node_paths[[i]]
      list(state = p[length(p)], value = values[i])
    }),
    first = lapply(seq_along(values), function(i) {
      list(state = node_paths[[i]][1L], value = values[i])
    }),
    all = lapply(seq_along(values), function(i) {
      p <- node_paths[[i]]
      list(state = p, value = rep(values[i], length(p)))
    })
  )
  st <- unlist(lapply(contrib, `[[`, "state"), use.names = FALSE)
  va <- unlist(lapply(contrib, `[[`, "value"), use.names = FALSE)
  keep <- st %in% states
  agg <- rowsum(va[keep], st[keep], reorder = FALSE)
  out <- stats::setNames(rep(0, length(states)), states)
  out[rownames(agg)] <- as.numeric(agg[, 1L])
  out
}

#' First-order distance matrix implied by a higher-order topology
#'
#' pathpy's projection (distance_matrix for HigherOrderNetwork) with the
#' per-node order generalization described at the top of this file.
#'
#' @param dist_ho Higher-order hop-distance matrix (nodes x nodes).
#' @param node_paths List of first-order paths, one per node.
#' @param states Character vector of first-order states.
#' @return Numeric matrix of first-order distances (states x states),
#'   zero on the diagonal, `Inf` where unreachable.
#' @noRd
.hoc_dist_first <- function(dist_ho, node_paths, states) {
  ns <- length(states)
  first_state <- vapply(node_paths, `[[`, character(1L), 1L)
  last_state <- vapply(node_paths, function(p) p[length(p)], character(1L))
  order_minus <- vapply(node_paths, length, integer(1L)) - 1L
  i_of <- match(first_state, states)
  j_of <- match(last_state, states)

  # Every finite higher-order hop count v -> w realizes a first-order
  # path from the first state of v to the last state of w, of length
  # dist_ho[v, w] + (order of v - 1). The diagonal of dist_ho is 0, so
  # this also covers what a node realizes on its own.
  idx <- which(is.finite(dist_ho), arr.ind = TRUE)
  i <- i_of[idx[, 1L]]
  j <- j_of[idx[, 2L]]
  d <- dist_ho[idx] + order_minus[idx[, 1L]]
  keep <- !is.na(i) & !is.na(j)

  out <- matrix(Inf, ns, ns, dimnames = list(states, states))
  if (any(keep)) {
    # aggregate per cell in one pass: assigning into duplicated matrix
    # indices would keep the last value written, not the smallest
    cell <- (j[keep] - 1L) * ns + i[keep]
    best <- tapply(d[keep], cell, min)
    out[as.integer(names(best))] <- as.numeric(best)
  }
  diag(out) <- 0
  out
}

# ---------------------------------------------------------------------------
# Betweenness on first-order states (pathpy's higher-order algorithm)
# ---------------------------------------------------------------------------

#' First-order betweenness from higher-order shortest paths
#'
#' Enumerates the shortest paths of the higher-order topology, maps each
#' to its first-order path, keeps the first-order-shortest ones per
#' (source, target) state pair (distinct paths only, as pathpy's set
#' semantics require), and credits intermediate states with the inverse
#' of the number of such paths.
#'
#' @param succ Successor list of the higher-order topology.
#' @param node_paths List of first-order paths, one per node.
#' @param states Character vector of first-order states.
#' @param max_paths Integer cap passed to the enumeration.
#' @return Named numeric vector, one entry per state.
#' @noRd
.hoc_betweenness_first <- function(succ, node_paths, states, max_paths) {
  n <- length(node_paths)
  best <- new.env(hash = TRUE, parent = emptyenv())

  register <- function(p1) {
    s1 <- p1[1L]
    d1 <- p1[length(p1)]
    key <- paste(s1, d1, sep = "\x02")
    l <- length(p1) - 1L
    cur <- best[[key]]
    path_key <- paste(p1, collapse = "\x01")
    if (is.null(cur) || l < cur$len) {
      best[[key]] <- list(len = l, paths = path_key)
    } else if (l == cur$len && !(path_key %in% cur$paths)) {
      cur$paths <- c(cur$paths, path_key)
      best[[key]] <- cur
    }
    NULL
  }

  lapply(seq_len(n), function(s) {
    # one BFS per source: sequential graph traversal
    paths <- .hoc_shortest_paths_from(succ, s, n, max_paths)
    lapply(paths, function(pl) {
      if (is.null(pl)) return(NULL)
      lapply(pl, function(p) {
        # higher_order_path_to_first_order: full first node, then the
        # last state of every subsequent node
        p1 <- c(node_paths[[p[1L]]],
                vapply(p[-1L], function(v) {
                  q <- node_paths[[v]]
                  q[length(q)]
                }, character(1L)))
        register(p1)
      })
      NULL
    })
    NULL
  })

  out <- stats::setNames(rep(0, length(states)), states)
  lapply(ls(best), function(key) {
    entry <- best[[key]]
    ends <- strsplit(key, "\x02", fixed = TRUE)[[1L]]
    share <- 1 / length(entry$paths)
    lapply(entry$paths, function(pk) {
      p1 <- strsplit(pk, "\x01", fixed = TRUE)[[1L]]
      if (length(p1) <= 2L) return(NULL)
      mid <- p1[-c(1L, length(p1))]
      mid <- mid[mid != ends[1L] & mid != ends[2L]]
      mid <- mid[mid %in% states]
      if (length(mid) == 0L) return(NULL)
      tab <- table(mid)
      out[names(tab)] <<- out[names(tab)] + as.integer(tab) * share
      NULL
    })
    NULL
  })
  out
}

# ---------------------------------------------------------------------------
# Public verb
# ---------------------------------------------------------------------------

#' Higher-order centralities with first-order projection
#'
#' Computes centralities on the higher-order topology of a
#' [build_hon()] network and, by default, projects them back onto the
#' first-order states, following Scholtes, Wider & Garas (2016). Because
#' a higher-order node carries the memory of how a state was reached,
#' these centralities can rank states differently from the same measures
#' on the first-order transition network - and that difference is the
#' evidence that memory matters for flow, not just for prediction.
#'
#' Semantics match pathpy 2.2.0 (the reference implementation by the
#' method's author), generalized from fixed-order to the variable-order
#' networks `build_hon()` produces: a higher-order node's order is the
#' number of states in the path it represents, and the first-order
#' distance implied by a higher-order hop count adds that node's order
#' minus one. With a uniform order the generalization reduces exactly to
#' pathpy's formulas.
#'
#' @param hon A `net_hon` object from [build_hon()].
#' @param type Character vector, any subset of
#'   `c("pagerank", "betweenness", "closeness")`. Default computes all
#'   three.
#' @param project Logical. `TRUE` (default) returns one row per
#'   first-order state; `FALSE` returns one row per higher-order node
#'   with the ordinary directed-graph centralities of the higher-order
#'   topology (useful to see which *contexts*, not states, dominate).
#' @param projection How a higher-order node's PageRank is distributed
#'   over the states of its path when `project = TRUE`: `"scaled"`
#'   (default, pathpy's default - each of the k states receives a k-th
#'   of the value, so the projected PageRank still sums to 1), `"last"`
#'   (all of it to the state actually occupied), `"first"`, or `"all"`
#'   (every state on the path receives the full value; the result no
#'   longer sums to 1). Ignored by `betweenness` and `closeness`, which
#'   are first-order by construction.
#' @param damping PageRank damping factor in (0, 1). Default `0.85`.
#' @param weighted Logical. `FALSE` (default, pathpy's convention) runs
#'   PageRank on the binary higher-order topology; `TRUE` uses the edge
#'   weights of `hon$matrix`.
#' @param max_iter,tol Power-iteration controls for PageRank.
#' @param max_paths Integer cap on the shortest paths enumerated for
#'   betweenness. Exceeding it raises a `hypernets_too_many_paths` error
#'   rather than silently truncating. Default `1e6`.
#' @param sort_by `NULL` (alphabetical by state or node, the default) or
#'   one of the requested `type`s - sort the table by that centrality,
#'   largest first, with ties broken by state/node so the order stays
#'   deterministic.
#'
#' @return A `data.frame`, one row per first-order state (or per
#'   higher-order node when `project = FALSE`), with the requested
#'   centralities as columns:
#'   \describe{
#'     \item{state}{First-order state (when `project = TRUE`).}
#'     \item{node, order}{Higher-order node and its order (when
#'       `project = FALSE`).}
#'     \item{pagerank}{Projected PageRank of the higher-order topology.}
#'     \item{betweenness}{Share of first-order shortest paths (implied by
#'       the higher-order topology) passing through the state.}
#'     \item{closeness}{Harmonic closeness on the implied first-order
#'       distances: the sum of inverse distances *into* the state.}
#'   }
#'   Rows are ordered by state (or by node) so the result is
#'   deterministic.
#'
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#'
#' @references
#' Scholtes, I., Wider, N., & Garas, A. (2016). Higher-order aggregate
#' networks in the analysis of temporal networks: path structures and
#' centralities. \emph{The European Physical Journal B} 89, 61.
#' \doi{10.1140/epjb/e2016-60663-0}
#'
#' Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
#' higher-order dependencies in networks. \emph{Science Advances} 2(5),
#' e1600028. \doi{10.1126/sciadv.1600028}
#'
#' @examples
#' seqs <- c(
#'   replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE),
#'   replicate(6, rep(c("x", "b", "d"), 4), simplify = FALSE)
#' )
#' hon <- build_hon(seqs, max_order = 2)
#' hon_centrality(hon)
#'
#' # Which contexts, rather than which states, carry the flow?
#' hon_centrality(hon, type = "pagerank", project = FALSE)
#'
#' @seealso [build_hon()], [bootstrap_hon()], [path_dependence()]
#'
#' @export
hon_centrality <- function(hon,
                           type = c("pagerank", "betweenness", "closeness"),
                           project = TRUE,
                           projection = c("scaled", "last", "first", "all"),
                           damping = 0.85,
                           weighted = FALSE,
                           max_iter = 1000L,
                           tol = 1e-12,
                           max_paths = 1e6,
                           sort_by = NULL,
                           top = NULL) {
  stopifnot(
    "`hon` must be a net_hon object from build_hon()" =
      inherits(hon, "net_hon"),
    "`project` must be TRUE or FALSE" =
      is.logical(project) && length(project) == 1L && !is.na(project),
    "`damping` must be a single number strictly between 0 and 1" =
      is.numeric(damping) && length(damping) == 1L && is.finite(damping) &&
      damping > 0 && damping < 1,
    "`weighted` must be TRUE or FALSE" =
      is.logical(weighted) && length(weighted) == 1L && !is.na(weighted),
    "`max_iter` must be a single positive number" =
      is.numeric(max_iter) && length(max_iter) == 1L && max_iter >= 1,
    "`tol` must be a single positive number" =
      is.numeric(tol) && length(tol) == 1L && tol > 0,
    "`max_paths` must be a single positive number" =
      is.numeric(max_paths) && length(max_paths) == 1L && max_paths >= 1
  )
  type <- match.arg(type, several.ok = TRUE)
  projection <- match.arg(projection)

  mat <- hon$matrix
  nodes <- rownames(mat)
  node_paths <- .hoc_node_paths(nodes)
  states <- hon$first_order_states
  n <- length(nodes)
  succ <- .hoc_succ(mat)

  out <- if (project) {
    data.frame(state = states, stringsAsFactors = FALSE)
  } else {
    data.frame(node = nodes,
               order = vapply(node_paths, length, integer(1L)),
               stringsAsFactors = FALSE)
  }

  if ("pagerank" %in% type) {
    pr <- .hoc_pagerank(mat, damping = damping, weighted = weighted,
                        max_iter = as.integer(max_iter), tol = tol)
    out$pagerank <- if (project) {
      as.numeric(.hoc_project(pr, node_paths, states, projection))
    } else {
      as.numeric(pr)
    }
  }

  if ("betweenness" %in% type) {
    out$betweenness <- if (project) {
      as.numeric(.hoc_betweenness_first(succ, node_paths, states, max_paths))
    } else {
      # ordinary directed-graph betweenness of the higher-order topology
      bw <- stats::setNames(rep(0, n), nodes)
      lapply(seq_len(n), function(s) {
        paths <- .hoc_shortest_paths_from(succ, s, n, max_paths)
        lapply(seq_len(n), function(t) {
          pl <- paths[[t]]
          if (is.null(pl) || t == s) return(NULL)
          share <- 1 / length(pl)
          lapply(pl, function(p) {
            if (length(p) <= 2L) return(NULL)
            mid <- p[-c(1L, length(p))]
            bw[mid] <<- bw[mid] + share
            NULL
          })
          NULL
        })
        NULL
      })
      as.numeric(bw)
    }
  }

  if ("closeness" %in% type) {
    dist_ho <- do.call(rbind, lapply(seq_len(n), function(s) {
      .hoc_bfs_dist(succ, s, n)
    }))
    out$closeness <- if (project) {
      d1 <- .hoc_dist_first(dist_ho, node_paths, states)
      # harmonic closeness over incoming distances (pathpy's formula)
      vapply(seq_along(states), function(j) {
        d <- d1[, j]
        d <- d[-j]
        sum(1 / d[is.finite(d) & d > 0])
      }, numeric(1L))
    } else {
      vapply(seq_len(n), function(j) {
        d <- dist_ho[, j]
        d <- d[-j]
        sum(1 / d[is.finite(d) & d > 0])
      }, numeric(1L))
    }
  }

  key <- if (project) out$state else out$node
  out <- if (is.null(sort_by)) {
    out[order(key), , drop = FALSE]
  } else {
    sort_by <- match.arg(sort_by, type)
    out[order(-out[[sort_by]], key), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
