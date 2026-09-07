# ---- HYPA: Hypothesis Testing for Path Anomalies ----
#
# Implements HYPA (LaRock et al. 2020) for detecting anomalous paths in
# sequential data. Uses a multi-hypergeometric null model on k-th order
# De Bruijn graphs to identify over/under-represented paths.

# ---------------------------------------------------------------------------
# Internal: Fit Xi matrix (iterative proportional fitting)
# ---------------------------------------------------------------------------

#' Compute propensity matrix Xi from node strengths
#'
#' Xi_{vw} = s_out(v) * s_in(w) for edges present in the De Bruijn graph.
#' The product of strengths gives N = sum(Xi) >> m = sum(adj), ensuring a
#' non-degenerate hypergeometric null model. This follows the HYPA null
#' model where edge propensity is proportional to the product of endpoint
#' weighted degrees.
#'
#' @param adj Square adjacency matrix of the De Bruijn graph.
#' @return Matrix Xi with same dimensions as adj (sum(Xi) >> sum(adj)).
#' @noRd
.hypa_fit_xi <- function(adj) {
  out_strength <- rowSums(adj)
  in_strength <- colSums(adj)
  mask <- adj > 0

  # Xi_{vw} = s_out(v) * s_in(w) where edges exist in the De Bruijn graph
  outer(out_strength, in_strength) * mask
}

# ---------------------------------------------------------------------------
# Internal: Compute HYPA scores
# ---------------------------------------------------------------------------

#' Compute hypergeometric p-values for each edge
#'
#' For each edge (v,w) with observed weight f, computes:
#'   \code{p_under = P(X <= f)} and \code{p_over = P(X >= f)}
#'   where \code{X ~ Hypergeometric(N, K, n)},
#'   \code{N = round(sum(Xi))}, \code{K = round(Xi[v,w])}, \code{n = sum(adj)}
#'
#' @param adj Adjacency matrix (edge weights = path frequencies).
#' @param xi Fitted propensity matrix.
#' @return Data frame with from, to, observed, expected, ratio, p_value,
#'   p_under, p_over, and anomaly.
#' @noRd
.hypa_compute_scores <- function(adj, xi) {
  n <- nrow(adj)
  nodes <- rownames(adj)

  # Total pool
  N_total <- round(sum(xi))
  n_draws <- sum(adj)

  # Collect edges
  edge_idx <- which(adj > 0, arr.ind = TRUE)
  if (nrow(edge_idx) == 0L) {
    return(data.frame(path = character(0L), from = character(0L),
                      to = character(0L), observed = integer(0L),
                      expected = numeric(0L), ratio = numeric(0L),
                      p_value = numeric(0L), p_under = numeric(0L),
                      p_over = numeric(0L), anomaly = character(0L),
                      stringsAsFactors = FALSE))
  }

  results <- lapply(seq_len(nrow(edge_idx)), function(idx) {
    i <- edge_idx[idx, 1L]
    j <- edge_idx[idx, 2L]
    f_obs <- adj[i, j]
    K <- round(xi[i, j])

    # Clamp parameters to valid range
    K <- min(K, N_total)
    K <- max(K, 0L)
    n_clamp <- min(n_draws, N_total)

    # Inclusive hypergeometric tails.
    p_under <- stats::phyper(f_obs, K, N_total - K, n_clamp)
    p_over <- stats::phyper(f_obs - 1, K, N_total - K, n_clamp,
                            lower.tail = FALSE)

    # Expected value: n * K / N
    expected <- if (N_total > 0) n_draws * K / N_total else 0

    # Reconstruct full path; from = context, to = next state
    from_parts <- strsplit(nodes[i], .HON_SEP, fixed = TRUE)[[1L]]
    to_parts <- strsplit(nodes[j], .HON_SEP, fixed = TRUE)[[1L]]
    next_state <- to_parts[length(to_parts)]
    path <- paste(c(from_parts, next_state), collapse = " -> ")

    ratio <- if (expected > 0) f_obs / expected else Inf

    data.frame(
      path = path,
      from = paste(from_parts, collapse = " -> "),
      to = next_state,
      observed = as.integer(f_obs),
      expected = expected,
      ratio = ratio,
      p_value = p_under,
      p_under = p_under,
      p_over = p_over,
      stringsAsFactors = FALSE
    )
  })

  result <- do.call(rbind, results)
  # Classify anomalies
  result$anomaly <- ifelse(result$p_under < 0.05, "under",
                           ifelse(result$p_over < 0.05, "over", "normal"))

  result
}

# ---------------------------------------------------------------------------
# Main function: build_hypa
# ---------------------------------------------------------------------------

#' HYPA scores for a single De Bruijn order
#'
#' Internal helper. Builds the order-\code{ord} De Bruijn graph from already
#' parsed \code{trajectories}, fits the propensity matrix, computes
#' hypergeometric tail probabilities, applies multiple-testing correction and
#' classifies each edge as over-/under-represented or normal.
#'
#' @param trajectories List of character vectors (parsed trajectories), as
#'   produced by \code{.hon_parse_input()}.
#' @param ord Integer scalar. Order of the De Bruijn graph for this layer.
#' @param alpha Numeric. Significance threshold for anomaly classification.
#' @param min_count Integer. Minimum observed count for an edge to be eligible
#'   for an anomaly classification.
#' @param p_adjust Character. Multiple-testing correction method (a value from
#'   \code{stats::p.adjust.methods} or \code{"none"}).
#' @return A plain \code{list} (no S3 class) with elements \code{scores},
#'   \code{over}, \code{under}, \code{adjacency}, \code{weights}, \code{xi},
#'   \code{edges}, \code{nodes}, \code{meta}, \code{n_over}, \code{n_under},
#'   \code{n_edges}, \code{order}; or \code{NULL} when the order produces no
#'   edges.
#' @noRd
.hypa_one_order <- function(trajectories, ord, alpha, min_count, p_adjust) {
  kg <- .mogen_count_kgrams(trajectories, ord)
  if (nrow(kg$edges) == 0L) {
    return(NULL)
  }
  nodes <- kg$nodes
  n <- length(nodes)
  adj <- matrix(0, nrow = n, ncol = n, dimnames = list(nodes, nodes))
  idx <- cbind(match(kg$edges$from, nodes), match(kg$edges$to, nodes))
  adj[idx] <- kg$edges$weight

  xi <- .hypa_fit_xi(adj)
  scores <- .hypa_compute_scores(adj, xi)

  if (p_adjust != "none" && nrow(scores) > 0L) {
    scores$p_adjusted_under <- stats::p.adjust(scores$p_under, method = p_adjust)
    scores$p_adjusted_over  <- stats::p.adjust(scores$p_over,  method = p_adjust)
  } else {
    scores$p_adjusted_under <- scores$p_under
    scores$p_adjusted_over  <- scores$p_over
  }

  scores$anomaly <- ifelse(
    scores$observed < min_count, "normal",
    ifelse(scores$p_adjusted_under < alpha, "under",
           ifelse(scores$p_adjusted_over < alpha, "over", "normal"))
  )
  scores$order <- ord

  anom_over  <- scores[scores$anomaly == "over",  , drop = FALSE]
  anom_under <- scores[scores$anomaly == "under", , drop = FALSE]
  normal     <- scores[scores$anomaly == "normal", , drop = FALSE]
  anom_over  <- anom_over[order(-anom_over$ratio),  , drop = FALSE]
  anom_under <- anom_under[order(anom_under$ratio), , drop = FALSE]
  scores <- rbind(anom_over, anom_under, normal)
  rownames(scores) <- NULL

  cg <- .ho_cograph_fields(adj, nodes, method = "hypa")
  list(
    scores = scores, over = anom_over, under = anom_under,
    adjacency = adj, weights = cg$weights, xi = xi,
    edges = cg$edges, nodes = cg$nodes, meta = cg$meta,
    n_over = nrow(anom_over), n_under = nrow(anom_under),
    n_edges = nrow(scores), order = ord
  )
}

#' Detect Path Anomalies via HYPA
#'
#' Constructs a k-th order De Bruijn graph from sequential trajectory data and
#' uses a hypergeometric null model to detect paths with anomalous frequencies.
#' Paths occurring more or less often than expected under the null model are
#' flagged as over- or under-represented.
#'
#' @param data A data.frame (rows = trajectories), list of character vectors,
#'   \code{tna} object, or \code{netobject} with sequence data. For
#'   \code{tna}/\code{netobject}, numeric state IDs are automatically
#'   converted to label names.
#' @param order Integer scalar or integer vector. Order(s) of the De Bruijn
#'   graph (default \code{2L}). An order of \code{k} detects anomalies in
#'   paths of length \code{k}. When a vector is supplied, one De Bruijn layer
#'   is built per order and the per-order results are stored in
#'   \code{$by_order} (named by order). The orders are sorted ascending
#'   internally, so the \code{cograph_network} slots
#'   (\code{$weights}, \code{$edges}, \code{$adjacency}, \code{$xi},
#'   \code{$nodes}, \code{$meta}) always describe the network of the
#'   \emph{lowest order that produced a layer} (a requested order with no
#'   edges is dropped from \code{$by_order}, \code{$order} and \code{$k}),
#'   regardless of the order in which the vector is given; \code{$scores}
#'   aggregates every built order.
#' @param alpha Numeric. Significance threshold for anomaly classification
#'   (default 0.05). Paths with HYPA score < alpha are under-represented;
#'   paths with score > 1-alpha are over-represented.
#' @param min_count Integer. Minimum observed count for a path to be
#'   classified as anomalous (default 5). Paths with fewer observations
#'   are always classified as \code{"normal"} regardless of their
#'   HYPA score, since rare occurrences are unreliable.
#' @param p_adjust Character. Method for multiple testing correction of
#'   p-values. Default \code{"BH"} (Benjamini-Hochberg FDR control).
#'   Accepts any method from \code{\link[stats]{p.adjust.methods}} or
#'   \code{"none"} to skip correction. Under- and over-representation
#'   p-values are adjusted separately (two-sided testing).
#' @param k Deprecated. Former name of \code{order}; if supplied it overrides
#'   \code{order} and emits a deprecation message. Use \code{order} instead.
#' @return An object of class \code{c("net_hypa", "cograph_network")} with
#'   components:
#'   \describe{
#'     \item{scores}{Data frame with path, from, to, observed, expected,
#'       ratio, p_value, p_under, p_over, p_adjusted_under,
#'       p_adjusted_over, anomaly, order
#'       columns (one block of rows per requested order). The \code{path}
#'       column shows the full state sequence
#'       (e.g., "A -> B -> C"); \code{from} is the context (conditioning
#'       states); \code{to} is the next state; \code{ratio} is
#'       observed / expected; \code{p_value} is retained as an alias for
#'       \code{p_under}, the raw lower-tail hypergeometric CDF value;
#'       \code{p_over} is the inclusive upper-tail probability
#'       \code{P(X >= observed)}; \code{p_adjusted_under} and
#'       \code{p_adjusted_over} are the corrected p-values for under- and
#'       over-representation tests respectively.}
#'     \item{ho_edges}{Alias for \code{scores} (all orders, arrow notation).}
#'     \item{over}{Subset of \code{scores} classified as over-represented.}
#'     \item{under}{Subset of \code{scores} classified as under-represented.}
#'     \item{adjacency}{Weighted adjacency matrix of the lowest-order De
#'       Bruijn graph.}
#'     \item{weights}{cograph weight matrix of the lowest-order graph.}
#'     \item{xi}{Fitted propensity matrix of the lowest-order graph.}
#'     \item{edges}{cograph edge data.frame of the lowest-order graph.}
#'     \item{by_order}{Named list of per-order result lists.}
#'     \item{order}{Integer vector of orders actually built (sorted ascending).}
#'     \item{k}{Back-compatibility alias for \code{order}.}
#'     \item{alpha}{Significance threshold used.}
#'     \item{p_adjust}{Multiple testing correction method used.}
#'     \item{n_anomalous}{Number of anomalous paths detected (all orders).}
#'     \item{n_over}{Number of over-represented paths (all orders).}
#'     \item{n_under}{Number of under-represented paths (all orders).}
#'     \item{n_edges}{Total number of edges (all orders).}
#'     \item{nodes}{data.frame (\code{id}, \code{label}, \code{name}) of the
#'       lowest-order De Bruijn graph nodes (arrow notation).}
#'     \item{directed}{Logical. Always \code{TRUE}.}
#'     \item{meta}{cograph meta list of the lowest-order graph.}
#'     \item{node_groups}{Always \code{NULL}.}
#'   }
#'
#' @references
#' LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
#' & Schweitzer, F. (2020). HYPA: Efficient Detection of Path Anomalies in
#' Time Series Data on Networks. \emph{SDM 2020}, 460-468.
#'
#' @examples
#' seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
#' hyp <- build_hypa(seqs, order = 2)
#'
#' \donttest{
#' trajs <- list(c("A","B","C"), c("A","B","C"), c("A","B","C"),
#'               c("A","B","D"), c("C","B","D"), c("C","B","A"))
#' h <- build_hypa(trajs, order = 2)
#' print(h)
#' }
#'
#' @export
build_hypa <- function(data, order = 2L, alpha = 0.05, min_count = 5L,
                       p_adjust = "BH", k = NULL) {
  if (!is.null(k)) {
    .Deprecated(msg = "'k' is deprecated; use 'order' instead.")
    order <- k
  }
  data <- .coerce_sequence_input(data)

  valid_methods <- c(stats::p.adjust.methods, "none")
  if (!is.character(p_adjust) || length(p_adjust) != 1L ||
      !p_adjust %in% valid_methods) {
    stop(sprintf("'p_adjust' must be one of: %s",
                 paste(valid_methods, collapse = ", ")), call. = FALSE)
  }
  stopifnot(
    "'data' must be a data.frame or list" =
      is.data.frame(data) || is.list(data),
    "'order' must contain only whole numbers" =
      is.numeric(order) && all(order == round(order)),
    "'order' must contain only integers >= 1" = all(order >= 1L),
    "'alpha' must be in (0, 0.5)" = alpha > 0 && alpha < 0.5,
    "'min_count' must be >= 1" = min_count >= 1L
  )
  # Sort ascending so the primary cograph slot (built from per_order[[1L]])
  # is deterministically the lowest order that produces a layer, independent
  # of how the `order` vector is arranged (see @param order).
  order <- sort(as.integer(order))
  min_count <- as.integer(min_count)

  trajectories <- .hon_parse_input(data, collapse_repeats = FALSE)
  if (length(trajectories) == 0L) {
    stop("No valid trajectories (each must have at least 2 states)")
  }

  per_order <- lapply(order, function(o) {
    .hypa_one_order(trajectories, o, alpha, min_count, p_adjust)
  })
  names(per_order) <- as.character(order)
  per_order <- Filter(Negate(is.null), per_order)
  if (length(per_order) == 0L) {
    stop("No edges at any requested order (paths too short or too few)")
  }

  scores_all <- do.call(rbind, lapply(per_order, `[[`, "scores"))
  over_all   <- do.call(rbind, lapply(per_order, `[[`, "over"))
  under_all  <- do.call(rbind, lapply(per_order, `[[`, "under"))
  rownames(scores_all) <- NULL
  if (!is.null(over_all))  rownames(over_all)  <- NULL
  if (!is.null(under_all)) rownames(under_all) <- NULL

  # Orders that actually produced a layer. Filter() above dropped any
  # requested order with no edges (e.g. an order higher than the sequence
  # lengths allow), so $order/$k must report only the built layers to stay
  # consistent with $by_order and $scores$order (Codex review).
  built_order <- as.integer(names(per_order))
  # Primary cograph slot = network of the lowest built order.
  # `order` was sorted ascending above, so per_order[[1L]] is the lowest.
  primary <- per_order[[1L]]
  result <- list(
    scores = scores_all,
    ho_edges = scores_all,
    edges = primary$edges,
    over = over_all,
    under = under_all,
    adjacency = primary$adjacency,
    weights = primary$weights,
    xi = primary$xi,
    by_order = per_order,
    order = built_order,
    k = built_order,  # back-compat alias; prefer $order
    alpha = alpha,
    p_adjust = p_adjust,
    n_anomalous = sum(vapply(per_order, function(x) x$n_over + x$n_under, integer(1))),
    n_over  = sum(vapply(per_order, `[[`, integer(1), "n_over")),
    n_under = sum(vapply(per_order, `[[`, integer(1), "n_under")),
    n_edges = sum(vapply(per_order, `[[`, integer(1), "n_edges")),
    nodes = primary$nodes,
    directed = TRUE,
    meta = primary$meta,
    node_groups = NULL
  )
  class(result) <- c("net_hypa", "cograph_network")
  result
}

# ---------------------------------------------------------------------------
# S3 methods
# ---------------------------------------------------------------------------

#' Print Method for net_hypa
#'
#' @param x A \code{net_hypa} object.
#' @param ... Additional arguments (ignored).
#'
#' @return The input object, invisibly.
#'
#' @examples
#' seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
#' hyp <- build_hypa(seqs, k = 2)
#' print(hyp)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B","C","A","B","C","A"),
#'   V2 = c("B","C","A","B","C","A","B","C","A","B"),
#'   V3 = c("C","A","B","C","A","B","C","A","B","C"),
#'   V4 = c("A","B","C","A","B","C","A","B","C","A")
#' )
#' hypa <- build_hypa(seqs, k = 2L)
#' print(hypa)
#' }
#'
#' @export
print.net_hypa <- function(x, ...) {
  ord_str <- paste(x$order, collapse = ", ")
  cat("HYPA: Path Anomaly Detection\n")
  cat(sprintf("  Order(s):     %s\n", ord_str))
  cat(sprintf("  Edges:        %d\n", x$n_edges))
  cat(sprintf("  Anomalous:    %d (alpha=%.2f, p_adjust=%s)\n",
              x$n_anomalous, x$alpha, x$p_adjust %||% "none"))
  cat(sprintf("    Over-repr:  %d\n", x$n_over))
  cat(sprintf("    Under-repr: %d\n", x$n_under))
  if (length(x$order) > 1L) {
    per <- vapply(x$by_order, function(po)
      sprintf("order %d: %d edges (%d over, %d under)",
              po$order, po$n_edges, po$n_over, po$n_under),
      character(1))
    cat("  Per-order:\n")
    for (line in per) cat(sprintf("    %s\n", line))
  }
  invisible(x)
}

#' Summary Method for net_hypa
#'
#' @param object A \code{net_hypa} object.
#' @param n Integer. Maximum number of paths to display per category
#'   (default: 10).
#' @param type Character. Which anomalies to show: \code{"all"} (default),
#'   \code{"over"}, or \code{"under"}.
#' @param order_by Character. Ranking used within each anomaly direction:
#'   \code{"sig"} ranks by the active tail probability, \code{"freq"} by
#'   observed count, \code{"ratio"} by observed/expected ratio, and
#'   \code{"path"} alphabetically.
#' @param ... Additional arguments (ignored).
#'
#' @return A data frame with path, observed, expected, ratio, p_tail, and
#'   direction columns.
#'
#'   Returned **invisibly**: `summary(x)` prints the summary and nothing
#'   else; assign the result to keep the table.
#' @examples
#' seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
#' hyp <- build_hypa(seqs, k = 2)
#' summary(hyp)
#'
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B","C","A","B","C","A"),
#'   V2 = c("B","C","A","B","C","A","B","C","A","B"),
#'   V3 = c("C","A","B","C","A","B","C","A","B","C"),
#'   V4 = c("A","B","C","A","B","C","A","B","C","A")
#' )
#' hypa <- build_hypa(seqs, k = 2L)
#' summary(hypa)
#' summary(hypa, type = "over", n = 5)
#' }
#'
#' @export
summary.net_hypa <- function(object, n = 10L,
                             type = c("all", "over", "under"),
                             order_by = c("sig", "freq", "frequency",
                                          "ratio", "path"),
                             ...) {
  type <- match.arg(type)
  order_by <- match.arg(order_by)
  if (order_by == "frequency") order_by <- "freq"

  reorder_df <- function(df, dir) {
    if (is.null(df) || nrow(df) == 0L) return(df)
    # p_value is the hypergeometric CDF P(X <= obs): close to 1 for
    # over-represented paths, close to 0 for under-represented. So "by
    # sig" picks the tail extremeness in each direction.
    idx <- switch(order_by,
      ratio = if (dir == "over") order(-df$ratio) else order(df$ratio),
      sig   = if (dir == "over") order(df$p_over) else order(df$p_under),
      freq  = order(-df$observed),
      path  = order(df$path)
    )
    df[idx, , drop = FALSE]
  }

  over_sorted  <- reorder_df(object$over,  "over")
  under_sorted <- reorder_df(object$under, "under")

  cat("HYPA Summary\n\n")
  cat(sprintf("  Order(s): %s | Edges: %d\n",
              paste(object$order, collapse = ", "), object$n_edges))
  cat(sprintf("  Alpha: %.2f | p_adjust: %s\n",
              object$alpha, object$p_adjust %||% "none"))
  cat(sprintf("  Anomalous: %d (over: %d, under: %d) | order_by: %s\n\n",
              object$n_anomalous, object$n_over, object$n_under, order_by))

  display_df <- function(df, dir) {
    if (is.null(df) || nrow(df) == 0L) return(df)
    out <- df
    out$p_tail <- if (dir == "over") out$p_over else out$p_under
    out
  }

  show_cols <- c("order", "path", "observed", "expected", "ratio", "p_tail")
  if (length(object$order) == 1L) show_cols <- setdiff(show_cols, "order")

  if (type %in% c("all", "over") && object$n_over > 0L) {
    cat("  Over-represented (top", min(n, object$n_over), "):\n")
    over_display <- display_df(over_sorted, "over")
    print(utils::head(over_display[, show_cols, drop = FALSE], n),
          row.names = FALSE)
    cat("\n")
  }

  if (type %in% c("all", "under") && object$n_under > 0L) {
    cat("  Under-represented (top", min(n, object$n_under), "):\n")
    under_display <- display_df(under_sorted, "under")
    print(utils::head(under_display[, show_cols, drop = FALSE], n),
          row.names = FALSE)
    cat("\n")
  }

  if (object$n_anomalous == 0L) {
    cat("  No anomalous paths detected.\n")
  }

  pick <- function(df, dir) {
    if (is.null(df) || nrow(df) == 0L) return(NULL)
    out <- display_df(df, dir)
    cols <- c("order", "path", "observed", "expected", "ratio", "p_tail")
    out <- out[, intersect(cols, names(out)), drop = FALSE]
    out$direction <- dir
    utils::head(out, n)
  }
  parts <- switch(type,
    all   = list(pick(over_sorted, "over"), pick(under_sorted, "under")),
    over  = list(pick(over_sorted, "over")),
    under = list(pick(under_sorted, "under"))
  )
  combined <- do.call(rbind, Filter(Negate(is.null), parts))
  if (is.null(combined)) {
    combined <- data.frame(path = character(0L),
                           observed = numeric(0L),
                           expected = numeric(0L),
                           ratio = numeric(0L),
                           p_tail = numeric(0L),
                           direction = character(0L),
                           stringsAsFactors = FALSE)
  } else {
    if (type == "all" && nrow(combined) > 0L) {
      idx <- switch(order_by,
        sig   = order(combined$p_tail),
        ratio = order(-abs(log(combined$ratio))),
        freq  = order(-combined$observed),
        path  = order(combined$path)
      )
      combined <- combined[idx, , drop = FALSE]
    }
    row.names(combined) <- NULL
  }
  invisible(combined)
}

#' Coerce a net_hypa to a tidy table
#'
#' @param x A `net_hypa` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"scores"` (default) for every scored path, `"over"` for the
#'   significantly over-represented paths only, or `"under"` for the
#'   significantly under-represented ones.
#' @param sort_by `NULL` (construction order, default), `"p"` (sorts by
#'   `p_value`, smallest first), `"ratio"` (largest observed/expected
#'   deviation first) or `"observed"` (heaviest paths first). `p_value` is
#'   the two-sided value; the direction-specific `p_adjusted_under` /
#'   `p_adjusted_over` columns are left for the caller to read, since a
#'   single ordering cannot serve both directions.
#' @return A data.frame, one row per path: `path`, `from`, `to`,
#'   `observed`, `expected`, `ratio`, the p-values (`p_value`, `p_under`,
#'   `p_over`, `p_adjusted_under`, `p_adjusted_over`), the `anomaly` label
#'   and the path `order`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_hypa <- function(x, row.names = NULL, optional = FALSE, ...,
                                   what = c("scores", "over", "under"),
                                   sort_by = NULL, top = NULL) {
  what <- match.arg(what)
  out <- switch(what, scores = x$scores, over = x$over, under = x$under)
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, c("p", "ratio", "observed"))
    col <- if (sort_by == "p") "p_value" else sort_by
    if (is.null(out[[col]])) {
      stop(errorCondition(
        sprintf("no `%s` column found on this net_hypa result", col),
        class = "hypernets_bad_input", call = NULL))
    }
    # p ascends (most significant first); the magnitude keys descend.
    key <- if (sort_by == "p") out[[col]] else -out[[col]]
    out <- out[order(key, out$from, out$to), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
