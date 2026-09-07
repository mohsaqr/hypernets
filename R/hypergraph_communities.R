# Repeated Infomap on the association projection, following the paper's
# representation-comparison workflow. hypernets owns the hypergraph projection;
# cograph owns the graph clustering and quality kernels.

.thg_partition_vector <- function(partition, nodes) {
  if (inherits(partition, "hg_communities")) partition <- partition$medoid
  if (is.data.frame(partition)) {
    label_col <- intersect(c("community", "cluster", "label", "predicted"),
                           names(partition))
    if (!"node" %in% names(partition) || !length(label_col)) {
      .thg_bad_input("a partition data.frame needs `node` and a label column")
    }
    partition <- stats::setNames(as.character(partition[[label_col[1L]]]),
                                  as.character(partition$node))
  }
  if (is.null(names(partition))) {
    if (length(partition) != length(nodes)) {
      .thg_bad_input("an unnamed partition must have one label per node")
    }
    names(partition) <- nodes
  }
  if (length(setdiff(nodes, names(partition)))) {
    .thg_bad_input("the partition does not cover every projected node")
  }
  as.character(partition[nodes])
}

#' Stable Infomap communities of a hypergraph projection
#'
#' Builds the normalized association graph of Coupette et al. (2024), runs
#' Infomap repeatedly, compares every pair of partitions with AMI, ARI and
#' NMI, and returns the run with the largest summed AMI as the medoid. The
#' paper uses 50 seeds and 100 Infomap trials per seed, which are the defaults.
#'
#' @param hg A static `net_hypergraph`.
#' @param n_runs Number of independent seeded Infomap runs (default 50).
#' @param trials Infomap trials within each run (default 100).
#' @param seeds Integer seeds. `NULL` uses `seq_len(n_runs)`.
#' @param method Which graph Infomap runs on: the `"association"` projection
#'   (default; the paper's hypergraph-derived representations `bh`, `bhs`,
#'   `mh`, `mhs`) or the `"citation"` projection (its classic graph
#'   representations `bg`, `mg`, and with `directed = FALSE` their undirected
#'   variants `bgu`, `mgu`). See [hg_project()].
#' @param duplicate_edges,self_association,edge_source Projection controls
#'   passed to [hg_project()]. Together these reproduce the paper's
#'   binary/multi and self-association representations.
#' @param directed For `method = "citation"`: run Infomap with directed flow
#'   on the source-to-member graph? Default `FALSE`.
#' @return An `hg_communities` object containing `medoid` (a tidy node/community
#'   table), all `partitions`, AMI/ARI/NMI similarity matrices, run metadata,
#'   community sizes, and the graph `projection` Infomap ran on.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#' hypergraphs. *Philosophical Transactions of the Royal Society A*,
#' 382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
#'   edge = rep(paste0("e", 1:4), each = 3)
#' )
#' h <- group_hypergraph(dat, "member", "edge")
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   fit <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
#'   as.data.frame(fit)
#' }
#' @export
hg_communities <- function(hg, n_runs = 50L, trials = 100L, seeds = NULL,
                           method = c("association", "citation"),
                           duplicate_edges = c("count", "collapse"),
                           self_association = FALSE, edge_source = NULL,
                           directed = FALSE) {
  .thg_check_hg(hg)
  method <- match.arg(method)
  duplicate_edges <- match.arg(duplicate_edges)
  if (!is.logical(directed) || length(directed) != 1L || is.na(directed)) {
    .thg_bad_input("`directed` must be TRUE or FALSE")
  }
  if (directed && !identical(method, "citation")) {
    .thg_bad_input("`directed` applies to `method = \"citation\"` only")
  }
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(errorCondition(
      "`hg_communities()` needs the suggested package `igraph` for Infomap",
      class = "hypernets_missing_dependency", call = NULL
    ))
  }
  whole <- function(x, name) {
    if (length(x) != 1L || !is.numeric(x) || !is.finite(x) || x < 1 ||
        abs(x - round(x)) > sqrt(.Machine$double.eps)) {
      .thg_bad_input(sprintf("`%s` must be one positive whole number", name))
    }
    as.integer(x)
  }
  n_runs <- whole(n_runs, "n_runs")
  trials <- whole(trials, "trials")
  seeds <- seeds %||% seq_len(n_runs)
  if (length(seeds) != n_runs || any(!is.finite(seeds)) ||
      any(abs(seeds - round(seeds)) > sqrt(.Machine$double.eps))) {
    .thg_bad_input("`seeds` must contain one whole-number seed per run")
  }
  seeds <- as.integer(seeds)

  projection <- if (identical(method, "association")) {
    hg_project(
      hg, method = "association", what = "matrix",
      duplicate_edges = duplicate_edges,
      self_association = self_association, edge_source = edge_source
    )
  } else {
    hg_project(
      hg, method = "citation", what = "matrix",
      duplicate_edges = duplicate_edges, edge_source = edge_source,
      directed = directed
    )
  }
  projection <- as.matrix(projection)
  nodes <- rownames(projection)
  partitions <- vector("list", n_runs)
  metadata <- vector("list", n_runs)

  if (length(nodes) == 0L) {
    .thg_bad_input("community detection requires at least one projected node")
  }
  edgeless <- !any(projection != 0)
  for (i in seq_len(n_runs)) {
    if (edgeless) {
      tab <- data.frame(node = nodes, community = seq_along(nodes),
                        stringsAsFactors = FALSE)
      codelength <- NA_real_
    } else {
      fit <- cograph::community_infomap(
        projection, nb.trials = trials, seed = seeds[i], directed = directed
      )
      tab <- data.frame(node = as.character(fit$node),
                        community = as.character(fit$community),
                        stringsAsFactors = FALSE)
      tab <- tab[match(nodes, tab$node), , drop = FALSE]
      codelength <- attr(fit, "igraph_result")$codelength %||% NA_real_
    }
    tab$run <- i
    tab$seed <- seeds[i]
    partitions[[i]] <- tab
    metadata[[i]] <- data.frame(
      run = i, seed = seeds[i], n_communities = length(unique(tab$community)),
      codelength = as.numeric(codelength), stringsAsFactors = FALSE
    )
  }

  similarities <- lapply(c("ami", "ari", "nmi"), function(metric) {
    m <- diag(1, n_runs, n_runs)
    dimnames(m) <- list(paste0("run_", seq_len(n_runs)),
                        paste0("run_", seq_len(n_runs)))
    if (n_runs > 1L) {
      for (i in seq_len(n_runs - 1L)) {
        for (j in (i + 1L):n_runs) {
          a <- partitions[[i]]$community
          b <- partitions[[j]]$community
          m[i, j] <- m[j, i] <- switch(metric,
            ami = .thg_ami(a, b), ari = .thg_ari(a, b), nmi = .thg_nmi(a, b)
          )
        }
      }
    }
    m
  })
  names(similarities) <- c("ami", "ari", "nmi")
  medoid_run <- unname(which.max(rowSums(similarities$ami))[1L])
  medoid <- partitions[[medoid_run]][, c("node", "community"), drop = FALSE]
  sizes <- as.data.frame(table(medoid$community), stringsAsFactors = FALSE)
  names(sizes) <- c("community", "n_nodes")
  sizes <- sizes[order(-sizes$n_nodes, sizes$community), , drop = FALSE]
  rownames(sizes) <- NULL

  structure(list(
    medoid = medoid,
    partitions = do.call(rbind, partitions),
    similarity = similarities,
    runs = do.call(rbind, metadata),
    sizes = sizes,
    projection = projection,
    medoid_run = medoid_run,
    params = list(n_runs = n_runs, trials = trials, seeds = seeds,
                  method = method, duplicate_edges = duplicate_edges,
                  self_association = self_association, directed = directed)
  ), class = "hg_communities")
}

#' Compare community structure across representations
#'
#' Sets several [hg_communities()] fits of the same data side by side, as the
#' paper does for its eight GFCC representations (Coupette et al. 2024,
#' Figure 8): the cluster-size distribution of each AMI medoid, the pairwise
#' AMI, ARI and NMI between medoids, and per-fit summaries (number of
#' communities, singletons, the balance between the two largest clusters).
#'
#' @param ... Named `hg_communities` fits, or one named list of them. The
#'   names label the representations (`bh`, `mhs`, `bgu`, ...).
#' @param hg Optional: the static `net_hypergraph` the fits were computed on.
#'   When given, every medoid is scored with [hg_community_quality()] on the
#'   projection its own fit used, and the scores are available as
#'   `what = "quality"`.
#' @param edge_source Hyperedge sources for the citation and self-association
#'   projections when scoring, as in [hg_project()].
#' @return A `hypernets_community_comparison` object. `as.data.frame()` returns
#'   its `"summary"` (default; one row per fit with `model`, `medoid_seed`,
#'   `n_communities`, `n_singletons`, `n_nontrivial`, `largest`, `second`
#'   and `balance` = second / largest), `"similarity"` (one row per pair of
#'   fits with `model_a`, `model_b`, `ami`, `ari`, `nmi`, on the nodes the
#'   two medoids share), `"sizes"` (one row per community of every medoid
#'   with `model`, `rank`, `n_nodes`) or, when `hg` was given, `"quality"`
#'   (one row per fit with the columns of [hg_community_quality()]).
#'   `plot()` draws the cluster-size distributions (`what = "sizes"`, the
#'   number of communities at least as large as each size, on logarithmic
#'   axes) or the similarity matrix (`what = "similarity"`, AMI below and
#'   ARI above the diagonal, values printed in the cells).
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
#'   edge = rep(paste0("e", 1:4), each = 3)
#' )
#' h <- group_hypergraph(dat, "member", "edge")
#' if (requireNamespace("igraph", quietly = TRUE)) {
#'   multi <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
#'   binary <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2,
#'                            duplicate_edges = "collapse")
#'   comparison <- hg_compare_communities(mh = multi, bh = binary)
#'   as.data.frame(comparison)
#'   as.data.frame(comparison, what = "similarity")
#' }
#' @export
hg_compare_communities <- function(..., hg = NULL, edge_source = NULL) {
  fits <- list(...)
  if (length(fits) == 1L && is.list(fits[[1L]]) &&
      !inherits(fits[[1L]], "hg_communities")) {
    fits <- fits[[1L]]
  }
  if (!is.null(hg)) .thg_check_hg(hg)
  if (length(fits) < 2L || is.null(names(fits)) || any(!nzchar(names(fits))) ||
      anyDuplicated(names(fits)) ||
      !all(vapply(fits, inherits, logical(1L), "hg_communities"))) {
    .thg_bad_input("supply at least two distinctly named hg_communities() fits")
  }
  models <- names(fits)
  summary_rows <- lapply(models, function(model) {
    fit <- fits[[model]]
    sizes <- fit$sizes$n_nodes
    data.frame(
      model = model, medoid_seed = fit$runs$seed[fit$medoid_run],
      n_runs = nrow(fit$runs), n_communities = length(sizes),
      n_singletons = sum(sizes == 1L), n_nontrivial = sum(sizes > 1L),
      largest = sizes[[1L]],
      second = if (length(sizes) > 1L) sizes[[2L]] else NA_integer_,
      balance = if (length(sizes) > 1L) sizes[[2L]] / sizes[[1L]] else NA_real_,
      stringsAsFactors = FALSE
    )
  })
  size_rows <- lapply(models, function(model) {
    sizes <- fits[[model]]$sizes$n_nodes
    data.frame(model = model, rank = seq_along(sizes), n_nodes = sizes,
               stringsAsFactors = FALSE)
  })
  pairs <- utils::combn(models, 2L)
  similarity_rows <- lapply(seq_len(ncol(pairs)), function(k) {
    a <- fits[[pairs[1L, k]]]$medoid
    b <- fits[[pairs[2L, k]]]$medoid
    joined <- merge(a, b, by = "node", suffixes = c("_a", "_b"))
    data.frame(
      model_a = pairs[1L, k], model_b = pairs[2L, k], n_nodes = nrow(joined),
      ami = .thg_ami(joined$community_a, joined$community_b),
      ari = .thg_ari(joined$community_a, joined$community_b),
      nmi = .thg_nmi(joined$community_a, joined$community_b),
      stringsAsFactors = FALSE
    )
  })
  quality <- if (is.null(hg)) NULL else do.call(rbind, lapply(models, function(model) {
    p <- fits[[model]]$params
    q <- hg_community_quality(
      hg, fits[[model]], method = p$method %||% "association",
      duplicate_edges = p$duplicate_edges,
      self_association = isTRUE(p$self_association), edge_source = edge_source
    )
    data.frame(model = model, q, stringsAsFactors = FALSE)
  }))
  structure(list(
    models = models,
    summary = do.call(rbind, summary_rows),
    sizes = do.call(rbind, size_rows),
    similarity = do.call(rbind, similarity_rows),
    quality = quality
  ), class = "hypernets_community_comparison")
}

#' @rdname hg_compare_communities
#' @param x A `hypernets_community_comparison` object.
#' @export
print.hypernets_community_comparison <- function(x, ...) {
  cat(sprintf("Community comparison across %d representations: %s\n",
              length(x$models), paste(x$models, collapse = ", ")))
  print(x$summary, row.names = FALSE)
  invisible(x)
}

#' @rdname hg_compare_communities
#' @param row.names,optional Unused; present for the base S3 contract.
#' @param what Which table: `"summary"` (default), `"similarity"`, `"sizes"`,
#'   `"quality"`, or `"matrix"` for the similarity as a square matrix with
#'   AMI below and ARI above the diagonal (the layout of `plot()`).
#' @export
as.data.frame.hypernets_community_comparison <- function(x, row.names = NULL,
                                                      optional = FALSE,
                                                      what = c("summary",
                                                               "similarity",
                                                               "sizes",
                                                               "quality",
                                                               "matrix"), ...) {
  what <- match.arg(what)
  if (identical(what, "matrix")) return(.thg_similarity_matrix(x))
  out <- x[[what]]
  if (is.null(out)) {
    .thg_bad_input("quality scores need the hypergraph: hg_compare_communities(..., hg = )")
  }
  rownames(out) <- NULL
  out
}

#' @rdname hg_compare_communities
#' @param object A `hypernets_community_comparison` object.
#' @export
summary.hypernets_community_comparison <- function(object, ...) {
  as.data.frame(object, what = "summary")
}

#' @rdname hg_compare_communities
#' @param ... For `plot`, unused.
#' @return For `plot`, a ggplot object.
#' @export
plot.hypernets_community_comparison <- function(x, what = c("sizes", "similarity"),
                                             ...) {
  what <- match.arg(what)
  if (identical(what, "sizes")) {
    sizes <- as.data.frame(x, what = "sizes")
    curve <- do.call(rbind, lapply(split(sizes, sizes$model), function(d) {
      s <- sort(d$n_nodes)
      distinct <- unique(s)
      data.frame(model = d$model[[1L]], size = distinct,
                 n_at_least = vapply(distinct, function(v) sum(s >= v), numeric(1L)),
                 stringsAsFactors = FALSE)
    }))
    curve$model <- factor(curve$model, levels = x$models)
    n_models <- length(x$models)
    return(
      ggplot2::ggplot(curve, ggplot2::aes(x = .data$size, y = .data$n_at_least,
                                          colour = .data$model,
                                          linetype = .data$model)) +
        ggplot2::geom_step(direction = "vh", linewidth = 0.7) +
        ggplot2::scale_x_log10() + ggplot2::scale_y_log10() +
        ggplot2::scale_colour_manual(values = rep_len(.thg_okabe_ito, n_models),
                                     name = NULL) +
        ggplot2::scale_linetype_manual(
          values = rep_len(c("solid", "dashed", "dotted", "dotdash",
                             "longdash", "twodash"), n_models), name = NULL) +
        ggplot2::labs(x = "cluster size", y = "number of clusters at least this large") +
        ggplot2::theme_minimal(base_size = 12) +
        ggplot2::theme(legend.position = "right")
    )
  }
  # AMI below the diagonal, ARI above it, as in the paper's Figure 8c;
  # cograph draws the matrix.
  cograph::plot_heatmap(
    as.data.frame(x, what = "matrix"), show_values = TRUE, limits = c(0, 1),
    colors = .thg_okabe_ito_ramp(), na_color = "white", show_diagonal = FALSE,
    legend_title = "AMI (below)\nARI (above)", axis_text_angle = 0,
    value_size = 3
  )
}

# The similarity matrix with AMI below and ARI above the diagonal.
.thg_similarity_matrix <- function(x) {
  models <- x$models
  sim <- x$similarity
  m <- matrix(NA_real_, length(models), length(models),
              dimnames = list(models, models))
  m[cbind(match(sim$model_b, models), match(sim$model_a, models))] <- sim$ami
  m[cbind(match(sim$model_a, models), match(sim$model_b, models))] <- sim$ari
  m
}

#' @rdname hg_communities
#' @export
print.hg_communities <- function(x, ...) {
  cat(sprintf("Hypergraph Infomap ensemble: %d nodes, %d runs\n",
              nrow(x$medoid), nrow(x$runs)))
  cat(sprintf("AMI medoid: run %d (seed %d), %d communities\n",
              x$medoid_run, x$runs$seed[x$medoid_run], nrow(x$sizes)))
  invisible(x)
}

#' @rdname hg_communities
#' @param x An `hg_communities` object.
#' @param row.names,optional Unused; present for the base S3 contract.
#' @param what Component to return: `"medoid"`, `"partitions"`, `"runs"`,
#'   `"sizes"`, `"ami"`, `"ari"`, or `"nmi"`.
#' @export
as.data.frame.hg_communities <- function(x, row.names = NULL, optional = FALSE,
                                         what = c("medoid", "partitions", "runs",
                                                  "sizes", "ami", "ari", "nmi"),
                                         ...) {
  what <- match.arg(what)
  if (what %in% c("ami", "ari", "nmi")) {
    return(as.data.frame(as.table(x$similarity[[what]]),
                         stringsAsFactors = FALSE,
                         responseName = what))
  }
  x[[what]]
}

#' @rdname hg_communities
#' @param ... Additional arguments passed to [cograph::splot()].
#' @export
plot.hg_communities <- function(x, ...) {
  membership <- x$medoid$community[match(rownames(x$projection), x$medoid$node)]
  cograph::splot(x$projection, groups = membership,
                 directed = isTRUE(x$params$directed), ...)
}

#' Quality of a projected-hypergraph partition
#'
#' Reports graph-partition diagnostics on the association projection:
#' unweighted coverage, weighted coverage, performance, weighted modularity,
#' and conductance. Conductance is the maximum (worst) community conductance
#' \eqn{cut(S, \bar S) / min(vol(S), vol(\bar S))}; lower is better.
#'
#' @param hg A static `net_hypergraph`.
#' @param partition An [hg_communities()] result, a tidy node/label table, or a
#'   named label vector.
#' @param method The projection the partition is scored on: the
#'   `"association"` graph (default) or the undirected `"citation"` graph.
#' @inheritParams hg_communities
#' @return A one-row data frame.
#' @export
hg_community_quality <- function(hg, partition,
                                 method = c("association", "citation"),
                                 duplicate_edges = c("count", "collapse"),
                                 self_association = FALSE, edge_source = NULL) {
  .thg_check_hg(hg)
  method <- match.arg(method)
  duplicate_edges <- match.arg(duplicate_edges)
  projection <- if (identical(method, "association")) {
    hg_project(
      hg, method = "association", what = "matrix",
      duplicate_edges = duplicate_edges,
      self_association = self_association, edge_source = edge_source
    )
  } else {
    hg_project(
      hg, method = "citation", what = "matrix",
      duplicate_edges = duplicate_edges, edge_source = edge_source
    )
  }
  projection <- as.matrix(projection)
  nodes <- rownames(projection)
  labels <- .thg_partition_vector(partition, nodes)
  binary <- (projection != 0) * 1
  upper <- upper.tri(binary)
  same <- outer(labels, labels, `==`)
  m <- sum(binary[upper])
  internal <- sum(binary[upper & same])
  total_weight <- sum(projection[upper])
  internal_weight <- sum(projection[upper & same])
  possible <- choose(length(nodes), 2)
  inter_nonedges <- sum(binary[upper & !same] == 0)
  weighted_quality <- cograph::cluster_quality(
    projection, labels, weighted = TRUE, directed = FALSE
  )
  # An undirected self-loop contributes twice to weighted degree. It never
  # crosses a cut, but it does increase the volume on its side.
  strength <- rowSums(projection) + diag(projection)
  total_volume <- sum(strength)
  community_conductance <- vapply(unique(labels), function(group) {
    inside <- labels == group
    volume <- sum(strength[inside])
    denominator <- min(volume, total_volume - volume)
    if (denominator <= 0) return(NA_real_)
    cut <- sum(projection[inside, !inside, drop = FALSE])
    cut / denominator
  }, numeric(1))
  conductance <- if (all(is.na(community_conductance))) NA_real_ else
    max(community_conductance, na.rm = TRUE)
  data.frame(
    coverage = if (m > 0) internal / m else NA_real_,
    weighted_coverage = if (total_weight > 0) internal_weight / total_weight else NA_real_,
    performance = if (possible > 0) (internal + inter_nonedges) / possible else NA_real_,
    modularity = weighted_quality$global$modularity,
    conductance = conductance,
    n_communities = length(unique(labels)), row.names = NULL
  )
}

#' @rdname hg_communities
#' @export
hypergraph_communities <- hg_communities

#' @rdname hg_community_quality
#' @export
hypergraph_community_quality <- hg_community_quality
