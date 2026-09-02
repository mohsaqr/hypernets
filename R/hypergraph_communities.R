# Repeated Infomap on the association projection, following the paper's
# representation-comparison workflow. honets owns the hypergraph projection;
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
#' @param duplicate_edges,self_association,edge_source Association-projection
#'   controls passed to [hg_project()]. Together these reproduce the paper's
#'   binary/multi-hypergraph and self-association representations.
#' @return An `hg_communities` object containing `medoid` (a tidy node/community
#'   table), all `partitions`, AMI/ARI/NMI similarity matrices, run metadata,
#'   community sizes, and the association `projection`.
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
                           duplicate_edges = c("count", "collapse"),
                           self_association = FALSE, edge_source = NULL) {
  .thg_check_hg(hg)
  duplicate_edges <- match.arg(duplicate_edges)
  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop(errorCondition(
      "`hg_communities()` needs the suggested package `igraph` for Infomap",
      class = "honets_missing_dependency", call = NULL
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

  projection <- hg_project(
    hg, method = "association", what = "matrix",
    duplicate_edges = duplicate_edges,
    self_association = self_association, edge_source = edge_source
  )
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
        projection, nb.trials = trials, seed = seeds[i]
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
                  duplicate_edges = duplicate_edges,
                  self_association = self_association)
  ), class = "hg_communities")
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
  cograph::splot(x$projection, groups = membership, directed = FALSE, ...)
}

#' Quality of a projected-hypergraph partition
#'
#' Reports the four graph-partition diagnostics used in the paper: unweighted
#' coverage, weighted coverage, performance, and weighted modularity.
#'
#' @param hg A static `net_hypergraph`.
#' @param partition An [hg_communities()] result, a tidy node/label table, or a
#'   named label vector.
#' @inheritParams hg_communities
#' @return A one-row data frame.
#' @export
hg_community_quality <- function(hg, partition,
                                 duplicate_edges = c("count", "collapse"),
                                 self_association = FALSE, edge_source = NULL) {
  .thg_check_hg(hg)
  duplicate_edges <- match.arg(duplicate_edges)
  projection <- as.matrix(hg_project(
    hg, method = "association", what = "matrix",
    duplicate_edges = duplicate_edges,
    self_association = self_association, edge_source = edge_source
  ))
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
  data.frame(
    coverage = if (m > 0) internal / m else NA_real_,
    weighted_coverage = if (total_weight > 0) internal_weight / total_weight else NA_real_,
    performance = if (possible > 0) (internal + inter_nonedges) / possible else NA_real_,
    modularity = weighted_quality$global$modularity,
    n_communities = length(unique(labels)), row.names = NULL
  )
}

#' @rdname hg_communities
#' @export
hypergraph_communities <- hg_communities

#' @rdname hg_community_quality
#' @export
hypergraph_community_quality <- hg_community_quality
