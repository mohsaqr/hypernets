# Induced Y/T/O motifs for 3-uniform hypergraphs, with the detailed
# configuration MCMC used by HypergraphX 1.5 in the paper's oracle code.

.thg_unique_edges <- function(hg) {
  b <- hg$incidence != 0
  edges <- lapply(seq_len(ncol(b)), function(j) sort(which(b[, j])))
  signatures <- vapply(edges, paste, collapse = ",", character(1L))
  edges[!duplicated(signatures)]
}

.thg_yto_counts <- function(edges) {
  if (length(edges) < 2L) return(c(Y = 0L, T = 0L, O = 0L))
  four_sets <- character()
  for (i in seq_len(length(edges) - 1L)) {
    for (j in (i + 1L):length(edges)) {
      if (length(intersect(edges[[i]], edges[[j]])) == 2L) {
        union_nodes <- sort(union(edges[[i]], edges[[j]]))
        if (length(union_nodes) == 4L) {
          four_sets <- c(four_sets, paste(union_nodes, collapse = ","))
        }
      }
    }
  }
  if (!length(four_sets)) return(c(Y = 0L, T = 0L, O = 0L))
  pair_count <- as.integer(table(four_sets))
  # A four-node set containing k distinct 3-edges generates choose(k, 2)
  # overlapping pairs: 1, 3, and 6 identify the induced Y, T, and O motifs.
  c(Y = sum(pair_count == 1L),
    T = sum(pair_count == 3L),
    O = sum(pair_count == 6L))
}

.thg_configuration_mcmc <- function(edges, steps = 10L * length(edges),
                                    collapse = TRUE) {
  if (length(edges) < 2L || steps < 1L) return(edges)
  current <- lapply(edges, sort)
  sizes <- lengths(current)
  for (step in seq_len(steps)) {
    # HypergraphX draws two edge positions independently, so a self-pair is
    # permitted (and is simply a no-op after reshuffling).
    ij <- sample.int(length(current), 2L, replace = TRUE)
    while (sizes[ij[1L]] != sizes[ij[2L]]) {
      ij <- sample.int(length(current), 2L, replace = TRUE)
    }
    f1 <- current[[ij[1L]]]
    f2 <- current[[ij[2L]]]
    common <- intersect(f1, f2)
    pool <- c(setdiff(f1, common), setdiff(f2, common))
    g1 <- common
    g2 <- common
    # Match HypergraphX 1.5's capacity-constrained sequence of fair coins,
    # including its original pool order, rather than sampling a subset.
    for (v in pool) {
      if (length(g1) < length(f1) && length(g2) < length(f2)) {
        if (stats::runif(1L) < 0.5) g1 <- c(g1, v) else g2 <- c(g2, v)
      } else if (length(g1) < length(f1)) {
        g1 <- c(g1, v)
      } else {
        g2 <- c(g2, v)
      }
    }
    g1 <- sort(g1)
    g2 <- sort(g2)
    current[[ij[1L]]] <- g1
    current[[ij[2L]]] <- g2
  }
  if (!collapse) return(current)
  signatures <- vapply(current, paste, collapse = ",", character(1L))
  current[!duplicated(signatures)]
}

.thg_restore_seed <- function(seed) {
  if (is.null(seed)) return(invisible(NULL))
  old <- if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    get(".Random.seed", envir = globalenv(), inherits = FALSE)
  } else NULL
  set.seed(seed)
  old
}

#' Y, T and O motif census with a hypergraph configuration null
#'
#' Counts the three possible induced four-node motifs in a simple 3-uniform
#' hypergraph: Y contains two of the four possible triples, T contains three,
#' and O contains all four. Repeated identical hyperedges are collapsed for
#' the motif census, matching HypergraphX 1.5. The null repeatedly reshuffles
#' pairs of equal-size hyperedges while preserving every node degree and edge
#' cardinality before duplicate-edge collapse.
#'
#' @param hg A 3-uniform `net_hypergraph` or a [temporal_hypergraph()].
#' @param n Number of configuration-model draws. The paper uses 1000.
#' @param seed Optional reproducibility seed; the caller's RNG state is
#'   restored on exit.
#' @param what `"test"` for observed/null summaries, `"counts"` for observed
#'   counts only, or `"draws"` for every null count.
#' @param alternative Empirical permutation-test direction.
#' @param at,snapshot_mode,multiedges Temporal snapshot arguments passed to
#'   [hypergraph_snapshots()] when `hg` is temporal.
#' @return A tidy data frame. Test output includes observed count, null mean
#'   and standard deviation, z-score, empirical p-value, relative abundance
#'   `delta`, and the normalized motif profile used by HypergraphX.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#' hypergraphs. *Philosophical Transactions of the Royal Society A*,
#' 382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' dat <- data.frame(
#'   member = c(1, 2, 3, 1, 2, 4), edge = rep(c("e1", "e2"), each = 3)
#' )
#' h <- group_hypergraph(dat, "member", "edge")
#' hg_motifs(h, what = "counts")
#' @export
hg_motifs <- function(hg, n = 1000L, seed = NULL,
                      what = c("test", "counts", "draws"),
                      alternative = c("two_sided", "greater", "less"),
                      at = NULL,
                      snapshot_mode = c("active", "cumulative", "all"),
                      multiedges = TRUE) {
  what <- match.arg(what)
  alternative <- match.arg(alternative)
  snapshot_mode <- match.arg(snapshot_mode)
  if (inherits(hg, "net_temporal_hypergraph")) {
    snaps <- hypergraph_snapshots(hg, at = at, mode = snapshot_mode,
                                  multiedges = multiedges)
    rows <- lapply(seq_along(snaps), function(i) {
      ans <- hg_motifs(
        snaps[[i]], n = n,
        seed = if (is.null(seed)) NULL else as.integer(seed) + i - 1L,
        what = what, alternative = alternative
      )
      data.frame(time = names(snaps)[i], ans, row.names = NULL,
                 stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    rownames(out) <- NULL
    return(out)
  }
  .thg_check_hg(hg)
  edges <- .thg_unique_edges(hg)
  if (length(edges) && any(lengths(edges) != 3L)) {
    .thg_bad_input("Y/T/O motifs require a 3-uniform hypergraph")
  }
  observed <- .thg_yto_counts(edges)
  if (identical(what, "counts")) {
    return(data.frame(motif = names(observed), count = as.integer(observed),
                      row.names = NULL))
  }
  if (length(n) != 1L || !is.finite(n) || n < 2L ||
      abs(n - round(n)) > sqrt(.Machine$double.eps)) {
    .thg_bad_input("`n` must be a whole number of at least 2")
  }
  n <- as.integer(n)

  old_seed <- .thg_restore_seed(seed)
  if (!is.null(seed)) on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv(), inherits = FALSE))
        rm(".Random.seed", envir = globalenv())
    } else assign(".Random.seed", old_seed, envir = globalenv())
  }, add = TRUE)

  null <- vapply(seq_len(n), function(i) {
    .thg_yto_counts(.thg_configuration_mcmc(edges))
  }, numeric(3L))
  rownames(null) <- names(observed)
  if (identical(what, "draws")) {
    return(data.frame(
      run = rep(seq_len(n), each = 3L), motif = rep(names(observed), n),
      count = as.numeric(null), row.names = NULL
    ))
  }

  null_mean <- rowMeans(null)
  null_sd <- sqrt(rowMeans((null - null_mean)^2))
  z <- ifelse(null_sd > 0, (observed - null_mean) / null_sd, NA_real_)
  extreme <- vapply(seq_along(observed), function(i) {
    switch(alternative,
      two_sided = sum(abs(null[i, ] - null_mean[i]) >=
                        abs(observed[i] - null_mean[i]) - sqrt(.Machine$double.eps)),
      greater = sum(null[i, ] >= observed[i]),
      less = sum(null[i, ] <= observed[i])
    )
  }, integer(1L))
  delta <- (observed - null_mean) / (observed + null_mean + 4)
  delta_norm <- sqrt(sum(delta^2))
  normalized_delta <- if (delta_norm > 0) delta / delta_norm else rep(0, 3L)
  data.frame(
    motif = names(observed), observed = as.integer(observed),
    null_mean = as.numeric(null_mean), null_sd = as.numeric(null_sd),
    z = as.numeric(z), p_value = (1 + extreme) / (n + 1),
    delta = as.numeric(delta), normalized_delta = as.numeric(normalized_delta),
    n = n, method = "configuration_mcmc", row.names = NULL
  )
}

#' @rdname hg_motifs
#' @export
hypergraph_motifs <- hg_motifs
