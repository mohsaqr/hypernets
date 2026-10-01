# Partition agreement, clustering stability, and seed-label construction.

# Adjusted Rand index (Hubert & Arabie 1985) from two label vectors.
.thg_ari <- function(a, b) {
  tab <- table(a, b)
  n <- sum(tab)
  if (n < 2L) return(1)
  sum_ij <- sum(choose(tab, 2))
  sum_a <- sum(choose(rowSums(tab), 2))
  sum_b <- sum(choose(colSums(tab), 2))
  expected <- sum_a * sum_b / choose(n, 2)
  denom <- (sum_a + sum_b) / 2 - expected
  # both partitions trivial (single cluster): identical by construction
  if (abs(denom) < sqrt(.Machine$double.eps)) return(1)
  (sum_ij - expected) / denom
}

# Information-theoretic partition agreement. The adjusted variant follows
# Vinh et al. (2010) and scikit-learn's arithmetic-mean normalization.
.thg_entropy <- function(margin) {
  p <- margin[margin > 0] / sum(margin)
  -sum(p * log(p))
}

.thg_mutual_information <- function(tab) {
  n <- sum(tab)
  rows <- rowSums(tab)
  cols <- colSums(tab)
  nz <- which(tab > 0, arr.ind = TRUE)
  if (!nrow(nz)) return(0)
  sum(vapply(seq_len(nrow(nz)), function(k) {
    i <- nz[k, 1L]
    j <- nz[k, 2L]
    nij <- tab[i, j]
    (nij / n) * log((nij * n) / (rows[i] * cols[j]))
  }, numeric(1L)))
}

.thg_expected_mi <- function(tab) {
  n <- sum(tab)
  if (n < 2L) return(0)
  rows <- rowSums(tab)
  cols <- colSums(tab)
  emi <- 0
  for (ai in rows) {
    for (bj in cols) {
      lo <- max(1, ai + bj - n)
      hi <- min(ai, bj)
      if (lo > hi) next
      nij <- seq.int(lo, hi)
      probability <- stats::dhyper(nij, ai, n - ai, bj)
      emi <- emi + sum(probability * (nij / n) *
                         log((nij * n) / (ai * bj)))
    }
  }
  emi
}

.thg_ami <- function(a, b) {
  tab <- table(a, b)
  mi <- .thg_mutual_information(tab)
  emi <- .thg_expected_mi(tab)
  normalizer <- (.thg_entropy(rowSums(tab)) +
                   .thg_entropy(colSums(tab))) / 2
  denominator <- normalizer - emi
  if (abs(denominator) < .Machine$double.eps) return(1)
  (mi - emi) / denominator
}

.thg_nmi <- function(a, b) {
  tab <- table(a, b)
  normalizer <- (.thg_entropy(rowSums(tab)) +
                   .thg_entropy(colSums(tab))) / 2
  if (normalizer < .Machine$double.eps) return(1)
  .thg_mutual_information(tab) / normalizer
}

# Extract the label column from a tidy labeling: `predicted`
# (classification results) first, else `cluster` / `community`
# (clustering results: hg_cluster(); hg_communities(), hg_mmsbm()).
.thg_labeling <- function(x, arg, node = "node", label = NULL) {
  stopifnot(
    "labelings must be data.frames" = is.data.frame(x),
    "`node` must be a single column name" =
      is.character(node) && length(node) == 1L && !is.na(node),
    "`label` must be NULL or a single column name" =
      is.null(label) || (is.character(label) && length(label) == 1L && !is.na(label))
  )
  if (!node %in% names(x)) {
    stop(errorCondition(
      sprintf("`%s` has no `%s` column", arg, node),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  column <- if (is.null(label)) {
    intersect(c("predicted", "cluster", "community", "label"), names(x))
  } else {
    intersect(label, names(x))
  }
  if (length(column) == 0L) {
    stop(errorCondition(
      if (is.null(label)) {
        sprintf(paste0("`%s` has no `predicted`, `cluster`, `community` or `label` ",
                       "column; pass a result from hg_cluster(), ",
                       "hg_classify(), hg_neural() or hg_hypergat(), a ",
                       "table of known labels, or name the column with `label`"), arg)
      } else {
        sprintf("`%s` has no `%s` column", arg, label)
      },
      class = "hypernets_bad_input", call = NULL
    ))
  }
  data.frame(node = as.character(x[[node]]),
             label = as.character(x[[column[[1]]]]),
             stringsAsFactors = FALSE)
}

# A column selector for `hg_agreement()`: one name applies to both
# labelings, two names apply to `x` and `y` in turn.
.thg_pair_selector <- function(value, arg) {
  if (is.null(value)) return(list(NULL, NULL))
  if (!is.character(value) || !length(value) %in% c(1L, 2L) || anyNA(value)) {
    .thg_bad_input(sprintf("`%s` must be one column name, or two (for `x` and `y`)", arg))
  }
  if (length(value) == 1L) list(value, value) else list(value[[1L]], value[[2L]])
}

# Coerce a `labels` argument to the named-character-vector contract the
# classifiers use internally. Accepts a named vector as-is, or a tidy
# data.frame with `node` plus a `label`, `cluster` or `predicted` column
# (first match wins), so a clustering or a subset of a corpus table can
# be passed directly.
.thg_labels_input <- function(labels) {
  if (!is.data.frame(labels)) return(labels)
  stopifnot(
    "a `labels` data.frame needs a `node` column and a `label`, `cluster`, `community` or `predicted` column" =
      "node" %in% names(labels) &&
        length(intersect(c("label", "cluster", "community", "predicted"),
                         names(labels))) > 0L
  )
  column <- intersect(c("label", "cluster", "community", "predicted"), names(labels))
  stats::setNames(as.character(labels[[column[[1]]]]),
                  as.character(labels$node))
}

#' Agreement between two labelings of the same nodes
#'
#' Compares any two tidy labelings -- partitions from [hg_cluster()],
#' predictions from [hg_classify()], [hg_neural()] or [hg_hypergat()] --
#' joined on their shared `node` column. `agreement` is the share of
#' nodes with literally equal labels (meaningful when both labelings use
#' the same label set, e.g. a classifier scored against the clustering
#' that produced its seeds); `ari` is the adjusted Rand index (Hubert &
#' Arabie 1985), which is label-permutation invariant and the right
#' statistic when the two label sets are arbitrary (e.g. two independent
#' clusterings). Adjusted mutual information (`ami`) and normalized mutual
#' information (`nmi`) are also available; the legal-hypergraphs workflow uses
#' AMI to select the medoid of repeated Infomap partitions.
#'
#' @param x,y Tidy labelings: data.frames with a `node` column and a
#'   `predicted`, `cluster` or `label` column (first match in that
#'   order wins). Nodes are matched by name; nodes present in only one
#'   labeling are dropped.
#' @param node,label Column names, one name for both labelings or two
#'   names for `x` and `y` in turn, that override the defaults above --
#'   `hg_agreement(predictions, corpus, node = c("node", "doc"), label =
#'   c("predicted", "year"))` scores a classifier against a column of the
#'   corpus table without reshaping it. `label = NULL` (default) keeps the
#'   `predicted` / `cluster` / `label` lookup.
#' @param what `"summary"` (default) for the one-row comparison,
#'   `"table"` for the tidy contingency table of the joined labels, or
#'   `"mapping"` for one row per label of `x` naming the label of `y` that
#'   holds most of its nodes -- how a partition survives a reweighting or a
#'   change of method, topic by topic.
#' @param method One or more label-permutation-invariant measures: `"ari"`
#'   (default), `"ami"`, or `"nmi"`. Ignored for `what = "table"`.
#' @return A base `data.frame`. For `what = "summary"`: one row with
#'   columns `n` (nodes compared), `agreement` (share of equal labels),
#'   `aligned` (nodes that stay with the majority of their `x` label in
#'   `y` -- the sum of `overlap` over `what = "mapping"`, so a
#'   label-name-free count of how many nodes a re-fit keeps together)
#'   and the requested measure columns. For `what = "table"`: one row per
#'   label pair with columns `label_x`, `label_y` and `n`. For
#'   `what = "mapping"`: one row per label of `x` with columns `label_x`,
#'   `n` (its nodes), `label_y` (the label of `y` holding most of them),
#'   `overlap` (how many) and `share` (`overlap / n`), in the natural order
#'   of `label_x`.
#' @references Hubert, L., & Arabie, P. (1985). Comparing partitions.
#'   *Journal of Classification*, 2, 193--218.
#'
#'   Vinh, N. X., Epps, J., & Bailey, J. (2010). Information theoretic
#'   measures for clusterings comparison. *Journal of Machine Learning
#'   Research*, 11, 2837--2854.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' fit <- hg_classify(hg, labels = c(cooking_1 = "Cluster 1",
#'                                   space_1 = "Cluster 2"))
#' hg_agreement(fit, topics)
#' hg_agreement(fit, topics, what = "table")
#' hg_agreement(fit, topics, what = "mapping")
#' # a labeling read from any table: name its node and label columns
#' known <- data.frame(doc = c("cooking_1", "space_2"), theme = c("cooking", "space"))
#' hg_agreement(topics, known, node = c("node", "doc"), label = c("cluster", "theme"),
#'              what = "table")
#' @export
hg_agreement <- function(x, y, what = c("summary", "table", "mapping"),
                         method = "ari", node = "node", label = NULL) {
  what <- match.arg(what)
  method <- match.arg(method, c("ari", "ami", "nmi"), several.ok = TRUE)
  node <- .thg_pair_selector(node, "node")
  label <- .thg_pair_selector(label, "label")
  joined <- merge(.thg_labeling(x, "x", node[[1L]], label[[1L]]),
                  .thg_labeling(y, "y", node[[2L]], label[[2L]]),
                  by = "node", suffixes = c("_x", "_y"))
  unlabelled <- is.na(joined$label_x) | is.na(joined$label_y)
  if (any(unlabelled)) {
    # table()-based indices drop NA labels; drop them for every column alike
    warning(warningCondition(sprintf(
      "%d node(s) without a label in `x` or `y` are left out of the comparison",
      sum(unlabelled)), class = "hypernets_missing_labels", call = NULL))
    joined <- joined[!unlabelled, , drop = FALSE]
  }
  if (nrow(joined) == 0L) {
    stop(errorCondition(
      "`x` and `y` share no node names; nothing to compare",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  if (identical(what, "table")) {
    counts <- stats::aggregate(node ~ label_x + label_y, data = joined,
                               length)
    names(counts) <- c("label_x", "label_y", "n")
    ord <- counts[order(counts$label_x, counts$label_y), , drop = FALSE]
    rownames(ord) <- NULL
    return(ord)
  }
  if (identical(what, "mapping")) {
    cross <- table(joined$label_x, joined$label_y)
    best <- apply(cross, 1L, which.max)
    out <- data.frame(
      label_x = rownames(cross),
      n = as.integer(rowSums(cross)),
      label_y = colnames(cross)[best],
      overlap = as.integer(apply(cross, 1L, max)),
      stringsAsFactors = FALSE
    )
    out$share <- out$overlap / out$n
    out <- out[order(match(out$label_x, .thg_kw_natural(out$label_x))), ,
               drop = FALSE]
    rownames(out) <- NULL
    return(out)
  }
  cross <- table(joined$label_x, joined$label_y)
  out <- data.frame(
    n = nrow(joined),
    agreement = mean(joined$label_x == joined$label_y),
    aligned = as.integer(sum(apply(cross, 1L, max)))
  )
  if ("ari" %in% method) out$ari <- .thg_ari(joined$label_x, joined$label_y)
  if ("ami" %in% method) out$ami <- .thg_ami(joined$label_x, joined$label_y)
  if ("nmi" %in% method) out$nmi <- .thg_nmi(joined$label_x, joined$label_y)
  out
}

#' Stability of a hypergraph clustering across resolutions
#'
#' Asks whether the partition [hg_cluster()] finds is a property of the data
#' or of the particular sample, one resolution `k` at a time.
#'
#' **`resample = "subset"` (default)** is cluster-wise subsampling stability
#' (Hennig 2007; the subsetting scheme of Ben-Hur et al. 2002). The nodes
#' (on a document hypergraph, the documents) are clustered once in full.
#' Then, `n_boot` times, a subsample of `floor(fraction * n)` nodes is drawn
#' without replacement, the sub-hypergraph is formed by restricting every
#' hyperedge to the drawn nodes (hyperedges left empty are dropped; the
#' incidence weights are kept, so a tf-idf corpus keeps the full corpus's
#' idf), and it is clustered again with the same `k`, `type`, `seed` and
#' `nstart`. Each original cluster, restricted to the drawn nodes, is matched
#' to its most similar subsample cluster by the Jaccard coefficient
#' \eqn{|A \cap B| / |A \cup B|} (0 when the cluster has no drawn member). A
#' cluster's stability is its mean Jaccard over the subsamples. Hennig
#' (2007) reads a mean of 0.5 or below as a dissolved cluster and 0.75 or
#' above as a stable one. The subsamples are drawn once, up front, after
#' `set.seed(seed)`, so they are the ones `fpc::clusterboot(bootmethod =
#' "subset", subtuning = floor(fraction * n), seed = seed)` draws, and the
#' same subsamples serve every `k`. The caller's random number stream is
#' restored on exit.
#'
#' A subsample can disconnect the hypergraph (a word that tied two groups of
#' documents may not be drawn); [hg_cluster()] cannot cut a disconnected
#' hypergraph, so that run is not scored, is counted in `n_failed`, and a
#' warning of class `hypernets_hypergraph_disconnected` reports the count.
#'
#' The `eigengap` column is the gap \eqn{\lambda_{k+1} - \lambda_k} between
#' consecutive eigenvalues of the full hypergraph's Laplacian (the `gap`
#' column of `hg_cluster(what = "eigenvalues")`); the eigengap heuristic
#' (von Luxburg 2007, Sec. 8.1) prefers the `k` where it is large relative
#' to the gaps before it. Read it beside the stability: a `k` with a large
#' gap and stable clusters is supported twice.
#'
#' **`resample = "seeds"`** keeps the check this verb made before 0.6.0: it
#' fits the full hypergraph twice per `k` with the two k-means `seeds` and
#' reports whether the partitions agree. It measures whether the solver is
#' deterministic on this embedding, not whether the structure survives a
#' change of sample, and it never resamples the data.
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hg`).
#' @param k Vector of cluster counts to test, each at least 2.
#' @param type `"zhou"` or `"random_walk"`, the Laplacian of [hg_cluster()].
#' @param resample `"subset"` (default) for subsampling stability, or
#'   `"seeds"` for the two-seed solver check.
#' @param n_boot Number of subsamples (default `100L`).
#' @param fraction Share of the nodes in each subsample, strictly between 0
#'   and 1 (default `0.5`, as in `fpc::clusterboot()`).
#' @param seed Seed for the subsample draws and the k-means starts of every
#'   fit (default `1L`).
#' @param what For `resample = "subset"`: `"k"` (default) for one row per
#'   resolution, or `"clusters"` for one row per resolution and cluster.
#' @param seeds For `resample = "seeds"`: two distinct k-means seeds
#'   (default `c(1L, 99L)`).
#' @param nstart Number of k-means starts per fit (default `25L`).
#' @return A base `data.frame`.
#'
#'   `resample = "subset"`, `what = "k"`: one row per resolution with `k`,
#'   `n_runs` (subsamples scored), `n_failed` (subsamples that disconnected
#'   the hypergraph), `mean_jaccard` and `min_jaccard` (the mean and the
#'   smallest of the clusters' stabilities), `n_stable` (clusters with
#'   stability at least 0.75), `n_dissolved` (clusters with stability 0.5 or
#'   below) and `eigengap`.
#'
#'   `resample = "subset"`, `what = "clusters"`: one row per resolution and
#'   cluster of the full partition, with `k`, `cluster`, `size` (its nodes),
#'   `jaccard` (its stability, `fpc::clusterboot()`'s `subsetmean`),
#'   `n_dissolved` (subsamples with Jaccard 0.5 or below, `subsetbrd`),
#'   `n_recovered` (subsamples with Jaccard above 0.75, `subsetrecover`) and
#'   `n_runs`.
#'
#'   `resample = "seeds"`: one row per resolution with `k`,
#'   `identical_partition` (are the two partitions literally identical,
#'   labels included) and `ari` (their adjusted Rand index).
#'
#'   Raises `hypernets_bad_input` for a `fraction` outside (0, 1), a
#'   subsample too small to hold `k` clusters, a non-positive `n_boot`, or
#'   `what = "clusters"` with `resample = "seeds"`.
#' @references
#' Hennig, C. (2007). Cluster-wise assessment of cluster stability.
#' *Computational Statistics & Data Analysis*, 52(1), 258--271.
#' \doi{10.1016/j.csda.2006.11.025}
#'
#' Ben-Hur, A., Elisseeff, A., & Guyon, I. (2002). A stability based method
#' for discovering structure in clustered data. *Pacific Symposium on
#' Biocomputing*, 7, 6--17. \doi{10.1142/9789812799623_0002}
#'
#' von Luxburg, U. (2007). A tutorial on spectral clustering. *Statistics and
#' Computing*, 17(4), 395--416. \doi{10.1007/s11222-007-9033-z}
#'
#' Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
#' Classification*, 2, 193--218.
#' @examples
#' hg <- text_hypergraph(head(covid_abstracts, 40), column = "abstract",
#'                       id = "doc", stop_words = stop_words_en(),
#'                       min_count = 3)
#' hg_stability(hg, k = 2:4, n_boot = 20)
#' hg_stability(hg, k = 3, n_boot = 20, what = "clusters")
#' hg_stability(hg, k = 2, resample = "seeds")
#' @export
hg_stability <- function(hg, k, type = c("zhou", "random_walk"),
                         resample = c("subset", "seeds"), n_boot = 100L,
                         fraction = 0.5, seed = 1L,
                         what = c("k", "clusters"),
                         seeds = c(1L, 99L), nstart = 25L) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  resample <- match.arg(resample)
  what <- match.arg(what)
  stopifnot(
    "`k` must be a vector of cluster counts, each at least 2" =
      is.numeric(k) && length(k) >= 1L && all(is.finite(k)) && all(k >= 2)
  )
  if (identical(resample, "seeds")) {
    stopifnot(
      "`seeds` must be two distinct seeds" =
        is.numeric(seeds) && length(seeds) == 2L &&
          !isTRUE(all.equal(seeds[[1]], seeds[[2]]))
    )
    if (identical(what, "clusters")) {
      .thg_bad_input(paste0("`what = \"clusters\"` needs `resample = ",
                            "\"subset\"`; the seed check has one row per k"))
    }
    return(.thg_seed_stability(hg, k, type, seeds, nstart))
  }
  if (!is.numeric(fraction) || length(fraction) != 1L || is.na(fraction) ||
      fraction <= 0 || fraction >= 1) {
    .thg_bad_input("`fraction` must be a single number strictly between 0 and 1")
  }
  if (!is.numeric(n_boot) || length(n_boot) != 1L || is.na(n_boot) ||
      n_boot < 1) {
    .thg_bad_input("`n_boot` must be a single number >= 1")
  }
  if (!is.numeric(seed) || length(seed) != 1L || is.na(seed)) {
    .thg_bad_input("`seed` must be a single number")
  }
  n <- hg$n_nodes
  m <- floor(fraction * n)
  if (m <= max(k)) {
    .thg_bad_input(sprintf(
      "a subsample of floor(%s * %d) = %d nodes cannot hold k = %d clusters; raise `fraction`",
      format(fraction), n, m, max(k)))
  }

  # the subsamples, drawn once from `seed` and shared by every k; the
  # caller's stream is restored on exit (hg_cluster() also calls set.seed())
  had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  saved_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
  on.exit(.thg_rng_restore(had_seed, saved_seed), add = TRUE)
  set.seed(seed)
  draws <- lapply(seq_len(n_boot), \(b) sample(n, m))

  full <- lapply(k, \(kk) {
    fit <- hg_cluster(hg, k = kk, type = type, seed = seed, nstart = nstart)
    stats::setNames(fit$cluster, fit$node)[hg$nodes]
  })
  spectrum <- hg_cluster(hg, k = max(k), type = type, seed = seed,
                         nstart = nstart, what = "eigenvalues")

  # one list per subsample: for every k, the Jaccard of each full cluster
  # with its best-matching subsample cluster, or NULL when disconnected
  runs <- lapply(draws, \(draw) {
    sub <- .thg_restrict_nodes(hg, draw)
    lapply(seq_along(k), \(i) {
      fit <- tryCatch(
        hg_cluster(sub, k = k[[i]], type = type, seed = seed,
                   nstart = nstart),
        hypernets_hypergraph_disconnected = function(e) NULL
      )
      if (is.null(fit)) return(NULL)
      .thg_best_jaccard(full[[i]][sub$nodes],
                        stats::setNames(fit$cluster, fit$node)[sub$nodes],
                        .thg_kw_natural(unique(full[[i]])))
    })
  })

  per_k <- lapply(seq_along(k), \(i) {
    labels <- .thg_kw_natural(unique(full[[i]]))
    scored <- Filter(Negate(is.null), lapply(runs, `[[`, i))
    jaccard <- if (length(scored) > 0L) {
      matrix(unlist(scored), nrow = length(labels))
    } else {
      matrix(numeric(0), nrow = length(labels), ncol = 0L)
    }
    data.frame(
      k = k[[i]], cluster = labels,
      size = as.integer(table(factor(full[[i]], levels = labels))),
      jaccard = if (ncol(jaccard) > 0L) rowMeans(jaccard) else NA_real_,
      n_dissolved = as.integer(rowSums(jaccard <= 0.5)),
      n_recovered = as.integer(rowSums(jaccard > 0.75)),
      n_runs = ncol(jaccard),
      stringsAsFactors = FALSE
    )
  })
  n_failed <- vapply(per_k, \(d) n_boot - d$n_runs[[1L]], numeric(1))
  if (any(n_failed > 0)) {
    warning(warningCondition(
      sprintf(paste0("%s of %d subsamples disconnected the hypergraph and ",
                     "were not scored (k = %s)"),
              paste(unique(n_failed), collapse = "/"), as.integer(n_boot),
              paste(k[n_failed > 0], collapse = ", ")),
      class = "hypernets_hypergraph_disconnected", call = NULL
    ))
  }
  if (identical(what, "clusters")) {
    out <- do.call(rbind, per_k)
    rownames(out) <- NULL
    return(out)
  }
  out <- data.frame(
    k = k,
    n_runs = vapply(per_k, \(d) d$n_runs[[1L]], integer(1)),
    n_failed = as.integer(n_failed),
    mean_jaccard = vapply(per_k, \(d) mean(d$jaccard), numeric(1)),
    min_jaccard = vapply(per_k, \(d) min(d$jaccard), numeric(1)),
    n_stable = vapply(per_k, \(d) sum(d$jaccard >= 0.75), integer(1)),
    n_dissolved = vapply(per_k, \(d) sum(d$jaccard <= 0.5), integer(1)),
    eigengap = spectrum$gap[match(k, spectrum$index)]
  )
  rownames(out) <- NULL
  out
}

# the pre-0.6.0 check: two k-means seeds on the full hypergraph
.thg_seed_stability <- function(hg, k, type, seeds, nstart) {
  rows <- lapply(k, \(kk) {
    a <- hg_cluster(hg, k = kk, type = type, seed = seeds[[1]],
                    nstart = nstart)
    b <- hg_cluster(hg, k = kk, type = type, seed = seeds[[2]],
                    nstart = nstart)
    joined <- merge(a, b, by = "node", suffixes = c("_a", "_b"))
    data.frame(k = kk,
               identical_partition = identical(a$cluster, b$cluster),
               ari = .thg_ari(joined$cluster_a, joined$cluster_b))
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  out
}

# Put back the random number stream saved as (had_seed, saved_seed).
.thg_rng_restore <- function(had_seed, saved_seed) {
  if (had_seed) {
    assign(".Random.seed", saved_seed, envir = globalenv())
  } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
    rm(".Random.seed", envir = globalenv())
  }
}

# The sub-hypergraph on the nodes `index` (in that order): every hyperedge
# restricted to them, hyperedges left empty dropped, weights kept.
.thg_restrict_nodes <- function(hg, index) {
  incidence <- hg$incidence[index, , drop = FALSE]
  keep_edge <- as.numeric(Matrix::colSums(incidence != 0)) > 0
  .thg_rebuild(hg, incidence[, keep_edge, drop = FALSE], keep_edge)
}

# For each cluster label of `full` (named by node, restricted to the
# subsample), the largest Jaccard coefficient with a cluster of `sub`
# (Hennig 2007); 0 when the cluster has no member in the subsample.
.thg_best_jaccard <- function(full, sub, labels) {
  sub_labels <- unique(sub)
  vapply(labels, \(j) {
    a <- full == j
    max(vapply(sub_labels, \(l) {
      b <- sub == l
      union <- sum(a | b)
      if (union == 0L) 0 else sum(a & b) / union
    }, numeric(1)), 0)
  }, numeric(1))
}

#' Seed labels from a clustering, for spreading or training
#'
#' Turns an unsupervised partition into the labeled examples that
#' [hg_classify()], [hg_neural()] and [hg_hypergat()] take: from each
#' cluster, the `n` nodes with the highest stationary probability `pi`
#' (the cluster's most representative members under the random walk).
#' Ties are broken alphabetically by node name so the selection is
#' deterministic.
#'
#' @param embedding The data.frame from
#'   `hg_cluster(what = "embedding")` -- columns `node`, `cluster`, `pi`.
#' @param n Seeds per cluster (default `5L`); clusters smaller than `n`
#'   contribute all their nodes.
#' @return A named character vector -- names are node identifiers, values
#'   their cluster labels -- ready to pass as the `labels` argument of
#'   [hg_classify()], [hg_neural()] or [hg_hypergat()].
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' pool <- hg_cluster(hg, k = 2, seed = 1, what = "embedding")
#' hg_seeds(pool, n = 1)
#' @export
hg_seeds <- function(embedding, n = 5L) {
  stopifnot(
    "`embedding` must be the data.frame from hg_cluster(what = \"embedding\") (columns `node`, `cluster`, `pi`)" =
      is.data.frame(embedding) &&
        all(c("node", "cluster", "pi") %in% names(embedding)),
    "`n` must be a single positive integer" =
      length(n) == 1L && is.finite(n) && n >= 1
  )
  chosen <- do.call(rbind, lapply(split(embedding, embedding$cluster), \(g)
    utils::head(g[order(-g$pi, g$node), , drop = FALSE], n)))
  stats::setNames(as.character(chosen$cluster), chosen$node)
}
