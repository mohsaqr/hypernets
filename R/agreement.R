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
# (classification results) first, else `cluster` (clustering results).
.thg_labeling <- function(x, arg) {
  stopifnot(
    "labelings must be data.frames with a `node` column" =
      is.data.frame(x) && "node" %in% names(x)
  )
  column <- intersect(c("predicted", "cluster", "label"), names(x))
  if (length(column) == 0L) {
    stop(errorCondition(
      sprintf(paste0("`%s` has no `predicted`, `cluster` or `label` ",
                     "column; pass a result from hg_cluster(), ",
                     "hg_classify(), hg_neural() or hg_hypergat(), or a ",
                     "table of known labels"), arg),
      class = "honets_bad_input", call = NULL
    ))
  }
  data.frame(node = x$node, label = as.character(x[[column[[1]]]]),
             stringsAsFactors = FALSE)
}

# Coerce a `labels` argument to the named-character-vector contract the
# classifiers use internally. Accepts a named vector as-is, or a tidy
# data.frame with `node` plus a `label`, `cluster` or `predicted` column
# (first match wins), so a clustering or a subset of a corpus table can
# be passed directly.
.thg_labels_input <- function(labels) {
  if (!is.data.frame(labels)) return(labels)
  stopifnot(
    "a `labels` data.frame needs a `node` column and a `label`, `cluster` or `predicted` column" =
      "node" %in% names(labels) &&
        length(intersect(c("label", "cluster", "predicted"),
                         names(labels))) > 0L
  )
  column <- intersect(c("label", "cluster", "predicted"), names(labels))
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
#' @param what `"summary"` (default) for the one-row comparison, or
#'   `"table"` for the tidy contingency table of the joined labels.
#' @param method One or more label-permutation-invariant measures: `"ari"`
#'   (default), `"ami"`, or `"nmi"`. Ignored for `what = "table"`.
#' @return A base `data.frame`. For `what = "summary"`: one row with
#'   columns `n` (nodes compared), `agreement` (share of equal labels)
#'   and the requested measure columns. For `what = "table"`: one row per label pair with
#'   columns `label_x`, `label_y` and `n`.
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
#' @export
hg_agreement <- function(x, y, what = c("summary", "table"), method = "ari") {
  what <- match.arg(what)
  method <- match.arg(method, c("ari", "ami", "nmi"), several.ok = TRUE)
  joined <- merge(.thg_labeling(x, "x"), .thg_labeling(y, "y"),
                  by = "node", suffixes = c("_x", "_y"))
  if (nrow(joined) == 0L) {
    stop(errorCondition(
      "`x` and `y` share no node names; nothing to compare",
      class = "honets_bad_input", call = NULL
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
  out <- data.frame(
    n = nrow(joined),
    agreement = mean(joined$label_x == joined$label_y)
  )
  if ("ari" %in% method) out$ari <- .thg_ari(joined$label_x, joined$label_y)
  if ("ami" %in% method) out$ami <- .thg_ami(joined$label_x, joined$label_y)
  if ("nmi" %in% method) out$nmi <- .thg_nmi(joined$label_x, joined$label_y)
  out
}

#' Seed stability of a hypergraph clustering across resolutions
#'
#' Fits [hg_cluster()] twice per requested `k` with different k-means
#' seeds and reports whether the partitions agree -- the reproducibility
#' criterion for choosing a resolution: a partition that changes with
#' the seed is not estimable from the data, whatever its eigengap looks
#' like.
#'
#' @param hg A [text_hypergraph()] (or any honets `net_hypergraph`).
#' @param k Vector of cluster counts to test, each at least 2.
#' @param type `"zhou"` or `"random_walk"`, as in [hg_cluster()].
#' @param seeds Two distinct k-means seeds (default `c(1L, 99L)`).
#' @param nstart Number of k-means starts per fit (default `25L`).
#' @return A base `data.frame`, one row per resolution, with columns `k`,
#'   `identical_partition` (are the two partitions literally identical,
#'   labels included) and `ari` (their adjusted Rand index; 1 means the
#'   same partition up to label names).
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' hg_stability(hg, k = 2, type = "random_walk")
#' @export
hg_stability <- function(hg, k, type = c("zhou", "random_walk"),
                         seeds = c(1L, 99L), nstart = 25L) {
  .thg_check_hg(hg)
  type <- match.arg(type)
  stopifnot(
    "`k` must be a vector of cluster counts, each at least 2" =
      is.numeric(k) && length(k) >= 1L && all(is.finite(k)) && all(k >= 2),
    "`seeds` must be two distinct seeds" =
      is.numeric(seeds) && length(seeds) == 2L &&
        !isTRUE(all.equal(seeds[[1]], seeds[[2]]))
  )
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

# Long-form aliases. Keep these as direct bindings so both public names have
# identical formals, bodies, and behavior without maintaining wrappers.
#' @rdname hg_agreement
#' @export
hypergraph_agreement <- hg_agreement

#' @rdname hg_stability
#' @export
hypergraph_stability <- hg_stability

#' @rdname hg_seeds
#' @export
hypergraph_seeds <- hg_seeds
