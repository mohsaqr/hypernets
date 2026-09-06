# Topic-level summaries of a clustered document hypergraph: sizes, quality
# (coherence and exclusivity of the owned vocabulary) and soft membership
# from the spectral embedding. Each verb returns a base data.frame with a
# class that carries a plot method; the data.frame is the data.

.thg_okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00",
                    "#CC79A7", "#999999", "#000000", "#F0E442")

# document -> cluster assignment as a factor in natural order, over hg$nodes
.thg_topic_groups <- function(hg, clusters) {
  .thg_check_hg(hg)
  assignment <- .thg_labels_input(clusters)
  stopifnot(
    "`clusters` must be a data.frame or a named vector" =
      !is.null(names(assignment))
  )
  unknown <- setdiff(names(assignment), hg$nodes)
  if (length(unknown) > 0L) {
    stop(errorCondition(
      paste0("Unknown node names in `clusters`: ",
             paste(unknown, collapse = ", ")),
      class = "honets_bad_input", call = NULL
    ))
  }
  groups <- factor(as.character(assignment)[match(hg$nodes,
                                                  names(assignment))])
  factor(groups, levels = .thg_kw_natural(levels(groups)))
}

#' Topic sizes
#'
#' Counts the documents of each cluster of a clustered hypergraph, as a
#' number and as a share of the clustered documents, and, when a weight per
#' document is given (a repeat count, a citation count), the weighted size
#' and share as well.
#'
#' @param hg The document hypergraph the clustering was computed on.
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels.
#' @param weights `NULL` (default), or a numeric vector named by document
#'   giving each document's weight.
#' @return A base `data.frame` of class `honets_topic_sizes`, one row per
#'   topic in natural order: `topic`, `n`, `share`, and with `weights`
#'   also `weighted_n` and `weighted_share`. `plot()` draws the shares as
#'   horizontal bars, weighted beside unweighted when both exist. Raises
#'   `honets_bad_input` for unknown node names or weights that do not name
#'   every clustered document.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_topic_sizes(hg, topics)
#' hg_topic_sizes(hg, topics, weights = c(cooking_1 = 3, cooking_2 = 1,
#'                                        space_1 = 1, space_2 = 1))
#' @export
hg_topic_sizes <- function(hg, clusters, weights = NULL) {
  groups <- .thg_topic_groups(hg, clusters)
  counts <- table(groups)
  out <- data.frame(
    topic = names(counts), n = as.integer(counts),
    share = as.numeric(counts) / sum(counts),
    stringsAsFactors = FALSE
  )
  if (!is.null(weights)) {
    assigned <- hg$nodes[!is.na(groups)]
    ok <- is.numeric(weights) && !is.null(names(weights)) &&
      all(assigned %in% names(weights)) && !anyNA(weights[assigned])
    if (!ok) {
      stop(errorCondition(
        "`weights` must be a numeric vector named by every clustered document",
        class = "honets_bad_input", call = NULL
      ))
    }
    weighted <- tapply(as.numeric(weights[assigned]), groups[!is.na(groups)],
                       sum)
    out$weighted_n <- as.numeric(weighted[out$topic])
    out$weighted_share <- out$weighted_n / sum(out$weighted_n)
  }
  rownames(out) <- NULL
  class(out) <- c("honets_topic_sizes", "data.frame")
  out
}

#' @rdname hg_topic_sizes
#' @param x A table returned by the verb.
#' @param ... Unused; for S3 consistency.
#' @export
plot.honets_topic_sizes <- function(x, ...) {
  d <- as.data.frame(x)
  long <- data.frame(topic = d$topic, measure = "share", value = d$share,
                     stringsAsFactors = FALSE)
  if ("weighted_share" %in% names(d)) {
    long <- rbind(long, data.frame(topic = d$topic, measure = "weighted share",
                                   value = d$weighted_share,
                                   stringsAsFactors = FALSE))
  }
  long$topic <- factor(long$topic, levels = rev(d$topic))
  long$measure <- factor(long$measure, levels = unique(long$measure))
  long$label <- sprintf("%.1f%%", 100 * long$value)
  ggplot2::ggplot(long, ggplot2::aes(x = .data$value, y = .data$topic,
                                     fill = .data$measure)) +
    ggplot2::geom_col(position = ggplot2::position_dodge(width = 0.8),
                      width = 0.75) +
    ggplot2::geom_text(ggplot2::aes(label = .data$label),
                       position = ggplot2::position_dodge(width = 0.8),
                       hjust = -0.15, size = 3) +
    ggplot2::scale_fill_manual(values = .thg_okabe_ito[c(4L, 1L)],
                               name = NULL) +
    ggplot2::scale_x_continuous(labels = \(v) sprintf("%.0f%%", 100 * v),
                                expand = ggplot2::expansion(mult = c(0, .2))) +
    ggplot2::labs(x = "share of documents", y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = if (nlevels(long$measure) > 1L) "top"
                   else "none")
}

#' Topic coherence and exclusivity
#'
#' Scores each cluster's owned vocabulary on two axes. Coherence is the
#' mean normalised pointwise mutual information (NPMI) between pairs of the
#' cluster's top words, computed on document co-occurrence within the
#' cluster: 1 when two words always occur together, 0 at independence,
#' negative when they avoid each other. Exclusivity is the mean share of
#' the same words, the fraction of each word's corpus count that falls in
#' the cluster: 1 when the cluster owns its vocabulary outright. The top
#' words are the ones [hg_keywords()] ranks by share with the given support
#' floor, so the two scores describe the vocabulary a reader is shown.
#'
#' @param hg The document hypergraph the clustering was computed on
#'   (bag-of-words, documents as nodes).
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels.
#' @param n Top words per cluster to score (default `10`).
#' @param min_docs Support floor passed to [hg_keywords()] (default `1`).
#' @return A base `data.frame` of class `honets_topic_quality`, one row per
#'   topic in natural order: `topic`, `size`, `n_words` (top words scored),
#'   `coherence`, `exclusivity`. `plot()` draws exclusivity against
#'   coherence with one labelled point per topic. Raises `honets_bad_input`
#'   for unknown node names.
#' @references
#' Bouma, G. (2009). Normalized (pointwise) mutual information in
#' collocation extraction. *Proceedings of GSCL*, 31--40.
#'
#' Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
#' for structural topic models. *Journal of Statistical Software*, 91(2).
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_topic_quality(hg, topics, n = 3)
#' @export
hg_topic_quality <- function(hg, clusters, n = 10L, min_docs = 1L) {
  groups <- .thg_topic_groups(hg, clusters)
  words <- hg_keywords(hg, clusters, n = n, type = "frequency",
                       sort_by = "share", min_docs = min_docs)
  words <- as.data.frame(words)
  present <- (hg$incidence != 0) * 1
  rows <- lapply(levels(groups), \(cl) {
    top <- words$word[words$cluster == cl]
    docs <- hg$nodes[!is.na(groups) & groups == cl]
    coherence <- .thg_npmi(present[docs, top, drop = FALSE])
    data.frame(
      topic = cl, size = length(docs), n_words = length(top),
      coherence = coherence,
      exclusivity = if (length(top) > 0L) {
        mean(words$share[words$cluster == cl])
      } else {
        NA_real_
      },
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("honets_topic_quality", "data.frame")
  out
}

# mean pairwise NPMI of the columns of a document x word presence matrix
.thg_npmi <- function(present) {
  present <- as.matrix(present)
  n_docs <- nrow(present)
  n_words <- ncol(present)
  if (n_words < 2L || n_docs == 0L) {
    return(NA_real_)
  }
  p <- colSums(present) / n_docs
  joint <- crossprod(present) / n_docs
  pairs <- which(upper.tri(joint), arr.ind = TRUE)
  p_ij <- joint[pairs]
  p_i <- p[pairs[, 1L]]
  p_j <- p[pairs[, 2L]]
  npmi <- ifelse(p_ij > 0,
                 log(p_ij / (p_i * p_j)) / (-log(p_ij)),
                 -1)
  # a pair present in every document has p_ij = 1 and NPMI 1 by convention
  npmi[p_ij >= 1] <- 1
  mean(npmi)
}

#' @rdname hg_topic_quality
#' @param x A table returned by the verb.
#' @param ... Unused; for S3 consistency.
#' @export
plot.honets_topic_quality <- function(x, ...) {
  d <- as.data.frame(x)
  ggplot2::ggplot(d, ggplot2::aes(x = .data$coherence, y = .data$exclusivity)) +
    ggplot2::geom_point(ggplot2::aes(size = .data$size),
                        colour = .thg_okabe_ito[4L], alpha = 0.8) +
    ggplot2::geom_text(ggplot2::aes(label = .data$topic), vjust = -0.9,
                       size = 3) +
    ggplot2::scale_size_area(max_size = 10, name = "documents") +
    ggplot2::labs(x = "coherence (mean NPMI of top words)",
                  y = "exclusivity (mean share of top words)") +
    ggplot2::theme_minimal(base_size = 12)
}

#' Soft topic membership from the spectral embedding
#'
#' Turns a hard partition into a membership of every document in every
#' topic. The documents are placed in the spectral embedding of the
#' hypergraph Laplacian that [hg_cluster()] cuts (the same `type` and
#' `edge_weights`), each topic's centre is the mean position of its
#' documents, and the membership of a document in a topic is the fuzzy
#' c-means weight with fuzziness 2: the inverse squared distance to that
#' centre, normalised over the topics. A document at a centre has
#' membership 1 there; a document halfway between two centres has 0.5 in
#' each. The values sum to one over the topics, so their level depends on
#' the number of topics `k`: the uniform value is `1 / k`, and with many
#' topics even a clearly assigned document has a modest membership. Read
#' them as ratios between topics (a document with 0.15 and 0.11 in two
#' topics is 1.4 times closer to the first), not as probabilities of
#' belonging.
#'
#' @param hg The document hypergraph the clustering was computed on.
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels.
#' @param type,edge_weights Passed to [hg_cluster()] to reproduce the
#'   embedding the partition was cut in (defaults as there).
#' @return A base `data.frame` of class `honets_membership`, one row per
#'   document and topic: `node`, `cluster` (the hard label), `topic`,
#'   `membership` (rows of one document sum to one). `plot()` draws, per
#'   topic, the distribution of its documents' membership in it: a topic
#'   whose documents sit near 1 is compact, one whose documents spread
#'   towards 0.5 overlaps its neighbours. Raises `honets_bad_input` for
#'   unknown node names.
#' @references
#' Bezdek, J. C. (1981). *Pattern Recognition with Fuzzy Objective Function
#' Algorithms*. Plenum.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_membership(hg, topics)
#' @export
hg_membership <- function(hg, clusters, type = c("zhou", "random_walk"),
                          edge_weights = NULL) {
  type <- match.arg(type)
  groups <- .thg_topic_groups(hg, clusters)
  k <- nlevels(groups)
  embedding <- hg_cluster(hg, k = k, type = type, edge_weights = edge_weights,
                          seed = 1L, what = "embedding")
  dims <- grep("^dim", names(embedding), value = TRUE)
  x <- as.matrix(embedding[, dims, drop = FALSE])
  rownames(x) <- embedding$node
  group_of <- groups[match(embedding$node, hg$nodes)]
  keep <- !is.na(group_of)
  x <- x[keep, , drop = FALSE]
  group_of <- droplevels(group_of[keep])
  centres <- do.call(rbind, lapply(levels(group_of), \(cl) {
    colMeans(x[group_of == cl, , drop = FALSE])
  }))
  # squared distances of every document to every centre, then the fuzzy
  # c-means membership with fuzziness 2: u_ik = 1 / sum_j (d_ik / d_jk)^2
  d2 <- outer(rowSums(x^2), rowSums(centres^2), `+`) - 2 * x %*% t(centres)
  d2 <- pmax(d2, 0)
  at_centre <- d2 < .Machine$double.eps
  inv <- 1 / pmax(d2, .Machine$double.eps)
  membership <- inv / rowSums(inv)
  membership[rowSums(at_centre) > 0, ] <- at_centre[rowSums(at_centre) > 0, ] * 1
  out <- data.frame(
    node = rep(rownames(x), times = ncol(membership)),
    cluster = rep(as.character(group_of), times = ncol(membership)),
    topic = rep(levels(group_of), each = nrow(membership)),
    membership = as.numeric(membership),
    stringsAsFactors = FALSE
  )
  out <- out[order(match(out$node, rownames(x)),
                   match(out$topic, levels(group_of))), , drop = FALSE]
  rownames(out) <- NULL
  class(out) <- c("honets_membership", "data.frame")
  out
}

#' @rdname hg_membership
#' @param x A table returned by the verb.
#' @param ... Unused; for S3 consistency.
#' @export
plot.honets_membership <- function(x, ...) {
  d <- as.data.frame(x)
  own <- d[d$cluster == d$topic, , drop = FALSE]
  own$topic <- factor(own$topic, levels = rev(.thg_kw_natural(unique(own$topic))))
  ggplot2::ggplot(own, ggplot2::aes(x = .data$membership, y = .data$topic)) +
    ggplot2::geom_boxplot(fill = .thg_okabe_ito[2L], colour = "grey30",
                          outlier.size = 0.8, width = 0.6) +
    ggplot2::scale_x_continuous(limits = c(0, 1)) +
    ggplot2::labs(x = "membership of a topic's documents in that topic",
                  y = NULL) +
    ggplot2::theme_minimal(base_size = 12)
}

#' @rdname hg_topic_sizes
#' @export
hypergraph_topic_sizes <- hg_topic_sizes

#' @rdname hg_topic_quality
#' @export
hypergraph_topic_quality <- hg_topic_quality

#' @rdname hg_membership
#' @export
hypergraph_membership <- hg_membership
