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
      class = "hypernets_bad_input", call = NULL
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
#' @return A base `data.frame` of class `hypernets_topic_sizes`, one row per
#'   topic in natural order: `topic`, `n`, `share`, and with `weights`
#'   also `weighted_n` and `weighted_share`. `plot()` draws the shares as
#'   horizontal bars, weighted beside unweighted when both exist. Raises
#'   `hypernets_bad_input` for unknown node names or weights that do not name
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
        class = "hypernets_bad_input", call = NULL
      ))
    }
    weighted <- tapply(as.numeric(weights[assigned]), groups[!is.na(groups)],
                       sum)
    out$weighted_n <- as.numeric(weighted[out$topic])
    out$weighted_share <- out$weighted_n / sum(out$weighted_n)
  }
  rownames(out) <- NULL
  class(out) <- c("hypernets_topic_sizes", "data.frame")
  out
}

#' @rdname hg_topic_sizes
#' @param x A table returned by the verb.
#' @param ... Unused; for S3 consistency.
#' @export
plot.hypernets_topic_sizes <- function(x, ...) {
  d <- .ho_plain(x)
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
#' Scores the top words of each topic on two axes, with the published
#' definitions that topic-model software reports, so hypergraph topics and
#' topics from another model (LDA, STM) can be compared on the same corpus.
#' The corpus is the document hypergraph: its documents are the reference
#' collection and a word's documents are the documents it occurs in.
#'
#' **Topics.** A topic is either a cluster of documents (`clusters`, from
#' [hg_cluster()] or any labelling) or a list of words (`words`, e.g. the
#' top terms of a fitted LDA or STM model), or both. For a cluster, the word
#' distribution is the plug-in estimate
#' \eqn{\beta_{kw} = c_{kw} / \sum_v c_{kv}}, where \eqn{c_{kw}} is the
#' token count of word \eqn{w} in the cluster's documents, over the words
#' that occur in some clustered document; its top words are the `n` most
#' probable (`sort_by = "probability"`, ties in vocabulary order, as
#' `stm` ranks \eqn{\beta}) or the `n` with the largest share of their
#' corpus count in the cluster (`sort_by = "share"`, the ranking of
#' [hg_keywords()] with `sort_by = "share"`), keeping words found in at least
#' `min_docs` of the cluster's documents. Given `words`, those words are
#' scored in their given order.
#'
#' **Coherence** (`coherence =`), with \eqn{D(v)} the number of documents
#' containing \eqn{v}, \eqn{D(v, u)} the number containing both, \eqn{N}
#' the number of documents and \eqn{v_1, \dots, v_M} the top words in rank
#' order:
#' * `"umass"` (default): the UMass coherence of Mimno et al. (2011),
#'   \eqn{\sum_{m=2}^{M} \sum_{l=1}^{m-1} \log\{(D(v_m, v_l) + \epsilon) /
#'   (D(v_l) + \epsilon)\}}, a sum over the \eqn{M(M-1)/2} pairs that
#'   conditions on the higher-ranked word, counted on the whole corpus, with
#'   \eqn{\epsilon = 0.01} in numerator and denominator as in
#'   `stm::semanticCoherence()` (Roberts et al. 2019); Mimno et al. add 1 to
#'   the numerator only. At most 0; closer to 0 is more coherent.
#' * `"npmi"`: the mean normalised pointwise mutual information over the
#'   word pairs (Bouma 2009; the topic-coherence measure of Lau et al.
#'   2014), with document co-occurrence probabilities \eqn{p(v) = D(v)/N}
#'   and \eqn{p(v, u) = D(v, u)/N} on the whole corpus:
#'   \eqn{\log\{p(v,u) / (p(v)p(u))\} / -\log p(v,u)}; -1 for a pair that
#'   never co-occurs (the limit), 1 for a pair present in every document.
#' * `"npmi_cluster"`: the same NPMI counted only on the cluster's own
#'   documents (the definition this verb used before 0.6.0; needs
#'   `clusters`). It rewards words that co-occur inside the topic whatever
#'   they do elsewhere, so it is not comparable with published values.
#'
#' **Exclusivity** (`exclusivity =`), which needs the cluster word
#' distributions and so `clusters`:
#' * `"frex"` (default): the sum over the top words of FREX (Bischof &
#'   Airoldi 2012), the weighted harmonic mean of a word's rank in the
#'   topic's distribution and its rank in exclusivity
#'   \eqn{\beta_{kw} / \sum_j \beta_{jw}}, both as empirical CDFs over the
#'   vocabulary: \eqn{1 / (w / \mathrm{ex} + (1 - w) / \mathrm{fr})} with
#'   weight `frexw` on exclusivity, exactly as `stm::exclusivity()` computes
#'   it. Between 0 and `n`; larger is more exclusive.
#' * `"share"`: the mean over the top words of the share of the word's
#'   corpus count that falls in the cluster, \eqn{c_{kw} / \sum_j c_{jw}}
#'   (the definition used before 0.6.0); 1 when the cluster owns its words.
#' * `"none"`: no exclusivity (`NA`); the choice when only `words` are
#'   given.
#'
#' The measures before 0.6.0 are
#' `coherence = "npmi_cluster", exclusivity = "share", sort_by = "share"`.
#'
#' @param hg The document hypergraph the topics describe: a bag-of-words
#'   [text_hypergraph()] with documents as nodes (the default construction).
#' @param clusters The tidy table returned by [hg_cluster()] (columns
#'   `node`, `cluster`), or a named vector of cluster labels. `NULL` when
#'   only `words` are scored.
#' @param topics A mixed-membership topic model of `hg` fitted by
#'   [hg_topics()]. Its topics are scored on their most probable words, and
#'   FREX is computed from its word distributions \eqn{P(w \mid z)}, as
#'   `stm` computes it from an STM's. `size` is then the topic's expected
#'   number of documents (its summed shares). Give `topics` instead of
#'   `clusters` and `words`.
#' @param words `NULL` (default); a data.frame of topic words with columns
#'   `word` and `topic` (or `cluster`, so [hg_keywords()] output is passed
#'   as it is) and optionally `rank` (the order within a topic; row order
#'   otherwise); or a character matrix with one column per topic and rows in
#'   rank order, as `topicmodels::terms()` returns it (column names are the
#'   topic labels). Every word must be in the vocabulary of `hg`. With
#'   `clusters` too, every topic must be a cluster label.
#' @param n Top words per topic to score (default `10`).
#' @param min_docs Support floor: with `clusters` and no `words`, keep only
#'   words found in at least this many of the cluster's documents (default
#'   `1`).
#' @param coherence `"umass"` (default), `"npmi"` or `"npmi_cluster"`; see
#'   Details.
#' @param exclusivity `"frex"` (default), `"share"` or `"none"`; see Details.
#' @param sort_by How a cluster's top words are chosen when `words` is not
#'   given: `"probability"` (default) or `"share"`; see Details.
#' @param frexw Weight on exclusivity in FREX, in (0, 1) (default `0.7`, as
#'   in `stm`).
#' @return A base `data.frame` of class `hypernets_topic_quality`, one row per
#'   topic (cluster labels in natural order, or topics of `words` in natural
#'   order): `topic`, `size` (the cluster's documents; `NA` without
#'   `clusters`), `n_words` (top words scored), `coherence`, `exclusivity`
#'   (`NA` for fewer than two words, or `exclusivity = "none"`),
#'   `coherence_type` and `exclusivity_type`. `plot()` draws exclusivity
#'   against coherence with one labelled point per topic. Raises
#'   `hypernets_bad_input` for unknown node names, a malformed `words`
#'   table, words outside the vocabulary, topics of `words` that are not
#'   cluster labels, a measure that needs `clusters` without them, or a
#'   hypergraph whose nodes are not documents.
#' @references
#' Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A.
#' (2011). Optimizing semantic coherence in topic models. *Proceedings of
#' the 2011 Conference on Empirical Methods in Natural Language Processing
#' (EMNLP)*, 262--272. <https://aclanthology.org/D11-1024/>
#'
#' Bouma, G. (2009). Normalized (pointwise) mutual information in
#' collocation extraction. *Proceedings of GSCL*, 31--40.
#'
#' Lau, J. H., Newman, D., & Baldwin, T. (2014). Machine reading tea leaves:
#' Automatically evaluating topic coherence and topic model quality.
#' *Proceedings of the 14th Conference of the European Chapter of the
#' Association for Computational Linguistics (EACL)*, 530--539.
#' \doi{10.3115/v1/E14-1056}
#'
#' Bischof, J. M., & Airoldi, E. M. (2012). Summarizing topical content with
#' word frequency and exclusivity. *Proceedings of the 29th International
#' Conference on Machine Learning (ICML)*. arXiv:1206.4631.
#'
#' Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
#' for structural topic models. *Journal of Statistical Software*, 91(2),
#' 1--40. \doi{10.18637/jss.v091.i02}
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' topics <- hg_cluster(hg, k = 2, seed = 1)
#' hg_topic_quality(hg, topics, n = 3)
#' hg_topic_quality(hg, topics, n = 3, coherence = "npmi")
#' # topic words from another model, scored on the same corpus
#' lda_words <- data.frame(topic = c("t1", "t1", "t2", "t2"),
#'                         word = c("soup", "salt", "telescope", "stars"))
#' hg_topic_quality(hg, words = lda_words, exclusivity = "none")
#' @export
hg_topic_quality <- function(hg, clusters = NULL, words = NULL, n = 10L,
                             min_docs = 1L,
                             coherence = c("umass", "npmi", "npmi_cluster"),
                             exclusivity = c("frex", "share", "none"),
                             sort_by = c("probability", "share"),
                             frexw = 0.7, topics = NULL) {
  .thg_check_hg(hg)
  layer <- hg$text
  if (!is.null(layer) && !(identical(layer$construction, "bag") &&
                           identical(layer$nodes, "doc"))) {
    .thg_bad_input(paste0("hg_topic_quality() needs a document hypergraph: ",
                          "text_hypergraph(construction = \"bag\", ",
                          "nodes = \"doc\")"))
  }
  coherence <- match.arg(coherence)
  exclusivity <- match.arg(exclusivity)
  sort_by <- match.arg(sort_by)
  stopifnot(
    "`n` must be a single number >= 1" =
      length(n) == 1L && is.numeric(n) && !is.na(n) && n >= 1,
    "`min_docs` must be a single number >= 1" =
      length(min_docs) == 1L && is.numeric(min_docs) && !is.na(min_docs) &&
        min_docs >= 1
  )
  if (!is.numeric(frexw) || length(frexw) != 1L || is.na(frexw) ||
      frexw <= 0 || frexw >= 1) {
    .thg_bad_input("`frexw` must be a single number strictly between 0 and 1")
  }
  if (!is.null(topics)) {
    return(.thg_topic_model_quality(hg, topics, clusters, words, n, min_docs,
                                    coherence, exclusivity, sort_by, frexw))
  }
  if (is.null(clusters) && is.null(words)) {
    .thg_bad_input("give `clusters`, `words`, or both")
  }
  if (is.null(clusters) && (identical(coherence, "npmi_cluster") ||
                            !identical(exclusivity, "none"))) {
    .thg_bad_input(paste0(
      "`coherence = \"npmi_cluster\"` and every `exclusivity` except ",
      "\"none\" need the cluster word distributions: give `clusters`, or ",
      "use exclusivity = \"none\""))
  }
  present <- .thg_presence(hg$incidence)
  vocabulary <- colnames(present)

  groups <- if (is.null(clusters)) NULL else .thg_topic_groups(hg, clusters)
  counts <- if (is.null(groups)) NULL else .thg_topic_counts(hg, groups)

  top <- if (!is.null(words)) {
    .thg_topic_words_input(words, vocabulary, n,
                           levels = if (is.null(groups)) NULL else levels(groups))
  } else if (identical(sort_by, "share")) {
    kw <- .ho_plain(hg_keywords(hg, clusters, n = n, type = "frequency",
                                    sort_by = "share", min_docs = min_docs))
    lapply(stats::setNames(levels(groups), levels(groups)),
           \(cl) kw$word[kw$cluster == cl])
  } else {
    .thg_topic_top_probability(counts, present, groups, n, min_docs)
  }

  rows <- lapply(names(top), \(topic) {
    top_words <- top[[topic]]
    docs <- if (is.null(groups)) NULL else
      hg$nodes[!is.na(groups) & groups == topic]
    coherence_value <- if (length(top_words) < 2L) {
      NA_real_
    } else {
      switch(coherence,
             umass = .thg_umass(present[, top_words, drop = FALSE]),
             npmi = .thg_npmi(present[, top_words, drop = FALSE]),
             npmi_cluster = .thg_npmi(present[docs, top_words, drop = FALSE]))
    }
    exclusivity_value <- if (length(top_words) == 0L) {
      NA_real_
    } else {
      switch(exclusivity,
             frex = .thg_frex(counts, topic, top_words, frexw),
             share = .thg_word_share(counts, topic, top_words),
             none = NA_real_)
    }
    data.frame(
      topic = topic,
      size = if (is.null(groups)) NA_integer_ else length(docs),
      n_words = length(top_words),
      coherence = coherence_value, exclusivity = exclusivity_value,
      coherence_type = coherence, exclusivity_type = exclusivity,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("hypernets_topic_quality", "data.frame")
  out
}

# Quality of a fitted hg_topics() model: its top words by probability, FREX
# from its topic x word distributions (the matrix stm passes as beta).
.thg_topic_model_quality <- function(hg, topics, clusters, words, n, min_docs,
                                     coherence, exclusivity, sort_by, frexw) {
  if (!inherits(topics, "net_hg_topics")) {
    .thg_bad_input("`topics` must be a topic model fitted by hg_topics()")
  }
  if (!is.null(clusters) || !is.null(words)) {
    .thg_bad_input("give `topics` alone, without `clusters` or `words`")
  }
  if (identical(coherence, "npmi_cluster") ||
      !identical(sort_by, "probability") || !isTRUE(min_docs == 1)) {
    .thg_bad_input(paste0(
      "a topic model has no documents of its own and ranks its words by ",
      "probability: `coherence = \"npmi_cluster\"`, `sort_by = \"share\"` ",
      "and `min_docs` apply to `clusters` only"))
  }
  present <- .thg_presence(hg$incidence)
  table <- topics$words
  unknown <- setdiff(unique(table$word), colnames(present))
  if (length(unknown) > 0L || !identical(nrow(present), topics$n_documents)) {
    .thg_bad_input("`topics` was not fitted on `hg`")
  }
  labels <- topics$topics$topic
  vocabulary <- unique(table$word)
  probability <- matrix(0, length(labels), length(vocabulary),
                        dimnames = list(labels, vocabulary))
  probability[cbind(match(table$topic, labels),
                    match(table$word, vocabulary))] <- table$probability
  probability <- probability[, colSums(probability) > 0, drop = FALSE]
  top <- lapply(stats::setNames(labels, labels), \(t) {
    ranked <- table[table$topic == t & table$probability > 0, , drop = FALSE]
    utils::head(ranked$word[order(ranked$rank)], n)
  })
  rows <- lapply(labels, \(topic) {
    top_words <- top[[topic]]
    coherence_value <- if (length(top_words) < 2L) {
      NA_real_
    } else {
      switch(coherence,
             umass = .thg_umass(present[, top_words, drop = FALSE]),
             npmi = .thg_npmi(present[, top_words, drop = FALSE]))
    }
    exclusivity_value <- switch(
      exclusivity,
      frex = .thg_frex(probability, topic, top_words, frexw),
      share = .thg_word_share(probability, topic, top_words),
      none = NA_real_)
    data.frame(
      topic = topic,
      size = topics$topics$documents[topics$topics$topic == topic],
      n_words = length(top_words),
      coherence = coherence_value, exclusivity = exclusivity_value,
      coherence_type = coherence, exclusivity_type = exclusivity,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  rownames(out) <- NULL
  class(out) <- c("hypernets_topic_quality", "data.frame")
  out
}

# document x word presence (0/1) of an incidence matrix, dense
.thg_presence <- function(incidence) {
  present <- as.matrix(incidence != 0) * 1
  dimnames(present) <- dimnames(incidence)
  present
}

# cluster x word token counts over the words that occur in some clustered
# document (a word no cluster uses has no exclusivity: 0 / 0)
.thg_topic_counts <- function(hg, groups) {
  scope <- .thg_kw_scope(hg)
  counts <- .thg_kw_aggregate(groups, .thg_kw_counts(hg, scope))
  counts[, colSums(counts) > 0, drop = FALSE]
}

# top words by the cluster's word distribution: count descending, ties in
# vocabulary order (what order(decreasing = TRUE) does on stm's beta),
# keeping words found in at least `min_docs` of the cluster's documents
.thg_topic_top_probability <- function(counts, present, groups, n, min_docs) {
  support <- .thg_kw_aggregate(groups, present[, colnames(counts),
                                               drop = FALSE])
  lapply(stats::setNames(rownames(counts), rownames(counts)), \(cl) {
    row <- counts[cl, ]
    ord <- order(-row, seq_along(row))
    eligible <- row[ord] > 0 & support[cl, ord] >= min_docs
    names(row)[utils::head(ord[eligible], n)]
  })
}

# a (topic, word[, rank]) table -> named list of word vectors in rank order
.thg_topic_words_input <- function(words, vocabulary, n, levels = NULL) {
  # topicmodels::terms(): one column per topic, rows in rank order
  if (is.matrix(words) && is.character(words)) {
    labels <- colnames(words) %||% paste("Topic", seq_len(ncol(words)))
    words <- data.frame(topic = rep(labels, each = nrow(words)),
                        word = as.vector(words),
                        rank = rep(seq_len(nrow(words)), times = ncol(words)),
                        stringsAsFactors = FALSE)
  }
  # hg_keywords() output labels its topics `cluster`
  if (is.data.frame(words) && !"topic" %in% names(words) &&
      "cluster" %in% names(words)) {
    words$topic <- words$cluster
  }
  if (!is.data.frame(words) || !all(c("topic", "word") %in% names(words)) ||
      nrow(words) == 0L) {
    .thg_bad_input(paste0("`words` must be a non-empty data.frame with ",
                          "columns `word` and `topic` (or `cluster`), or a ",
                          "character matrix with one column per topic"))
  }
  topic <- as.character(words$topic)
  word <- as.character(words$word)
  rank <- if ("rank" %in% names(words)) words$rank else seq_along(word)
  if (anyNA(topic) || anyNA(word) || !is.numeric(rank) || anyNA(rank)) {
    .thg_bad_input("`words` has missing topics, words or ranks")
  }
  if (anyDuplicated(paste(topic, word, sep = "\r")) > 0L) {
    .thg_bad_input("`words` lists a word twice within one topic")
  }
  unknown <- setdiff(word, vocabulary)
  if (length(unknown) > 0L) {
    .thg_bad_input(sprintf("words not in the vocabulary of `hg`: %s",
                           paste(utils::head(unknown, 10L), collapse = ", ")))
  }
  topics <- .thg_kw_natural(unique(topic))
  if (!is.null(levels)) {
    stray <- setdiff(topics, levels)
    if (length(stray) > 0L) {
      .thg_bad_input(sprintf("topics of `words` that are not cluster labels: %s",
                             paste(stray, collapse = ", ")))
    }
    topics <- levels
  }
  lapply(stats::setNames(topics, topics), \(t) {
    in_topic <- topic == t
    utils::head(word[in_topic][order(rank[in_topic])], n)
  })
}

# UMass coherence (Mimno et al. 2011) in stm's form: the columns of
# `present` are the top words in rank order; sum over pairs m > l of
# log((D(m, l) + 0.01) / (D(l) + 0.01)) on the whole corpus.
.thg_umass <- function(present) {
  co <- crossprod(present)
  pairs <- which(lower.tri(co), arr.ind = TRUE)  # row m > column l
  sum(log(co[pairs] + 0.01) - log(diag(co)[pairs[, 2L]] + 0.01))
}

# FREX (Bischof & Airoldi 2012) summed over `top_words`, exactly as
# stm::exclusivity(): empirical-CDF ranks of exclusivity and frequency over
# the vocabulary, combined by a weighted harmonic mean.
.thg_frex <- function(counts, topic, top_words, frexw) {
  beta <- counts / rowSums(counts)
  share <- sweep(beta, 2L, colSums(beta), `/`)
  n_vocab <- ncol(beta)
  ex <- rank(share[topic, ]) / n_vocab
  fr <- rank(beta[topic, ]) / n_vocab
  frex <- 1 / (frexw / ex + (1 - frexw) / fr)
  sum(frex[top_words])
}

# mean share of each top word's corpus count that falls in the cluster
.thg_word_share <- function(counts, topic, top_words) {
  total <- colSums(counts)
  mean(counts[topic, top_words] / total[top_words])
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
plot.hypernets_topic_quality <- function(x, ...) {
  d <- .ho_plain(x)
  coherence_label <- switch(d$coherence_type[[1L]] %||% "npmi_cluster",
                            umass = "coherence (UMass, sum over word pairs)",
                            npmi = "coherence (mean NPMI of top words)",
                            "coherence (mean in-topic NPMI of top words)")
  exclusivity_label <- switch(d$exclusivity_type[[1L]] %||% "share",
                              frex = "exclusivity (FREX, sum over top words)",
                              none = "exclusivity (not computed)",
                              "exclusivity (mean share of top words)")
  sized <- !all(is.na(d$size))
  d$exclusivity_shown <- ifelse(is.na(d$exclusivity), 0, d$exclusivity)
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$coherence,
                                       y = .data$exclusivity_shown))
  p <- if (sized) {
    p + ggplot2::geom_point(ggplot2::aes(size = .data$size),
                            colour = .thg_okabe_ito[4L], alpha = 0.8) +
      ggplot2::scale_size_area(max_size = 10, name = "documents")
  } else {
    p + ggplot2::geom_point(colour = .thg_okabe_ito[4L], size = 3)
  }
  p +
    ggplot2::geom_text(ggplot2::aes(label = .data$topic), vjust = -0.9,
                       size = 3) +
    ggplot2::labs(x = coherence_label, y = exclusivity_label) +
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
#' @return A base `data.frame` of class `hypernets_membership`, one row per
#'   document and topic: `node`, `cluster` (the hard label), `topic`,
#'   `membership` (rows of one document sum to one). `plot()` draws, per
#'   topic, the distribution of its documents' membership in it: a topic
#'   whose documents sit near 1 is compact, one whose documents spread
#'   towards 0.5 overlaps its neighbours. Raises `hypernets_bad_input` for
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
  class(out) <- c("hypernets_membership", "data.frame")
  out
}

#' @rdname hg_membership
#' @param x A table returned by the verb.
#' @param ... Unused; for S3 consistency.
#' @export
plot.hypernets_membership <- function(x, ...) {
  d <- .ho_plain(x)
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

#' Co-clustering of nodes and hyperedges (documents and words)
#'
#' Clusters the nodes and the hyperedges of a hypergraph together, so that
#' each cluster is a pair: a set of documents and the set of words that
#' characterises them, from one spectral cut of the document-word bipartite
#' graph (Dhillon 2001). With \eqn{A} the incidence matrix (nodes by
#' hyperedges, the stored weights: counts or tf-idf on a
#' [text_hypergraph()]), \eqn{D_1} and \eqn{D_2} its row and column sums,
#' the scaled matrix \eqn{A_n = D_1^{-1/2} A D_2^{-1/2}} is decomposed by
#' SVD; its leading singular pair is trivial, the next
#' \eqn{\ell = \lceil \log_2 k \rceil} left and right singular vectors
#' \eqn{U, V} give the embedding
#' \eqn{Z = [D_1^{-1/2} U;\ D_2^{-1/2} V]}, one row per node and per
#' hyperedge, and k-means on the rows of \eqn{Z} gives the `k` co-clusters
#' (Dhillon 2001, Algorithm Multipartition).
#'
#' For a binary incidence with unit hyperedge weights, the node half of the
#' embedding spans the same space as the Zhou et al. (2006) hypergraph
#' Laplacian eigenvectors that [hg_cluster()] cuts (the squared singular
#' values are one minus the Laplacian eigenvalues), so co-clustering reads
#' the same cut on both sides of the bipartite graph; with other weights the
#' two differ, because the Zhou Laplacian reads only the binary pattern.
#' A cluster may hold nodes and no hyperedges, or the reverse.
#'
#' @param hg A [text_hypergraph()] (or any hypernets `net_hg`) with no empty
#'   node or hyperedge.
#' @param k Number of co-clusters, at least 2.
#' @param seed Seed for the k-means starts; the caller's random number stream
#'   is restored on exit. `NULL` (default) uses the current stream.
#' @param nstart Number of k-means starts (default `25L`).
#' @param what `"clusters"` (default) for the assignment, or `"embedding"` for
#'   the assignment plus the embedding coordinates `dim1..dimL`.
#' @param role Which rows to return: `"all"` (default), `"node"` (on a
#'   document hypergraph, the documents) or `"hyperedge"` (the words).
#' @return A base `data.frame`, one row per node and then per hyperedge
#'   (filtered by `role`), with columns `node` (the node or hyperedge name),
#'   `role` (`"node"` or `"hyperedge"`) and `cluster` (`"Cluster 1"`, ...,
#'   numbered by first appearance in that row order); with
#'   `what = "embedding"` also `dim1..dimL`. Raises `hypernets_bad_input` for
#'   a `k` below 2 or above the number of rows, or an empty node or
#'   hyperedge (zero degree, where the scaling is undefined).
#' @references
#' Dhillon, I. S. (2001). Co-clustering documents and words using bipartite
#' spectral graph partitioning. *Proceedings of the Seventh ACM SIGKDD
#' International Conference on Knowledge Discovery and Data Mining*,
#' 269--274. \doi{10.1145/502512.502550}
#'
#' Zhou, D., Huang, J., & Schölkopf, B. (2006). Learning with hypergraphs:
#' Clustering, classification, and embedding. *Advances in Neural Information
#' Processing Systems 19*, 1601--1608.
#' @examples
#' hg <- text_hypergraph(c(
#'   cooking_1 = "simmer the soup with onions and carrots",
#'   cooking_2 = "this soup recipe needs salt on a cold night",
#'   space_1 = "the telescope revealed a distant galaxy and stars",
#'   space_2 = "astronomers aimed the telescope at the stars all night"
#' ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
#' hg_cocluster(hg, k = 2, seed = 1)
#' hg_cocluster(hg, k = 2, seed = 1, role = "hyperedge")
#' @export
hg_cocluster <- function(hg, k, seed = NULL, nstart = 25L,
                         what = c("clusters", "embedding"),
                         role = c("all", "node", "hyperedge")) {
  .thg_check_hg(hg)
  what <- match.arg(what)
  role <- match.arg(role)
  incidence <- hg$incidence
  n_rows <- nrow(incidence) + ncol(incidence)
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k < 2 || k > n_rows) {
    .thg_bad_input(sprintf(
      "`k` must be a single number between 2 and %d (nodes plus hyperedges)",
      n_rows))
  }
  stopifnot("`nstart` must be a single number >= 1" =
              length(nstart) == 1L && is.numeric(nstart) && nstart >= 1)
  row_sum <- as.numeric(Matrix::rowSums(incidence))
  col_sum <- as.numeric(Matrix::colSums(incidence))
  if (any(row_sum <= 0) || any(col_sum <= 0)) {
    .thg_bad_input(paste0("co-clustering needs every node and hyperedge to ",
                          "carry positive weight; drop empty ones first"))
  }
  n_vectors <- as.integer(ceiling(log2(k)))
  scaled <- Matrix::Diagonal(x = 1 / sqrt(row_sum)) %*% incidence %*%
    Matrix::Diagonal(x = 1 / sqrt(col_sum))
  decomposition <- if (methods::is(incidence, "sparseMatrix") &&
                       n_vectors + 1L < min(dim(incidence))) {
    RSpectra::svds(scaled, k = n_vectors + 1L)
  } else {
    svd(as.matrix(scaled), nu = n_vectors + 1L, nv = n_vectors + 1L)
  }
  keep <- seq_len(n_vectors) + 1L   # drop the trivial leading pair
  u <- decomposition$u[, keep, drop = FALSE]
  v <- decomposition$v[, keep, drop = FALSE]
  # fix each singular pair's sign (largest |u| entry positive) so the
  # embedding is reproducible across SVD backends
  flip <- apply(u, 2L, \(col) sign(col[which.max(abs(col))]))
  u <- sweep(u, 2L, flip, `*`)
  v <- sweep(v, 2L, flip, `*`)
  z <- rbind(u / sqrt(row_sum), v / sqrt(col_sum))
  colnames(z) <- paste0("dim", seq_len(n_vectors))
  names_out <- c(hg$nodes,
                 colnames(incidence) %||% paste0("e", seq_len(ncol(incidence))))

  if (!is.null(seed)) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    saved_seed <- if (had_seed) get(".Random.seed", envir = globalenv())
    on.exit(.thg_rng_restore(had_seed, saved_seed), add = TRUE)
    set.seed(seed)
  }
  fit <- stats::kmeans(z, centers = k, nstart = nstart, iter.max = 100L)
  # deterministic labels: "Cluster 1" is the first row's cluster, and so on
  relabel <- match(fit$cluster, unique(fit$cluster))
  out <- data.frame(
    node = names_out,
    role = rep(c("node", "hyperedge"), c(nrow(incidence), ncol(incidence))),
    cluster = paste("Cluster", relabel),
    stringsAsFactors = FALSE
  )
  if (identical(what, "embedding")) {
    out <- cbind(out, as.data.frame(z, row.names = NULL))
  }
  if (!identical(role, "all")) {
    out <- out[out$role == role, , drop = FALSE]
  }
  rownames(out) <- NULL
  out
}
