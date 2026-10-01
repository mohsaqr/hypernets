# ---- Mixed-membership topic model: KL non-negative matrix factorization ----
#
# hg_topics() factorizes the document-word counts of a hypergraph, X ~ W H,
# with non-negative W (documents x topics) and H (topics x words), by the
# multiplicative updates that minimise the generalised Kullback-Leibler
# divergence (Lee & Seung 1999, 2001). With that divergence the fit is the
# maximum-likelihood solution of probabilistic latent semantic analysis
# (Hofmann 1999; Gaussier & Goutte 2005): after normalisation a row of W is
# a document's distribution over topics and a row of H a topic's
# distribution over words. The updates follow scikit-learn's multiplicative
# solver for beta = 1 (W then H in each step, the model's values floored at
# the float32 machine epsilon where the data are non-zero, and word weights
# below the double epsilon set to zero), so the two agree step for step;
# the local oracle checks that.
#
# The fit depends on its random start. Every start is kept, the best one
# (lowest divergence) is reported, and the others are compared with it by
# the average Jaccard agreement of the topics' ranked top words, topics
# matched one to one by the Hungarian method (Greene, O'Callaghan &
# Cunningham 2014).

# The float32 machine epsilon, the floor scikit-learn applies.
.TM_EPSILON <- 1.1920928955078125e-07

# Document-word counts of a hypergraph as a sparse documents x words matrix.
# A text hypergraph carries its counts in the text layer, whatever its
# orientation or weighting; any other hypergraph uses its incidence, nodes
# as documents and hyperedges as words.
.tm_counts <- function(hg) {
  layer <- hg$text
  if (!is.null(layer) && is.data.frame(layer$weights) &&
      all(c("doc", "word", "count") %in% names(layer$weights))) {
    w <- layer$weights
    docs <- sort(unique(as.character(w$doc)))
    words <- sort(unique(as.character(w$word)))
    return(Matrix::sparseMatrix(
      i = match(w$doc, docs), j = match(w$word, words),
      x = as.numeric(w$count), dims = c(length(docs), length(words)),
      dimnames = list(docs, words)))
  }
  inc <- methods::as(hg$incidence, "CsparseMatrix")
  dimnames(inc) <- list(hg$nodes, colnames(hg$incidence) %||%
                          sprintf("h%d", seq_len(ncol(inc))))
  inc
}

# The model's value W H at the non-zero cells of X, floored as scikit-learn
# does. The factors are stored transposed (topics x documents, topics x
# words) so that gathering the cells reads contiguous columns.
.tm_fitted <- function(Wt, Htt, i, j) {
  pmax(colSums(Wt[, i, drop = FALSE] * Htt[, j, drop = FALSE]), .TM_EPSILON)
}

# Generalised KL divergence D(X || W H) over every cell: the non-zero cells
# contribute x log(x / wh) - x, and every cell contributes its wh.
.tm_divergence <- function(x, wh, Wt, Htt) {
  sum(x * log(x / wh)) - sum(x) + sum(rowSums(Wt) * rowSums(Htt))
}

# One fit from given starting factors (Wt: k x documents, Htt: k x words),
# each step updating W and then H as scikit-learn's multiplicative solver
# does for beta = 1. The divergence is checked every `check_every` steps.
# Multiplicative updates are sequential by construction: each step reads
# the factors the previous step produced, so the loop cannot be vectorised.
.tm_fit <- function(X, Wt, Htt, max_iter, tol, check_every = 10L) {
  X <- methods::as(X, "CsparseMatrix")
  i <- X@i + 1L
  j <- rep.int(seq_len(ncol(X)), diff(X@p))
  x <- X@x
  ratio <- X
  floor_zero <- function(v) {
    v[v == 0] <- .TM_EPSILON
    v
  }
  previous <- .tm_divergence(x, .tm_fitted(Wt, Htt, i, j), Wt, Htt)
  current <- previous
  converged <- FALSE
  iterations <- 0L
  for (step in seq_len(max_iter)) {
    iterations <- step
    ratio@x <- x / .tm_fitted(Wt, Htt, i, j)
    Wt <- Wt * as.matrix(Matrix::tcrossprod(Htt, ratio)) /
      floor_zero(rowSums(Htt))
    ratio@x <- x / .tm_fitted(Wt, Htt, i, j)
    Htt <- Htt * as.matrix(Wt %*% ratio) / floor_zero(rowSums(Wt))
    # word weights that underflow the double epsilon are set to zero, the
    # stabilisation scikit-learn applies to H for this divergence
    Htt[Htt < .Machine$double.eps] <- 0
    if (step %% check_every == 0L || step == max_iter) {
      current <- .tm_divergence(x, .tm_fitted(Wt, Htt, i, j), Wt, Htt)
      if (tol > 0 && abs(previous - current) <= tol * abs(previous)) {
        converged <- TRUE
        break
      }
      previous <- current
    }
  }
  list(W = t(Wt), Ht = t(Htt), divergence = current,
       iterations = iterations, converged = converged)
}

# Starting factors: uniform on (0.5, 1.5) scaled so that W H has the mean of
# X, drawn from the caller's (already seeded) stream.
.tm_start <- function(X, k) {
  scale <- sqrt(sum(X) / (prod(dim(X)) * k))
  list(Wt = matrix(stats::runif(k * nrow(X), 0.5, 1.5) * scale, k, nrow(X)),
       Htt = matrix(stats::runif(k * ncol(X), 0.5, 1.5) * scale, k, ncol(X)))
}

# Average Jaccard of two ranked word lists (Greene et al. 2014): the mean,
# over depths 1..t, of the Jaccard coefficient of the two top-d sets.
.tm_average_jaccard <- function(a, b) {
  depth <- min(length(a), length(b))
  mean(vapply(seq_len(depth), \(d) {
    top_a <- a[seq_len(d)]
    top_b <- b[seq_len(d)]
    length(intersect(top_a, top_b)) / length(union(top_a, top_b))
  }, numeric(1L)))
}

# The Hungarian method (Kuhn 1955; Munkres 1957) for a square cost matrix:
# the permutation p minimising sum(cost[i, p[i]]). Shortest augmenting
# paths with row and column potentials; each augmentation reads the
# potentials the previous one left, so the loops are sequential.
.tm_hungarian <- function(cost) {
  n <- nrow(cost)
  stopifnot("`cost` must be a square numeric matrix" =
              is.matrix(cost) && is.numeric(cost) && ncol(cost) == n)
  u <- numeric(n + 1L)
  v <- numeric(n + 1L)
  match_col <- integer(n + 1L)   # row matched to each column (0 = none)
  way <- integer(n + 1L)
  for (row in seq_len(n)) {
    match_col[1L] <- row
    col0 <- 1L
    minv <- rep(Inf, n + 1L)
    used <- rep(FALSE, n + 1L)
    repeat {
      used[col0] <- TRUE
      row0 <- match_col[col0]
      delta <- Inf
      col1 <- 0L
      for (col in 2:(n + 1L)) {
        if (!used[col]) {
          reduced <- cost[row0, col - 1L] - u[row0 + 1L] - v[col]
          if (reduced < minv[col]) {
            minv[col] <- reduced
            way[col] <- col0
          }
          if (minv[col] < delta) {
            delta <- minv[col]
            col1 <- col
          }
        }
      }
      for (col in seq_len(n + 1L)) {
        if (used[col]) {
          u[match_col[col] + 1L] <- u[match_col[col] + 1L] + delta
          v[col] <- v[col] - delta
        } else {
          minv[col] <- minv[col] - delta
        }
      }
      col0 <- col1
      if (match_col[col0] == 0L) break
    }
    repeat {
      col1 <- way[col0]
      match_col[col0] <- match_col[col1]
      col0 <- col1
      if (col0 == 1L) break
    }
  }
  assignment <- integer(n)
  assignment[match_col[2:(n + 1L)]] <- seq_len(n)
  assignment
}

# Ranked top words of every topic of a fit (columns of Ht).
.tm_top_words <- function(Ht, words, depth) {
  lapply(seq_len(ncol(Ht)), \(t) words[order(-Ht[, t], words)][seq_len(depth)])
}

# Agreement of a fit with the reference: topics matched one to one by the
# Hungarian method on their average Jaccard, then averaged. Returns the
# mean agreement and the matched agreement of every reference topic.
.tm_agreement <- function(reference, other) {
  k <- length(reference)
  similarity <- outer(seq_len(k), seq_len(k), Vectorize(\(a, b) {
    .tm_average_jaccard(reference[[a]], other[[b]])
  }))
  matched <- .tm_hungarian(max(similarity) - similarity)
  per_topic <- similarity[cbind(seq_len(k), matched)]
  list(mean = mean(per_topic), per_topic = per_topic)
}

#' Mixed-membership topic model of a hypergraph
#'
#' Factorizes the document-word counts of a hypergraph into topics, so that
#' every document is a mixture of topics and every topic a distribution over
#' words. The counts \eqn{X} are approximated by \eqn{W H}, with
#' non-negative \eqn{W} (documents by topics) and \eqn{H} (topics by
#' words), by the multiplicative updates that minimise the generalised
#' Kullback-Leibler divergence (Lee & Seung 1999, 2001). With this
#' divergence the fit is the maximum-likelihood solution of probabilistic
#' latent semantic analysis (Hofmann 1999; Gaussier & Goutte 2005): the
#' normalised rows of \eqn{W} are the documents' topic shares
#' \eqn{P(z \mid d)} and the normalised rows of \eqn{H} the topics' word
#' probabilities \eqn{P(w \mid z)}. LDA is the same model with Dirichlet
#' priors on these distributions.
#'
#' The fit depends on its random start, so it is run from `nstart` starts
#' and the start with the lowest divergence is kept. Every other start is
#' compared with it by the agreement of the topics' top words (Greene,
#' O'Callaghan & Cunningham 2014): the topics of the two fits are matched
#' one to one by the Hungarian method (Kuhn 1955) on the average Jaccard
#' coefficient of their ranked top `depth` words, and the matched
#' coefficients are averaged. A topic that the other starts find again has a
#' high `agreement`.
#'
#' A text hypergraph from [text_hypergraph()] is factorized on its word
#' counts, whatever its node orientation and incidence weighting; any other
#' hypergraph on its incidence, with its nodes as documents and its
#' hyperedges as words.
#'
#' @param hg A `net_hg` or `text_hypergraph`.
#' @param k Number of topics, a whole number of at least 2 and smaller than
#'   the numbers of documents and words.
#' @param nstart Number of random starts (default `10`).
#' @param max_iter Maximum multiplicative updates per start (default
#'   `1000`).
#' @param tol Relative change of the divergence between two checks, ten
#'   updates apart, below which a start has converged (default `1e-5`); `0`
#'   runs every start for `max_iter` updates.
#' @param depth Number of top words compared between starts (default `10`).
#' @param seed Seed of the first start; start `r` uses `seed + r - 1`. The
#'   caller's random-number stream is restored on exit. Default `1`.
#' @param parallel Run the starts with `parallel::mclapply()` (not on
#'   Windows). Every start seeds itself, so the result equals the serial
#'   one. Default `FALSE`.
#' @param n_cores Cores when `parallel = TRUE` (default `2`).
#' @return A `net_hg_topics` object. Read it with [hg_get()]: one row per
#'   topic (`what = "topics"`, the default: `topic`, `prevalence` (its mean
#'   share over the documents), `documents` (the summed shares, the
#'   expected number of documents), `agreement` (its mean matched average
#'   Jaccard over the other starts) and `top_words`); one row per document
#'   and topic (`"shares"`: `node`, `topic`, `share`, summing to one over a
#'   document's topics); one row per topic and word (`"words"`: `topic`,
#'   `rank`, `word`, `probability`, `n` words per topic); one row per
#'   document with its largest share (`"documents"`: `node`, `topic`,
#'   `share`); or one row per start (`"restarts"`: `run`, `divergence`,
#'   `iterations`, `converged`, `best`, `agreement`). `print()` shows the
#'   topics, `summary()` every table, and `plot()` the top words of every
#'   topic. A start that reaches `max_iter` without converging raises the
#'   warning `hypernets_no_converge`.
#' @section Conditions:
#' `hypernets_bad_input` for an input that is not a hypergraph, a `k`
#' outside its range, or invalid control arguments.
#' @references
#' Lee, D. D., & Seung, H. S. (1999). Learning the parts of objects by
#' non-negative matrix factorization. \emph{Nature}, 401, 788-791.
#' \doi{10.1038/44565}
#'
#' Lee, D. D., & Seung, H. S. (2001). Algorithms for non-negative matrix
#' factorization. \emph{Advances in Neural Information Processing Systems},
#' 13, 556-562.
#'
#' Hofmann, T. (1999). Probabilistic latent semantic indexing.
#' \emph{Proceedings of SIGIR 1999}, 50-57. \doi{10.1145/312624.312649}
#'
#' Gaussier, E., & Goutte, C. (2005). Relation between PLSA and NMF and
#' implications. \emph{Proceedings of SIGIR 2005}, 601-602.
#' \doi{10.1145/1076034.1076148}
#'
#' Greene, D., O'Callaghan, D., & Cunningham, P. (2014). How many topics?
#' Stability analysis for topic models. \emph{Machine Learning and Knowledge
#' Discovery in Databases (ECML PKDD 2014)}, 498-513.
#' \doi{10.1007/978-3-662-44848-9_32}
#'
#' Kuhn, H. W. (1955). The Hungarian method for the assignment problem.
#' \emph{Naval Research Logistics Quarterly}, 2, 83-97.
#' \doi{10.1002/nav.3800020109}
#' @seealso [hg_topic_quality()] scores the topics' words,
#'   [hg_cluster()] for a partition of the documents.
#' @examples
#' corpus <- c(
#'   a = "soup salt onion soup broth", b = "salt soup broth onion",
#'   c = "stars sky moon night", d = "sky stars night moon moon",
#'   e = "soup stars salt sky night broth")
#' fit <- hg_topics(text_hypergraph(corpus), k = 2, nstart = 3)
#' fit
#' hg_get(fit, what = "shares")
#' hg_get(fit, what = "words", n = 4)
#' @export
hg_topics <- function(hg, k, nstart = 10L, max_iter = 1000L, tol = 1e-5,
                      depth = 10L, seed = 1L, parallel = FALSE,
                      n_cores = 2L) {
  .thg_check_hg(hg)
  whole <- function(v, low) {
    is.numeric(v) && length(v) == 1L && is.finite(v) && v >= low &&
      v == round(v)
  }
  X <- .tm_counts(hg)
  if (!whole(k, 2) || k >= min(dim(X))) {
    .thg_bad_input(sprintf(paste0(
      "`k` must be a whole number from 2 to %d (smaller than the numbers of ",
      "documents and words)"), min(dim(X)) - 1L))
  }
  if (!whole(nstart, 1) || !whole(max_iter, 1) || !whole(depth, 1) ||
      !whole(seed, -.Machine$integer.max)) {
    .thg_bad_input(paste0("`nstart`, `max_iter` and `depth` must be whole ",
                          "numbers of at least 1, and `seed` a whole number"))
  }
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    .thg_bad_input("`tol` must be a single non-negative number")
  }
  k <- as.integer(k)
  depth <- as.integer(min(depth, ncol(X)))
  old_seed <- if (exists(".Random.seed", envir = globalenv())) {
    get(".Random.seed", envir = globalenv())
  }
  on.exit({
    if (is.null(old_seed)) {
      if (exists(".Random.seed", envir = globalenv())) {
        rm(".Random.seed", envir = globalenv())
      }
    } else {
      assign(".Random.seed", old_seed, envir = globalenv())
    }
  }, add = TRUE)
  # every start seeds itself, so the parallel run reproduces the serial one
  one_start <- \(r) {
    set.seed(as.integer(seed) + r - 1L)
    start <- .tm_start(X, k)
    .tm_fit(X, start$Wt, start$Htt, as.integer(max_iter), tol)
  }
  runs <- seq_len(as.integer(nstart))
  fits <- if (isTRUE(parallel) && .Platform$OS.type != "windows") {
    parallel::mclapply(runs, one_start, mc.cores = as.integer(n_cores))
  } else {
    lapply(runs, one_start)
  }
  converged <- vapply(fits, \(f) f$converged, logical(1L))
  divergence <- vapply(fits, \(f) f$divergence, numeric(1L))
  best <- which.min(divergence)
  if (tol > 0 && !converged[best]) {
    warning(warningCondition(sprintf(paste0(
      "the kept start did not converge in %d updates (tol = %g); raise ",
      "`max_iter`"), max_iter, tol),
      class = "hypernets_no_converge", call = NULL))
  }
  words <- colnames(X)
  tops <- lapply(fits, \(f) .tm_top_words(f$Ht, words, depth))
  agreements <- lapply(seq_along(fits), \(r) .tm_agreement(tops[[best]],
                                                          tops[[r]]))
  others <- setdiff(seq_along(fits), best)
  topic_agreement <- if (length(others)) {
    rowMeans(matrix(unlist(lapply(agreements[others], `[[`, "per_topic")),
                    nrow = k))
  } else {
    rep(NA_real_, k)
  }
  .tm_result(fits[[best]], X, k, topic_agreement, fits, divergence, best,
             converged, agreements, others,
             list(k = k, nstart = as.integer(nstart),
                  max_iter = as.integer(max_iter), tol = tol, depth = depth,
                  seed = as.integer(seed)))
}

# Assemble the result: topics ordered by prevalence, shares and word
# probabilities normalised, tables built once.
.tm_result <- function(fit, X, k, topic_agreement, fits, divergence, best,
                       converged, agreements, others, params) {
  W <- fit$W
  Ht <- fit$Ht
  word_totals <- colSums(Ht)
  # P(z | d): W scaled by each topic's word mass, normalised per document
  mass <- W * matrix(word_totals, nrow(W), k, byrow = TRUE)
  doc_totals <- rowSums(mass)
  shares <- mass / ifelse(doc_totals > 0, doc_totals, 1)
  # P(w | z)
  probabilities <- Ht / matrix(ifelse(word_totals > 0, word_totals, 1),
                               nrow(Ht), k, byrow = TRUE)
  order_topics <- order(-colMeans(shares), seq_len(k))
  shares <- shares[, order_topics, drop = FALSE]
  probabilities <- probabilities[, order_topics, drop = FALSE]
  topic_agreement <- topic_agreement[order_topics]
  labels <- sprintf("Topic %d", seq_len(k))
  docs <- rownames(X)
  words <- colnames(X)
  share_table <- data.frame(
    node = rep(docs, times = k),
    topic = rep(labels, each = length(docs)),
    share = as.vector(shares), stringsAsFactors = FALSE)
  share_table <- share_table[order(match(share_table$node, docs),
                                   match(share_table$topic, labels)), ]
  rownames(share_table) <- NULL
  word_table <- do.call(rbind, lapply(seq_len(k), \(t) {
    ranked <- order(-probabilities[, t], words)
    data.frame(topic = labels[t], rank = seq_along(ranked),
               word = words[ranked], probability = probabilities[ranked, t],
               stringsAsFactors = FALSE)
  }))
  rownames(word_table) <- NULL
  top_words <- vapply(seq_len(k), \(t) {
    topic_words <- word_table$word[word_table$topic == labels[t]]
    paste(utils::head(topic_words, 5L), collapse = ", ")
  }, character(1L))
  topic_table <- data.frame(
    topic = labels, prevalence = colMeans(shares),
    documents = colSums(shares), agreement = topic_agreement,
    top_words = top_words, stringsAsFactors = FALSE)
  dominant <- max.col(shares, ties.method = "first")
  document_table <- data.frame(
    node = docs, topic = labels[dominant],
    share = shares[cbind(seq_along(docs), dominant)],
    stringsAsFactors = FALSE)
  restarts <- data.frame(
    run = seq_along(fits), divergence = divergence,
    iterations = vapply(fits, \(f) f$iterations, integer(1L)),
    converged = converged, best = seq_along(fits) == best,
    agreement = vapply(agreements, `[[`, numeric(1L), "mean"))
  structure(list(topics = topic_table, shares = share_table,
                 words = word_table, documents = document_table,
                 restarts = restarts, n_documents = length(docs),
                 n_words = length(words), params = params),
            class = "net_hg_topics")
}

#' @rdname hg_topics
#' @param x A `net_hg_topics` object.
#' @param what Which table: `"topics"` (default), `"shares"`, `"words"`,
#'   `"documents"` or `"restarts"`.
#' @param n For `what = "words"`: the number of words per topic (default
#'   `10`). For `print()`: rows of the topic table shown.
#' @param topic For `what = "shares"` and `"words"`: keep only these topics
#'   (labels such as `"Topic 3"`).
#' @param top Keep only the first `top` rows of the returned table.
#' @param ... Unused.
#' @export
hg_get.net_hg_topics <- function(x, what = c("topics", "shares", "words",
                                             "documents", "restarts"), ...,
                                 n = 10L, topic = NULL, top = NULL) {
  what <- match.arg(what)
  out <- x[[what]]
  if (!is.null(topic) && what %in% c("shares", "words")) {
    if (!all(topic %in% x$topics$topic)) {
      .thg_bad_input(sprintf("`topic` must be among %s",
                             paste(x$topics$topic, collapse = ", ")))
    }
    out <- out[out$topic %in% topic, , drop = FALSE]
  }
  if (identical(what, "words")) {
    if (!is.numeric(n) || length(n) != 1L || n < 1) {
      .thg_bad_input("`n` must be a single number of at least 1")
    }
    out <- out[out$rank <= n, , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}

#' @rdname hg_topics
#' @export
print.net_hg_topics <- function(x, n = 10L, ...) {
  kept <- x$restarts[x$restarts$best, , drop = FALSE]
  cat(sprintf(paste0("Topic model (KL factorization): %d topics, %d ",
                     "documents, %d words; best of %d starts (divergence ",
                     "%.6g, %s)\n"),
              x$params$k, x$n_documents, x$n_words, nrow(x$restarts),
              kept$divergence,
              if (isTRUE(kept$converged)) "converged" else "not converged"))
  .ho_print_table(x, n)
  invisible(x)
}

#' @rdname result-summary
#' @export
summary.net_hg_topics <- function(object, ...) .ho_summary(object)

#' @rdname hg_topics
#' @param y Unused.
#' @return `plot()` returns a ggplot of the `n` most probable words of every
#'   topic, one panel per topic.
#' @export
plot.net_hg_topics <- function(x, y, n = 8L, ...) {
  words <- hg_get(x, what = "words", n = n)
  words$label <- factor(paste(words$topic, words$word, sep = "\r"),
                        levels = rev(paste(words$topic, words$word,
                                           sep = "\r")))
  words$topic <- factor(words$topic, levels = x$topics$topic)
  ggplot2::ggplot(words, ggplot2::aes(x = .data$probability,
                                      y = .data$label)) +
    ggplot2::geom_col(fill = "#0072B2", width = 0.7) +
    ggplot2::facet_wrap(ggplot2::vars(.data$topic), scales = "free",
                        ncol = 4) +
    ggplot2::scale_y_discrete(labels = \(l) sub(".*\r", "", l)) +
    ggplot2::labs(x = "word probability", y = NULL) +
    ggplot2::theme_minimal(base_size = 11)
}

#' Choose the number of topics by coherence and exclusivity
#'
#' Fits [hg_topics()] for every number of topics in `k` and scores each
#' model by the mean semantic coherence (UMass; Mimno et al. 2011) and the
#' mean exclusivity (FREX; Bischof & Airoldi 2012) of its topics, computed
#' by [hg_topic_quality()] on the `n` most probable words. This is the
#' diagnostic of `stm::searchK()` (Roberts, Stewart & Tingley 2019). Adding
#' topics usually makes them more exclusive and less coherent, so no single
#' number is best; Roberts et al. (2014) choose among the models on the
#' frontier of the two, those that no other model beats on both. A model
#' is on the `frontier` when no other number of topics has both a higher
#' mean coherence and a higher mean exclusivity.
#'
#' @inheritParams hg_topics
#' @param k Numbers of topics to compare, at least two whole numbers.
#' @param n Top words per topic scored (default `10`).
#' @param frexw Weight on exclusivity in FREX, in (0, 1) (default `0.7`).
#' @param parallel Fit the numbers of topics with `parallel::mclapply()`
#'   (not on Windows); every fit seeds itself, so the result equals the
#'   serial one. Default `FALSE`.
#' @param n_cores Cores when `parallel = TRUE` (default `2`).
#' @return A base `data.frame` of class `hypernets_topic_search`, one row
#'   per number of topics: `k`, `coherence` and `exclusivity` (means over
#'   the topics), `divergence` (of the kept start), `agreement` (mean
#'   agreement of the topics across the starts; `NA` with `nstart = 1`),
#'   `converged` and `frontier`. `plot()` draws exclusivity against
#'   coherence, one labelled point per number of topics, the frontier
#'   joined by a line and drawn as filled points. Conditions are those of
#'   [hg_topics()]; a start that does not converge raises
#'   `hypernets_no_converge`.
#' @references
#' Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
#' for structural topic models. \emph{Journal of Statistical Software},
#' 91(2), 1-40. \doi{10.18637/jss.v091.i02}
#'
#' Roberts, M. E., Stewart, B. M., Tingley, D., Lucas, C., Leder-Luis, J.,
#' Gadarian, S. K., Albertson, B., & Rand, D. G. (2014). Structural topic
#' models for open-ended survey responses. \emph{American Journal of
#' Political Science}, 58(4), 1064-1082. \doi{10.1111/ajps.12103}
#'
#' Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A.
#' (2011). Optimizing semantic coherence in topic models. \emph{Proceedings
#' of EMNLP 2011}, 262-272.
#'
#' Bischof, J. M., & Airoldi, E. M. (2012). Summarizing topical content with
#' word frequency and exclusivity. \emph{Proceedings of ICML 2012}, 201-208.
#' @seealso [hg_topics()], [hg_topic_quality()].
#' @examples
#' corpus <- c(
#'   a = "soup salt onion soup broth", b = "salt soup broth onion",
#'   c = "stars sky moon night", d = "sky stars night moon moon",
#'   e = "soup stars salt sky night broth", f = "onion salt stars broth")
#' search <- hg_topic_search(text_hypergraph(corpus), k = 2:3, nstart = 2,
#'                           n = 4)
#' search
#' @export
hg_topic_search <- function(hg, k, nstart = 1L, n = 10L, frexw = 0.7,
                            max_iter = 1000L, tol = 1e-5, seed = 1L,
                            parallel = FALSE, n_cores = 2L) {
  .thg_check_hg(hg)
  if (!is.numeric(k) || length(k) < 2L || anyNA(k) || any(k != round(k)) ||
      anyDuplicated(k) > 0L) {
    .thg_bad_input("`k` must hold at least two distinct whole numbers")
  }
  k <- sort(as.integer(k))
  one_k <- \(topics) {
    fit <- hg_topics(hg, k = topics, nstart = nstart, max_iter = max_iter,
                     tol = tol, depth = n, seed = seed)
    quality <- hg_topic_quality(hg, topics = fit, n = n, frexw = frexw)
    kept <- fit$restarts[fit$restarts$best, , drop = FALSE]
    data.frame(k = topics, coherence = mean(quality$coherence),
               exclusivity = mean(quality$exclusivity),
               divergence = kept$divergence,
               agreement = mean(fit$topics$agreement),
               converged = kept$converged)
  }
  rows <- if (isTRUE(parallel) && .Platform$OS.type != "windows") {
    parallel::mclapply(k, one_k, mc.cores = as.integer(n_cores))
  } else {
    lapply(k, one_k)
  }
  # mclapply() returns a worker's error as a try-error; re-raise its
  # condition so the class survives
  failed <- vapply(rows, inherits, logical(1L), "try-error")
  if (any(failed)) {
    stop(attr(rows[[which(failed)[1L]]], "condition"))
  }
  out <- do.call(rbind, rows)
  dominated <- vapply(seq_len(nrow(out)), \(r) {
    any(out$coherence > out$coherence[r] &
          out$exclusivity > out$exclusivity[r])
  }, logical(1L))
  out$frontier <- !dominated
  rownames(out) <- NULL
  class(out) <- c("hypernets_topic_search", "data.frame")
  out
}

#' @rdname hg_topic_search
#' @param x A table returned by the verb.
#' @param ... Unused.
#' @export
plot.hypernets_topic_search <- function(x, ...) {
  d <- .ho_plain(x)
  d$model <- ifelse(d$frontier, "on the frontier", "dominated")
  frontier <- d[d$frontier, , drop = FALSE]
  frontier <- frontier[order(frontier$coherence), , drop = FALSE]
  ggplot2::ggplot(d, ggplot2::aes(x = .data$coherence,
                                  y = .data$exclusivity)) +
    ggplot2::geom_line(data = frontier, colour = "#0072B2") +
    ggplot2::geom_point(ggplot2::aes(shape = .data$model), size = 3,
                        colour = "#0072B2", fill = "#0072B2") +
    ggplot2::scale_shape_manual(values = c("on the frontier" = 21,
                                           dominated = 1), name = NULL) +
    ggplot2::geom_text(ggplot2::aes(label = .data$k), vjust = -1,
                       size = 3.5) +
    ggplot2::scale_y_continuous(
      expand = ggplot2::expansion(mult = c(0.05, 0.12))) +
    ggplot2::labs(x = "mean coherence (UMass)",
                  y = "mean exclusivity (FREX)") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(legend.position = "bottom")
}
