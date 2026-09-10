# Constructor: corpus -> weighted text hypergraph, in four constructions.
.thg_sparse_threshold <- 1e6

# Storage rule for the bag and sentence constructions: sparse once the
# incidence would exceed .thg_sparse_threshold cells. Both counts arrive as
# integers, so the product is taken in double -- a corpus large enough to need
# the sparse path is exactly the one whose integer product overflows to NA,
# and NA would fall through `isTRUE()` and send the largest corpora down the
# dense path.
.thg_choose_sparse <- function(n_docs, n_vocab, construction) {
  stopifnot(
    "`n_docs` must be a single count" =
      length(n_docs) == 1L && is.finite(n_docs) && n_docs >= 0,
    "`n_vocab` must be a single count" =
      length(n_vocab) == 1L && is.finite(n_vocab) && n_vocab >= 0
  )
  construction %in% c("bag", "sentence") &&
    as.double(n_docs) * as.double(n_vocab) >= .thg_sparse_threshold
}
#
# "bag": the bipartite document-word incidence (Hayashi et al. 2020 analyze
# documents as vertices with words as hyperedges; HyperGAT-style uses the
# reverse orientation). "window": sliding/tumbling token windows as
# hyperedges over words (the HyperGAT Sec. 3.2 sequential construction, our
# weighted extension). "knn": documents plus their k nearest embedding
# neighbors as hyperedges (the tier-3 -> tier-2 bridge, via knn_hypergraph()).
# All spectral and structural machinery lives in the hypergraph family; this
# file owns tokenization, weighting, and the tidy text-facing surface.

# Tokenize a character vector into a list of word vectors: lowercase
# (optionally), normalize curly apostrophes (U+2019, written escaped so the
# source stays ASCII) so possessives stay one token, split on non-letter
# runs, trim edge apostrophes, drop empties.
.thg_tokenize <- function(text, lowercase) {
  if (isTRUE(lowercase)) {
    text <- tolower(text)
  }
  text <- gsub("\u2019", "'", text)
  tokens <- strsplit(text, "[^[:alpha:]']+")
  lapply(tokens, \(x) {
    x <- gsub("^'+|'+$", "", x)
    x[nzchar(x)]
  })
}

# Drop tokens shorter than `min_chars`. `min_chars = 1` is a no-op, so the
# default costs nothing.
.thg_drop_short <- function(tokens, min_chars) {
  if (min_chars <= 1L) return(tokens)
  lapply(tokens, \(x) x[nchar(x) >= min_chars])
}

# Window one token sequence: sliding (step 1, full windows; a sequence
# shorter than the window is one whole-sequence window) or tumbling
# (consecutive chunks of `window`, trailing partial chunk included). Each window becomes its sorted set of
# distinct words.
.thg_window_sets <- function(toks, window, mode) {
  n <- length(toks)
  wins <- if (identical(mode, "tumbling")) {
    unname(split(toks, ceiling(seq_along(toks) / window)))
  } else if (n <= window) {
    list(toks)
  } else {
    lapply(seq_len(n - window + 1L), \(i) toks[seq.int(i, length.out = window)])
  }
  lapply(wins, \(win) sort(unique(win)))
}
# Which words survive the vocabulary filters. Ranking is by decreasing corpus
# count with ties broken alphabetically, so the kept set is deterministic and
# every filter is a prefix of the same ranking: min_count trims the rare tail,
# max_words caps the head, coverage keeps the shortest prefix carrying that
# share of all tokens.
.thg_keep_vocabulary <- function(total, min_count = 1L, max_words = Inf,
                                 coverage = 1) {
  stopifnot(
    "internal: `total` must be a named count vector" =
      length(total) > 0L && !is.null(names(total))
  )
  counts <- as.numeric(total)
  ranked <- order(-counts, names(total))
  counts <- counts[ranked]
  words <- names(total)[ranked]

  keep <- counts >= min_count
  if (is.finite(max_words)) {
    keep <- keep & seq_along(counts) <= max_words
  }
  # cumulative token share; the tolerance keeps coverage = 1 from dropping the
  # last word to floating-point error in cumsum()/sum()
  share <- cumsum(counts) / sum(counts)
  reached <- which(share >= coverage - sqrt(.Machine$double.eps))
  if (length(reached) > 0L) {
    keep <- keep & seq_along(counts) <= reached[[1L]]
  }
  sort(words[keep])
}


#' Build a weighted hypergraph from a text corpus
#'
#' Tokenizes a corpus (base R, deterministic) and builds a weighted
#' hypergraph (a [group_hypergraph()] object with a text layer), in one of
#' three constructions:
#'
#' * `construction = "bag"` (default): the document-word incidence. With
#'   `nodes = "doc"`, documents are vertices and each word is a hyperedge
#'   over the documents containing it (the orientation for document
#'   clustering and transductive classification, Hayashi et al. 2020); with
#'   `nodes = "word"` the orientation reverses (the HyperGAT
#'   sentence-as-hyperedge setup generalized to whole documents).
#' * `construction = "window"`: words are vertices and every token window is
#'   a hyperedge -- the sliding-window sequential construction of Ding et al.
#'   (2020, Sec. 3.2), extended with weights: each hyperedge is a distinct
#'   window content (the sorted set of words co-occurring in a window) and
#'   its incidence weight is the number of windows with that content.
#'   Windows never cross document boundaries. `window_mode = "sliding"`
#'   moves one token at a time (a document shorter than `window` forms one
#'   whole-document window); `"tumbling"` uses consecutive chunks, trailing
#'   partial chunk included. This is the text front end of
#'   [window_hypergraph()]: on token sequences with no repeated word inside
#'   a window, the two constructions produce the same incidence matrix
#'   (tested), with the hyperedges here named by their word content. They
#'   differ only at the edges of the definition -- a document shorter than
#'   `window`, a trailing partial chunk, or a word repeated inside one
#'   window, which a set-valued hyperedge collapses.
#' * `construction = "sentence"`: words are vertices and every sentence is
#'   a hyperedge -- HyperGAT's sentence-as-hyperedge construction (Ding et
#'   al. 2020, Sec. 3.2) over the whole corpus. Documents are split at
#'   `.`, `!`, `?` and `;`; each sentence with at least one surviving token
#'   is one hyperedge named `<doc>#<sentence>`, and its incidence weight
#'   is the word's occurrence count in that sentence. Two words are then
#'   related when they share a sentence, not merely a document -- the
#'   scope [hg_keywords()] uses for `type = "sentence_centrality"`.
#' * `construction = "knn"`: documents are vertices and each document plus
#'   its `k` nearest neighbors in an embedding space is one hyperedge,
#'   weighted by cosine similarity (see [knn_hypergraph()]). Pass a
#'   precomputed `embeddings` matrix, or leave it `NULL` to encode the text
#'   with the `sbert` package (if installed; models download only on
#'   explicit user confirmation, per sbert's policy).
#'
#' Tokens are maximal runs of letters (with internal apostrophes); digits
#' and punctuation are separators. `weight = "tfidf"` (bag construction
#' only) uses the smoothed formula `tf * (log((1 + N) / (1 + df)) + 1)` with
#' `tf` the raw count, `N` the number of documents, and `df` the word's
#' document frequency -- the `smooth_idf` variant of Manning, Raghavan and
#' Schutze (2008) as popularized by scikit-learn, which is never zero and so
#' never silently disconnects a vertex. `stop_words` and `min_count`
#' filtering happen before windowing, so removing a token closes the gap it
#' leaves in the sequence.
#'
#' @param x A character vector of documents, or a `data.frame` containing a
#'   text column.
#' @param column Name of the text column when `x` is a `data.frame`.
#' @param id Optional name of an ID column when `x` is a `data.frame`; its
#'   values (unique, non-missing) become the document identifiers. Defaults
#'   to the names of `x` when it is a named character vector, otherwise
#'   `"doc_1"`, `"doc_2"`, ...
#' @param construction `"bag"` (default), `"window"`, `"sentence"` or
#'   `"knn"` -- see
#'   Details.
#' @param nodes Which entity is the vertex set for the bag construction:
#'   `"doc"` (default) or `"word"`. Ignored by `"window"` (vertices are
#'   words) and `"knn"` (vertices are documents).
#' @param weight Term weighting for the bag construction: `"n"` (raw count,
#'   default) or `"tfidf"`. The window construction always weights by
#'   window counts, the sentence construction by in-sentence counts, the
#'   knn construction by cosine similarity.
#' @param stop_words Optional character vector of words to drop after
#'   tokenization (compared after lowercasing when `lowercase = TRUE`); see
#'   [stop_words_en()]. Not applicable to `"knn"`.
#' @param min_count Minimum total corpus count for a word to be kept
#'   (default `1L`, keep everything). Not applicable to `"knn"`.
#' @param max_words Keep at most this many words, the most frequent first
#'   (ties broken alphabetically). `Inf`, the default, keeps every word.
#'   Pruning the long tail is the usual way to make a large corpus tractable:
#'   the rare words carry little signal but dominate the vocabulary, and the
#'   incidence has one hyperedge per word.
#' @param coverage Keep the fewest most-frequent words whose combined corpus
#'   count reaches this share of all tokens, a number in `(0, 1]`. `1`, the
#'   default, keeps every word; `0.99` keeps the words carrying 99% of the
#'   tokens. Applied together with `max_words` and `min_count`, the strictest
#'   wins. Documents left with no tokens are dropped with a
#'   `hypernets_dropped_documents` warning.
#' @param min_chars Minimum number of characters for a word to be kept
#'   (default `1L`, keep everything). Raise it to drop the single letters
#'   and short fragments that initials, enumerations and hyphenated
#'   line breaks leave behind. Not applicable to `"knn"`.
#' @param lowercase Lowercase the text before tokenization (default `TRUE`).
#' @param window Window size in tokens for `construction = "window"`
#'   (default `3L`).
#' @param window_mode `"sliding"` (default) or `"tumbling"`, for
#'   `construction = "window"`.
#' @param k Number of nearest neighbors per hyperedge for
#'   `construction = "knn"` (default `10L`).
#' @param embeddings Numeric matrix of document embeddings for
#'   `construction = "knn"`: one row per document, either rownames matching
#'   the document IDs or rows in document order. `NULL` (default) encodes
#'   the text with `sbert`.
#' @param model Passed to `sbert::encode()` when embeddings are computed
#'   (`NULL` = sbert's default model).
#' @param sparse Store the incidence as a `Matrix::dgCMatrix` (bag and
#'   sentence constructions). Default `NULL` chooses: sparse when the
#'   document-by-word incidence would have a million cells or more, dense
#'   otherwise. Sparse hypergraphs scale to tens
#'   of thousands of documents; [hg_cluster()], [hg_classify()],
#'   [hg_pagerank()], and [hg_measures()] use sparse operator paths that
#'   agree with the dense engines (tested), while tensor centralities and
#'   the null test currently require the dense representation.
#'
#' @return An object of class `c("text_hypergraph", "net_hypergraph")` -- a
#'   [group_hypergraph()] hypergraph accepted by every hypernets
#'   hypergraph verb and by [hg_measures()], [hg_centrality()],
#'   [hg_cluster()], and [hg_classify()] -- with a `text` field recording the
#'   corpus tables. Use [as.data.frame.text_hypergraph()] for the tidy
#'   weight table, and its `what` argument for the document and vocabulary
#'   tables.
#'
#' @section Conditions: Raises `hypernets_bad_input` (broken argument contract,
#'   including bag-only arguments passed to other constructions),
#'   `hypernets_empty_corpus` (no document survives tokenization and filtering),
#'   `hypernets_missing_embeddings` (`construction = "knn"` with neither
#'   `embeddings` nor the sbert package), and warns with
#'   `hypernets_dropped_documents` when some documents end up empty.
#'
#' @references
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
#' Hypergraph attention networks for inductive text classification.
#' *EMNLP 2020*. \doi{10.18653/v1/2020.emnlp-main.399}
#'
#' Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
#' random walks, Laplacians, and clustering. *CIKM 2020*.
#' \doi{10.1145/3340531.3412034}
#'
#' Manning, C. D., Raghavan, P., & Schutze, H. (2008). *Introduction to
#' Information Retrieval*. Cambridge University Press.
#'
#' @examples
#' corpus <- c(
#'   cooking_1 = "Simmer the soup with onions and carrots",
#'   cooking_2 = "This soup recipe needs a pinch of salt",
#'   space_1   = "The telescope revealed a distant galaxy",
#'   space_2   = "Astronomers aimed the telescope at the night sky"
#' )
#' hg <- text_hypergraph(corpus, weight = "tfidf")
#' hg
#' as.data.frame(hg)
#'
#' win <- text_hypergraph(corpus, construction = "window", window = 3)
#' win
#'
#' # knn from precomputed embeddings (sbert-free, offline)
#' emb <- matrix(c(1, 0,  0.9, 0.1,  0, 1,  0.1, 0.9),
#'               nrow = 4, byrow = TRUE,
#'               dimnames = list(names(corpus), NULL))
#' text_hypergraph(corpus, construction = "knn", k = 1, embeddings = emb)
#' @export
text_hypergraph <- function(x, column = NULL, id = NULL,
                            construction = c("bag", "window", "sentence",
                                             "knn"),
                            nodes = c("doc", "word"),
                            weight = c("n", "tfidf"),
                            stop_words = NULL,
                            min_count = 1L,
                            max_words = Inf,
                            coverage = 1,
                            min_chars = 1L,
                            lowercase = TRUE,
                            window = 3L,
                            window_mode = c("sliding", "tumbling"),
                            k = 10L,
                            embeddings = NULL,
                            model = NULL,
                            sparse = NULL) {
  construction <- match.arg(construction)
  nodes <- match.arg(nodes)
  stopifnot("`sparse` must be NULL, TRUE or FALSE" =
              is.null(sparse) || isTRUE(sparse) || isFALSE(sparse))
  if (isTRUE(sparse) && !construction %in% c("bag", "sentence")) {
    stop(errorCondition(
      "`sparse = TRUE` currently supports the bag and sentence constructions only",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  weight <- match.arg(weight)
  window_mode <- match.arg(window_mode)
  stopifnot(
    "`x` must be a character vector or a data.frame" =
      is.character(x) || is.data.frame(x),
    "`stop_words` must be NULL or a character vector" =
      is.null(stop_words) || is.character(stop_words),
    "`min_count` must be a single count >= 1" =
      length(min_count) == 1L && is.finite(min_count) && min_count >= 1,
    "`max_words` must be a single count >= 1, or Inf" =
      length(max_words) == 1L && is.numeric(max_words) && !is.na(max_words) &&
        max_words >= 1,
    "`coverage` must be a single share in (0, 1]" =
      length(coverage) == 1L && is.numeric(coverage) && is.finite(coverage) &&
        coverage > 0 && coverage <= 1,
    "`min_chars` must be a single count >= 1" =
      length(min_chars) == 1L && is.finite(min_chars) && min_chars >= 1,
    "`lowercase` must be TRUE or FALSE" =
      isTRUE(lowercase) || isFALSE(lowercase),
    "`window` must be a single count >= 2" =
      length(window) == 1L && is.finite(window) && window >= 2
  )

  if (is.data.frame(x)) {
    if (is.null(column) || !is.character(column) || length(column) != 1L ||
        !column %in% names(x)) {
      stop(errorCondition(
        "when `x` is a data.frame, `column` must name one of its columns",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    text <- as.character(x[[column]])
    if (is.null(id)) {
      doc_id <- sprintf("doc_%d", seq_len(nrow(x)))
    } else {
      if (!is.character(id) || length(id) != 1L || !id %in% names(x)) {
        stop(errorCondition(
          "`id` must name a column of `x`",
          class = "hypernets_bad_input", call = NULL
        ))
      }
      doc_id <- as.character(x[[id]])
      if (anyNA(doc_id) || anyDuplicated(doc_id) > 0L) {
        stop(errorCondition(
          sprintf("`%s` must hold unique, non-missing document IDs", id),
          class = "hypernets_bad_input", call = NULL
        ))
      }
    }
    meta <- x[setdiff(names(x), c(column, id))]
  } else {
    text <- x
    doc_id <- names(x) %||% sprintf("doc_%d", seq_along(x))
    if (anyNA(doc_id) || anyDuplicated(doc_id) > 0L || !all(nzchar(doc_id))) {
      stop(errorCondition(
        "names of `x` must be unique, non-empty document IDs",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    meta <- NULL
  }
  text[is.na(text)] <- ""

  if (identical(construction, "knn")) {
    return(.thg_knn_text(text, doc_id, meta, k = k, embeddings = embeddings,
                         model = model, stop_words = stop_words,
                         min_count = min_count, min_chars = min_chars,
                         max_words = max_words, coverage = coverage,
                         weight = weight))
  }

  sentence_tokens <- NULL
  if (identical(construction, "sentence")) {
    # split first, tokenize each sentence, drop empty sentences; the
    # document's tokens are the concatenation so the shared count /
    # min_count / vocabulary pipeline below applies unchanged
    sentence_tokens <- lapply(strsplit(text, "[.!?;]+"), \(pieces) {
      toks <- .thg_tokenize(pieces, lowercase = lowercase)
      if (!is.null(stop_words)) {
        toks <- lapply(toks, \(x) x[!x %in% stop_words])
      }
      toks <- .thg_drop_short(toks, min_chars)
      toks[lengths(toks) > 0L]
    })
    tokens <- lapply(sentence_tokens,
                     \(s) unlist(s, use.names = FALSE) %||% character(0))
  } else {
    tokens <- .thg_tokenize(text, lowercase = lowercase)
    if (!is.null(stop_words)) {
      tokens <- lapply(tokens, \(x) x[!x %in% stop_words])
    }
    tokens <- .thg_drop_short(tokens, min_chars)
  }

  words <- unlist(tokens, use.names = FALSE) %||% character(0)
  if (length(words) == 0L) {
    stop(errorCondition(
      "no document contains any token after tokenization and filtering",
      class = "hypernets_empty_corpus", call = NULL
    ))
  }
  long <- data.frame(
    doc = rep(doc_id, lengths(tokens)),
    word = words,
    n = 1L
  )
  counts <- stats::aggregate(n ~ doc + word, data = long, FUN = sum)

  # min_count, max_words and coverage are one filter over one ranking
  total <- tapply(counts$n, counts$word, sum)
  keep <- .thg_keep_vocabulary(total, min_count = min_count,
                               max_words = max_words, coverage = coverage)
  if (length(keep) == 0L) {
    stop(errorCondition(
      sprintf(
        "no word survives `min_count = %d`, `max_words = %s` and `coverage = %g`",
        as.integer(min_count), format(max_words), coverage
      ),
      class = "hypernets_empty_corpus", call = NULL
    ))
  }
  vocabulary_full <- length(total)
  token_share <- sum(as.numeric(total[keep])) / sum(as.numeric(total))
  if (length(keep) < vocabulary_full) {
    counts <- counts[counts$word %in% keep, , drop = FALSE]
    tokens <- lapply(tokens, \(t) t[t %in% keep])
    if (!is.null(sentence_tokens)) {
      sentence_tokens <- lapply(sentence_tokens, \(doc) {
        doc <- lapply(doc, \(t) t[t %in% keep])
        doc[lengths(doc) > 0L]
      })
    }
    message(sprintf(
      "vocabulary pruned to %d of %d words, retaining %.2f%% of tokens",
      length(keep), vocabulary_full, 100 * token_share
    ))
  }

  kept_docs <- doc_id[doc_id %in% counts$doc]
  dropped <- setdiff(doc_id, kept_docs)
  if (length(dropped) > 0L) {
    warning(warningCondition(
      sprintf(
        "%d document(s) had no remaining tokens and were dropped: %s",
        length(dropped), paste(dropped, collapse = ", ")
      ),
      class = "hypernets_dropped_documents"
    ))
  }

  n_docs <- length(kept_docs)
  doc_freq <- tapply(counts$doc, counts$word, \(d) length(unique(d)))
  vocabulary <- data.frame(
    word = names(doc_freq),
    n = as.integer(tapply(counts$n, counts$word, sum)),
    doc_freq = as.integer(doc_freq)
  )
  vocabulary <- vocabulary[order(vocabulary$word), , drop = FALSE]
  rownames(vocabulary) <- NULL

  # storage: sparse when the incidence would exceed a million cells, dense
  # otherwise (the dense engines cover every algorithm; the sparse ones the
  # spectral / transductive core)
  if (is.null(sparse)) {
    sparse <- .thg_choose_sparse(n_docs, nrow(vocabulary), construction)
  }

  n_windows <- NULL
  sentences <- NULL
  if (identical(construction, "sentence")) {
    if (identical(weight, "tfidf")) {
      stop(errorCondition(
        "`weight = \"tfidf\"` applies to the bag construction only; sentence hyperedges are weighted by in-sentence counts",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    names(sentence_tokens) <- doc_id
    sentence_tokens <- sentence_tokens[kept_docs]
    edge_doc <- rep(kept_docs, lengths(sentence_tokens))
    edge_index <- unlist(lapply(lengths(sentence_tokens), seq_len),
                         use.names = FALSE)
    edge_id <- paste0(edge_doc, "#", edge_index)
    flat <- unlist(sentence_tokens, recursive = FALSE, use.names = FALSE)
    sent_long <- data.frame(
      edge = rep(edge_id, lengths(flat)),
      word = unlist(flat, use.names = FALSE),
      n = 1L
    )
    sent_counts <- stats::aggregate(n ~ edge + word, data = sent_long,
                                    FUN = sum)
    sent_counts$w <- as.numeric(sent_counts$n)
    builder <- if (isTRUE(sparse)) .thg_sparse_bipartite else
      group_hypergraph
    hg <- builder(sent_counts, actor = "word", group = "edge", weight = "w")
    weights <- data.frame(edge = sent_counts$edge, word = sent_counts$word,
                          weight = sent_counts$w)
    weights <- weights[order(weights$edge, weights$word), , drop = FALSE]
    rownames(weights) <- NULL
    sentences <- data.frame(edge = edge_id, doc = edge_doc,
                            sentence = edge_index,
                            n_tokens = lengths(flat),
                            stringsAsFactors = FALSE)
    nodes <- "word"
  } else if (identical(construction, "window")) {
    if (identical(weight, "tfidf")) {
      stop(errorCondition(
        "`weight = \"tfidf\"` applies to the bag construction only; windowed hyperedges are weighted by window counts",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    doc_tokens <- tokens[lengths(tokens) > 0L]
    sets <- unlist(lapply(doc_tokens, .thg_window_sets, window = window,
                          mode = window_mode),
                   recursive = FALSE)
    n_windows <- length(sets)
    keys <- vapply(sets, paste, character(1), collapse = "+")
    win_long <- data.frame(
      edge = rep(keys, lengths(sets)),
      word = unlist(sets, use.names = FALSE),
      n = 1L
    )
    win_counts <- stats::aggregate(n ~ edge + word, data = win_long,
                                   FUN = sum)
    win_counts$w <- as.numeric(win_counts$n)
    hg <- group_hypergraph(win_counts, actor = "word",
                                      group = "edge", weight = "w")
    weights <- data.frame(edge = win_counts$edge, word = win_counts$word,
                          weight = win_counts$w)
    weights <- weights[order(weights$edge, weights$word), , drop = FALSE]
    rownames(weights) <- NULL
    nodes <- "word"
  } else {
    if (identical(weight, "tfidf")) {
      idf <- log((1 + n_docs) / (1 + doc_freq)) + 1
      vocabulary$idf <- as.numeric(idf[vocabulary$word])
      counts$w <- counts$n * as.numeric(idf[counts$word])
    } else {
      counts$w <- as.numeric(counts$n)
    }
    stopifnot("internal: non-positive weights produced" = all(counts$w > 0))
    builder <- if (isTRUE(sparse)) .thg_sparse_bipartite else
      group_hypergraph
    hg <- if (identical(nodes, "doc")) {
      builder(counts, actor = "doc", group = "word", weight = "w")
    } else {
      builder(counts, actor = "word", group = "doc", weight = "w")
    }
    weights <- data.frame(doc = counts$doc, word = counts$word,
                          n = counts$n, weight = counts$w)
    weights <- weights[order(weights$doc, weights$word), , drop = FALSE]
    rownames(weights) <- NULL
  }

  n_tokens <- tapply(counts$n, counts$doc, sum)
  n_types <- tapply(counts$word, counts$doc, \(w) length(unique(w)))
  documents <- data.frame(
    doc = kept_docs,
    n_tokens = as.integer(n_tokens[kept_docs]),
    n_types = as.integer(n_types[kept_docs])
  )
  if (!is.null(meta) && ncol(meta) > 0L) {
    meta_kept <- meta[match(kept_docs, doc_id), , drop = FALSE]
    rownames(meta_kept) <- NULL
    documents <- cbind(documents, meta_kept)
  }
  rownames(documents) <- NULL

  hg$text <- list(
    documents = documents,
    vocabulary = vocabulary,
    weights = weights,
    construction = construction,
    nodes = nodes,
    sentences = sentences,
    weighting = switch(construction, window = "window_count",
                       sentence = "sentence_count", weight),
    window = if (identical(construction, "window")) as.integer(window) else NULL,
    window_mode = if (identical(construction, "window")) window_mode else NULL,
    n_windows = n_windows,
    n_dropped = length(dropped),
    min_count = as.integer(min_count),
    max_words = max_words,
    coverage = coverage,
    n_vocabulary_full = vocabulary_full,
    token_share = token_share
  )
  class(hg) <- c("text_hypergraph", class(hg))
  hg
}

# The knn construction: documents as vertices, each document plus its k
# nearest embedding neighbors as one cosine-weighted hyperedge.
.thg_knn_text <- function(text, doc_id, meta, k, embeddings, model,
                          stop_words, min_count, min_chars, max_words,
                          coverage, weight) {
  if (!is.null(stop_words) || min_count > 1L || min_chars > 1L ||
      is.finite(max_words) || coverage < 1 ||
      identical(weight, "tfidf")) {
    stop(errorCondition(
      "`stop_words`, `min_count`, `min_chars`, `max_words`, `coverage`, and `weight` apply to token-based constructions, not construction = \"knn\"",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  if (is.null(embeddings)) {
    if (!requireNamespace("sbert", quietly = TRUE)) {
      stop(errorCondition(
        "construction = \"knn\" needs an `embeddings` matrix, or the sbert package installed to compute one",
        class = "hypernets_missing_embeddings", call = NULL
      ))
    }
    embeddings <- sbert::encode(text, model = model)
    rownames(embeddings) <- doc_id
  } else {
    stopifnot(
      "`embeddings` must be a numeric matrix" =
        is.matrix(embeddings) && is.numeric(embeddings)
    )
    if (is.null(rownames(embeddings))) {
      if (nrow(embeddings) != length(doc_id)) {
        stop(errorCondition(
          sprintf(
            "`embeddings` has %d rows but the corpus has %d documents",
            nrow(embeddings), length(doc_id)
          ),
          class = "hypernets_bad_input", call = NULL
        ))
      }
      rownames(embeddings) <- doc_id
    } else {
      if (!setequal(rownames(embeddings), doc_id)) {
        stop(errorCondition(
          "rownames of `embeddings` must match the document IDs",
          class = "hypernets_bad_input", call = NULL
        ))
      }
      embeddings <- embeddings[doc_id, , drop = FALSE]
    }
  }

  hg <- knn_hypergraph(embeddings, k = k, weight = "cosine")

  documents <- data.frame(doc = doc_id)
  if (!is.null(meta) && ncol(meta) > 0L) {
    meta_reset <- meta
    rownames(meta_reset) <- NULL
    documents <- cbind(documents, meta_reset)
  }

  nz <- which(hg$incidence != 0, arr.ind = TRUE)
  weights <- data.frame(
    doc = rownames(hg$incidence)[nz[, "row"]],
    edge = colnames(hg$incidence)[nz[, "col"]],
    weight = as.numeric(hg$incidence[nz])
  )
  weights <- weights[order(weights$edge, weights$doc), , drop = FALSE]
  rownames(weights) <- NULL

  hg$text <- list(
    documents = documents,
    vocabulary = data.frame(word = character(0), n = integer(0),
                            doc_freq = integer(0)),
    weights = weights,
    construction = "knn",
    nodes = "doc",
    weighting = "cosine",
    k = as.integer(k),
    n_dropped = 0L
  )
  class(hg) <- c("text_hypergraph", class(hg))
  hg
}

#' Print a text hypergraph
#'
#' @param x A [text_hypergraph()] object.
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.text_hypergraph <- function(x, ...) {
  head_line <- switch(x$text$construction,
    bag = sprintf(
      "Text hypergraph: %d documents, %d words (%s as nodes, weight = %s)",
      nrow(x$text$documents), nrow(x$text$vocabulary),
      if (identical(x$text$nodes, "doc")) "documents" else "words",
      x$text$weighting
    ),
    window = sprintf(
      "Text hypergraph: %d documents, %d words (windowed hyperedges: w = %d, %s, %d windows)",
      nrow(x$text$documents), nrow(x$text$vocabulary),
      x$text$window, x$text$window_mode, x$text$n_windows
    ),
    sentence = sprintf(
      "Text hypergraph: %d documents, %d words (sentence hyperedges: %d sentences)",
      nrow(x$text$documents), nrow(x$text$vocabulary),
      nrow(x$text$sentences)
    ),
    knn = sprintf(
      "Text hypergraph: %d documents (kNN embedding hyperedges: k = %d, cosine)",
      nrow(x$text$documents), x$text$k
    )
  )
  cat(head_line, "\n", sep = "")
  size <- as.integer(sub("size_", "", names(x$size_distribution)))
  sizes <- rep(size, x$size_distribution)
  cat(sprintf(
    "Hyperedges: %d (%s); sizes %d-%d, median %g\n",
    x$n_hyperedges,
    switch(x$text$construction,
      bag = if (identical(x$text$nodes, "doc")) "words" else "documents",
      window = "distinct windows",
      sentence = "sentences",
      knn = "kNN neighborhoods"
    ),
    min(sizes), max(sizes), stats::median(sizes)
  ))
  invisible(x)
}

#' Tidy tables of a text hypergraph
#'
#' @param x A [text_hypergraph()] object.
#' @param row.names,optional Ignored; present for S3 consistency.
#' @param what Which table: `"weights"` (default) -- for the bag construction
#'   one row per document-word pair (`doc`, `word`, `n`, `weight`); for the
#'   window construction one row per window-content/word membership
#'   (`edge`, `word`, `weight` = window count); for the sentence
#'   construction one row per sentence/word membership (`edge`, `word`,
#'   `weight` = in-sentence count); for the knn construction one
#'   row per hyperedge membership (`doc`, `edge`, `weight` = cosine
#'   similarity). `"sentences"` (sentence construction only) gives one row
#'   per sentence hyperedge (`edge`, `doc`, `sentence`, `n_tokens`).
#'   `"documents"` gives one row per document (with
#'   `n_tokens`/`n_types` for token-based constructions, plus any metadata
#'   columns carried from the input). `"vocabulary"` gives one row per word
#'   (`word`, `n`, `doc_freq`, and `idf` under tf-idf weighting; empty for
#'   the knn construction, which has no token layer).
#' @param ... Unused.
#' @return A base `data.frame` as described under `what`.
#' @examples
#' hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
#' as.data.frame(hg)
#' as.data.frame(hg, what = "documents")
#' @export
as.data.frame.text_hypergraph <- function(x, row.names = NULL,
                                          optional = FALSE,
                                          what = c("weights", "documents",
                                                   "vocabulary",
                                                   "sentences"),
                                          ...) {
  what <- match.arg(what)
  if (identical(what, "sentences") && is.null(x$text$sentences)) {
    stop(errorCondition(
      "`what = \"sentences\"` needs construction = \"sentence\"",
      class = "hypernets_bad_input", call = NULL
    ))
  }
  x$text[[what]]
}
