# v0.5 benchmark harness -------------------------------------------------
#
# R8 / R52 / MR / Ohsumed / 20NG on the TextGCN files (Yao 2019), the same
# cleaned corpora and train/test splits HyperGAT's preprocessing reads --
# any other cleaning would make accuracy incomparable with the published
# tables. Data files are fetched by benchmarks/download.sh into
# benchmarks/data/ (gitignored).
#
# Methods benchmarked:
#   - bench_transduction(): transductive label spreading on the document-
#     word hypergraph (Zhou 2006 "zhou" / Hayashi 2020 EDVW "random_walk"),
#     via text_hypergraph(sparse = TRUE) + hg_classify().
#   - bench_centroid(): tf-idf nearest-centroid (Rocchio) baseline on the
#     same incidence weights.
#
# Every bench_* verb returns one tidy row: accuracy with a percentile
# bootstrap CI, macro-F1, sizes, and wall-clock timings.

# Read one dataset: doc id, text, split ("train"/"test"), label.
bench_load <- function(name, dir = file.path("benchmarks", "data")) {
  stopifnot(
    "`name` must be a single dataset name" =
      is.character(name) && length(name) == 1L
  )
  split_file <- file.path(dir, paste0(name, ".txt"))
  corpus_file <- file.path(dir, "corpus", paste0(name, ".clean.txt"))
  stopifnot(
    "split file missing -- run benchmarks/download.sh" =
      file.exists(split_file),
    "corpus file missing -- run benchmarks/download.sh" =
      file.exists(corpus_file)
  )
  meta <- utils::read.delim(split_file, header = FALSE,
                            col.names = c("row", "split", "label"),
                            colClasses = "character", quote = "")
  docs <- readLines(corpus_file, encoding = "UTF-8", warn = FALSE)
  stopifnot(
    "corpus and split files must have the same number of documents" =
      length(docs) == nrow(meta)
  )
  data.frame(
    id = sprintf("%s_%05d", name, seq_along(docs)),
    text = docs,
    split = ifelse(grepl("train", meta$split, fixed = TRUE),
                   "train", "test"),
    label = meta$label,
    stringsAsFactors = FALSE
  )
}

# Accuracy + percentile bootstrap CI + macro-F1 for one prediction vector.
# NA predictions (documents the model could not score) count as errors.
.bench_score <- function(truth, predicted, n_boot = 1000L, seed = 1L) {
  stopifnot(
    "`truth` and `predicted` must have equal length" =
      length(truth) == length(predicted)
  )
  correct <- !is.na(predicted) & predicted == truth
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
            add = TRUE, after = FALSE)
  }
  set.seed(seed)
  boot <- vapply(
    seq_len(n_boot),
    \(i) mean(sample(correct, replace = TRUE)),
    numeric(1)
  )
  ci <- stats::quantile(boot, c(0.025, 0.975), names = FALSE)
  classes <- sort(unique(truth))
  f1 <- vapply(classes, \(cl) {
    tp <- sum(correct & truth == cl)
    fp <- sum(!is.na(predicted) & predicted == cl & truth != cl)
    fn <- sum(truth == cl) - tp
    if (tp == 0L) return(0)  # precision and recall both zero
    2 * tp / (2 * tp + fp + fn)
  }, numeric(1))
  data.frame(
    accuracy = mean(correct),
    ci_low = ci[1L],
    ci_high = ci[2L],
    macro_f1 = mean(f1),
    n_unscored = sum(is.na(predicted))
  )
}

# Largest connected component of a sparse document-word hypergraph, by
# alternating frontier propagation (docs -> words -> docs), never
# materializing the doc-doc adjacency. Returns the node names inside it.
# Label spreading cannot reach documents outside the seeded component, so
# the transduction benchmark runs on the giant component and scores the
# excluded test documents as errors (the honest treatment).
.bench_giant_component <- function(hg) {
  membership <- (hg$incidence > 0) * 1
  n <- nrow(membership)
  unseen <- rep(TRUE, n)
  best <- logical(n)
  # while loops: sequential component discovery -- each pass consumes one
  # component and each BFS wave depends on the previous frontier
  while (any(unseen) && sum(best) < sum(unseen)) {
    x <- numeric(n)
    x[which(unseen)[1L]] <- 1
    repeat {
      words <- as.numeric(Matrix::crossprod(membership, x) > 0)
      x_new <- as.numeric((membership %*% words) > 0)
      if (sum(x_new) == sum(x)) break
      x <- x_new
    }
    comp <- x > 0
    if (sum(comp) > sum(best)) best <- comp
    unseen <- unseen & !comp
  }
  hg$nodes[best]
}

# Transductive hypergraph classification on one dataset.
bench_transduction <- function(name, xi = 0.99,
                               type = c("zhou", "random_walk"),
                               normalization = c("class_mass", "none"),
                               weight = c("tfidf", "n"), min_count = 1L,
                               dir = file.path("benchmarks", "data")) {
  type <- match.arg(type)
  normalization <- match.arg(normalization)
  weight <- match.arg(weight)
  corpus <- bench_load(name, dir)
  docs <- corpus$text
  names(docs) <- corpus$id
  t_build <- system.time({
    hg <- hypernets::text_hypergraph(
      docs, weight = weight, stop_words = character(0),
      min_count = min_count, sparse = TRUE
    )
    keep <- .bench_giant_component(hg)
    if (length(keep) < hg$n_nodes) {
      # rebuild on the giant component only (see .bench_giant_component)
      hg <- hypernets::text_hypergraph(
        docs[names(docs) %in% keep], weight = weight,
        stop_words = character(0), min_count = min_count, sparse = TRUE
      )
    }
  })[["elapsed"]]
  # dropped or unreachable train documents cannot seed; dropped or
  # unreachable test documents are scored as errors via the NA rule
  train <- corpus[corpus$split == "train" & corpus$id %in% hg$nodes, ]
  seeds <- train$label
  names(seeds) <- train$id
  t_fit <- system.time(
    fit <- hypernets::hg_classify(hg, labels = seeds, xi = xi,
                                       type = type,
                                       normalization = normalization)
  )[["elapsed"]]
  test <- corpus[corpus$split == "test", ]
  predicted <- fit$predicted[match(test$id, fit$node)]
  cbind(
    data.frame(
      dataset = name,
      method = paste0("hypergraph_", type,
                      ifelse(identical(normalization, "class_mass"),
                             "_cmn", "")),
      weight = weight, xi = xi,
      n_train = sum(corpus$split == "train"), n_test = nrow(test),
      n_nodes = hg$n_nodes, n_edges = hg$n_hyperedges
    ),
    .bench_score(test$label, predicted),
    data.frame(build_s = t_build, fit_s = t_fit)
  )
}

# tf-idf nearest-centroid (Rocchio) baseline on the same incidence.
bench_centroid <- function(name, weight = c("tfidf", "n"), min_count = 1L,
                           dir = file.path("benchmarks", "data")) {
  weight <- match.arg(weight)
  corpus <- bench_load(name, dir)
  docs <- corpus$text
  names(docs) <- corpus$id
  t_build <- system.time(
    hg <- hypernets::text_hypergraph(
      docs, weight = weight, stop_words = character(0),
      min_count = min_count, sparse = TRUE
    )
  )[["elapsed"]]
  t_fit <- system.time({
    x <- hg$incidence  # documents x words, chosen weighting
    norms <- sqrt(Matrix::rowSums(x^2))
    norms[norms == 0] <- 1
    x <- Matrix::Diagonal(x = 1 / norms) %*% x
    train <- corpus[corpus$split == "train" & corpus$id %in% hg$nodes, ]
    rows <- match(train$id, hg$nodes)
    classes <- sort(unique(train$label))
    indicator <- Matrix::sparseMatrix(
      i = match(train$label, classes), j = seq_along(rows), x = 1,
      dims = c(length(classes), length(rows))
    )
    centroids <- indicator %*% x[rows, , drop = FALSE]
    c_norms <- sqrt(Matrix::rowSums(centroids^2))
    c_norms[c_norms == 0] <- 1
    centroids <- Matrix::Diagonal(x = 1 / c_norms) %*% centroids
    cosine <- as.matrix(Matrix::tcrossprod(x, centroids))
    predicted_all <- classes[max.col(cosine, ties.method = "first")]
  })[["elapsed"]]
  test <- corpus[corpus$split == "test", ]
  predicted <- predicted_all[match(test$id, hg$nodes)]
  cbind(
    data.frame(
      dataset = name, method = "tfidf_centroid", weight = weight,
      xi = NA_real_,
      n_train = sum(corpus$split == "train"), n_test = nrow(test),
      n_nodes = hg$n_nodes, n_edges = hg$n_hyperedges
    ),
    .bench_score(test$label, predicted),
    data.frame(build_s = t_build, fit_s = t_fit)
  )
}

# Low-label regime: subsample the training seeds to a fraction (stratified
# by class, at least one seed per class), run both the transductive
# classifier and the centroid baseline on the SAME seed draws, and report
# mean accuracy with the SD over draws. Label spreading is designed for
# exactly this regime; the full-split table above is its worst case.
bench_lowlabel <- function(name, fractions = c(0.01, 0.05, 0.1, 0.2),
                           n_draws = 5L, xi = 0.99,
                           weight = c("tfidf", "n"),
                           dir = file.path("benchmarks", "data"),
                           seed = 1L) {
  weight <- match.arg(weight)
  corpus <- bench_load(name, dir)
  docs <- corpus$text
  names(docs) <- corpus$id
  hg <- hypernets::text_hypergraph(
    docs, weight = weight, stop_words = character(0), sparse = TRUE
  )
  keep <- .bench_giant_component(hg)
  if (length(keep) < hg$n_nodes) {
    hg <- hypernets::text_hypergraph(
      docs[names(docs) %in% keep], weight = weight,
      stop_words = character(0), sparse = TRUE
    )
  }
  train <- corpus[corpus$split == "train" & corpus$id %in% hg$nodes, ]
  test <- corpus[corpus$split == "test", ]
  if (exists(".Random.seed", envir = globalenv())) {
    old_seed <- get(".Random.seed", envir = globalenv())
    on.exit(assign(".Random.seed", old_seed, envir = globalenv()),
            add = TRUE, after = FALSE)
  }
  set.seed(seed)

  x <- hg$incidence
  norms <- sqrt(Matrix::rowSums(x^2))
  norms[norms == 0] <- 1
  x <- Matrix::Diagonal(x = 1 / norms) %*% x

  draw_seeds <- function(fraction) {
    # stratified subsample: ceiling() guarantees >= 1 seed per class
    picked <- unlist(lapply(split(seq_len(nrow(train)), train$label), \(i)
      sample(i, max(1L, ceiling(fraction * length(i))))))
    seeds <- train$label[picked]
    names(seeds) <- train$id[picked]
    seeds
  }

  eval_draw <- function(fraction, draw) {
    seeds <- draw_seeds(fraction)
    fit <- hypernets::hg_classify(hg, labels = seeds, xi = xi,
                                       normalization = "class_mass")
    trans_pred <- fit$predicted[match(test$id, fit$node)]
    rows <- match(names(seeds), hg$nodes)
    classes <- sort(unique(seeds))
    indicator <- Matrix::sparseMatrix(
      i = match(seeds, classes), j = seq_along(rows), x = 1,
      dims = c(length(classes), length(rows))
    )
    centroids <- indicator %*% x[rows, , drop = FALSE]
    c_norms <- sqrt(Matrix::rowSums(centroids^2))
    c_norms[c_norms == 0] <- 1
    centroids <- Matrix::Diagonal(x = 1 / c_norms) %*% centroids
    cent_all <- classes[max.col(as.matrix(Matrix::tcrossprod(x, centroids)),
                                ties.method = "first")]
    cent_pred <- cent_all[match(test$id, hg$nodes)]
    data.frame(
      dataset = name, fraction = fraction, draw = draw,
      n_seeds = length(seeds),
      method = c("hypergraph_zhou_cmn", "tfidf_centroid"),
      accuracy = c(
        mean(!is.na(trans_pred) & trans_pred == test$label),
        mean(!is.na(cent_pred) & cent_pred == test$label)
      )
    )
  }

  grid <- expand.grid(fraction = fractions, draw = seq_len(n_draws))
  per_draw <- do.call(rbind, Map(eval_draw, grid$fraction, grid$draw))
  agg <- aggregate(accuracy ~ dataset + fraction + n_seeds + method,
                   data = per_draw, \(a) c(mean = mean(a), sd = sd(a)))
  data.frame(
    agg[c("dataset", "fraction", "n_seeds", "method")],
    accuracy = agg$accuracy[, "mean"],
    sd = agg$accuracy[, "sd"],
    n_draws = n_draws
  )
}

# HGNN (hg_neural) on one dataset: full published split, n_seeds runs
# (neural training is stochastic -- house rule: never a single seed),
# 90/10 stratified validation inside the training labels, accuracy and
# macro-F1 reported as mean and SD over seeds. No giant-component
# restriction: propagation is defined on disconnected hypergraphs.
bench_neural <- function(name, hidden = 128L, epochs = 600L, lr = 0.01,
                         n_seeds = 3L, weight = c("tfidf", "n"),
                         dir = file.path("benchmarks", "data")) {
  weight <- match.arg(weight)
  corpus <- bench_load(name, dir)
  docs <- corpus$text
  names(docs) <- corpus$id
  hg <- hypernets::text_hypergraph(
    docs, weight = weight, stop_words = character(0), sparse = TRUE
  )
  train <- corpus[corpus$split == "train" & corpus$id %in% hg$nodes, ]
  seeds <- train$label
  names(seeds) <- train$id
  test <- corpus[corpus$split == "test", ]
  one <- function(s) {
    t_fit <- system.time(
      fit <- hypernets::hg_neural(hg, labels = seeds, hidden = hidden,
                                       epochs = epochs, lr = lr, seed = s)
    )[["elapsed"]]
    predicted <- fit$predicted[match(test$id, fit$node)]
    cbind(.bench_score(test$label, predicted), data.frame(fit_s = t_fit))
  }
  runs <- do.call(rbind, lapply(seq_len(n_seeds), one))
  data.frame(
    dataset = name, method = "hgnn", weight = weight,
    hidden = hidden, epochs = epochs, lr = lr, n_seeds = n_seeds,
    n_train = sum(corpus$split == "train"), n_test = nrow(test),
    accuracy = mean(runs$accuracy), sd = sd(runs$accuracy),
    acc_min = min(runs$accuracy), acc_max = max(runs$accuracy),
    macro_f1 = mean(runs$macro_f1), fit_s = mean(runs$fit_s)
  )
}

# HyperGAT (hg_hypergat) on the official repo's own corpus files (they
# keep sentence structure; the TextGCN .clean.txt files have none). The
# official preprocessing keeps words with 5 <= freq < freq(10th most
# common word) -- mirrored via min_count and a top-9 frequency stop
# extension -- and lemmatizes with WordNet, which this package does not
# (documented difference). n_seeds runs, mean/SD/range reported.
bench_hypergat <- function(name, n_seeds = 3L, epochs = 10L,
                           dir = file.path("local_testing_and_equivalence",
                                           "HyperGAT_TextClassification",
                                           "data")) {
  corpus_file <- file.path(dir, paste0(name, "_corpus.txt"))
  labels_file <- file.path(dir, paste0(name, "_labels.txt"))
  stopifnot("official corpus files must exist" =
              file.exists(corpus_file) && file.exists(labels_file))
  docs <- readLines(corpus_file, encoding = "UTF-8", warn = FALSE)
  meta <- utils::read.delim(labels_file, header = FALSE,
                            col.names = c("row", "split", "label"),
                            colClasses = "character", quote = "")
  stopifnot("corpus and labels must align" = length(docs) == nrow(meta))
  names(docs) <- sprintf("%s_%05d", name, seq_along(docs))
  is_train <- grepl("train", meta$split, fixed = TRUE)
  seeds <- meta$label[is_train]
  names(seeds) <- names(docs)[is_train]
  # mirror the official highbar rule: drop the 9 most frequent tokens
  toks <- table(unlist(hypernets:::.thg_tokenize(docs, TRUE)))
  top_toks <- names(head(sort(toks, decreasing = TRUE), 9))
  stops <- union(stop_words_en(), top_toks)
  one <- function(s) {
    t_fit <- system.time(
      fit <- hypernets::hg_hypergat(
        docs, labels = seeds, stop_words = stops, min_count = 5L,
        epochs = epochs, seed = s
      )
    )[["elapsed"]]
    predicted <- fit$predicted[match(names(docs), fit$node)]
    test <- !is_train
    cbind(.bench_score(meta$label[test], predicted[test]),
          data.frame(fit_s = t_fit))
  }
  runs <- do.call(rbind, lapply(seq_len(n_seeds), one))
  data.frame(
    dataset = name, method = "hypergat", epochs = epochs,
    n_seeds = n_seeds, n_train = sum(is_train), n_test = sum(!is_train),
    accuracy = mean(runs$accuracy), sd = sd(runs$accuracy),
    acc_min = min(runs$accuracy), acc_max = max(runs$accuracy),
    macro_f1 = mean(runs$macro_f1), fit_s = mean(runs$fit_s)
  )
}
