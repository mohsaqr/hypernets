# Multi-seed benchmark of hypernets unsupervised hypergraph clustering against
# the actual BERTopic package (Grootendorst 2022) on R8's labelled topics.
# Run from the package root:
#   Rscript benchmarks/run_bertopic_benchmark.R
#
# Python side (benchmarks/bertopic_core.py) needs the local venv:
#   python3.13 -m venv benchmarks/.venv
#   benchmarks/.venv/bin/pip install bertopic sentence-transformers
#
# Outputs (benchmarks/results/, gitignored):
#   bertopic_python_r8.csv        BERTopic per-seed metrics (three variants)
#   bertopic_embeddings_r8.{npy,csv}  all-MiniLM-L6-v2 document embeddings
#   bertopic_r8.csv               per-seed metrics, hypernets arms + BERTopic
#   bertopic_effects.csv          paired effects with bootstrap CIs
#
# hypernets arms:
#   hypernets_hypergraph  tf-idf document-word hypergraph -> Zhou spectral
#                      embedding -> k-means (as in run_umap_hdbscan_benchmark.R)
#   hypernets_knn_sbert   kNN hypergraph on the identical MiniLM embeddings ->
#                      Zhou spectral embedding -> k-means; same vectors as
#                      BERTopic, so this arm isolates the clustering paradigm.

suppressPackageStartupMessages(library(hypernets))
source("benchmarks/harness.R")
source("R/agreement.R")

seeds <- 1:10
dataset <- "R8"
knn_k <- 10L
result_dir <- file.path("benchmarks", "results")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

python <- file.path("benchmarks", ".venv", "bin", "python")
stopifnot("BERTopic venv missing -- see header of this script" =
            file.exists(python))
python_out <- file.path(result_dir, "bertopic_python_r8.csv")
embedding_stem <- file.path(result_dir, "bertopic_embeddings_r8")
numba_cache <- file.path(tempdir(), "hypernets_numba")
dir.create(numba_cache, showWarnings = FALSE)
old_numba <- Sys.getenv("NUMBA_CACHE_DIR", unset = NA_character_)
Sys.setenv(NUMBA_CACHE_DIR = numba_cache, TOKENIZERS_PARALLELISM = "false")
on.exit(if (is.na(old_numba)) Sys.unsetenv("NUMBA_CACHE_DIR") else
  Sys.setenv(NUMBA_CACHE_DIR = old_numba), add = TRUE)

if (!file.exists(python_out)) {
  status <- system2(
    python,
    c("benchmarks/bertopic_core.py", "--dataset", dataset,
      "--data", file.path("benchmarks", "data"), "--output", python_out,
      "--embedding-cache", embedding_stem,
      "--seeds", paste(seeds, collapse = ",")),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(status, "status"))) {
    cat(status, sep = "\n")
    stop("BERTopic runner failed")
  }
}
baseline <- utils::read.csv(python_out, stringsAsFactors = FALSE)
comparators <- c("bertopic_default", "bertopic_k", "bertopic_k_assigned")
stopifnot(
  "cached BERTopic table does not contain the requested seeds" =
    identical(sort(unique(baseline$seed)), seeds),
  "cached BERTopic table does not contain all three variants" =
    setequal(unique(baseline$method), comparators),
  "one row per variant per seed expected" =
    nrow(baseline) == length(seeds) * length(comparators)
)

corpus <- bench_load(dataset)
docs <- corpus$text
names(docs) <- corpus$id
truth <- stats::setNames(corpus$label, corpus$id)
n_classes <- length(unique(truth))

embeddings <- as.matrix(utils::read.csv(paste0(embedding_stem, ".csv"),
                                        header = FALSE))
stopifnot("embedding rows must match the corpus" =
            nrow(embeddings) == length(docs))
dimnames(embeddings) <- list(corpus$id, NULL)

# One spectral embedding per arm (seeded eigen-solve), then k-means per seed;
# both arms use the same k as BERTopic's nr_topics.
spectral_arm <- function(hg) {
  started <- proc.time()[["elapsed"]]
  embedding <- hg_cluster(hg, k = n_classes, seed = seeds[1L], nstart = 25L,
                          what = "embedding")
  list(embedding = embedding,
       seconds = proc.time()[["elapsed"]] - started)
}

build_started <- proc.time()[["elapsed"]]
hg_bag <- text_hypergraph(docs, weight = "tfidf", stop_words = character(),
                          sparse = TRUE)
bag <- spectral_arm(hg_bag)
bag$seconds <- bag$seconds + proc.time()[["elapsed"]] - build_started

build_started <- proc.time()[["elapsed"]]
hg_knn <- text_hypergraph(docs, construction = "knn", k = knn_k,
                          embeddings = embeddings)
knn <- spectral_arm(hg_knn)
knn$seconds <- knn$seconds + proc.time()[["elapsed"]] - build_started

arms <- list(hypernets_hypergraph = bag, hypernets_knn_sbert = knn)

kmeans_rows <- function(method, arm) {
  dims <- grep("^dim", names(arm$embedding), value = TRUE)
  x <- as.matrix(arm$embedding[, dims, drop = FALSE])
  reference <- truth[arm$embedding$node]
  stopifnot("every clustered node needs a label" = !anyNA(reference))
  do.call(rbind, lapply(seeds, function(seed) {
    set.seed(seed)
    started <- proc.time()[["elapsed"]]
    fit <- stats::kmeans(x, centers = n_classes, nstart = 25L,
                         iter.max = 100L)
    data.frame(
      dataset = dataset, method = method, seed = seed,
      ari = .thg_ari(reference, fit$cluster),
      ami = .thg_ami(reference, fit$cluster),
      nmi = .thg_nmi(reference, fit$cluster),
      n_clusters = length(unique(fit$cluster)), outlier_fraction = 0,
      embedding_s = arm$seconds,
      fit_s = proc.time()[["elapsed"]] - started
    )
  }))
}
hypernets_rows <- do.call(rbind, Map(kmeans_rows, names(arms), arms))

keep <- c("dataset", "method", "seed", "ari", "ami", "nmi", "n_clusters",
          "outlier_fraction", "embedding_s", "fit_s")
per_seed <- rbind(hypernets_rows[, keep], baseline[, keep])
per_seed <- per_seed[order(per_seed$method, per_seed$seed), ]
rownames(per_seed) <- NULL

bootstrap_mean <- function(x, n_boot = 5000L, seed = 902L) {
  set.seed(seed)
  draws <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  unname(stats::quantile(draws, c(.025, .975)))
}

metrics <- c("ari", "ami", "nmi")
grid <- expand.grid(reference = names(arms), comparator = comparators,
                    metric = metrics, stringsAsFactors = FALSE)
effects <- do.call(rbind, lapply(seq_len(nrow(grid)), function(i) {
  reference <- grid$reference[i]
  comparator <- grid$comparator[i]
  metric <- grid$metric[i]
  h <- per_seed[per_seed$method == reference, metric]
  b <- per_seed[per_seed$method == comparator, metric]
  delta <- h - b
  ci_h <- bootstrap_mean(h, seed = 100L + i)
  ci_b <- bootstrap_mean(b, seed = 200L + i)
  ci_d <- bootstrap_mean(delta, seed = 300L + i)
  data.frame(
    dataset = dataset, reference = reference, comparator = comparator,
    metric = metric, n_seeds = length(seeds),
    hypernets_mean = mean(h), hypernets_ci_low = ci_h[1L], hypernets_ci_high = ci_h[2L],
    baseline_mean = mean(b), baseline_ci_low = ci_b[1L],
    baseline_ci_high = ci_b[2L], mean_difference = mean(delta),
    difference_ci_low = ci_d[1L], difference_ci_high = ci_d[2L],
    paired_dz = if (stats::sd(delta) > 0) mean(delta) / stats::sd(delta) else NA
  )
}))

utils::write.csv(per_seed, file.path(result_dir, "bertopic_r8.csv"),
                 row.names = FALSE)
utils::write.csv(effects, file.path(result_dir, "bertopic_effects.csv"),
                 row.names = FALSE)
cat("BERTopic", unique(baseline$bertopic_version), "with",
    unique(baseline$embedding_model), "\n")
print(aggregate(cbind(ari, ami, nmi, n_clusters, outlier_fraction) ~ method,
                data = per_seed, FUN = mean), digits = 4)
print(effects[, c("reference", "comparator", "metric", "mean_difference",
                  "difference_ci_low", "difference_ci_high", "paired_dz")],
      digits = 4)
