# Multi-seed UMAP/HDBSCAN benchmark against honets unsupervised
# hypergraph clustering on R8's labelled topics. Run from the package root:
#   Rscript benchmarks/run_umap_hdbscan_benchmark.R
#
# Outputs:
#   benchmarks/results/umap_hdbscan_r8.csv      per-seed metrics
#   benchmarks/results/umap_hdbscan_effects.csv paired effects

suppressPackageStartupMessages(library(honets))
source("benchmarks/harness.R")
source("R/agreement.R")

seeds <- 1:10
dataset <- "R8"
result_dir <- file.path("benchmarks", "results")
dir.create(result_dir, showWarnings = FALSE, recursive = TRUE)

python_out <- file.path(result_dir, "umap_hdbscan_python_r8.csv")
numba_cache <- file.path(tempdir(), "honets_numba")
dir.create(numba_cache, showWarnings = FALSE)
old_numba <- Sys.getenv("NUMBA_CACHE_DIR", unset = NA_character_)
Sys.setenv(NUMBA_CACHE_DIR = numba_cache)
on.exit(if (is.na(old_numba)) Sys.unsetenv("NUMBA_CACHE_DIR") else
  Sys.setenv(NUMBA_CACHE_DIR = old_numba), add = TRUE)

if (!file.exists(python_out)) {
  status <- system2(
    "python3",
    c("benchmarks/umap_hdbscan_core.py", "--dataset", dataset,
      "--data", file.path("benchmarks", "data"), "--output", python_out,
      "--seeds", paste(seeds, collapse = ",")),
    stdout = TRUE, stderr = TRUE
  )
  if (!is.null(attr(status, "status"))) {
    cat(status, sep = "\n")
    stop("UMAP/HDBSCAN assignment runner failed")
  }
}
baseline <- utils::read.csv(python_out, stringsAsFactors = FALSE)
if (!identical(sort(unique(baseline$seed)), seeds)) {
  stop("cached UMAP/HDBSCAN table does not contain the requested seeds")
}

corpus <- bench_load(dataset)
docs <- corpus$text
names(docs) <- corpus$id
truth <- stats::setNames(corpus$label, corpus$id)

build_started <- proc.time()[["elapsed"]]
hg <- text_hypergraph(
  docs, weight = "tfidf", stop_words = character(), sparse = TRUE
)
embedding <- hg_cluster(
  hg, k = length(unique(truth)), seed = seeds[1L], nstart = 25L,
  what = "embedding"
)
build_seconds <- proc.time()[["elapsed"]] - build_started
dims <- grep("^dim", names(embedding), value = TRUE)
x <- as.matrix(embedding[, dims, drop = FALSE])

honets_rows <- lapply(seeds, function(seed) {
  set.seed(seed)
  started <- proc.time()[["elapsed"]]
  fit <- stats::kmeans(x, centers = length(unique(truth)), nstart = 25L,
                       iter.max = 100L)
  elapsed <- proc.time()[["elapsed"]] - started
  reference <- truth[embedding$node]
  data.frame(
    dataset = dataset, method = "honets_hypergraph", seed = seed,
    ari = .thg_ari(reference, fit$cluster),
    ami = .thg_ami(reference, fit$cluster),
    nmi = .thg_nmi(reference, fit$cluster),
    n_clusters = length(unique(fit$cluster)), outlier_fraction = 0,
    embedding_s = build_seconds, fit_s = elapsed
  )
})
per_seed <- rbind(do.call(rbind, honets_rows), baseline)
per_seed <- per_seed[order(per_seed$method, per_seed$seed), ]
rownames(per_seed) <- NULL

bootstrap_mean <- function(x, n_boot = 5000L, seed = 902L) {
  set.seed(seed)
  draws <- replicate(n_boot, mean(sample(x, replace = TRUE)))
  unname(stats::quantile(draws, c(.025, .975)))
}

metrics <- c("ari", "ami", "nmi")
comparators <- c("umap_hdbscan", "umap_kmeans_k8")
effects <- do.call(rbind, lapply(seq_along(comparators), function(j) {
  comparator <- comparators[j]
  do.call(rbind, lapply(metrics, function(metric) {
    metric_i <- match(metric, metrics)
    h <- per_seed[per_seed$method == "honets_hypergraph", metric]
    b <- per_seed[per_seed$method == comparator, metric]
    delta <- h - b
    ci_h <- bootstrap_mean(h, seed = 100L + metric_i)
    ci_b <- bootstrap_mean(b, seed = 200L + 10L * j + metric_i)
    ci_d <- bootstrap_mean(delta, seed = 300L + 10L * j + metric_i)
    data.frame(
      dataset = dataset, comparator = comparator, metric = metric,
      n_seeds = length(seeds), honets_mean = mean(h),
      honets_ci_low = ci_h[1L], honets_ci_high = ci_h[2L],
      baseline_mean = mean(b), baseline_ci_low = ci_b[1L],
      baseline_ci_high = ci_b[2L], mean_difference = mean(delta),
      difference_ci_low = ci_d[1L], difference_ci_high = ci_d[2L],
      paired_dz = if (stats::sd(delta) > 0) mean(delta) / stats::sd(delta) else NA
    )
  }))
}))

utils::write.csv(per_seed,
                 file.path(result_dir, "umap_hdbscan_r8.csv"),
                 row.names = FALSE)
utils::write.csv(effects,
                 file.path(result_dir, "umap_hdbscan_effects.csv"),
                 row.names = FALSE)
print(effects, digits = 4)
