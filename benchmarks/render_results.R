# Render benchmarks/results/*.csv as markdown tables (stdout), so every
# number in RESULTS.md and the pkgdown article is generated, never typed.
# Usage: /usr/local/bin/Rscript benchmarks/render_results.R

results_dir <- file.path("benchmarks", "results")
files <- list.files(results_dir, pattern = "^(R8|R52|mr|ohsumed|20ng)\\.csv$",
                    full.names = TRUE)
stopifnot("no result files found -- run benchmarks/run_benchmarks.R" =
            length(files) > 0L)
results <- do.call(rbind, lapply(files, utils::read.csv))

fmt <- function(x, digits = 4) formatC(x, format = "f", digits = digits)

row_md <- function(r) {
  sprintf("| %s | %s | %s [%s, %s] | %s | %d | %.1f | %.2f |",
          r$dataset, r$method, fmt(r$accuracy), fmt(r$ci_low),
          fmt(r$ci_high), fmt(r$macro_f1), r$n_unscored, r$build_s, r$fit_s)
}

cat("| Dataset | Method | Accuracy [95% CI] | Macro-F1 | Unscored |",
    "Build (s) | Fit (s) |\n")
cat("|---|---|---|---|---|---|---|\n")
ord <- order(match(results$dataset, c("20ng", "R8", "R52", "ohsumed", "mr")),
             results$method)
invisible(lapply(seq_len(nrow(results)),
                 \(i) cat(row_md(results[ord, ][i, ]), "\n")))

lowlabel_file <- file.path(results_dir, "lowlabel.csv")
if (file.exists(lowlabel_file)) {
  low <- utils::read.csv(lowlabel_file)
  cat("\n| Dataset | Seed fraction | Seeds | Method | Accuracy (mean) |",
      "SD over draws |\n")
  cat("|---|---|---|---|---|---|\n")
  lord <- order(low$dataset, low$fraction, low$method)
  invisible(lapply(seq_len(nrow(low)), \(i) {
    r <- low[lord, ][i, ]
    cat(sprintf("| %s | %d%% | %d | %s | %s | %s |\n",
                r$dataset, round(100 * r$fraction), r$n_seeds, r$method,
                fmt(r$accuracy), fmt(r$sd)))
  }))
}

neural_files <- file.path(results_dir, c("neural.csv", "neural_default.csv"))
neural_files <- neural_files[file.exists(neural_files)]
if (length(neural_files) > 0L) {
  neural <- do.call(rbind, lapply(neural_files, utils::read.csv))
  cat("\n| Dataset | Method | Accuracy (mean of seeds) | SD | Range |",
      "Macro-F1 | Fit (s) |\n")
  cat("|---|---|---|---|---|---|---|\n")
  nord <- order(match(neural$dataset, c("20ng", "R8", "R52", "ohsumed", "mr")))
  invisible(lapply(seq_len(nrow(neural)), \(i) {
    r <- neural[nord, ][i, ]
    cat(sprintf("| %s | hgnn (%d ep, lr %.3g) | %s | %s | [%s, %s] | %s | %.0f |\n",
                r$dataset, r$epochs, r$lr, fmt(r$accuracy), fmt(r$sd),
                fmt(r$acc_min), fmt(r$acc_max), fmt(r$macro_f1), r$fit_s))
  }))
}

hypergat_file <- file.path(results_dir, "hypergat.csv")
if (file.exists(hypergat_file)) {
  hga <- utils::read.csv(hypergat_file)
  cat("\n| Dataset | Method | Accuracy (mean of seeds) | SD | Range |",
      "Macro-F1 | Fit (s) |\n")
  cat("|---|---|---|---|---|---|---|\n")
  invisible(lapply(seq_len(nrow(hga)), \(i) {
    r <- hga[i, ]
    cat(sprintf("| %s | hypergat (%d ep, semantic = %s) | %s | %s | [%s, %s] | %s | %.0f |\n",
                r$dataset, r$epochs, r$semantic, fmt(r$accuracy), fmt(r$sd),
                fmt(r$acc_min), fmt(r$acc_max), fmt(r$macro_f1), r$fit_s))
  }))

  # Compare like with like: sentence-only runs (semantic = "none") against
  # the paper's "w/o semantic" ablation (Ding et al. 2020, Table 4); runs
  # with semantic hyperedges against full HyperGAT (Table 2).
  published <- utils::read.csv(file.path("benchmarks", "published_table2.csv"))
  hga$comparator <- ifelse(hga$semantic == "none", "HyperGAT w/o semantic",
                           "HyperGAT")
  full <- published[published$model == "HyperGAT", c("dataset", "accuracy")]
  names(full) <- c("dataset", "full_accuracy")
  cmp <- merge(hga, published[, c("dataset", "model", "accuracy", "sd")],
               by.x = c("dataset", "comparator"), by.y = c("dataset", "model"),
               suffixes = c("", "_published"))
  cmp <- merge(cmp, full, by = "dataset")
  cat("\n| Dataset | Ours (semantic) | Published comparator | Published",
      "accuracy | Ours - comparator | Full HyperGAT (Table 2) |\n")
  cat("|---|---|---|---|---|---|\n")
  invisible(lapply(seq_len(nrow(cmp)), \(i) {
    r <- cmp[i, ]
    cat(sprintf("| %s | %s (%s) | %s (Table %s) | %s +/- %s | %+.4f | %s |\n",
                r$dataset, fmt(r$accuracy), r$semantic, r$comparator,
                if (r$comparator == "HyperGAT") "2" else "4",
                fmt(r$accuracy_published), fmt(r$sd_published),
                r$accuracy - r$accuracy_published, fmt(r$full_accuracy)))
  }))
}
