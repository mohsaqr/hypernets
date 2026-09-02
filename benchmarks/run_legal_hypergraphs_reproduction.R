#!/usr/bin/env Rscript

# Reproduce the numerical inputs to Figures 6 and 8 of Coupette, Hartung &
# Katz (2024), "Legal hypergraphs", from the authors' Zenodo archive 8081507.
#
# Usage:
#   Rscript benchmarks/run_legal_hypergraphs_reproduction.R \
#     /path/to/legalhypergraphs-8081507 [output-directory]

suppressPackageStartupMessages(library(honets))

args <- commandArgs(trailingOnly = TRUE)
archive <- if (length(args)) args[[1L]] else
  Sys.getenv("LEGAL_HYPERGRAPHS_ARCHIVE")
output <- if (length(args) >= 2L) args[[2L]] else
  file.path("benchmarks", "results", "legal-hypergraphs")
if (!nzchar(archive) || !dir.exists(archive)) {
  stop("Pass the unpacked Zenodo 8081507 archive as the first argument.",
       call. = FALSE)
}
if (!requireNamespace("jsonlite", quietly = TRUE)) {
  stop("This reproduction needs the suggested package `jsonlite`.",
       call. = FALSE)
}
dir.create(output, recursive = TRUE, showWarnings = FALSE)

data_path <- file.path(archive, "data_new")
results_path <- file.path(archive, "results")
write_result <- function(x, name) {
  utils::write.csv(x, file.path(output, name), row.names = FALSE,
                   na = "")
}

# ---- GFCC association graphs (Figure 8) ---------------------------------

gfcc_nodes <- utils::read.csv(
  file.path(data_path, "bverfge_nodes_volumes-001-160.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
gfcc_edges <- utils::read.csv(
  file.path(data_path, "bverfge_edges_volumes-001-160.csv"),
  stringsAsFactors = FALSE, check.names = FALSE
)
decision_date <- stats::setNames(as.Date(gfcc_nodes$date), gfcc_nodes$key)
backward <- decision_date[gfcc_edges$source] > decision_date[gfcc_edges$target]
gfcc_edges <- gfcc_edges[!is.na(backward) & backward, , drop = FALSE]

# Sparse incidence is essential here: the dense 3618 x 46165 matrix would
# store more than 160 million cells.  Explicit `nodes` retains the decisions
# that are isolated in a particular projection, exactly as the paper does.
gfcc <- group_hypergraph(
  gfcc_edges, member = "target", group = "chunk_id",
  nodes = gfcc_nodes$key, sparse = TRUE
)
edge_source <- vapply(split(gfcc_edges$source, gfcc_edges$chunk_id),
                      function(x) unique(x)[[1L]], character(1L))

representations <- data.frame(
  model = c("bh", "bhs", "mh", "mhs"),
  duplicate_edges = c("collapse", "collapse", "count", "count"),
  self_association = c(FALSE, TRUE, FALSE, TRUE),
  stringsAsFactors = FALSE
)
projection_rows <- lapply(seq_len(nrow(representations)), function(i) {
  spec <- representations[i, ]
  projection <- hg_project(
    gfcc, method = "association", what = "matrix",
    duplicate_edges = spec$duplicate_edges,
    self_association = spec$self_association,
    edge_source = if (spec$self_association) edge_source else NULL
  )
  upper <- upper.tri(projection)
  data.frame(
    model = spec$model,
    n_nodes = nrow(projection),
    n_edges = sum(projection[upper] != 0),
    total_weight = sum(projection) / 2,
    n_isolates = sum(Matrix::rowSums(projection != 0) == 0),
    stringsAsFactors = FALSE
  )
})
projection_summary <- do.call(rbind, projection_rows)
expected_projection <- data.frame(
  model = c("bh", "bhs", "mh", "mhs"),
  n_nodes = 3618L,
  n_edges = c(25114L, 51494L, 25114L, 51494L),
  total_weight = c(16978, 20205, 23196, 26423),
  n_isolates = c(664L, 216L, 664L, 216L)
)
stopifnot(
  gfcc$n_nodes == 3618L,
  gfcc$n_hyperedges == 46165L,
  Matrix::nnzero(gfcc$incidence) == 77187L,
  isTRUE(all.equal(projection_summary, expected_projection,
                   check.attributes = FALSE))
)
write_result(projection_summary, "gfcc-association-summary.csv")

# The archive records the AMI-selected medoid in the first column of its
# evaluation table.  Read the corresponding raw JSON partitions rather than
# relying on the old pandas pickle, which is not portable across pandas 1/2.
evaluation <- utils::read.csv(
  file.path(results_path, "gfcc_clusterings_evaluation.csv"),
  stringsAsFactors = FALSE
)
medoid_ids <- unique(evaluation$medoid)
model_order <- c("bb", "bbu", "mb", "mbu", "bh", "bhs", "mh", "mhs")
medoid_ids <- medoid_ids[match(model_order,
                               sub("[0-9]+$", "", medoid_ids))]

read_partition <- function(identifier) {
  model <- sub("[0-9]+$", "", identifier)
  seed <- as.integer(sub("^[a-z]+", "", identifier))
  path <- file.path(results_path, "clusterings", sprintf(
    "gfcc_%s_s-%02d.json.gz", model, seed
  ))
  connection <- gzfile(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  parsed <- jsonlite::fromJSON(paste(readLines(connection, warn = FALSE),
                                     collapse = ""))
  data.frame(node = as.character(parsed$key),
             label = as.character(parsed$module),
             stringsAsFactors = FALSE)
}
medoids <- stats::setNames(lapply(medoid_ids, read_partition), model_order)

size_rows <- lapply(seq_along(medoids), function(i) {
  sizes <- sort(table(medoids[[i]]$label), decreasing = TRUE)
  data.frame(
    model = names(medoids)[i], rank = seq_along(sizes),
    n_nodes = as.integer(sizes), stringsAsFactors = FALSE
  )
})
medoid_sizes <- do.call(rbind, size_rows)
medoid_summary <- do.call(rbind, lapply(seq_along(medoids), function(i) {
  sizes <- medoid_sizes$n_nodes[medoid_sizes$model == names(medoids)[i]]
  data.frame(
    model = names(medoids)[i],
    medoid_seed = as.integer(sub("^[a-z]+", "", medoid_ids[[i]])),
    n_communities = length(sizes), n_singletons = sum(sizes == 1L),
    n_nontrivial = sum(sizes > 1L), largest = sizes[[1L]],
    second = sizes[[2L]], balance = sizes[[2L]] / sizes[[1L]],
    stringsAsFactors = FALSE
  )
}))
write_result(medoid_summary, "gfcc-medoid-summary.csv")
write_result(medoid_sizes, "gfcc-medoid-sizes.csv")

similarity_rows <- list()
cursor <- 1L
for (i in seq_len(length(medoids) - 1L)) {
  for (j in (i + 1L):length(medoids)) {
    score <- hg_agreement(medoids[[i]], medoids[[j]],
                          method = c("ami", "ari", "nmi"))
    similarity_rows[[cursor]] <- data.frame(
      model_a = names(medoids)[i], model_b = names(medoids)[j],
      ami = score$ami, ari = score$ari, nmi = score$nmi,
      stringsAsFactors = FALSE
    )
    cursor <- cursor + 1L
  }
}
medoid_similarity <- do.call(rbind, similarity_rows)
write_result(medoid_similarity, "gfcc-medoid-similarity.csv")

# Re-evaluate each association medoid on its own honets projection and verify
# the four diagnostics reported in the authors' archived evaluation table.
quality_rows <- lapply(seq_len(nrow(representations)), function(i) {
  spec <- representations[i, ]
  quality <- hg_community_quality(
    gfcc, medoids[[spec$model]], duplicate_edges = spec$duplicate_edges,
    self_association = spec$self_association,
    edge_source = if (spec$self_association) edge_source else NULL
  )
  data.frame(model = spec$model, quality, stringsAsFactors = FALSE)
})
association_quality <- do.call(rbind, quality_rows)
for (i in seq_len(nrow(association_quality))) {
  model <- association_quality$model[[i]]
  archived <- evaluation[
    evaluation$medoid == medoid_ids[match(model, model_order)] &
      evaluation$evaluated_on == model, , drop = FALSE
  ]
  stopifnot(
    abs(association_quality$weighted_coverage[[i]] -
          archived$weighted_coverage[[1L]]) < 1e-12,
    abs(association_quality$coverage[[i]] - archived$coverage[[1L]]) < 1e-12,
    abs(association_quality$performance[[i]] -
          archived$performance[[1L]]) < 1e-12,
    abs(association_quality$modularity[[i]] -
          archived$modularity[[1L]]) < 1e-12
  )
}
write_result(association_quality, "gfcc-association-medoid-quality.csv")

# ---- ICSID temporal hyperedge centrality (Figure 6) ---------------------

icsid <- utils::read.csv(
  file.path(data_path, "icsid_normalized.csv"), stringsAsFactors = FALSE,
  check.names = FALSE
)
icsid$date_of_constitution_of_tribunal <-
  as.Date(icsid$date_of_constitution_of_tribunal)
icsid$date_concluded <- as.Date(icsid$date_concluded)
cutoff <- as.Date("2023-06-15")
icsid$date_concluded[is.na(icsid$date_concluded)] <- cutoff
icsid_temporal <- temporal_hypergraph(
  icsid,
  member = c("president_name", "arbitrator_1_name", "arbitrator_2_name"),
  edge = "caseno", start = "date_of_constitution_of_tribunal",
  end = "date_concluded", evolution = "interval"
)

centrality_rows <- list()
for (representation in c("active", "aggregate")) {
  snapshot <- hypergraph_snapshot(
    icsid_temporal, at = cutoff,
    mode = if (representation == "active") "active" else "all",
    multiedges = FALSE
  )
  score <- hg_edge_centrality(snapshot, s = 1, normalized = TRUE)
  score$rank <- ave(-score$value, score$measure,
                    FUN = function(x) rank(x, ties.method = "min"))
  score$representation <- representation
  score$tribunal <- vapply(score$edge, function(edge) {
    index <- match(edge, colnames(snapshot$incidence))
    paste(sort(snapshot$nodes[snapshot$hyperedges[[index]]]), collapse = " | ")
  }, character(1L))
  centrality_rows[[representation]] <- score[
    score$rank <= 10L,
    c("representation", "measure", "rank", "edge", "tribunal", "value")
  ]
}
icsid_centrality <- do.call(rbind, centrality_rows)
rownames(icsid_centrality) <- NULL
stopifnot(
  nrow(icsid) == 742L,
  centrality_rows$active$edge[
    centrality_rows$active$measure == "betweenness" &
      centrality_rows$active$rank == 1L
  ] == "ARB/17/21",
  abs(centrality_rows$active$value[
    centrality_rows$active$measure == "betweenness" &
      centrality_rows$active$rank == 1L
  ] - 0.026540129226878584) < 1e-14,
  abs(centrality_rows$aggregate$value[
    centrality_rows$aggregate$measure == "closeness" &
      centrality_rows$aggregate$rank == 1L
  ] - 0.5480570264418387) < 1e-14
)
write_result(icsid_centrality, "icsid-figure6-centrality.csv")

# Compact, dependency-free versions of the three result figures.  They are
# intended as numerical reproductions; the original hand-tuned label layout
# remains available as PDFs in the authors' archive.
grDevices::pdf(file.path(output, "gfcc-figure8b-cluster-sizes.pdf"),
               width = 8, height = 6)
plot(NA, xlim = range(medoid_sizes$n_nodes), ylim = c(1, 120), log = "x",
     xlab = "Cluster size", ylab = "Number of clusters at least this large")
palette <- grDevices::hcl.colors(length(medoids), "Dark 3")
for (i in seq_along(medoids)) {
  sizes <- sort(medoid_sizes$n_nodes[medoid_sizes$model == names(medoids)[i]])
  lines(sizes, rev(seq_along(sizes)), col = palette[[i]], lwd = 2)
}
legend("topright", legend = names(medoids), col = palette, lwd = 2,
       bty = "n", ncol = 2)
grDevices::dev.off()

grDevices::pdf(file.path(output, "gfcc-figure8c-medoid-similarity.pdf"),
               width = 9, height = 5)
old_par <- par(mfrow = c(1, 2), mar = c(5, 5, 3, 1))
for (metric in c("ami", "ari")) {
  matrix_score <- diag(1, length(model_order))
  dimnames(matrix_score) <- list(model_order, model_order)
  for (i in seq_len(nrow(medoid_similarity))) {
    a <- medoid_similarity$model_a[[i]]
    b <- medoid_similarity$model_b[[i]]
    matrix_score[a, b] <- matrix_score[b, a] <- medoid_similarity[[metric]][[i]]
  }
  image(seq_along(model_order), seq_along(model_order), t(matrix_score),
        zlim = c(0, 1), col = grDevices::hcl.colors(32, "Inferno"),
        axes = FALSE, xlab = "", ylab = "", main = toupper(metric))
  axis(1, seq_along(model_order), model_order)
  axis(2, seq_along(model_order), model_order, las = 1)
}
par(old_par)
grDevices::dev.off()

grDevices::pdf(file.path(output, "icsid-figure6-centrality-ranks.pdf"),
               width = 10, height = 8)
old_par <- par(mfrow = c(2, 2), mar = c(5, 9, 3, 1))
for (representation in c("active", "aggregate")) {
  for (metric in c("betweenness", "closeness")) {
    panel <- icsid_centrality[
      icsid_centrality$representation == representation &
        icsid_centrality$measure == metric, , drop = FALSE
    ]
    panel <- panel[order(panel$value), , drop = FALSE]
    graphics::barplot(panel$value, names.arg = panel$edge, horiz = TRUE,
                      las = 1, border = NA, col = "#377eb8",
                      main = paste(representation, metric), xlab = "centrality")
  }
}
par(old_par)
grDevices::dev.off()

cat("PASS: Legal hypergraphs full-data reproduction\n")
cat("Outputs:", normalizePath(output), "\n")
