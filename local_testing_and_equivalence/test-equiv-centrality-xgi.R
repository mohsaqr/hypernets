# Independent XGI oracle for clique, Z and H eigenvector centralities.
# XGI is deliberately not a package dependency. Run in an environment with
# Python xgi, numpy and pandas installed:
#   XGI_PYTHON=/path/to/python \
#     Rscript local_testing_and_equivalence/test-equiv-centrality-xgi.R

source("R/utils.R")
source("R/accessors.R")
source("R/hypergraph_groups.R")
source("R/hypergraph_centrality.R")

python <- Sys.getenv("XGI_PYTHON", unset = "python3")
probe <- suppressWarnings(system2(
  python, c("-c", shQuote("import xgi, numpy, pandas")),
  stdout = TRUE, stderr = TRUE
))
if (!is.null(attr(probe, "status"))) {
  stop("XGI oracle unavailable; set XGI_PYTHON to an environment with xgi")
}

long <- data.frame(
  node = c("a", "b", "c", "a", "b", "d", "a", "c", "e",
           "b", "d", "e"),
  edge = rep(paste0("e", 1:4), each = 3),
  stringsAsFactors = FALSE
)
hg <- group_hypergraph(long, "node", "edge")
ours <- hypergraph_centrality(
  hg, type = c("clique", "Z", "H"), max_iter = 10000L, tol = 1e-12
)

exchange <- tempfile("xgi_centrality_")
dir.create(exchange)
utils::write.csv(long, file.path(exchange, "memberships.csv"), row.names = FALSE)
oracle <- file.path(exchange, "oracle.py")
writeLines(c(
  "import sys, pandas as pd, xgi",
  "root = sys.argv[1]",
  "d = pd.read_csv(f'{root}/memberships.csv')",
  "edges = {e: list(g.node) for e, g in d.groupby('edge', sort=True)}",
  "H = xgi.Hypergraph(edges)",
  "clique = xgi.clique_eigenvector_centrality(H)",
  "zec = xgi.z_eigenvector_centrality(H)",
  "try:",
  "    hec = xgi.h_eigenvector_centrality(H)",
  "except (AttributeError, TypeError):",
  "    hec = xgi.uniform_h_eigenvector_centrality(H)",
  "nodes = sorted(str(n) for n in H.nodes)",
  "out = pd.DataFrame({'node': nodes,",
  "                    'clique': [clique[n] for n in nodes],",
  "                    'Z': [zec[n] for n in nodes],",
  "                    'H': [hec[n] for n in nodes]})",
  "out.to_csv(f'{root}/xgi.csv', index=False)"
), oracle)
status <- system2(python, c(oracle, exchange), stdout = TRUE, stderr = TRUE)
if (!is.null(attr(status, "status"))) {
  cat(status, sep = "\n")
  stop("XGI centrality oracle failed")
}
theirs <- utils::read.csv(file.path(exchange, "xgi.csv"),
                          stringsAsFactors = FALSE)
theirs <- theirs[match(ours$node, theirs$node), ]

cosine <- function(x, y) sum(x * y) / sqrt(sum(x^2) * sum(y^2))
agreement <- vapply(c("clique", "Z", "H"), function(type) {
  cosine(ours[[type]], theirs[[type]])
}, numeric(1))
print(agreement)
stopifnot("clique/Z/H directions must match XGI" = all(agreement > 1 - 1e-7))
cat("PASS: clique, Z and H eigenvector centralities match XGI directions\n")
