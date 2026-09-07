# Independent assignment oracle for wasserstein_distance(). The R code builds
# persistence-diagram costs; SciPy, not hypernets' Hungarian implementation,
# solves the augmented assignment. Run manually from the package root:
#   Rscript local_testing_and_equivalence/test-equiv-wasserstein-scipy.R

source("R/simplicial_distances.R")

set.seed(20260902)
d1 <- data.frame(
  dimension = 0L,
  birth = runif(12, 0, 2),
  death = runif(12, 3, 6)
)
d2 <- data.frame(
  dimension = 0L,
  birth = runif(9, 0, 2),
  death = runif(9, 3, 6)
)

exchange <- tempfile("wasserstein_equiv_")
dir.create(exchange)
utils::write.csv(d1[, c("birth", "death")],
                 file.path(exchange, "d1.csv"), row.names = FALSE)
utils::write.csv(d2[, c("birth", "death")],
                 file.path(exchange, "d2.csv"), row.names = FALSE)

oracle <- file.path(exchange, "oracle.py")
writeLines(c(
  "import sys, numpy as np, pandas as pd",
  "from scipy.optimize import linear_sum_assignment",
  "root = sys.argv[1]",
  "a = pd.read_csv(f'{root}/d1.csv').to_numpy()",
  "b = pd.read_csv(f'{root}/d2.csv').to_numpy()",
  "n, m = len(a), len(b)",
  "pair = np.max(np.abs(a[:, None, :] - b[None, :, :]), axis=2)",
  "da = np.abs(a[:, 1] - a[:, 0]) / 2",
  "db = np.abs(b[:, 1] - b[:, 0]) / 2",
  "cost = np.zeros((n + m, n + m))",
  "cost[:n, :m] = pair ** 2",
  "cost[:n, m:] = np.repeat((da ** 2)[:, None], n, axis=1)",
  "cost[n:, :m] = np.repeat((db ** 2)[None, :], m, axis=0)",
  "rows, cols = linear_sum_assignment(cost)",
  "np.savetxt(f'{root}/distance.txt', [np.sqrt(cost[rows, cols].sum())])"
), oracle)

status <- system2("python3", c(oracle, exchange), stdout = TRUE, stderr = TRUE)
stopifnot("SciPy oracle must run" = is.null(attr(status, "status")))
theirs <- scan(file.path(exchange, "distance.txt"), quiet = TRUE)
ours <- unname(wasserstein_distance(d1, d2, order = 2))
delta <- abs(ours - theirs)
stopifnot("Wasserstein must match SciPy assignment" = delta < 1e-12)
cat(sprintf("PASS: Wasserstein matches scipy.optimize; abs diff = %.2e\n", delta))
