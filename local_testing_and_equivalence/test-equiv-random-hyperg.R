# Exact seeded equivalence for the three directly shared HyperG samplers.
# The SBM implementations share the published graph-augmentation definition,
# but use independent native graph-SBM engines and are tested by invariants.
# Run manually from the package root:
#   Rscript local_testing_and_equivalence/test-equiv-random-hyperg.R

stopifnot("HyperG must be installed" = requireNamespace("HyperG", quietly = TRUE))
source("R/utils.R")
source("R/hypergraph_random.R")

as_node_edge <- function(h) (t(as.matrix(HyperG::incidence_matrix(h))) > 0) * 1L

set.seed(71)
ref_gnp <- HyperG::sample_gnp_hypergraph(n = 18, m = 11, p = .25)
our_gnp <- hg_sample_gnp(n = 18, m = 11, p = .25, seed = 71)
stopifnot("G(n,p) seeded incidence parity" =
            identical(unname(our_gnp$incidence), unname(as_node_edge(ref_gnp))))

set.seed(72)
ref_uniform <- HyperG::sample_k_uniform_hypergraph(n = 18, m = 11, k = 4)
our_uniform <- hg_sample_uniform(n = 18, m = 11, k = 4, seed = 72)
stopifnot("k-uniform seeded incidence parity" =
            identical(unname(our_uniform$incidence),
                      unname(as_node_edge(ref_uniform))))

set.seed(73)
ref_regular <- HyperG::sample_k_regular_hypergraph(n = 18, m = 11, k = 3)
our_regular <- hg_sample_regular(n = 18, m = 11, k = 3, seed = 73)
stopifnot("k-regular seeded incidence parity" =
            identical(unname(our_regular$incidence),
                      unname(as_node_edge(ref_regular))))

cat("PASS: G(n,p), k-uniform and k-regular seeded incidence match HyperG\n")
