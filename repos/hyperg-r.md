# HyperG — CRAN (R)

- **What**: the only general hypergraph toolkit on CRAN. v1.0.0, David J.
  Marchette, single release 2021-03-04, dormant since. Deps: igraph, mclust,
  proxy, RSpectra, gtools, Matrix.
- **Verified 2026-08-24** (refman): ~50 functions — constructors/converters
  (incidence, edgelist, literal, membership, dual, line graph, complement,
  add/delete), structure predicates (Helly, conformal, linear, hypertree),
  spectral (`hypergraph_laplacian_matrix`, ASE/LSE embeddings, spectrum,
  entropy), `cluster_spectral` (mclust on the embedding), `kCores`, seven
  random samplers (gnp, SBM, k-uniform/regular, geometric, epsilon, knn),
  plotting. **Undirected, UNWEIGHTED only** (stated in its own docs); no
  tensor centralities; no sequence/text input.
- **Role for us**: local oracle for the unweighted Laplacian, kNN and dual
  cases already shipped in honets, and for the open random-generator work.
  It is never a declared dependency. Attaching both packages masks
  `dual_hypergraph()` and `knn_hypergraph()`; oracle scripts must namespace
  both sides explicitly.
- **Links**: https://cran.r-project.org/package=HyperG
