# SimplicialComplex — CRAN (R)

- **What**: geometric point-cloud TDA package, v0.1.2, ChiChien Wang
  (TDA-R org), MIT. Deps: Matrix, gtools, igraph, ggplot2, geometry, RANN,
  clue.
- **Verified 2026-08-24** (refman): Alpha/Cech/Delaunay/VR/Witness complexes,
  cubical + "flood" filtrations (landmark-based, sparse GF(2) reduction,
  scales to ~1e5 simplices), persistence pairs, Betti/Euler, boundary
  operators, `bottleneck_distance` + `wasserstein_distance` (both keep
  essential death = Inf), `persistence_landscape`, matching visualization.
  Inputs are Euclidean point clouds only — no graphs, no weighted networks,
  no sequences, no q-analysis.
- **Role for us**: oracle for the open honets `wasserstein_distance()`
  (and a second oracle for the shipped `bottleneck_distance`). NOT a
  delegation target (third-party, 0.1.x, different conventions).
- **Collision caution**: its expanding TDA surface overlaps honets names;
  oracle scripts must use explicit namespaces rather than attach both.
- **Links**: https://cran.r-project.org/package=SimplicialComplex ·
  https://github.com/TDA-R/SimplicialComplex
