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
- **Role for us**: *candidate* independent oracle for the shipped
  `wasserstein_distance()` and `bottleneck_distance()`. **Not used yet**
  (checked 2026-09-29: no test in `tests/` or
  `local_testing_and_equivalence/` calls it). The oracles actually in use
  are TDAstats (ripser) for persistent homology and SciPy
  `linear_sum_assignment` for the Wasserstein assignment. NOT a delegation
  target (third-party, 0.1.x, different conventions).
- **Collision caution**: its expanding TDA surface overlaps hypernets names;
  oracle scripts must use explicit namespaces rather than attach both.
- **Links**: https://cran.r-project.org/package=SimplicialComplex ·
  https://github.com/TDA-R/SimplicialComplex
