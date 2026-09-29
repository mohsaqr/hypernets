# TDA in Python — gudhi / ripser.py / giotto-tda

- **What** (one note for the family): GUDHI (INRIA; the reference C++/Python
  TDA library — complexes, persistence, bottleneck/Wasserstein), ripser.py
  (fast VR persistence), giotto-tda (scikit-learn-compatible TDA pipelines).
- **Verified**: established primary project documentation inspected for the
  persistence/diagram-distance scope. GUDHI 3.13.0 is current on PyPI
  (checked 2026-09-29).
- **Status (2026-09-29): none of these is used as an oracle yet.** The
  oracles actually in use are the R package TDAstats (ripser-backed;
  `local_testing_and_equivalence/test-equiv-persistent_homology.R`) and SciPy
  `linear_sum_assignment` for `wasserstein_distance()`
  (`test-equiv-wasserstein-scipy.R`). `bottleneck_distance()` and
  `persistence_landscape()` have no external oracle.
- **Role for us**: planned oracles (ROADMAP.md Phase 1) — GUDHI SimplexTree
  Betti numbers, `bottleneck_distance`, `wasserstein_distance` (the field
  reference) and `representations.Landscape`. Nothing else is needed from
  them — hypernets' TDA scope (clique/pathway complexes, q-analysis) is
  relational, not point-cloud.
- **Links**: https://gudhi.inria.fr · https://github.com/scikit-tda/ripser.py ·
  https://github.com/giotto-ai/giotto-tda
