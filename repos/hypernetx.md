# HyperNetX (HNX) — PNNL (Python)

- **What**: hypergraph analysis library from Pacific Northwest National
  Laboratory; data structures on sparse incidence matrices plus an algorithms
  layer.
- **Verified 2026-08-25** (GitHub source, `hypernetx/algorithms/`): modules =
  `clustering` (`hypergraph_modularity.py` — Kumar-style modularity;
  `laplacians_clustering.py`), `homology`, `generation`, `matching`,
  `metrics`, `temporal`, `concepts`, `embeddings`.
- **The key fact**: `laplacians_clustering.py` implements the **Hayashi,
  Aksoy, Park & Park (2020) EDVW pipeline exactly** — `prob_trans()`
  (edge-dependent-vertex-weight random-walk transition matrix), `get_pi()`
  (stationary distribution), `norm_lap()` (normalized Laplacian),
  `spec_clus()` (spectral clustering) — citing the paper in the docstrings.
  Sinan Aksoy is at PNNL, so this is the author-adjacent reference
  implementation.
- **Role for us**: primary local oracle for the shipped weighted Laplacian,
  RDC-Spec clustering and EDVW PageRank paths. Its homology module remains a
  possible second TDA oracle. HyperNetX does not supply the Hayashi
  SymNMF/JointNMF branches; hypernets now implements RDC-Sym (Algorithm 2),
  J-NMF (Eq. 18), and JS-NMF (Eq. 19), verified directly against their
  paper equations and objective invariants.
- **Links**: https://github.com/pnnl/HyperNetX ·
  https://hypernetx.readthedocs.io · `pip install hypernetx`
