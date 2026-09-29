# HyperNetX (HNX) — PNNL (Python)

- **What**: hypergraph analysis library from Pacific Northwest National
  Laboratory; data structures on sparse incidence matrices plus an algorithms
  layer.
- **Verified 2026-08-25** (GitHub source, `hypernetx/algorithms/`): modules =
  `clustering` (`hypergraph_modularity.py` — Kumar-style modularity;
  `laplacians_clustering.py`), `homology`, `generation`, `matching`,
  `metrics`, `temporal`, `concepts`, `embeddings`.
- **Verified 2026-09-29 on 2.4.3** (current on PyPI, released 2026-07-23;
  wheel source inspected): `metrics/s_centrality_measures.py` has
  `s_betweenness_centrality`, `s_closeness_centrality`,
  `s_harmonic_centrality`, `s_harmonic_closeness_centrality`,
  `s_eccentricity`; `clustering/hypergraph_modularity.py` has `modularity`,
  `kumar` (Kumar et al. 2020), `last_step`, `two_section`;
  `generation/generative_models.py` has `erdos_renyi_hypergraph`,
  `chung_lu_hypergraph`, `dcsbm_hypergraph`; `temporal/contagion.py` has
  `discrete_SIR`, `discrete_SIS`, `Gillespie_SIR`; `hif.py` has `to_hif` /
  `from_hif`; `homology` has `betti_numbers`, `homology_basis`. Planned
  oracles for ROADMAP.md Phases 2 and 4 (s-centralities, Kumar modularity,
  generators, HIF).
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
