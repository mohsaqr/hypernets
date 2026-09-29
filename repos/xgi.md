# XGI — CompleX Group Interactions (Python)

- **What**: the modern, actively maintained Python library for higher-order
  networks (hypergraphs + simplicial complexes): data structures, measures,
  algorithms, generators, drawing.
- **Verified 2026-08-25** (readthedocs API pages): algorithms modules =
  assortativity, centrality, clustering coefficients, connected components,
  shortest path, simpliciality. Centrality module includes
  `h_eigenvector_centrality` (Qi tensor H-eigenpair / Benson),
  `z_eigenvector_centrality`, `uniform_h_eigenvector_centrality`,
  `clique_eigenvector_centrality`, `katz_centrality`,
  `line_vector_centrality`.
- **Relevance**: the centrality set overlaps hypernets' shipped
  `hypergraph_centrality()` (CEC / Z / H) exactly.
- **Role for us**: **second independent oracle** for the shipped centralities.
  A subprocess-based cross-check for clique, Z and H directions now lives in
  `local_testing_and_equivalence/test-equiv-centrality-xgi.R`; XGI stays out
  of hypernets' runtime dependency graph.
  Also the reference for measure naming when new measures are considered.
- **Verified 2026-09-29 on 0.10.2** (current on PyPI, released 2026-05-15;
  wheel source inspected): `linalg` has `hodge_laplacian`,
  `multiorder_laplacian`, `laplacian` and `normalized_hypergraph_laplacian`;
  `algorithms.simpliciality` has `simplicial_fraction`,
  `edit_simpliciality`, `face_edit_simpliciality`,
  `simplicial_edit_distance`, `mean_face_edit_distance`; `readwrite.hif` has
  `read_hif` / `write_hif` (HIF); generators include `random_hypergraph`,
  `chung_lu_hypergraph`, `dcsbm_hypergraph`, `watts_strogatz_hypergraph`,
  `uniform_HSBM`, `uniform_HPPM`, `uniform_erdos_renyi_hypergraph`,
  `uniform_hypergraph_configuration_model`, `random_simplicial_complex`,
  `flag_complex`; also `katz_centrality`, `local_clustering_coefficient`,
  `degree_assortativity`. These are the planned oracles for ROADMAP.md
  Phase 2 (Hodge Laplacians, simpliciality, clustering, Katz, HIF).
- **Links**: https://xgi.readthedocs.io · https://github.com/xgi-org/xgi ·
  `pip install xgi`
