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
- **Relevance**: the centrality set overlaps honets' shipped
  `hypergraph_centrality()` (CEC / Z / H) exactly.
- **Role for us**: **second independent oracle** for the shipped centralities.
  A subprocess-based cross-check for clique, Z and H directions now lives in
  `local_testing_and_equivalence/test-equiv-centrality-xgi.R`; XGI stays out
  of honets' runtime dependency graph.
  Also the reference for measure naming when new measures are considered.
- **Not verified this session**: its generator and stats submodules' exact
  surface.
- **Links**: https://xgi.readthedocs.io · https://github.com/xgi-org/xgi ·
  `pip install xgi`
