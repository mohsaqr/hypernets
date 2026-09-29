# HypergraphX (HGX) — Battiston group (Python)

- **What**: Python library for higher-order network analysis from the group
  behind the higher-order-networks physics reviews (Battiston et al.):
  measures, community detection, motif analysis, filtering, visualization.
- **Version**: 1.8.0 is current on PyPI (released 2026-05-18, checked
  2026-09-29). The API conventions were inspected on 1.5 (2026-09-02):
  static/directed/temporal/multiplex structures, motifs, generators, random
  walks and multiple community algorithms.
- **How it is used**: *not* as a live test dependency — no test imports
  hypergraphx. Its conventions enter through the Legal Hypergraphs
  reproduction: frozen formula fixtures for the motif (Y/T/O),
  configuration-null and s-centrality results, and the Zenodo 8081507
  notebooks (`local_testing_and_equivalence/test-legal-hypergraphs-oracle.R`).
- **Role for us**: best source for future broader motif (`compute_motifs`,
  order 3/4, directed) and generator coverage (ROADMAP.md Phases 1-2).
- **Links**: https://github.com/HGX-Team/hypergraphx · `pip install hypergraphx`
