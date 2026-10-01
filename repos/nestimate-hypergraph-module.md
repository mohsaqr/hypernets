# Nestimate hypergraph module — historical source/oracle

This note records the pre-consolidation implementation. As of 2026-09-29
hypernets' DESCRIPTION does not list Nestimate, and the retained local
identity tests treat it as a frozen historical oracle. **Planned change**
(ROADMAP.md Phase 0b, decided 2026-09-29): hypernets will `Import`
Nestimate (>= 0.8.5) and re-export its memory and simplicial verbs, while
the hypergraph and text families stay hypernets' own; Nestimate itself does
not change. Update this note when that lands. The original five exports
(Nestimate CLAUDE.md, `R/hypergraph*.R`, `R/bipartite_groups.R`,
`R/clique_expansion.R`) were:

- `build_hypergraph()` — hyperedges from ESTIMATED WEIGHTED networks
  (netobject / cograph_network / simplicial_complex input; clique promotion
  with threshold binarisation, Bernoulli `p` sampling, `max_size`,
  pairwise inclusion).
- `bipartite_groups()` — hypergraph from long-format affiliation data
  (player/group columns, optional summed weights). Already ingests any
  document-word table: `bipartite_groups(long, player = "word",
  group = "doc", weight = "n")`.
- `hypergraph_measures()` — tidy suite: hyperdegree, node strength, max edge
  size, co-degree matrix, edge sizes, pairwise overlap / overlap-coefficient
  / Jaccard matrices, density, avg edge size, size distribution,
  intersection profile.
- `hypergraph_centrality()` — clique-expansion eigenvector (CEC) + tensor
  **Z- and H-eigenvector** centralities (power iteration, normalize option).
- `clique_expansion()` — weighted back-projection to a graph.

Base R + BLAS only, zero added dependencies. Equivalence-tested in Nestimate
`local_testing_and_equivalence/test-equiv-hypergraph.R` (40 configs; CEC vs
igraph, Z/H vs clean-room tensor iteration; TOL 1e-10, cosine 1e-6 for
sign/rotation-ambiguous eigenvectors).

All listed capabilities plus Laplacians, clustering, transduction, window and
kNN construction, PageRank and duals now live natively in hypernets, as does
the Wasserstein distance between persistence diagrams (`wasserstein_distance()`,
`R/simplicial_distances.R`, shipped 2026-09-02; SciPy assignment oracle in
`local_testing_and_equivalence/test-equiv-wasserstein-scipy.R`).
