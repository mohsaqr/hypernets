# Nestimate hypergraph module — historical source/oracle

This note records the pre-consolidation implementation. hypernets now owns the
methods and does not depend on Nestimate. No current hypernets work should edit
or call Nestimate; retained local identity tests treat it only as a frozen
historical oracle. The original five exports were:
CLAUDE.md and `R/hypergraph*.R`, `R/bipartite_groups.R`,
`R/clique_expansion.R`):

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
kNN construction, PageRank and duals now live natively in hypernets. Wasserstein
distance remains open in hypernets.
