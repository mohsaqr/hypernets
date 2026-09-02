# Legal hypergraphs — paper and reproducibility archive

- **Citation**: Corinna Coupette, Dirk Hartung and Daniel Martin Katz (2024),
  “Legal hypergraphs”, *Philosophical Transactions of the Royal Society A*
  382(2270), 20230141. DOI: [10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141).
- **Local paper**:
  [`../papers/2024-PhilTrans-LegalHypergraphs-Coupette.pdf`](../papers/2024-PhilTrans-LegalHypergraphs-Coupette.pdf)
  (19 pages; version of record; CC BY 4.0).
- **What**: a temporal-hypergraph treatment of legal systems, demonstrated on
  more than 70 years of GFCC citation and ICSID collaboration data. Its
  analysis moves from individual vertices and hyperedges to paths,
  projections, motifs and communities.
- **Methods implemented in honets**: temporal growing/interval snapshots;
  hyperedge degree/cardinality/neighbourhood distributions; s-line graphs
  with s-betweenness and s-closeness; log-subhypergraph centrality;
  configuration-model nulls and induced Y/T/O motifs; an association projection in
  which an edge of size |e| contributes `1 / (|e| - 1)` to each incident
  pair; repeated seeded Infomap with AMI-medoid selection; ARI/AMI/NMI
  comparisons; and partition coverage, performance and modularity.
- **API mapping**: `temporal_hypergraph()` / `hypergraph_snapshots()`,
  `hypergraph_edges()`, `hypergraph_edge_centrality()`,
  `hypergraph_centrality(type = "subhypergraph")`,
  `hypergraph_project(method = "association")`, `hypergraph_motifs()`,
  `hypergraph_communities()`, `hypergraph_agreement(method = "ami")`, and
  `hypergraph_community_quality()`.
- **Verified 2026-09-02**: bibliographic metadata, open-access license, the
  19-page paper, and the complete Zenodo source were inspected. Formula-level
  tests match HypergraphX 1.5 and scikit-learn conventions. Rebuilding the
  441-node, 742-case ICSID aggregate from the authors' released CSV reproduces
  their observed motif census exactly: **Y = 478, T = 7, O = 0**.
- **Links**: [publisher / DOI](https://doi.org/10.1098/rsta.2023.0141) ·
  [PubMed Central](https://pmc.ncbi.nlm.nih.gov/articles/PMC10894694/) ·
  [reproducibility archive (Zenodo 8081507)](https://doi.org/10.5281/zenodo.8081507)
