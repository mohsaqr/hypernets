# honets text and hypergraph backlog

This is the maintained backlog after the texthypergraph merge. Historical
Nestimate-era tasks are recorded in `ROADMAP-text.md`; honets does not depend
on Nestimate and no work in this list requires modifying it.

## Paper-method completions

- [x] Hayashi et al. (2020): RDC-SymNMF, J-NMF and JS-NMF shipped beside
  RDC-Spec, with direct paper-equation tests.
- [x] Ding et al. (2020): official LDA semantic-hyperedge path shipped.
- [ ] Run full semantic HyperGAT on 20NG, MR and Ohsumed as well as R8/R52.
- [x] Linmei et al. (2019): heterogeneous multi-feature-space HGAT with
  type- and node-level attention.
- [x] Yadati et al. (2019): dynamic, fast and one-edge mediator HyperGCN.
- [x] Dong et al. (2020): HNHN with node/edge normalization exponents.
- [x] Chien et al. (2022): AllDeepSets and AllSetTransformer.

## Statistical and ecosystem completions

- [x] Add `wasserstein_distance()` for persistence diagrams, using the
  SimplicialComplex/GUDHI finite-order and ground-metric conventions and an
  exact native assignment solver.
- [x] Add random hypergraph generators: Bernoulli-incidence G(n,p), SBM,
  k-uniform and k-regular, with `hg_*` and `hypergraph_*` names.
- [x] Add an independent XGI oracle runner for clique, Z and H centralities;
  XGI remains a local-validation dependency, never a runtime dependency.
- [x] Run a ten-seed UMAP/HDBSCAN baseline on R8, including a matched
  eight-cluster control, paired effects, and bootstrap confidence intervals.
- [ ] Run the promised benchmark through the actual Python `BERTopic` class;
  do not label the modular UMAP/HDBSCAN baseline as BERTopic.
- [x] Add worst-community weighted conductance to hypergraph
  community-quality output.
- [x] Ship one dedicated embedding wrapper under both `hg_embed()` and
  `hypergraph_embed()`, delegating to the existing spectral/NMF engine.

## Legal Hypergraphs reproduction

- [x] Implement temporal snapshots, hyperedge distributions,
  s-centralities, subhypergraph centrality, association projections, Y/T/O
  motifs, configuration nulls, Infomap ensembles and partition quality.
- [x] Reproduce the ICSID aggregate motif census: Y = 478, T = 7, O = 0.
- [x] Reproduce the full GFCC association/community analysis and the
  paper's temporal centrality figures from Zenodo 8081507. The runner checks
  all four association graphs, all eight archived AMI medoids, Figure 8
  similarities and quality scores, and the Figure 6 ICSID rankings/values.

## Documentation and release

- [ ] Keep `papers/README.md` as the complete method bibliography, including
  sources not stored locally.
- [ ] Keep every `repos/*.md` note synchronized with shipped code and the
  current upstream release.
- [ ] CRAN submission and flagship honets paper.
