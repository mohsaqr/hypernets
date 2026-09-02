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
- [ ] Add an independent XGI oracle for clique, Z and H centralities.
- [ ] Run the promised multi-seed BERTopic benchmark with effect sizes and
  confidence intervals.
- [ ] Add conductance to hypergraph community-quality output.
- [ ] Decide whether a dedicated `hypergraph_embed()` wrapper improves the
  API over `hypergraph_cluster(..., what = "embedding")`; do not document a
  nonexistent `hg_embed()` in the meantime.

## Legal Hypergraphs reproduction

- [x] Implement temporal snapshots, hyperedge distributions,
  s-centralities, subhypergraph centrality, association projections, Y/T/O
  motifs, configuration nulls, Infomap ensembles and partition quality.
- [x] Reproduce the ICSID aggregate motif census: Y = 478, T = 7, O = 0.
- [ ] Reproduce the full GFCC association/community analysis and the
  paper's temporal centrality figures from the Zenodo archives.

## Documentation and release

- [ ] Keep `papers/README.md` as the complete method bibliography, including
  sources not stored locally.
- [ ] Keep every `repos/*.md` note synchronized with shipped code and the
  current upstream release.
- [ ] CRAN submission and flagship honets paper.
