# honets text and hypergraph backlog

This is the maintained backlog after the texthypergraph merge.
Historical Nestimate-era tasks are recorded in `ROADMAP-text.md`; honets
does not depend on Nestimate and no work in this list requires modifying
it.

## Paper-method completions

Hayashi et al. (2020): RDC-SymNMF, J-NMF and JS-NMF shipped beside
RDC-Spec, with direct paper-equation tests.

Ding et al. (2020): official LDA semantic-hyperedge path shipped.

Run full semantic HyperGAT on 20NG, MR and Ohsumed as well as R8/R52.

Linmei et al. (2019): heterogeneous multi-feature-space HGAT with type-
and node-level attention.

Yadati et al. (2019): dynamic, fast and one-edge mediator HyperGCN.

Dong et al. (2020): HNHN with node/edge normalization exponents.

Chien et al. (2022): AllDeepSets and AllSetTransformer.

## Statistical and ecosystem completions

Add
[`wasserstein_distance()`](https://mohsaqr.github.io/hypernets/reference/wasserstein_distance.md)
for persistence diagrams, using the SimplicialComplex/GUDHI finite-order
and ground-metric conventions and an exact native assignment solver.

Add random hypergraph generators: Bernoulli-incidence G(n,p), SBM,
k-uniform and k-regular, with `hg_*` and `hypergraph_*` names.

Add an independent XGI oracle runner for clique, Z and H centralities;
XGI remains a local-validation dependency, never a runtime dependency.

Run a ten-seed UMAP/HDBSCAN baseline on R8, including a matched
eight-cluster control, paired effects, and bootstrap confidence
intervals.

Run the promised benchmark through the actual Python `BERTopic` class
(0.17.4, `all-MiniLM-L6-v2`, ten seeds, three variants, two honets arms;
`benchmarks/run_bertopic_benchmark.R`, results in
`benchmarks/RESULTS.md` and the benchmarks article). The modular
UMAP/HDBSCAN baseline keeps its own label.

Add worst-community weighted conductance to hypergraph community-quality
output.

Ship one dedicated embedding wrapper under both
[`hg_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md)
and
[`hypergraph_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md),
delegating to the existing spectral/NMF engine.

## Legal Hypergraphs reproduction

Implement temporal snapshots, hyperedge distributions, s-centralities,
subhypergraph centrality, association projections, Y/T/O motifs,
configuration nulls, Infomap ensembles and partition quality.

Reproduce the ICSID aggregate motif census: Y = 478, T = 7, O = 0.

Reproduce the full GFCC association/community analysis and the paper’s
temporal centrality figures from Zenodo 8081507. The runner checks all
four association graphs, all eight archived AMI medoids, Figure 8
similarities and quality scores, and the Figure 6 ICSID rankings/values.

## Documentation and release

Keep `papers/README.md` as the complete method bibliography, including
sources not stored locally.

Keep every `repos/*.md` note synchronized with shipped code and the
current upstream release.

CRAN submission and flagship honets paper.
