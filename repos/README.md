# Package and repository map

Maintained 2026-09-02. Runtime dependencies are deliberately separated from
equivalence oracles and research references. Oracles are local-test tools,
never hypernets runtime dependencies.

## Runtime and peer packages

| Package | Relationship to hypernets | Note |
|---|---|---|
| cograph | Imported graph/plot engine after a hypergraph is projected; no export collisions | [`cograph.md`](cograph.md) |
| Dynet | Separate temporal-network peer; never a dependency | [`dynet.md`](dynet.md) |
| sbert | Suggested native embedding frontend for kNN text hypergraphs | [`sbert.md`](sbert.md) |
| Nestimate | Historical source and identity oracle only; hypernets has no runtime dependency | [`nestimate-hypergraph-module.md`](nestimate-hypergraph-module.md) |

## Equivalence oracles

| Repository/package | Oracle role | Current status |
|---|---|---|
| HyperNetX | Zhou/Hayashi Laplacians, EDVW transition/PageRank, s-line graph | Shipped methods have local parity tests; note needs no “parked” interpretation |
| XGI | Clique/Z/H centralities and measure naming | Independent Python equivalence runner shipped for all three centralities |
| HyperG | Unweighted kNN/dual constructions and random generators | kNN/dual plus G(n,p), SBM, uniform and regular models shipped; three samplers have exact seeded parity |
| SimplicialComplex | Persistence diagrams, landscapes and diagram distances | Bottleneck and Wasserstein shipped |
| GUDHI / ripser.py / giotto-tda | Independent persistence and Wasserstein conventions | Wasserstein conventions implemented; SciPy assignment oracle shipped |
| pathpy / pyHON / HYPA / HONEM | Memory-network construction, MOGen and anomaly methods | Integrated in `local_testing_and_equivalence/` |
| HypergraphX | Legal motifs, s-centralities, temporal and community conventions | Legal formula fixtures shipped; broader oracle coverage remains useful |
| Legal Hypergraphs archive | Published GFCC/ICSID workflow | Full GFCC Figure 8 and ICSID Figures 6--7 reproduction passes against Zenodo 8081507 |

## Neural references and baselines

| Repository/package | Role | Current status |
|---|---|---|
| HyperGAT_TextClassification | Official Ding et al. implementation | Sentence and LDA semantic paths shipped; 20NG/MR/Ohsumed runs open |
| DHG / DeepHypergraph | Official-adjacent HGNN and packaged neural-model zoo | HGNN parity plus native HyperGCN/HNHN shipped |
| AllSet | Official AllDeepSets/AllSetTransformer implementation | Both native models shipped with multiset-invariance tests |
| BERTopic | Practical embedding/topic baseline | Actual-package multi-seed benchmark remains open; UMAP/HDBSCAN baseline is reported separately |
| text | Historical reticulate/Hugging Face route | Superseded by sbert for the native pipeline |

Detailed notes: [`hypernetx.md`](hypernetx.md), [`xgi.md`](xgi.md),
[`hyperg-r.md`](hyperg-r.md), [`simplicialcomplex-r.md`](simplicialcomplex-r.md),
[`tda-python.md`](tda-python.md), [`pathpy-hon-python.md`](pathpy-hon-python.md),
[`hypergraphx.md`](hypergraphx.md), [`legal-hypergraphs.md`](legal-hypergraphs.md),
[`hypergat-textclassification.md`](hypergat-textclassification.md),
[`deephypergraph.md`](deephypergraph.md), [`allset.md`](allset.md),
[`bertopic.md`](bertopic.md), and [`text-r.md`](text-r.md).
