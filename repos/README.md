# Package and repository map

Maintained 2026-09-02; synced with the code 2026-09-29 (every "In use"
below was re-checked by grepping `tests/` and `local_testing_and_equivalence/`
for a call to that tool; upstream versions re-checked on PyPI/CRAN/GitHub). Runtime dependencies are deliberately separated from
equivalence oracles and research references. Oracles are local-test tools,
never hypernets runtime dependencies.

## Runtime and peer packages

| Package | Relationship to hypernets | Note |
|---|---|---|
| cograph | Imported graph/plot engine after a hypergraph is projected; no export collisions | [`cograph.md`](cograph.md) |
| Dynet | Separate temporal-network peer; never a dependency (installed 0.4.10; source `../Dyna` and r-universe 0.5.1; not on CRAN) | [`dynet.md`](dynet.md) |
| sbert | Suggested native embedding frontend for kNN text hypergraphs | [`sbert.md`](sbert.md) |
| Nestimate | Historical source and identity oracle; no runtime dependency in DESCRIPTION as of 2026-09-29. ROADMAP.md Phase 0b (decided 2026-09-29) plans `Imports: Nestimate (>= 0.8.5)` | [`nestimate-hypergraph-module.md`](nestimate-hypergraph-module.md) |

## Equivalence oracles

| Repository/package | Oracle role | Current status |
|---|---|---|
| HyperNetX 2.4.3 | Zhou/Hayashi Laplacians, EDVW transition/PageRank, s-line graph | Shipped methods have local parity tests; note needs no “parked” interpretation |
| XGI 0.10.2 | Clique/Z/H centralities and measure naming; planned Hodge/simpliciality/HIF oracle | Independent Python equivalence runner shipped for all three centralities |
| HyperG | Unweighted kNN/dual constructions and random generators | kNN/dual plus G(n,p), SBM, uniform and regular models shipped; three samplers have exact seeded parity |
| SimplicialComplex | Candidate oracle for diagram distances and landscapes | **Not used yet**; bottleneck and Wasserstein ship without it |
| GUDHI / ripser.py / giotto-tda | Planned oracle for Betti numbers, bottleneck, Wasserstein, landscapes (ROADMAP.md Phase 1) | **Not used yet** (GUDHI 3.13.0 current) |
| TDAstats (R, ripser) | Vietoris-Rips persistent homology | In use: `test-equiv-persistent_homology.R` |
| SciPy | `linear_sum_assignment` for the Wasserstein matching | In use: `test-equiv-wasserstein-scipy.R` |
| pathpy 2.2.0 / pyHON / pyMOGen | Memory-network construction, MOGen, HON centralities | In use in `local_testing_and_equivalence/`; pathpyG (v0.2.0-alpha, GitHub only) not used |
| HYPA reference code | `build_hypa()` | **Not used yet** (ROADMAP.md Phase 1) |
| HONEM (Saebi et al.) | `build_honem()` | No external code: checked only against a clean-room R reimplementation (`test-equiv-honem.R`) |
| BiasedUrn (R) + Monte Carlo | Wallenius null and configuration-model semantics for `build_hypa()` | In use: `test-equiv-hypa.R` |
| markovchain (R) | `verifyMarkovProperty()` chi-square for `markov_order_test()` | In use: `test-equiv-markov_order.R` |
| infomap 2.15.1 | Map-equation flow and codelength for `hon_communities()` | In use: `test-equiv-hon-communities-infomap.R` ([`infomap.md`](infomap.md)) |
| scikit-learn 1.8.0 | ARI / AMI / NMI for `hg_agreement()` | In use: frozen reference values in `tests/testthat/test-agreement-sklearn.R` |
| HypergraphX 1.8.0 | Legal motifs, s-centralities, temporal and community conventions | Not a live test: used via frozen formula fixtures and the Zenodo notebooks; `compute_motifs` planned |
| Legal Hypergraphs archive | Published GFCC/ICSID workflow | Full GFCC Figure 8 and ICSID Figures 6--7 reproduction passes against Zenodo 8081507 |

## Neural references and baselines

| Repository/package | Role | Current status |
|---|---|---|
| HyperGAT_TextClassification | Official Ding et al. implementation | Sentence and LDA semantic paths shipped; 20NG/MR/Ohsumed runs open |
| DHG / DeepHypergraph | Official-adjacent HGNN and packaged neural-model zoo | HGNN parity plus native HyperGCN/HNHN shipped |
| AllSet | Official AllDeepSets/AllSetTransformer implementation | Both native models shipped with multiset-invariance tests |
| BERTopic 0.17.4 | Practical embedding/topic baseline | Actual-package ten-seed R8 benchmark run: `benchmarks/RESULTS.md`; UMAP/HDBSCAN baseline reported separately |
| text 1.9 | Historical reticulate/Hugging Face route | Superseded by sbert for the native pipeline |

Detailed notes: [`hypernetx.md`](hypernetx.md), [`xgi.md`](xgi.md),
[`hyperg-r.md`](hyperg-r.md), [`simplicialcomplex-r.md`](simplicialcomplex-r.md),
[`tda-python.md`](tda-python.md), [`pathpy-hon-python.md`](pathpy-hon-python.md),
[`hypergraphx.md`](hypergraphx.md), [`legal-hypergraphs.md`](legal-hypergraphs.md),
[`hypergat-textclassification.md`](hypergat-textclassification.md),
[`deephypergraph.md`](deephypergraph.md), [`allset.md`](allset.md),
[`bertopic.md`](bertopic.md), [`text-r.md`](text-r.md), and [`infomap.md`](infomap.md).
