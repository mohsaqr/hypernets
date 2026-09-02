# Research paper library

Local copies of papers used to design or verify honets. This directory is
kept in the source repository for research and reproducibility, but excluded
from the built R package by `.Rbuildignore`.

| Year | Paper / topic | Local copy |
|---|---|---|
| 2006 | Learning with hypergraphs (Zhou et al.) | [`2006-NeurIPS-LearningWithHypergraphs-Zhou.pdf`](2006-NeurIPS-LearningWithHypergraphs-Zhou.pdf) |
| 2019 | HGNN (Feng et al.) | [`2019-AAAI-HGNN-Feng.pdf`](2019-AAAI-HGNN-Feng.pdf) |
| 2019 | Hypergraph attention networks (Linmei et al.) | [`2019-EMNLP-HGAT-Linmei.pdf`](2019-EMNLP-HGAT-Linmei.pdf) |
| 2019 | Edge-dependent vertex-weight random walks (Chitra & Raphael) | [`2019-ICML-HypergraphRW-Chitra.pdf`](2019-ICML-HypergraphRW-Chitra.pdf) |
| 2019 | HyperGCN (Yadati et al.) | [`2019-NeurIPS-HyperGCN-Yadati.pdf`](2019-NeurIPS-HyperGCN-Yadati.pdf) |
| 2019 | XLNet (Yang et al.) | [`2019-NeurIPS-XLNet-Yang.pdf`](2019-NeurIPS-XLNet-Yang.pdf) |
| 2020 | Hypergraph clustering (Hayashi et al.) | [`2020-CIKM-HypergraphClustering-Hayashi.pdf`](2020-CIKM-HypergraphClustering-Hayashi.pdf) |
| 2020 | HyperGAT (Ding et al.) | [`2020-EMNLP-HyperGAT-Ding.pdf`](2020-EMNLP-HyperGAT-Ding.pdf) |
| 2020 | HNHN (Dong et al.) | [`2020-HNHN-Dong.pdf`](2020-HNHN-Dong.pdf) |
| 2022 | AllSet (Chien et al.) | [`2022-ICLR-AllSet-Chien.pdf`](2022-ICLR-AllSet-Chien.pdf) |
| 2022 | Text-classification survey (Li et al.) | [`2022-TIST-TextClassificationSurvey-Li.pdf`](2022-TIST-TextClassificationSurvey-Li.pdf) |
| 2024 | **Legal hypergraphs** (Coupette, Hartung & Katz) | [`2024-PhilTrans-LegalHypergraphs-Coupette.pdf`](2024-PhilTrans-LegalHypergraphs-Coupette.pdf) |

## Hypergraph random walks, Laplacians, and clustering

Hayashi, Aksoy, Park and Park (2020), CIKM, DOI
[10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034).

- `hypergraph_laplacian(type = "random_walk")`: representative-digraph
  transition, stationary distribution, and Chung normalized Laplacian.
- `hypergraph_cluster(algorithm = "spectral")`: RDC-Spec (Algorithm 1).
- `hypergraph_cluster(algorithm = "symnmf")`: RDC-Sym (Algorithm 2 and
  Eq. 16).
- `hypergraph_joint_cluster(method = "joint")`: J-NMF (Eq. 18), requiring
  an auxiliary vertex relation such as the patent citation matrix.
- `hypergraph_joint_cluster(method = "joint_symmetric")`: JS-NMF (Eq. 19).
- `hypergraph_transduction()`: the Zhou et al. semi-supervised method used
  in the surrounding normalized-hypergraph framework.

The local tests use direct transition/Laplacian formula fixtures,
non-negative-factor and assignment invariants, monotone RDC-Sym objectives,
and explicit evaluations of Eqs. 16, 18, and 19. HyperNetX remains the
author-adjacent external oracle for the shared RDC-Spec pipeline; it does not
implement the NMF branches.

## HyperGAT

Ding, Wang, Li, Li and Liu (2020), EMNLP, DOI
[10.18653/v1/2020.emnlp-main.399](https://doi.org/10.18653/v1/2020.emnlp-main.399).

- `hg_hypergat(semantic = "none")`: document-level word nodes, sentence
  hyperedges, two dual-attention layers, masked pooling and classification;
  this is the paper's “w/o semantic” ablation.
- `hg_hypergat(semantic = "lda")`: full semantic construction. LDA is fitted
  only to labeled training documents; the topic count defaults to the class
  count; the ten highest-probability words per topic define semantic edges
  inside each document.
- `lda_keywords =`: bypasses fitting with a saved topic-keyword list, which
  is the exact boundary used by the official `generate_lda.py` and
  `utils.py::get_slice()` pipeline.

The attention layer has float32 forward parity with the official PyTorch
implementation. Package tests independently check its equations, the exact
semantic-edge layout, deterministic LDA fitting, training-only input scope,
and both fitted and precomputed public paths.

## Remaining neural hypergraph models

- Linmei et al. (2019): `heterogeneous_hgat()` keeps independent document,
  topic and entity feature spaces and implements Eqs. 2--7, including both
  type-level and node-level attention. It is deliberately distinct from
  HyperGAT.
- Yadati et al. (2019): `hypergraph_hypergcn()` implements Algorithms 1--3.
  Tests pin the farthest-pair choice, mediator edge set, `2|e|-3` weights,
  self-loops and symmetric normalization.
- Dong et al. (2020): `hypergraph_hnhn()` implements Algorithm 1 and the
  Section 2.3 alpha/beta degree normalizers in both propagation directions.
- Chien et al. (2022): `hypergraph_allset()` implements AllDeepSets (Eq. 7)
  and AllSetTransformer (Eq. 8), with independent node-to-edge and
  edge-to-node multiset learners. Descriptive wrappers select either model.

## Legal hypergraphs

Corinna Coupette, Dirk Hartung and Daniel Martin Katz (2024), “Legal
hypergraphs”, *Philosophical Transactions of the Royal Society A* 382(2270),
20230141. DOI: [10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141).

- Open-access record: [PubMed Central PMC10894694](https://pmc.ncbi.nlm.nih.gov/articles/PMC10894694/)
- Reproducibility archive: [Zenodo 8081507](https://doi.org/10.5281/zenodo.8081507)
- honets method mapping: [`../repos/legal-hypergraphs.md`](../repos/legal-hypergraphs.md)
- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Method bibliography without a local PDF

These sources directly define exported honets methods but are indexed by DOI
or bibliographic record rather than copied into this repository.

### Memory networks

- Xu, Wickramarathne & Chawla (2016), BuildHON,
  [doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028).
- Saebi et al. (2020), HONEM,
  [doi:10.1089/big.2019.0169](https://doi.org/10.1089/big.2019.0169).
- LaRock et al. (2020), HYPA,
  [doi:10.1137/1.9781611976236.52](https://doi.org/10.1137/1.9781611976236.52).
- Scholtes (2017), multi-order generative models,
  [doi:10.1145/3097983.3098145](https://doi.org/10.1145/3097983.3098145).
- Scholtes, Wider & Garas (2016), higher-order centralities,
  [doi:10.1140/epjb/e2016-60663-0](https://doi.org/10.1140/epjb/e2016-60663-0).

### Simplicial topology

- Zomorodian & Carlsson (2005), persistent homology,
  [doi:10.1007/s00454-004-1146-y](https://doi.org/10.1007/s00454-004-1146-y).
- Bubenik (2015), persistence landscapes,
  [doi:10.1111/cgf.12447](https://doi.org/10.1111/cgf.12447).
- Atkin (1974), q-analysis, *Mathematical Structure in Human Affairs*.
- Edelsbrunner & Harer (2010), *Computational Topology: An Introduction*.

### Hypergraphs and statistical decisions

- Battiston et al. (2020), higher-order-network review,
  [doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004).
- Benson (2019), three hypergraph eigenvector centralities,
  [doi:10.1137/18M1203031](https://doi.org/10.1137/18M1203031).
- Aksoy et al. (2020), high-order hypergraph walks,
  [doi:10.1140/epjds/s13688-020-00231-0](https://doi.org/10.1140/epjds/s13688-020-00231-0).
- Chodrow (2020), configuration models,
  [doi:10.1093/comnet/cnaa018](https://doi.org/10.1093/comnet/cnaa018).
- Zhu, Ghahramani & Lafferty (2003), class-mass-normalized label spreading.
- Hubert & Arabie (1985), adjusted Rand index.
- Vinh, Epps & Bailey (2010), adjusted mutual information.
- Gotelli (2000) and Phipson & Smyth (2010), null-model and permutation
  inference conventions.
