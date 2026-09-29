# Research paper library

Local copies of papers used to design or verify hypernets. This directory is
kept in the source repository for research and reproducibility, but excluded
from the built R package by `.Rbuildignore`.

| Year | Paper / topic | Local copy |
|---|---|---|
| 2006 | Learning with hypergraphs (Zhou et al.) | [`2006-NeurIPS-LearningWithHypergraphs-Zhou.pdf`](2006-NeurIPS-LearningWithHypergraphs-Zhou.pdf) |
| 2019 | HGNN (Feng et al.) | [`2019-AAAI-HGNN-Feng.pdf`](2019-AAAI-HGNN-Feng.pdf) |
| 2019 | Hypergraph attention networks (Linmei et al.) | [`2019-EMNLP-HGAT-Linmei.pdf`](2019-EMNLP-HGAT-Linmei.pdf) |
| 2019 | Edge-dependent vertex-weight random walks (Chitra & Raphael) | [`2019-ICML-HypergraphRW-Chitra.pdf`](2019-ICML-HypergraphRW-Chitra.pdf) |
| 2019 | HyperGCN (Yadati et al.) | [`2019-NeurIPS-HyperGCN-Yadati.pdf`](2019-NeurIPS-HyperGCN-Yadati.pdf) |
| 2019 | XLNet (Yang et al.) — **not used by the package** (no reference in `R/`; kept as background) | [`2019-NeurIPS-XLNet-Yang.pdf`](2019-NeurIPS-XLNet-Yang.pdf) |
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
- hypernets method mapping: [`../repos/legal-hypergraphs.md`](../repos/legal-hypergraphs.md)
- License: [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

## Method bibliography without a local PDF

These sources directly define exported hypernets methods but are indexed by DOI
or bibliographic record rather than copied into this repository. Entries
added 2026-09-29 were checked against Crossref / publisher / arXiv / CRAN
records; the file each is cited in follows the entry. **Used, not yet
cited** marks a method used in `R/` without a roxygen `@references` entry.

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
- Saebi, Xu, Kaplan, Ribeiro & Chawla (2020), BuildHON+, *EPJ Data Science*
  9(1), 15,
  [doi:10.1140/epjds/s13688-020-00233-y](https://doi.org/10.1140/epjds/s13688-020-00233-y).
  `R/memory_hon.R`.
- Gote, Casiraghi, Schweitzer & Scholtes (2023), multi-order generative
  models for variable-length paths, *Applied Network Science* 8, 68,
  [doi:10.1007/s41109-023-00596-x](https://doi.org/10.1007/s41109-023-00596-x).
  `R/memory_mogen.R`, `R/hypergraph_null.R` — **the roxygen cites it as
  "Gote & Scholtes, 8, 62": two authors and the article number are wrong.**
- Rosvall, Esquivel, Lancichinetti, West & Lambiotte (2014), memory in network
  flows, *Nature Communications* 5, 4630,
  [doi:10.1038/ncomms5630](https://doi.org/10.1038/ncomms5630).
  `R/memory_communities.R`, `R/memory_markov_stability.R`.
- Rosvall & Bergstrom (2008), the map equation, *PNAS* 105(4), 1118-1123,
  [doi:10.1073/pnas.0706851105](https://doi.org/10.1073/pnas.0706851105).
  `R/memory_communities.R`.
- Edler, Bohlin & Rosvall (2017), Infomap for memory and multilayer networks,
  *Algorithms* 10(4), 112,
  [doi:10.3390/a10040112](https://doi.org/10.3390/a10040112).
  `R/memory_communities.R`.
- Lambiotte & Rosvall (2012), smart teleportation, *Physical Review E* 85(5),
  056107,
  [doi:10.1103/PhysRevE.85.056107](https://doi.org/10.1103/PhysRevE.85.056107).
  `R/memory_communities.R`, `R/memory_markov_stability.R`.
- Scholtes, Wider, Pfitzner, Garas, Tessone & Schweitzer (2014),
  causality-driven slow-down and speed-up of diffusion, *Nature
  Communications* 5, 5024,
  [doi:10.1038/ncomms6024](https://doi.org/10.1038/ncomms6024).
  `R/memory_markov_stability.R`.
- Kemeny & Snell (1976), *Finite Markov Chains*, Springer-Verlag
  (Undergraduate Texts in Mathematics; reprint of the 1960 edition), ISBN
  0-387-90192-2. `R/memory_markov_stability.R`.
- Anderson & Goodman (1957), statistical inference about Markov chains,
  *Annals of Mathematical Statistics* 28(1), 89-110,
  [doi:10.1214/aoms/1177707039](https://doi.org/10.1214/aoms/1177707039).
  `R/memory_markov_order.R`.
- Agresti (1992), exact inference for contingency tables, *Statistical
  Science* 7(1), 131-153,
  [doi:10.1214/ss/1177011454](https://doi.org/10.1214/ss/1177011454).
  `R/memory_markov_order.R`.
- Tong (1975), Markov-chain order by AIC, *Journal of Applied Probability*
  12(3), 488-497, [doi:10.2307/3212863](https://doi.org/10.2307/3212863).
  `R/memory_markov_order.R`.
- Katz (1981), criteria for estimating the order of a Markov chain,
  *Technometrics* 23(3), 243-249 (end page unverified: Crossref gives only the
  start page), [doi:10.2307/1267787](https://doi.org/10.2307/1267787).
  `R/memory_markov_order.R` (its roxygen gives the Taylor & Francis DOI
  10.1080/00401706.1981.10486293, not checked).

### Simplicial topology

- Zomorodian & Carlsson (2005), persistent homology,
  [doi:10.1007/s00454-004-1146-y](https://doi.org/10.1007/s00454-004-1146-y).
- Bubenik (2015), persistence landscapes,
  [doi:10.1111/cgf.12447](https://doi.org/10.1111/cgf.12447).
- Atkin (1974), q-analysis, *Mathematical Structure in Human Affairs*.
- Edelsbrunner & Harer (2010), *Computational Topology: An Introduction*.
- Edelsbrunner, Letscher & Zomorodian (2002), topological persistence and
  simplification, *Discrete & Computational Geometry* 28(4), 511-533,
  [doi:10.1007/s00454-002-2885-2](https://doi.org/10.1007/s00454-002-2885-2).
  The journal version is 2002; the FOCS conference version is 2000.
  `R/simplicial_homology.R` (corrected to 2002 on 2026-09-29); the code
  comment in `R/simplicial.R` still says 2000.
- Kerber, Morozov & Nigmetov (2017), geometry helps to compare persistence
  diagrams, *ACM Journal of Experimental Algorithmics* 22, Article 1.4, 1-20,
  [doi:10.1145/3064175](https://doi.org/10.1145/3064175) (article number from
  a search snippet; the ACM page was unreachable). `R/simplicial_distances.R`.
- Hatcher (2002), *Algebraic Topology*, Cambridge University Press, ISBN
  0-521-79540-0. `R/simplicial.R`, `R/simplicial_homology.R`.

### Hypergraphs and statistical decisions

- Battiston et al. (2020), higher-order-network review,
  [doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004).
- Benson (2019), three hypergraph eigenvector centralities,
  [doi:10.1137/18M1203031](https://doi.org/10.1137/18M1203031).
- Aksoy et al. (2020), high-order hypergraph walks,
  [doi:10.1140/epjds/s13688-020-00231-0](https://doi.org/10.1140/epjds/s13688-020-00231-0).
- Chodrow (2020), configuration models,
  [doi:10.1093/comnet/cnaa018](https://doi.org/10.1093/comnet/cnaa018).
- Bianconi (2021), *Higher-Order Networks*, Cambridge University Press
  (Elements in the Structure and Dynamics of Complex Networks),
  [doi:10.1017/9781108770996](https://doi.org/10.1017/9781108770996).
  `R/hypernets-package.R`.
- Estrada & Rodríguez-Velázquez (2005), complex networks as hypergraphs,
  arXiv:physics/0505137; journal version under a different title, "Subgraph
  centrality and clustering in complex hyper-networks", *Physica A* 364,
  581-594 (2006),
  [doi:10.1016/j.physa.2005.12.002](https://doi.org/10.1016/j.physa.2005.12.002).
  `R/hypergraph_centrality.R` (cites the preprint).
- Chung (2005), Laplacians and the Cheeger inequality for directed graphs,
  *Annals of Combinatorics* 9(1), 1-19,
  [doi:10.1007/s00026-005-0237-z](https://doi.org/10.1007/s00026-005-0237-z).
  `R/hypergraph_laplacian.R`.
- Burgio, Matamalas, Gómez & Arenas (2020), cooperation from networks to
  hypergraphs, *Entropy* 22(7), 744,
  [doi:10.3390/e22070744](https://doi.org/10.3390/e22070744). `R/hypergraph.R`.
- Do, Yoon, Hooi & Shin (2020), structural patterns and generative models of
  real-world hypergraphs, *KDD '20*, 176-186,
  [doi:10.1145/3394486.3403060](https://doi.org/10.1145/3394486.3403060).
  `R/hypergraph_measures.R` (cites the arXiv version).
- **Unverified / probably wrong:** `R/hypergraph_measures.R` cites "Lee, G.,
  Choe, M., & Shin, K. (2024). A survey on hypergraph representation,
  learning and mining. *Data Mining & Knowledge Discovery* 37, 1-39." No such
  paper could be found. Candidates: Lee, Bu, Eliassi-Rad & Shin (2025), "A
  survey on hypergraph mining: patterns, tools, and generators", *ACM
  Computing Surveys* 57(8), 1-36,
  [doi:10.1145/3719002](https://doi.org/10.1145/3719002) (arXiv:2401.08878);
  or Lee, Choe & Shin (2021), "How do hyperedges overlap in real-world
  hypergraphs? Patterns, measures, and generators", *WWW 2021*. The roxygen
  needs correcting.
- Marchette (2021), *HyperG: Hypergraphs in R*, R package 1.0.0 (CRAN
  2021-03-04). `R/hypergraph_random.R`.
- Page, Brin, Motwani & Winograd (1999), the PageRank citation ranking,
  Stanford InfoLab Technical Report 1999-66 (report number from search
  results; the InfoLab server was unreachable). `R/hypergraph_pagerank.R`.
- Fruchterman & Reingold (1991), force-directed placement, *Software:
  Practice and Experience* 21(11), 1129-1164,
  [doi:10.1002/spe.4380211102](https://doi.org/10.1002/spe.4380211102).
  `R/hypergraph_plot.R`.

### Text and neural models

- Grootendorst (2022), BERTopic: neural topic modeling with a class-based
  TF-IDF procedure, arXiv:2203.05794. `R/text_verbs.R` (`"ctfidf"`).
- Roberts, Stewart & Tingley (2019), stm, *Journal of Statistical Software*
  91(2), [doi:10.18637/jss.v091.i02](https://doi.org/10.18637/jss.v091.i02).
  `R/text_topics.R`.
- Hoffman, Blei & Bach (2010), online learning for LDA, *NIPS 23*, 856-864.
  Cited only in a code comment of `R/text_hypergat.R` (the `semantic = "lda"`
  path), not in roxygen.
- Blei, Ng & Jordan (2003), latent Dirichlet allocation, *JMLR* 3, 993-1022.
  **Used, not yet cited**: `hg_hypergat(semantic = "lda")`.
- Zaheer, Kottur, Ravanbakhsh, Póczos, Salakhutdinov & Smola (2017), Deep
  Sets, *NIPS 30*, 3391-3401 (arXiv:1703.06114). **Used, not yet cited**:
  `hypergraph_allset(model = "deepsets")`.
- Lee, Lee, Kim, Kosiorek, Choi & Teh (2019), Set Transformer, *ICML 2019*,
  PMLR 97, 3744-3753. **Used, not yet cited**: the PMA block of
  `hypergraph_allset(model = "transformer")`.
- Kipf & Welling (2017), graph convolutional networks, *ICLR 2017*
  (arXiv:1609.02907). **Used, not yet cited**: the graph-convolution basis of
  `hg_neural()`, `hg_hypergcn()`, `heterogeneous_hgat()`.
- Yao, Mao & Luo (2019), graph convolutional networks for text
  classification (TextGCN), *AAAI* 33(01), 7370-7377,
  [doi:10.1609/aaai.v33i01.33017370](https://doi.org/10.1609/aaai.v33i01.33017370).
  **Used, not yet cited**: the TextGCN corpora, splits and baseline numbers in
  `benchmarks/RESULTS.md` and the benchmarks article.
- Reimers & Gurevych (2019), Sentence-BERT, *EMNLP-IJCNLP 2019*, 3980-3990,
  [doi:10.18653/v1/D19-1410](https://doi.org/10.18653/v1/D19-1410).
  **Used, not yet cited**: `sbert` embeddings for `knn_hypergraph()` /
  `text_hypergraph(construction = "knn")`.
- Not used by the package: Yang et al. (2019), XLNet, *NeurIPS 32*
  (arXiv:1906.08237); the local PDF is kept.
- Zhu, Ghahramani & Lafferty (2003), class-mass-normalized label spreading.
- Hubert & Arabie (1985), adjusted Rand index.
- Vinh, Epps & Bailey (2010), adjusted mutual information.
- Gotelli (2000) and Phipson & Smyth (2010), null-model and permutation
  inference conventions.
