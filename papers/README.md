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

These sources define or support hypernets methods and are indexed by DOI or
bibliographic record rather than copied here. **Re-verified 2026-09-29**
against Crossref (`api.crossref.org/works/<doi>`, DOI content negotiation),
the publisher or proceedings page (ACL Anthology, JMLR, PMLR, NeurIPS
proceedings, the ACM version-of-record PDF), arXiv, or a library record.
Where a detail could only be confirmed from a search-result snippet, or not
at all, the entry says so. The source file citing each work follows it.

Flags: **Used, not cited** = the method is used in `R/` but no roxygen
`@references` names it. **Roxygen wrong** = the `R/` citation disagrees with
the verified record (fix in `R/`, not here).

### Local PDFs above: verified records

- Zhou, Huang & Schölkopf (2006), *NeurIPS 19* (MIT Press). The proceedings
  record gives no page numbers; the "1601-1608" in `R/hypergraph_project.R`
  and `R/text_sequences.R` is **unverified**.
- Feng, You, Zhang, Ji & Gao (2019), *AAAI* 33(01), 3558-3565,
  [doi:10.1609/aaai.v33i01.33013558](https://doi.org/10.1609/aaai.v33i01.33013558).
  `R/text_neural.R`.
- Linmei, Yang, Shi, Ji & Li (2019), *EMNLP-IJCNLP 2019*, 4821-4830 (ACL
  Anthology; Crossref lists 4820-4829),
  [doi:10.18653/v1/D19-1488](https://doi.org/10.18653/v1/D19-1488).
  `R/heterogeneous_hgat.R`.
- Chitra & Raphael (2019), *ICML 2019*, PMLR 97, 1172-1181.
  `R/hypergraph_centrality.R`, `R/hypergraph_laplacian.R`,
  `R/hypergraph_pagerank.R`, `R/hypergraph_window.R`.
- Yadati, Nimishakavi, Yadav, Nitin, Louis & Talukdar (2019), *NeurIPS 32*
  (no page numbers in the proceedings record). `R/hypergraph_hypergcn.R`.
- Hayashi, Aksoy, Park & Park (2020), *CIKM '20*, 495-504,
  [doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034).
- Ding, Wang, Li, Li & Liu (2020), *EMNLP 2020*, 4927-4936,
  [doi:10.18653/v1/2020.emnlp-main.399](https://doi.org/10.18653/v1/2020.emnlp-main.399).
- Dong, Sawin & Bengio (2020), HNHN, Graph Representation Learning and
  Beyond workshop at ICML 2020, arXiv:2006.12278. `R/hypergraph_hnhn.R`.
- Chien, Pan, Peng & Milenkovic (2022), You are AllSet, *ICLR 2022*,
  arXiv:2106.13264. `R/hypergraph_allset.R`.
- Li, Peng, Li, Xia, Yang, Sun, Yu & He (2022), a survey on text
  classification, *ACM TIST* 13(2), 1-41,
  [doi:10.1145/3495162](https://doi.org/10.1145/3495162). Background only;
  not cited in `R/`.
- Yang, Dai, Yang, Carbonell, Salakhutdinov & Le (2019), XLNet, *NeurIPS 32*,
  arXiv:1906.08237. **Not used by the package**: no occurrence of "XLNet" in
  `R/`, `tests/`, `benchmarks/` or `vignettes/` (grep, 2026-09-29). The local
  PDF is background only and could be removed.

### Memory networks

- Xu, Wickramarathne & Chawla (2016), representing higher-order dependencies
  in networks (BuildHON), *Science Advances* 2(5), e1600028,
  [doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028).
  `R/memory_hon.R`, `R/memory_centrality.R`, `R/memory_inference.R`,
  `R/text_sequences.R`.
- Saebi, Xu, Kaplan, Ribeiro & Chawla (2020), BuildHON+, *EPJ Data Science*
  9(1), 15,
  [doi:10.1140/epjds/s13688-020-00233-y](https://doi.org/10.1140/epjds/s13688-020-00233-y).
  `R/memory_hon.R`.
- Saebi, Ciampaglia, Kaplan & Chawla (2020), HONEM, *Big Data* 8(4),
  255-269, [doi:10.1089/big.2019.0169](https://doi.org/10.1089/big.2019.0169).
  `R/memory_honem.R`.
- LaRock, Nanumyan, Scholtes, Casiraghi, Eliassi-Rad & Schweitzer (2020),
  HYPA, *Proc. 2020 SIAM International Conference on Data Mining*, 460-468,
  [doi:10.1137/1.9781611976236.52](https://doi.org/10.1137/1.9781611976236.52).
  `R/memory_hypa.R`, `R/hypergraph_hypa.R`.
- Scholtes (2017), When is a network a network? Multi-order graphical model
  selection, *KDD '17*, 1037-1046,
  [doi:10.1145/3097983.3098145](https://doi.org/10.1145/3097983.3098145).
  `R/memory_mogen.R`, `R/memory_markov_order.R`.
- Gote, Casiraghi, Schweitzer & Scholtes (2023), predicting variable-length
  paths with multi-order generative models, *Applied Network Science* 8, 68,
  [doi:10.1007/s41109-023-00596-x](https://doi.org/10.1007/s41109-023-00596-x).
  `R/memory_mogen.R` (roxygen now correct).
- Scholtes, Wider & Garas (2016), higher-order aggregate networks: path
  structures and centralities, *European Physical Journal B* 89(3), 61,
  [doi:10.1140/epjb/e2016-60663-0](https://doi.org/10.1140/epjb/e2016-60663-0).
  `R/memory_centrality.R`.
- Scholtes, Wider, Pfitzner, Garas, Tessone & Schweitzer (2014),
  causality-driven slow-down and speed-up of diffusion, *Nature
  Communications* 5, 5024,
  [doi:10.1038/ncomms6024](https://doi.org/10.1038/ncomms6024).
  `R/memory_markov_stability.R`.
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
- Blondel, Guillaume, Lambiotte & Lefebvre (2008), fast unfolding of
  communities (Louvain), *J. Stat. Mech.* 2008(10), P10008,
  [doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008).
  **Used, not cited**: `hon_communities()` names "the Louvain machinery"
  (`R/memory_communities.R`).
- Kemeny & Snell (1976), *Finite Markov Chains*, Springer-Verlag
  (Undergraduate Texts in Mathematics; reprint of the 1960 edition), ISBN
  0-387-90192-2. `R/memory_markov_stability.R`.
- Anderson & Goodman (1957), statistical inference about Markov chains,
  *Annals of Mathematical Statistics* 28(1), 89-110,
  [doi:10.1214/aoms/1177707039](https://doi.org/10.1214/aoms/1177707039).
  `R/memory_markov_order.R`.
- Agresti (1992), a survey of exact inference for contingency tables,
  *Statistical Science* 7(1), 131-153,
  [doi:10.1214/ss/1177011454](https://doi.org/10.1214/ss/1177011454) (pages
  from Project Euclid via search; Crossref carries none).
  `R/memory_markov_order.R`.
- Besag & Mondal (2013), exact goodness-of-fit tests for Markov chains,
  *Biometrics* 69(2), 488-496,
  [doi:10.1111/biom.12009](https://doi.org/10.1111/biom.12009).
  `R/memory_markov_order.R`.
- Tong (1975), Markov-chain order by AIC, *Journal of Applied Probability*
  12(3), 488-497, [doi:10.2307/3212863](https://doi.org/10.2307/3212863).
  `R/memory_markov_order.R`.
- Katz (1981), on some criteria for estimating the order of a Markov chain,
  *Technometrics* 23(3), 243-249,
  [doi:10.2307/1267787](https://doi.org/10.2307/1267787). Pages verified from
  the 7-page article PDF in the Technometrics archive at stat.cmu.edu
  (Crossref gives only the start page). **Roxygen wrong**:
  `R/memory_markov_order.R` gives `10.1080/00401706.1981.10486293`, which
  doi.org does not resolve (HTTP 404) and Crossref does not know, though
  Taylor & Francis uses that string in its URL; cite `10.2307/1267787`.
- Efron & Tibshirani (1993), *An Introduction to the Bootstrap*, Chapman &
  Hall, [doi:10.1201/9780429246593](https://doi.org/10.1201/9780429246593)
  (the CRC reissue's DOI). `R/memory_inference.R`.
- Good (2005), *Permutation, Parametric and Bootstrap Tests of Hypotheses*,
  3rd ed., Springer,
  [doi:10.1007/b138696](https://doi.org/10.1007/b138696).
  `R/memory_inference.R`.
- Cover & Thomas (2006), *Elements of Information Theory*, 2nd ed., Wiley,
  [doi:10.1002/047174882X](https://doi.org/10.1002/047174882X).
  `R/memory_path_dependence.R`.

### Simplicial topology

- Zomorodian & Carlsson (2005), computing persistent homology, *Discrete &
  Computational Geometry* 33(2), 249-274,
  [doi:10.1007/s00454-004-1146-y](https://doi.org/10.1007/s00454-004-1146-y).
  Not currently cited in `R/`.
- Edelsbrunner, Letscher & Zomorodian (2002), topological persistence and
  simplification, *Discrete & Computational Geometry* 28(4), 511-533,
  [doi:10.1007/s00454-002-2885-2](https://doi.org/10.1007/s00454-002-2885-2).
  The journal version is 2002; the FOCS conference version is 2000.
  `R/simplicial_homology.R`; the code comment at `R/simplicial.R` (~l. 498)
  still says 2000.
- Edelsbrunner & Harer (2010), *Computational Topology: An Introduction*,
  American Mathematical Society,
  [doi:10.1090/mbk/069](https://doi.org/10.1090/mbk/069).
  `R/simplicial_distances.R`.
- Kerber, Morozov & Nigmetov (2017), geometry helps to compare persistence
  diagrams, *ACM Journal of Experimental Algorithmics* 22(1), Article 1.4,
  1-20, [doi:10.1145/3064175](https://doi.org/10.1145/3064175). Article
  number verified from the ACM version-of-record PDF (author-hosted).
  `R/simplicial_distances.R`.
- Kuhn (1955), the Hungarian method for the assignment problem, *Naval
  Research Logistics Quarterly* 2(1-2), 83-97,
  [doi:10.1002/nav.3800020109](https://doi.org/10.1002/nav.3800020109).
  **Used, not cited**: the native Hungarian solver in
  `wasserstein_distance()` (`R/simplicial_distances.R`).
- Bubenik (2015), statistical topological data analysis using persistence
  landscapes, *Journal of Machine Learning Research* 16(3), 77-102 (JMLR
  page; no DOI). The DOI 10.1111/cgf.12447 formerly listed here was wrong
  (it resolves to an unrelated 2014 *Computer Graphics Forum* paper).
  `R/simplicial_distances.R`.
- Atkin (1974), *Mathematical Structure in Human Affairs*, Heinemann
  Educational, London, ISBN 0-435-82025-7. `R/simplicial.R`.
- Hatcher (2002), *Algebraic Topology*, Cambridge University Press, ISBN
  0-521-79540-0. `R/simplicial.R`, `R/simplicial_homology.R`.

### Hypergraphs and statistical decisions

- Battiston, Cencetti, Iacopini, Latora, Lucas, Patania, Young & Petri
  (2020), networks beyond pairwise interactions, *Physics Reports* 874, 1-92,
  [doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004).
  `R/hypernets-package.R`.
- Bianconi (2021), *Higher-Order Networks*, Cambridge University Press
  (Elements in the Structure and Dynamics of Complex Networks),
  [doi:10.1017/9781108770996](https://doi.org/10.1017/9781108770996).
  `R/hypernets-package.R`.
- Benson (2019), three hypergraph eigenvector centralities, *SIAM Journal on
  Mathematics of Data Science* 1(2), 293-312,
  [doi:10.1137/18M1203031](https://doi.org/10.1137/18M1203031).
  `R/hypergraph_centrality.R`.
- Kolda & Mayo (2011), shifted power method for computing tensor eigenpairs,
  *SIAM Journal on Matrix Analysis and Applications* 32(4), 1095-1124,
  [doi:10.1137/100801482](https://doi.org/10.1137/100801482). Named only in a
  code comment of `R/hypergraph_centrality.R` (the SS-HOPM kernel).
- Chang, Pearson & Zhang (2009), on eigenvalue problems of real symmetric
  tensors, *Journal of Mathematical Analysis and Applications* 350(1),
  416-422,
  [doi:10.1016/j.jmaa.2008.09.067](https://doi.org/10.1016/j.jmaa.2008.09.067).
  Named only in a code comment of `R/hypergraph_centrality.R`. The same
  authors' Perron-Frobenius paper (*Communications in Mathematical Sciences*
  6(2), 507-520, 2008,
  [doi:10.4310/CMS.2008.v6.n2.a12](https://doi.org/10.4310/CMS.2008.v6.n2.a12))
  may be the one meant; the comment gives only "2009".
- Estrada & Rodríguez-Velázquez (2005), complex networks as hypergraphs,
  arXiv:physics/0505137 (v1, 16 pages). Journal version under a different
  title: "Subgraph centrality and clustering in complex hyper-networks",
  *Physica A* 364, 581-594 (2006),
  [doi:10.1016/j.physa.2005.12.002](https://doi.org/10.1016/j.physa.2005.12.002).
  Both records verified (arXiv API, Crossref). `R/hypergraph_centrality.R`
  cites the preprint; citing the Physica A version is preferable.
- Aksoy, Joslyn, Ortiz Marrero, Praggastis & Purvine (2020), hypernetwork
  science via high-order hypergraph walks, *EPJ Data Science* 9(1), 16,
  [doi:10.1140/epjds/s13688-020-00231-0](https://doi.org/10.1140/epjds/s13688-020-00231-0).
  `R/hypergraph_project.R`.
- Chodrow (2020), configuration models of random hypergraphs, *Journal of
  Complex Networks* 8(3), cnaa018,
  [doi:10.1093/comnet/cnaa018](https://doi.org/10.1093/comnet/cnaa018).
  `R/hypergraph_null.R`.
- Gotelli (2000), null model analysis of species co-occurrence patterns,
  *Ecology* 81(9), 2606-2621,
  [doi:10.1890/0012-9658(2000)081\[2606:NMAOSC\]2.0.CO;2](https://doi.org/10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2).
  `R/hypergraph_null.R`.
- Phipson & Smyth (2010), permutation p-values should never be zero,
  *Statistical Applications in Genetics and Molecular Biology* 9(1),
  [doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585).
  `R/hypergraph_null.R`.
- Chung (2005), Laplacians and the Cheeger inequality for directed graphs,
  *Annals of Combinatorics* 9(1), 1-19,
  [doi:10.1007/s00026-005-0237-z](https://doi.org/10.1007/s00026-005-0237-z).
  `R/hypergraph_laplacian.R`.
- Zhu, Ghahramani & Lafferty (2003), semi-supervised learning using Gaussian
  fields and harmonic functions, *ICML 2003*, 912-919 (pages from search
  results citing the proceedings; no DOI). `R/hypergraph_laplacian.R`.
- Page, Brin, Motwani & Winograd (1999), the PageRank citation ranking:
  bringing order to the web, Stanford InfoLab Technical Report 1999-66.
  **Report number unverified**: it appears only in search results; the
  InfoLab server (ilpubs.stanford.edu:8090/422/) returned errors on
  2026-09-29. `R/hypergraph_pagerank.R`.
- Burgio, Matamalas, Gómez & Arenas (2020), cooperation from networks to
  hypergraphs, *Entropy* 22(7), 744,
  [doi:10.3390/e22070744](https://doi.org/10.3390/e22070744). `R/hypergraph.R`.
- Agrawal & Srikant (1994), fast algorithms for mining association rules,
  *Proc. 20th VLDB*, 487-499, Morgan Kaufmann (pages from search results; no
  DOI). `R/hypergraph_groups.R`.
- Perc, Gómez-Gardeñes, Szolnoki, Floría & Moreno (2013), evolutionary
  dynamics of group interactions on structured populations, *Journal of the
  Royal Society Interface* 10(80), 20120997,
  [doi:10.1098/rsif.2012.0997](https://doi.org/10.1098/rsif.2012.0997).
  `R/hypergraph_groups.R`.
- Tian & Zafarani (2024), higher-order networks representation and learning:
  a survey, *ACM SIGKDD Explorations Newsletter* 26(1), 1-18,
  [doi:10.1145/3682112.3682114](https://doi.org/10.1145/3682112.3682114)
  (arXiv:2402.19414; its Section 5.1.5 is "Downgrading a hypergraph to a
  dyadic graph"). **Roxygen wrong** in `R/hypergraph_expansion.R`: first
  author is Hao Tian ("H.", not "Y."), and the title is not "Higher-order
  network analysis methods".
- Do, Yoon, Hooi & Shin (2020), structural patterns and generative models of
  real-world hypergraphs, *KDD '20*, 176-186,
  [doi:10.1145/3394486.3403060](https://doi.org/10.1145/3394486.3403060).
  `R/hypergraph_measures.R` (cites the arXiv version, 2006.07060).
- Lee, Bu, Eliassi-Rad & Shin (2025), a survey on hypergraph mining:
  patterns, tools, and generators, *ACM Computing Surveys* 57(8), 1-36,
  [doi:10.1145/3719002](https://doi.org/10.1145/3719002) (arXiv:2401.08878).
  `R/hypergraph_measures.R` (roxygen now corrected from the earlier
  non-existent "Lee, Choe & Shin 2024"). Its "Article 203" is from search
  results only; the ACM page returned 403.
- Hubert & Arabie (1985), comparing partitions, *Journal of Classification*
  2(1), 193-218,
  [doi:10.1007/BF01908075](https://doi.org/10.1007/BF01908075).
  `R/agreement.R`, `R/memory_communities.R`.
- Vinh, Epps & Bailey (2010), information theoretic measures for clusterings
  comparison, *Journal of Machine Learning Research* 11(95), 2837-2854 (JMLR
  page; no DOI). `R/agreement.R`.
- Marchette (2021), *HyperG: Hypergraphs in R*, R package 1.0.0 (CRAN
  2021-03-04). `R/hypergraph_random.R`.
- Fruchterman & Reingold (1991), force-directed placement, *Software:
  Practice and Experience* 21(11), 1129-1164,
  [doi:10.1002/spe.4380211102](https://doi.org/10.1002/spe.4380211102).
  `R/hypergraph_plot.R`.
- Wasserman & Faust (1994), *Social Network Analysis: Methods and
  Applications*, Cambridge University Press,
  [doi:10.1017/CBO9780511815478](https://doi.org/10.1017/CBO9780511815478).
  **Used, not cited**: the "Wasserman–Faust correction" for disconnected
  line graphs in `R/hypergraph_edge_centrality.R`.

### Text and neural models

- Manning, Raghavan & Schütze (2008), *Introduction to Information
  Retrieval*, Cambridge University Press,
  [doi:10.1017/CBO9780511809071](https://doi.org/10.1017/CBO9780511809071).
  `R/text_hypergraph.R`.
- van Eck & Waltman (2009), how to normalize cooccurrence data?, *JASIST*
  60(8), 1635-1651,
  [doi:10.1002/asi.21075](https://doi.org/10.1002/asi.21075).
  `R/text_verbs.R`.
- Bouma (2009), normalized (pointwise) mutual information in collocation
  extraction, *Proceedings of GSCL* 30, 31-40 (from search results; no DOI).
  `R/text_topics.R`.
- Bezdek (1981), *Pattern Recognition with Fuzzy Objective Function
  Algorithms*, Plenum,
  [doi:10.1007/978-1-4757-0450-1](https://doi.org/10.1007/978-1-4757-0450-1).
  `R/text_topics.R`.
- Roberts, Stewart & Tingley (2019), stm, *Journal of Statistical Software*
  91(2), [doi:10.18637/jss.v091.i02](https://doi.org/10.18637/jss.v091.i02).
  `R/text_topics.R`.
- Grootendorst (2022), BERTopic: neural topic modeling with a class-based
  TF-IDF procedure, arXiv:2203.05794. `R/text_verbs.R` (`"ctfidf"`).
- Blei, Ng & Jordan (2003), latent Dirichlet allocation, *Journal of Machine
  Learning Research* 3, 993-1022 (JMLR page; the MIT Press DOI
  10.1162/jmlr.2003.3.4-5.993 is listed by Crossref as deleted, so none is
  given). **Used, not cited**: `hg_hypergat(semantic = "lda")`
  (`R/text_hypergat.R`).
- Hoffman, Blei & Bach (2010), online learning for latent Dirichlet
  allocation, *NeurIPS 23* (the proceedings record gives no pages and lists
  the authors as Hoffman, Bach, Blei; the "856-864" formerly given here is
  unverified and dropped). Cited only in a code comment of
  `R/text_hypergat.R`, not in roxygen.
- Zaheer, Kottur, Ravanbakhsh, Póczos, Salakhutdinov & Smola (2017), Deep
  Sets, *NeurIPS 30* (arXiv:1703.06114; no pages in the proceedings record).
  **Used, not cited**: `hg_allset()` AllDeepSets (`R/hypergraph_allset.R`).
- Lee, Lee, Kim, Kosiorek, Choi & Teh (2019), Set Transformer, *ICML 2019*,
  PMLR 97, 3744-3753 (arXiv:1810.00825). **Used, not cited**: the PMA block
  of the AllSetTransformer (`R/hypergraph_allset.R`, `.thg_pma_module()`).
- Kipf & Welling (2017), semi-supervised classification with graph
  convolutional networks, *ICLR 2017*, arXiv:1609.02907. **Used, not cited**:
  the graph-convolution basis of `hg_neural()`, `hg_hypergcn()`,
  `heterogeneous_hgat()`.
- Kingma & Ba (2015), Adam: a method for stochastic optimization, *ICLR
  2015*, arXiv:1412.6980. **Used, not cited**: the optimizer of every torch
  model (`R/text_neural.R`, `R/text_hypergat.R`, `R/hypergraph_hypergcn.R`,
  `R/hypergraph_hnhn.R`, `R/hypergraph_allset.R`).
- Srivastava, Hinton, Krizhevsky, Sutskever & Salakhutdinov (2014), dropout,
  *Journal of Machine Learning Research* 15(56), 1929-1958 (JMLR page).
  **Used, not cited**: dropout in the same torch models.
- Yao, Mao & Luo (2019), graph convolutional networks for text
  classification (TextGCN), *AAAI* 33(01), 7370-7377,
  [doi:10.1609/aaai.v33i01.33017370](https://doi.org/10.1609/aaai.v33i01.33017370).
  **Used, not cited**: the TextGCN corpora, splits and baseline numbers in
  `benchmarks/RESULTS.md` and the benchmarks article.
- Reimers & Gurevych (2019), Sentence-BERT, *EMNLP-IJCNLP 2019*, 3982-3992
  (ACL Anthology; Crossref lists 3980-3990),
  [doi:10.18653/v1/D19-1410](https://doi.org/10.18653/v1/D19-1410).
  **Used, not cited**: `sbert` embeddings for `knn_hypergraph()` /
  `text_hypergraph(construction = "knn")`.
