# hypernets

Higher-order network analysis in R.

A network is *higher-order* when a relation reaches beyond a single pair of
nodes at a single moment — either it binds **more than two nodes at once**, or
it depends on **more than the current node**. The literature (Battiston et al.
2020; Bianconi 2021) organises that idea into three structure families.
hypernets implements all three behind one taxonomy, plus a text family that
reads a corpus as a hypergraph.

| Family | A relation is… | Built from |
|---|---|---|
| **Memory networks** | a state *plus the memory of how it was reached* | categorical sequences |
| **Simplicial complexes** | a set of nodes that are *all* mutually related, plus every subset | a network, a distance matrix, or a memory network |
| **Hypergraphs** | an arbitrary set of nodes bound as a unit | network cliques, sliding windows over sequences, group membership, or embedding neighbours |
| **Text hypergraphs** | a document bound to the words it contains (or a word to its documents) | a corpus: bag of words, token windows, or sentence embeddings |

Base R, ggplot2 and Matrix (with RSpectra for the sparse spectral paths).
No Python, no NLP dependencies; the neural tier needs `torch`, the sentence
embeddings `sbert`, both optional.

## Verbs

### Memory networks

| Verb | Method | Reference |
|---|---|---|
| `hon()` | Higher-order network with rule extraction (BuildHON+) | Xu, Wickramarathne & Chawla (2016) |
| `honem()` | Higher-order network embedding | Saebi, Ciampaglia, Kaplan & Chawla (2020) |
| `hypa()` | Hypergeometric path anomaly detection | LaRock et al. (2020) |
| `mogen()` | Multi-order generative model | Scholtes (2017) |
| `hg_centrality()` | PageRank / betweenness / closeness on the higher-order topology, projected to states | Scholtes, Wider & Garas (2016) |
| `hg_bootstrap()` | Bootstrap CIs + rule support | Efron & Tibshirani (1993) |
| `hg_compare()` | Two-sample permutation comparison of rules | Good (2005) |
| `markov_order()` | Permutation-based Markov order test | Anderson & Goodman (1957) |
| `hg_markov_stability()` | Persistence, stationary flow and first-passage times per state | Kemeny & Snell (1976) |
| `memory()` | Per-context order-k vs order-1 KL diagnostic | Cover & Thomas (2006) |

### Simplicial complexes

| Verb | Method | Reference |
|---|---|---|
| `simplicial()` | Clique, Vietoris-Rips, or pathway complex, verified against igraph's cliques and Euler-Poincaré on construction | Hatcher (2002) |
| `hg_betti()`, `hg_euler()` | Homology ranks and the Euler characteristic | Hatcher (2002) |
| `hg_homology()` | Betti curves and persistence diagrams over a filtration | Edelsbrunner & Harer (2010) |
| `hg_landscape()`, `hg_bottleneck()`, `hg_wasserstein()` | Comparable summaries of persistence diagrams | Bubenik (2015); Cohen-Steiner et al. (2007); Kerber et al. (2017) |
| `hg_degree()`, `hg_qanalysis()` | Higher-order degree and q-connectivity | Atkin (1974) |

### Hypergraphs

| Verb | Method | Reference |
|---|---|---|
| `network_hypergraph()` | Hyperedges from a network's cliques | Burgio et al. (2020) |
| `window_hypergraph()` | Hyperedges from sliding windows over sequences | Ding et al. (2020) |
| `group_hypergraph()` | Hyperedges from bipartite group membership | Perc et al. (2013) |
| `hg_centrality()` | Tensor eigenvector centralities (clique, Z, H) and EDVW PageRank | Benson (2019); Chitra & Raphael (2019) |
| `hg_measures()` | Structural measures (hyperdegree, overlap, density) | — |
| `hg_laplacian()`, `hg_cluster()`, `hg_joint_cluster()`, `hg_classify()` | Weighted normalised Laplacian; RDC-Spec, RDC-SymNMF, J-NMF and JS-NMF clustering; label spreading | Zhou, Huang & Schölkopf (2006); Hayashi et al. (2020) |
| `hg_pagerank()` | EDVW PageRank with personalization, sparse-capable | Chitra & Raphael (2019); Page et al. (1999) |
| `knn_hypergraph()`, `dual_hypergraph()` | Embedding nearest-neighbour hyperedges; vertex/hyperedge role swap | — |
| `hg_clique_expansion()`, `hg_project()`, `hg_line_graph()` | Projections: weighted pairwise network, association-weighted graph, s-line graph | [Coupette, Hartung & Katz (2024)](https://doi.org/10.1098/rsta.2023.0141); Aksoy et al. (2020) |
| `temporal_hypergraph()`, `hg_snapshot()`, `hg_snapshots()` | Interval and contact hypergraphs in Dynet's vocabulary (`actor`/`group`, alias detection, calendar clocks, observation bounds); active and cumulative snapshots on a `step`/`window` grid | Coupette, Hartung & Katz (2024) |
| `hg_edges()`, `hg_edge_centrality()` | Hyperedge distributions and s-betweenness/s-closeness | Coupette, Hartung & Katz (2024); Aksoy et al. (2020) |
| `hg_motifs()` | Induced Y/T/O census and configuration-model null profile | Coupette, Hartung & Katz (2024) |
| `hg_communities()`, `hg_community_quality()` | Repeated Infomap, AMI-medoid selection, coverage/performance/modularity/conductance (on a memory network: the map equation for memory networks) | Coupette, Hartung & Katz (2024); Rosvall et al. (2014) |
| `hg_embed()` | Spectral or symmetric-NMF node coordinates | Zhou et al. (2006); Hayashi et al. (2020) |
| `hg_null_test()` | Degree-preserving swap and configuration-model nulls | Chodrow (2020) |
| `hg_neural()`, `text_hypergat()`, `heterogeneous_hgat()`, `hg_hypergcn()`, `hg_hnhn()`, `hg_allset()` | HGNN, HyperGAT, heterogeneous HGAT, HyperGCN, HNHN, AllDeepSets and AllSetTransformer in native torch | Feng et al. (2019); Linmei et al. (2019); Yadati et al. (2019); Ding et al. (2020); Dong et al. (2020); Chien et al. (2022) |

### Text hypergraphs

| Verb | Method | Reference |
|---|---|---|
| `text_hypergraph()` | Corpus to weighted document-word hypergraph: bag of words with smoothed tf-idf, token windows, or embedding kNN | Hayashi et al. (2020); Ding et al. (2020); Manning, Raghavan & Schütze (2008) |
| `hg_measures()`, `hg_centrality()` | Structural measures and tensor centralities as tidy tables | Benson (2019) |
| `hg_cluster()`, `hg_keywords()` | Spectral or symmetric-NMF topic clustering and per-cluster keywords | Zhou, Huang & Schölkopf (2006); Hayashi et al. (2020) |
| `hg_stability()`, `hg_agreement()`, `hg_seeds()` | Subsampling cluster stability (per-cluster Jaccard) with the eigengap, partition agreement (ARI/AMI/NMI), seed selection | Hennig (2007); von Luxburg (2007); Hubert & Arabie (1985); Vinh et al. (2010) |
| `hg_topic_quality()` | Topic coherence (UMass, NPMI) and FREX exclusivity, as in stm; scores external topic word lists too | Mimno et al. (2011); Lau, Newman & Baldwin (2014); Roberts, Stewart & Tingley (2019) |
| `hg_cocluster()` | Spectral co-clustering of documents and words together | Dhillon (2001) |
| `hg_classify()` | Few-label transductive classification, with class-mass normalization | Zhou et al. (2006); Zhu, Ghahramani & Lafferty (2003) |

A `text_hypergraph` *is* a `net_hg`, so every hypergraph verb takes
it. Constructors are nouns (`hon()`, `mogen()`, `simplicial()`,
`<source>_hypergraph()`); every other verb is `hg_*()`.
Bundled data: `covid_abstracts` (165 abstracts)
and `covid_embeddings`.

### Research sources

The source repository includes an indexed [research paper library](papers/)
and [method/repository reference notes](repos/). For the projection tier, see
the local [*Legal hypergraphs* paper](papers/2024-PhilTrans-LegalHypergraphs-Coupette.pdf)
and its [method mapping and reproducibility links](repos/legal-hypergraphs.md).

## The taxonomy

The same naming rules hold in every family, so a verb from one reads like a
verb from another:

- **Constructors are nouns**: `hon()`, `honem()`, `mogen()`, `simplicial()`,
  `markov_order()`, `memory()`, `text_hypergraph()`, `window_hypergraph()`, ... Where a
  family admits several construction routes, they are selected with `type =`.
- **Everything else is `hg_*()`**. Measures (`hg_centrality()`,
  `hg_measures()`, `hg_degree()`) return a tidy `data.frame`, one row per node
  or per structure; inference (`hg_bootstrap()`, `hg_compare()`,
  `hg_null_test()`) returns a result carrying estimates, intervals and
  p-values. One verb names one idea and dispatches on its input:
  `hg_centrality()`, `hg_communities()` and `hypa()` take a memory network
  or a hypergraph.
- **`hg_get()` is the one reader.** Every object hands over its tables through
  `hg_get(x, what = )` — reaching into a result with `$` is never required.
  `what =` selects a secondary table; filters, `sort_by` and `top` are named
  arguments:

```r
hg_get(net)                                # the rules of a memory network
hg_get(net, what = "nodes")                # its higher-order nodes
hg_get(net, what = "pathways")             # its higher-order pathways
hg_get(mg,  what = "paths", k = 3)         # path counts of a MOGen model
hg_get(ph,  what = "betti")                # Betti curves, not the diagram
hg_get(hg,  what = "memberships")          # one row per node in a hyperedge
```

Every sequence-taking verb reads a long event table (`action =`, `actor =`,
`time =`, `session =`), a wide data frame, a list of sequences, or a network
model object carrying its sequences.

## Crossing between families

The families are entry points into one another, not islands. The same coded
sessions can be read all three ways:

```r
library(hypernets)

# memory network: what follows what, given how you got here
net <- hon(human_long, action = "code", actor = "session_id",
           time = "timestamp", max_order = 3)
hg_centrality(net, sort_by = "pagerank")

# simplicial complex: the same sessions as a pathway complex
sc <- simplicial(net, type = "pathway")
hg_betti(sc)

# hypergraph: each session as one multi-way interaction over its codes
hg <- group_hypergraph(human_long, actor = "code", group = "session_id")
hg_centrality(hg, type = "pagerank")

# text: a corpus as a document-word hypergraph, clustered into topics
thg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
                       weight = "tfidf", stop_words = stop_words_en(),
                       min_count = 3L)
topics <- hg_cluster(thg, k = 4, type = "random_walk", seed = 1)
hg_keywords(thg, topics, n = 5, collapse = TRUE)                    # tf-idf mass
hg_keywords(thg, topics, n = 5, type = "ctfidf", collapse = TRUE)  # BERTopic's c-TF-IDF
hg_keywords(thg, topics, n = 5, type = "centrality", collapse = TRUE)
```

`hg_clique_expansion()` projects a hypergraph back to a pairwise network that any
first-order tool accepts; `hg_get(net, what = "pathways")` hands
sequence-derived path strings from the memory family to the other two.

Bundled data: `human_long`, `ai_long` (coded human-AI pair-programming
sessions); `covid_abstracts`, `covid_embeddings` (COVID-19 education
research abstracts and their sentence embeddings).

## Vignettes

The package ships `vignette("text-hypergraphs")` and
`vignette("text-constructions")`, plus a benchmark article
(`vignettes/articles/benchmarks.Rmd`, R8/R52/MR/Ohsumed/20NG). Family-level
worked analyses, including the Legal hypergraphs paper workflow, are
hand-knit reports kept outside the repository. Longer tutorials live in
`Tutorial_docs/`.

## Provenance

hypernets is the home package for higher-order structure in its family of
network packages, alongside
[psychnets](https://github.com/mohsaqr/psychnets) (cross-sectional psychometric
networks) and [idiographic](https://github.com/mohsaqr/idiographic)
(person-specific temporal models). The memory-network and simplicial-complex
estimators are imported from a sibling estimation package and returned
unchanged under hypernets' names, so every number is that estimator's own;
identity is proven by exact `identical()` tests. The text family arrived in
0.3.0 (2026-09-01) from the retired `texthypergraph` package. See `NEWS.md`
for every rename.

The inherited external-equivalence suites (pyHON, pathpy multi-order models,
Wallenius / Monte-Carlo HYPA references, igraph, HyperNetX) live in
`local_testing_and_equivalence/`, which is build-ignored and gated:

```sh
NOT_CRAN=true HYPERNETS_EQUIV_TESTS=true Rscript -e \
  'library(hypernets); testthat::test_dir("local_testing_and_equivalence")'
```

## Installation

```r
# development version
devtools::install_github("mohsaqr/hypernets")
```

## References

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M., Patania, A.,
Young, J.-G., & Petri, G. (2020). Networks beyond pairwise interactions:
Structure and dynamics. *Physics Reports*, 874, 1-92.

Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University Press.
