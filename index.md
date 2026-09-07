# honets

Higher-order network analysis in R.

A network is *higher-order* when a relation reaches beyond a single pair
of nodes at a single moment — either it binds **more than two nodes at
once**, or it depends on **more than the current node**. The literature
(Battiston et al. 2020; Bianconi 2021) organises that idea into three
structure families. honets implements all three behind one taxonomy,
plus a text family that reads a corpus as a hypergraph.

| Family | A relation is… | Built from |
|----|----|----|
| **Memory networks** | a state *plus the memory of how it was reached* | categorical sequences |
| **Simplicial complexes** | a set of nodes that are *all* mutually related, plus every subset | a network, a distance matrix, or a memory network |
| **Hypergraphs** | an arbitrary set of nodes bound as a unit | network cliques, sliding windows over sequences, group membership, or embedding neighbours |
| **Text hypergraphs** | a document bound to the words it contains (or a word to its documents) | a corpus: bag of words, token windows, or sentence embeddings |

Base R, ggplot2 and Matrix (with RSpectra for the sparse spectral
paths). No Python, no NLP dependencies; the neural tier needs `torch`,
the sentence embeddings `sbert`, both optional.

## Verbs

### Memory networks

| Verb | Method | Reference |
|----|----|----|
| [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md) | Higher-order network with rule extraction (BuildHON+) | Xu, Wickramarathne & Chawla (2016) |
| [`build_honem()`](https://mohsaqr.github.io/hypernets/reference/build_honem.md) | Higher-order network embedding | Saebi, Ciampaglia, Kaplan & Chawla (2020) |
| [`build_hypa()`](https://mohsaqr.github.io/hypernets/reference/build_hypa.md) | Hypergeometric path anomaly detection | LaRock et al. (2020) |
| [`build_mogen()`](https://mohsaqr.github.io/hypernets/reference/build_mogen.md) | Multi-order generative model | Scholtes (2017) |
| [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md) | PageRank / betweenness / closeness on the higher-order topology, projected to states | Scholtes, Wider & Garas (2016) |
| [`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md) | Bootstrap CIs + rule support | Efron & Tibshirani (1993) |
| [`compare_hon()`](https://mohsaqr.github.io/hypernets/reference/compare_hon.md) | Two-sample permutation comparison of rules | Good (2005) |
| [`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md) | Permutation-based Markov order test | — |
| [`path_dependence()`](https://mohsaqr.github.io/hypernets/reference/path_dependence.md) | Per-context order-k vs order-1 KL diagnostic | Cover & Thomas (2006) |

### Simplicial complexes

| Verb | Method | Reference |
|----|----|----|
| [`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md) | Clique, Vietoris-Rips, or pathway complex | — |
| [`betti_numbers()`](https://mohsaqr.github.io/hypernets/reference/betti_numbers.md), [`euler_characteristic()`](https://mohsaqr.github.io/hypernets/reference/euler_characteristic.md) | Homology ranks and the Euler characteristic | — |
| [`persistent_homology()`](https://mohsaqr.github.io/hypernets/reference/persistent_homology.md) | Betti curves and persistence diagrams over a filtration | Edelsbrunner & Harer (2010) |
| [`persistence_landscape()`](https://mohsaqr.github.io/hypernets/reference/persistence_landscape.md), [`bottleneck_distance()`](https://mohsaqr.github.io/hypernets/reference/bottleneck_distance.md) | Comparable summaries of persistence diagrams | Bubenik (2015) |
| [`simplicial_degree()`](https://mohsaqr.github.io/hypernets/reference/simplicial_degree.md), [`q_analysis()`](https://mohsaqr.github.io/hypernets/reference/q_analysis.md) | Higher-order degree and q-connectivity | Atkin (1974) |
| [`verify_simplicial()`](https://mohsaqr.github.io/hypernets/reference/verify_simplicial.md) | Independent check of the clique enumeration | — |

### Hypergraphs

| Verb | Method | Reference |
|----|----|----|
| [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md) | Hyperedges from a network’s cliques | Burgio et al. (2020) |
| [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md) | Hyperedges from sliding windows over sequences | Ding et al. (2020) |
| [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md) | Hyperedges from bipartite group membership | Perc et al. (2013) |
| [`hypergraph_centrality()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_centrality.md) | Tensor eigenvector centralities (clique, Z, H) and EDVW PageRank | Benson (2019); Chitra & Raphael (2019) |
| [`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md) | Structural measures (hyperdegree, overlap, density) | — |
| [`hypergraph_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_laplacian.md), [`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md), [`hypergraph_joint_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_joint_cluster.md), [`hypergraph_transduction()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_transduction.md) | Weighted normalised Laplacian; RDC-Spec, RDC-SymNMF, J-NMF and JS-NMF clustering; label spreading | Zhou, Huang & Schölkopf (2006); Hayashi et al. (2020) |
| [`hypergraph_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md) | EDVW PageRank with personalization, sparse-capable | Chitra & Raphael (2019); Page et al. (1999) |
| [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md), [`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md) | Embedding nearest-neighbour hyperedges; vertex/hyperedge role swap | — |
| [`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md), [`hypergraph_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md), [`hypergraph_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md) | Projections: weighted pairwise network, association-weighted graph, s-line graph | [Coupette, Hartung & Katz (2024)](https://doi.org/10.1098/rsta.2023.0141); Aksoy et al. (2020) |
| [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md), [`hypergraph_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md), [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md) | Interval and contact hypergraphs in Dynet’s vocabulary (`actor`/`group`, alias detection, calendar clocks, observation bounds); active and cumulative snapshots on a `step`/`window` grid | Coupette, Hartung & Katz (2024) |
| [`hypergraph_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md), [`hypergraph_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md) | Hyperedge distributions and s-betweenness/s-closeness | Coupette, Hartung & Katz (2024); Aksoy et al. (2020) |
| [`hypergraph_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md) | Induced Y/T/O census and configuration-model null profile | Coupette, Hartung & Katz (2024) |
| [`hypergraph_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md), [`hypergraph_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md) | Repeated Infomap, AMI-medoid selection, coverage/performance/modularity/conductance | Coupette, Hartung & Katz (2024) |
| [`hg_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md), [`hypergraph_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md) | Spectral or symmetric-NMF node coordinates | Zhou et al. (2006); Hayashi et al. (2020) |
| [`hypergraph_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md) | Degree-preserving swap and configuration-model nulls | Chodrow (2020) |
| [`hypergraph_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md), [`text_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md), [`heterogeneous_hgat()`](https://mohsaqr.github.io/hypernets/reference/heterogeneous_hgat.md), [`hypergraph_hypergcn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hypergcn.md), [`hypergraph_hnhn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hnhn.md), [`hypergraph_allset()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_allset.md) | HGNN, HyperGAT, heterogeneous HGAT, HyperGCN, HNHN, AllDeepSets and AllSetTransformer in native torch | Feng et al. (2019); Linmei et al. (2019); Yadati et al. (2019); Ding et al. (2020); Dong et al. (2020); Chien et al. (2022) |

### Text hypergraphs

| Verb | Method | Reference |
|----|----|----|
| [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md) | Corpus to weighted document-word hypergraph: bag of words with smoothed tf-idf, token windows, or embedding kNN | Hayashi et al. (2020); Ding et al. (2020); Manning, Raghavan & Schütze (2008) |
| [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md), [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md) | Structural measures and tensor centralities as tidy tables | Benson (2019) |
| [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md), [`hypergraph_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md) | Spectral or symmetric-NMF topic clustering and per-cluster keywords | Zhou, Huang & Schölkopf (2006); Hayashi et al. (2020) |
| [`hypergraph_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_stability.md), [`hypergraph_agreement()`](https://mohsaqr.github.io/hypernets/reference/hg_agreement.md), [`hypergraph_seeds()`](https://mohsaqr.github.io/hypernets/reference/hg_seeds.md) | Seed stability, partition agreement (ARI/AMI/NMI), seed selection | Hubert & Arabie (1985); Vinh et al. (2010) |
| [`hypergraph_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md) | Few-label transductive classification, with class-mass normalization | Zhou et al. (2006); Zhu, Ghahramani & Lafferty (2003) |

A `text_hypergraph` *is* a `net_hypergraph`, so every hypergraph verb
takes it. Long-form verbs have direct compact aliases
([`hypergraph_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
and
[`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md),
for example): both names are the same function, not wrappers. The three
existing engine names
[`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md),
[`hypergraph_centrality()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_centrality.md)
and
[`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md)
remain distinct from their tidy `hg_*()` views. Bundled data:
`covid_abstracts` (165 abstracts) and `covid_embeddings`.

### Research sources

The source repository includes an indexed [research paper
library](https://mohsaqr.github.io/hypernets/papers/) and
[method/repository reference
notes](https://mohsaqr.github.io/hypernets/repos/). For the projection
tier, see the local [*Legal hypergraphs*
paper](https://mohsaqr.github.io/hypernets/papers/2024-PhilTrans-LegalHypergraphs-Coupette.pdf)
and its [method mapping and reproducibility
links](https://mohsaqr.github.io/hypernets/repos/legal-hypergraphs.md).

## The taxonomy

The same naming rules hold in every family, so a verb from one reads
like a verb from another:

- **`build_*()`** constructors always return a `net_*` object. Where a
  family admits several construction routes, they are selected with
  `type =`.
- **`*_centrality()`, `*_measures()`, `*_degree()`** always return a
  tidy `data.frame`, one row per node or per structure.
- **`bootstrap_*()`, `compare_*()`, `*_test()`** return a `net_*` result
  carrying estimates, intervals and p-values.
- **Every `net_*` result supports
  [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) and
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html).** The
  tidy accessor is the supported way in — reaching into a result with
  `$` is never required, and a secondary table is selected with
  `what =`:

``` r

as.data.frame(hon)                       # the rules
as.data.frame(hon, what = "nodes")       # the higher-order nodes
as.data.frame(ph,  what = "betti")       # Betti curves, not the diagram
as.data.frame(hm,  what = "global")      # whole-hypergraph summary
```

## Crossing between families

The families are entry points into one another, not islands. The same
coded sessions can be read all three ways:

``` r

library(honets)

# memory network: what follows what, given how you got here
hon <- build_hon(human_long, action = "code", actor = "session_id",
                 time = "timestamp", max_order = 3)
hon_centrality(hon, sort_by = "pagerank")

# simplicial complex: the same sessions as a pathway complex
sc <- build_simplicial(hon, type = "pathway")
betti_numbers(sc)

# hypergraph: each session as one multi-way interaction over its codes
hg <- group_hypergraph(human_long, actor = "code", group = "session_id")
hypergraph_centrality(hg, type = "pagerank")

# text: a corpus as a document-word hypergraph, clustered into topics
thg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
                       weight = "tfidf", stop_words = stop_words_en(),
                       min_count = 3L)
topics <- hg_cluster(thg, k = 4, type = "random_walk", seed = 1)
hg_keywords(thg, topics, n = 5, collapse = TRUE)                    # tf-idf mass
hg_keywords(thg, topics, n = 5, type = "ctfidf", collapse = TRUE)  # BERTopic's c-TF-IDF
hg_keywords(thg, topics, n = 5, type = "centrality", collapse = TRUE)
```

[`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md)
projects a hypergraph back to a pairwise network that any first-order
tool accepts;
[`pathways()`](https://mohsaqr.github.io/hypernets/reference/pathways.md)
hands sequence-derived path strings from the memory family to the other
two.

Bundled data: `human_long`, `ai_long` (coded human-AI pair-programming
sessions); `covid_abstracts`, `covid_embeddings` (COVID-19 education
research abstracts and their sentence embeddings).

## Vignettes

The package ships
[`vignette("legal-hypergraphs")`](https://mohsaqr.github.io/hypernets/articles/legal-hypergraphs.md),
[`vignette("text-hypergraphs")`](https://mohsaqr.github.io/hypernets/articles/text-hypergraphs.md)
and
[`vignette("text-constructions")`](https://mohsaqr.github.io/hypernets/articles/text-constructions.md),
plus a benchmark article (`vignettes/articles/benchmarks.Rmd`,
R8/R52/MR/Ohsumed/20NG). The other three families have one worked
document each in `docs/` — a worked analysis on real data, not a syntax
reference. Start with **`docs/index.html`** for the map, or
**`docs/overview.html`** to see the same dataset read all three ways.
Rebuild them with `Rscript docs/_knit_all.R`. Longer tutorials live in
`Tutorial_docs/`.

## Provenance

honets is the home package for higher-order structure in the
[Nestimate](https://github.com/mohsaqr/Nestimate) family, alongside
[psychnets](https://github.com/mohsaqr/psychnets) (cross-sectional
psychometric networks) and
[idiographic](https://github.com/mohsaqr/idiographic) (person-specific
temporal models).

The code was moved verbatim from Nestimate 0.9.0 in two steps — the
memory family in 0.1.0 (2026-08-24), the simplicial and hypergraph
families in 0.2.0 (2026-08-26, the latter by way of the short-lived
`hypernets` package, now folded in). The text family arrived in 0.3.0
(2026-09-01) from the `texthypergraph` package, which is retired; its
generic hypergraph methods (PageRank, projections, null models, sparse
engines, neural tier) joined the hypergraph family, and the merge
removed the Nestimate dependency. Cross-package identity is proven by
exact [`identical()`](https://rdrr.io/r/base/identical.html) tests in
`local_testing_and_equivalence/`, including the seeded permutation and
clique-sampling RNG paths. The 0.2.0 consolidation renamed part of the
*surface* to put all three families under one taxonomy — never a
computed value; see `NEWS.md` for the full mapping.

The inherited external-equivalence suites (pyHON, pathpy multi-order
models, Wallenius / Monte-Carlo HYPA references, igraph, HyperNetX) live
in `local_testing_and_equivalence/`, which is build-ignored and gated:

``` sh
NOT_CRAN=true HONETS_EQUIV_TESTS=true Rscript -e \
  'library(honets); testthat::test_dir("local_testing_and_equivalence")'
```

## Installation

``` r

# development version
devtools::install_github("mohsaqr/hypernets")
```

## References

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.

Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University
Press.
