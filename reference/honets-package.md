# honets: Higher-Order Network Analysis

A higher-order network is one in which a relation reaches beyond a
single pair of nodes at a single moment: either it *binds more than two
nodes at once*, or it *depends on more than the current node*. The
literature (Battiston et al. 2020; Bianconi 2021) organizes that idea
into three structure families, and honets implements all three behind
one taxonomy.

## The three families

- **Memory networks**:

  A node is a state *plus the memory of how it was reached*, so a
  relation depends on history rather than only on the present state.
  Built from categorical sequences. Constructors
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md),
  [`build_honem()`](https://mohsaqr.github.io/hypernets/reference/build_honem.md),
  [`build_hypa()`](https://mohsaqr.github.io/hypernets/reference/build_hypa.md),
  [`build_mogen()`](https://mohsaqr.github.io/hypernets/reference/build_mogen.md);
  diagnostics
  [`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md),
  [`path_dependence()`](https://mohsaqr.github.io/hypernets/reference/path_dependence.md);
  inference
  [`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md),
  [`compare_hon()`](https://mohsaqr.github.io/hypernets/reference/compare_hon.md);
  measures
  [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md).

- **Simplicial complexes**:

  A relation is a set of nodes that are *all* mutually related, together
  with every one of its subsets, which gives the object a geometry and
  hence a topology. Constructor
  [`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md);
  measures
  [`betti_numbers()`](https://mohsaqr.github.io/hypernets/reference/betti_numbers.md),
  [`euler_characteristic()`](https://mohsaqr.github.io/hypernets/reference/euler_characteristic.md),
  [`simplicial_degree()`](https://mohsaqr.github.io/hypernets/reference/simplicial_degree.md),
  [`q_analysis()`](https://mohsaqr.github.io/hypernets/reference/q_analysis.md);
  topology
  [`persistent_homology()`](https://mohsaqr.github.io/hypernets/reference/persistent_homology.md),
  [`persistence_landscape()`](https://mohsaqr.github.io/hypernets/reference/persistence_landscape.md),
  [`bottleneck_distance()`](https://mohsaqr.github.io/hypernets/reference/bottleneck_distance.md),
  [`wasserstein_distance()`](https://mohsaqr.github.io/hypernets/reference/wasserstein_distance.md).

- **Hypergraphs**:

  A relation is an arbitrary set of nodes bound as a unit, with no
  requirement that its subsets also be relations. Constructors
  [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md),
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md);
  random models
  [`hg_sample_gnp()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_gnp.md),
  [`hg_sample_sbm()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_sbm.md),
  [`hg_sample_uniform()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_uniform.md);
  measures
  [`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md),
  [`hypergraph_centrality()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_centrality.md);
  spectral methods
  [`hypergraph_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_laplacian.md),
  [`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md),
  [`hypergraph_transduction()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_transduction.md);
  PageRank
  [`hypergraph_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md);
  embeddings
  [`hypergraph_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md);
  projections
  [`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md),
  [`hypergraph_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md),
  [`hypergraph_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md),
  [`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md);
  temporal views
  [`hypergraph_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md),
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md);
  hyperedge tables
  [`hypergraph_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
  and s-centrality
  [`hypergraph_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md);
  communities
  [`hypergraph_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  and
  [`hypergraph_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md);
  motifs
  [`hypergraph_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md);
  null models
  [`hypergraph_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md);
  neural networks
  [`hypergraph_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md),
  [`text_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md),
  [`heterogeneous_hgat()`](https://mohsaqr.github.io/hypernets/reference/heterogeneous_hgat.md),
  [`hypergraph_hypergcn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hypergcn.md),
  [`hypergraph_hnhn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hnhn.md),
  [`hypergraph_allset()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_allset.md);
  embedding constructor
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md).

- **Text hypergraphs**:

  A corpus is a bipartite document-word structure, which is a hypergraph
  in either orientation: documents as nodes bound by shared words, or
  words as nodes bound by shared documents. Constructor
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (bag of words, token windows, or embedding nearest neighbours); tidy
  verbs
  [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
  [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
  [`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md),
  [`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md),
  [`hg_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_stability.md),
  [`hg_agreement()`](https://mohsaqr.github.io/hypernets/reference/hg_agreement.md),
  [`hg_seeds()`](https://mohsaqr.github.io/hypernets/reference/hg_seeds.md).
  Every verb also accepts any `net_hypergraph`.

Long-form hypergraph verbs also have direct compact `hg_*` aliases; for
example,
[`hypergraph_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
and
[`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
are the same function. Existing engine names
([`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md),
[`hypergraph_centrality()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_centrality.md),
[`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md))
remain distinct from their tidy `hg_*` views.

## Verb grammar

The same naming rules hold across all three families:

- `build_*()`:

  Constructors. Always return an object of class `net_*`. Where a family
  admits several construction routes, they are selected with `type =`,
  never with a differently named argument.

- `*_centrality()`, `*_measures()`, `*_degree()`:

  Measures. Always return a tidy `data.frame`, one row per node or per
  structure.

- `bootstrap_*()`, `compare_*()`, `*_test()`:

  Inference. Return a `net_*` result object carrying estimates,
  intervals and p-values.

- [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html),
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

  Every `net_*` result object supports all four.
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) is the
  tidy accessor – reaching into a result with `$` is never required, and
  a secondary table is selected with `what =`.

## Crossing between families

The families are entry points into one another, not islands. The text
family is the corpus front end of the hypergraph family: a
`text_hypergraph` *is* a `net_hypergraph`, so every hypergraph verb
takes it, and every `hg_*()` verb takes any `net_hypergraph` in return.
[`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md)
with `type = "pathway"` turns a memory network into a simplicial
complex;
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
turns the same sequences into a hypergraph;
[`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md)
projects a hypergraph back to a pairwise network that any of the
first-order tools accept.
[`pathways()`](https://mohsaqr.github.io/hypernets/reference/pathways.md)
is the shared accessor that hands sequence-derived path strings from the
memory family to the other two.

## References

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.

Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University
Press.

## See also

Useful links:

- <https://github.com/mohsaqr/hypernets>

- <https://mohsaqr.github.io/hypernets>

- Report bugs at <https://github.com/mohsaqr/hypernets/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me>
([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]
