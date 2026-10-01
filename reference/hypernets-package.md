# hypernets: Higher-Order Network Analysis

A higher-order network is one in which a relation reaches beyond a
single pair of nodes at a single moment: either it *binds more than two
nodes at once*, or it *depends on more than the current node*. The
literature (Battiston et al. 2020; Bianconi 2021) organizes that idea
into three structure families, and hypernets implements all three behind
one taxonomy.

## Ownership

The memory-network and simplicial-complex estimators are imported from a
sibling estimation package; hypernets wraps them under its own names
([`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
[`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md),
[`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md),
[`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md),
[`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
[`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md),
[`hg_markov_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_markov_stability.md),
[`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md),
[`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md),
[`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md),
[`hg_bottleneck()`](https://mohsaqr.github.io/hypernets/reference/hg_bottleneck.md),
[`hg_betti()`](https://mohsaqr.github.io/hypernets/reference/hg_betti.md),
[`hg_euler()`](https://mohsaqr.github.io/hypernets/reference/hg_euler.md),
[`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md),
[`hg_degree()`](https://mohsaqr.github.io/hypernets/reference/hg_degree.md))
and returns their result objects unchanged, so every number is the
estimator's own. The wrappers add one sequence-input contract (long,
wide, list or model input; see
[sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md)),
the
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
reader for every result class, and the verbs hypernets computes itself
on top of them
([`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md),
[`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md),
[`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md),
[`hg_wasserstein()`](https://mohsaqr.github.io/hypernets/reference/hg_wasserstein.md)).
The hypergraph and text families are hypernets' own.

## The three families

- **Memory networks**:

  A node is a state *plus the memory of how it was reached*, so a
  relation depends on history rather than only on the present state.
  Built from categorical sequences. Constructors
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
  [`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md),
  [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md),
  [`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md);
  diagnostics
  [`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
  [`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md);
  inference
  [`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md),
  [`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md);
  measures
  [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md).

- **Simplicial complexes**:

  A relation is a set of nodes that are *all* mutually related, together
  with every one of its subsets, which gives the object a geometry and
  hence a topology. Constructor
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md);
  measures
  [`hg_betti()`](https://mohsaqr.github.io/hypernets/reference/hg_betti.md),
  [`hg_euler()`](https://mohsaqr.github.io/hypernets/reference/hg_euler.md),
  [`hg_degree()`](https://mohsaqr.github.io/hypernets/reference/hg_degree.md),
  [`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md);
  topology
  [`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md),
  [`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md),
  [`hg_bottleneck()`](https://mohsaqr.github.io/hypernets/reference/hg_bottleneck.md),
  [`hg_wasserstein()`](https://mohsaqr.github.io/hypernets/reference/hg_wasserstein.md).

- **Hypergraphs**:

  A relation is an arbitrary set of nodes bound as a unit, with no
  requirement that its subsets also be relations. Constructors
  [`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md);
  random models
  [`hg_sample_gnp()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_gnp.md),
  [`hg_sample_sbm()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_sbm.md),
  [`hg_sample_uniform()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_uniform.md);
  measures
  [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
  [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md);
  spectral methods
  [`hg_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hg_laplacian.md),
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
  [`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md);
  PageRank
  [`hg_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md);
  embeddings
  [`hg_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md);
  projections
  [`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md),
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md),
  [`hg_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md),
  [`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md);
  temporal views
  [`hg_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hg_snapshot.md),
  [`hg_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hg_snapshots.md);
  hyperedge tables
  [`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
  and s-centrality
  [`hg_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md);
  communities
  [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  and
  [`hg_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md);
  motifs
  [`hg_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md);
  null models
  [`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md);
  neural networks
  [`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md),
  [`text_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md),
  [`heterogeneous_hgat()`](https://mohsaqr.github.io/hypernets/reference/heterogeneous_hgat.md),
  [`hg_hypergcn()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergcn.md),
  [`hg_hnhn()`](https://mohsaqr.github.io/hypernets/reference/hg_hnhn.md),
  [`hg_allset()`](https://mohsaqr.github.io/hypernets/reference/hg_allset.md);
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
  Every verb also accepts any `net_hg`.

## Verb grammar

The same naming rules hold across all families:

- Constructors are nouns:

  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
  [`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md),
  [`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md),
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md),
  [`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
  [`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md),
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  ... Where a family admits several construction routes, they are
  selected with `type =`.

- Everything else is `hg_*()`:

  Measures return a tidy `data.frame`, one row per node or per structure
  ([`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
  [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
  [`hg_degree()`](https://mohsaqr.github.io/hypernets/reference/hg_degree.md));
  inference returns a result object carrying estimates, intervals and
  p-values
  ([`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md),
  [`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md),
  [`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)).
  One verb names one idea and dispatches on its input:
  [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
  [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  and [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md)
  take a memory network or a hypergraph.

- [`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
  is the one reader:

  Every result object hands over its tables through `hg_get(x, what = )`
  – never through `$`. `what =` selects a secondary table; filters,
  `sort_by` and `top` are named arguments. hypernets' result classes
  also have [`print()`](https://rdrr.io/r/base/print.html),
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) methods.

## Crossing between families

The families are entry points into one another, not islands. The text
family is the corpus front end of the hypergraph family: a
`text_hypergraph` *is* a `net_hg`, so every hypergraph verb takes it,
and every `hg_*()` verb takes any `net_hg` in return.
[`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)
with `type = "pathway"` turns a memory network into a simplicial
complex;
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
turns the same sequences into a hypergraph;
[`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md)
projects a hypergraph back to a pairwise network that any of the
first-order tools accept. `hg_get(x, what = "pathways")` on a memory
network hands its sequence-derived path strings to the other two.

## References

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.

Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University
Press.

## See also

Useful links:

- <https://github.com/mohsaqr/hypernets>

- <https://mohsaqr.github.io/hypernets/>

- Report bugs at <https://github.com/mohsaqr/hypernets/issues>

## Author

**Maintainer**: Mohammed Saqr <saqr@saqr.me>
([ORCID](https://orcid.org/0000-0001-5881-3109)) \[copyright holder\]
