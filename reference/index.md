# Package index

## Package

- [`hypernets`](https://mohsaqr.github.io/hypernets/reference/hypernets-package.md)
  [`hypernets-package`](https://mohsaqr.github.io/hypernets/reference/hypernets-package.md)
  : hypernets: Higher-Order Network Analysis

## Memory networks

A node is a state plus the memory of how it was reached.

- [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
  : Build a Higher-Order Network (HON)
- [`build_honem()`](https://mohsaqr.github.io/hypernets/reference/build_honem.md)
  : Build HONEM Embeddings for Higher-Order Networks
- [`build_hypa()`](https://mohsaqr.github.io/hypernets/reference/build_hypa.md)
  : Detect Path Anomalies via HYPA
- [`build_mogen()`](https://mohsaqr.github.io/hypernets/reference/build_mogen.md)
  : Build Multi-Order Generative Model (MOGen)
- [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md)
  : Higher-order centralities with first-order projection
- [`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md)
  : Bootstrap inference for higher-order network rules
- [`compare_hon()`](https://mohsaqr.github.io/hypernets/reference/compare_hon.md)
  : Two-sample permutation comparison of higher-order network rules
- [`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md)
  : Test the Markov order of a sequential process
- [`path_dependence()`](https://mohsaqr.github.io/hypernets/reference/path_dependence.md)
  : Per-Context Path Dependence at Order k
- [`mogen_transitions()`](https://mohsaqr.github.io/hypernets/reference/mogen_transitions.md)
  : Extract Transition Table from a MOGen Model
- [`path_counts()`](https://mohsaqr.github.io/hypernets/reference/path_counts.md)
  : Count Path Frequencies in Trajectory Data
- [`pathways()`](https://mohsaqr.github.io/hypernets/reference/pathways.md)
  : Extract Pathways from Higher-Order Network Objects

## Simplicial complexes

A set of mutually related nodes together with every subset.

- [`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md)
  : Build a Simplicial Complex
- [`betti_numbers()`](https://mohsaqr.github.io/hypernets/reference/betti_numbers.md)
  : Betti Numbers
- [`euler_characteristic()`](https://mohsaqr.github.io/hypernets/reference/euler_characteristic.md)
  : Euler Characteristic
- [`simplicial_degree()`](https://mohsaqr.github.io/hypernets/reference/simplicial_degree.md)
  : Simplicial Degree
- [`q_analysis()`](https://mohsaqr.github.io/hypernets/reference/q_analysis.md)
  : Q-Analysis
- [`persistent_homology()`](https://mohsaqr.github.io/hypernets/reference/persistent_homology.md)
  : Persistent Homology
- [`persistence_landscape()`](https://mohsaqr.github.io/hypernets/reference/persistence_landscape.md)
  : Persistence Landscape
- [`bottleneck_distance()`](https://mohsaqr.github.io/hypernets/reference/bottleneck_distance.md)
  : Bottleneck Distance Between Persistence Diagrams
- [`wasserstein_distance()`](https://mohsaqr.github.io/hypernets/reference/wasserstein_distance.md)
  : Wasserstein Distance Between Persistence Diagrams
- [`verify_simplicial()`](https://mohsaqr.github.io/hypernets/reference/verify_simplicial.md)
  : Verify Simplicial Complex Against igraph

## Hypergraphs

An arbitrary set of nodes bound as a unit.

### Constructors

- [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  [`print(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  [`summary(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  : Higher-order hypergraph from a network's clique structure
- [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
  : Windowed Sequence Hyperedges
- [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
  : Hypergraph from co-occurrence data or an edge list
- [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md)
  : Build a k-nearest-neighbor hypergraph from embeddings
- [`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md)
  : Dual of a hypergraph
- [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md)
  : Temporal hypergraph from an edge list or co-occurrence data
- [`hypergraph_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md)
  [`hg_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md)
  : Extract one snapshot from a temporal hypergraph
- [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  [`hg_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  : Extract a sequence of temporal-hypergraph snapshots
- [`hg_subset()`](https://mohsaqr.github.io/hypernets/reference/hg_subset.md)
  [`hypergraph_subset()`](https://mohsaqr.github.io/hypernets/reference/hg_subset.md)
  : Sub-hypergraph by hyperedges, nodes or hyperedge attributes
- [`hg_sample_gnp()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_gnp.md)
  [`hypergraph_sample_gnp()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_gnp.md)
  : Sample Bernoulli-Incidence Random Hypergraphs
- [`hg_sample_sbm()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_sbm.md)
  [`hypergraph_sample_sbm()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_sbm.md)
  : Sample Stochastic-Block-Model Hypergraphs
- [`hg_sample_uniform()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_uniform.md)
  [`hg_sample_regular()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_uniform.md)
  [`hypergraph_sample_uniform()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_uniform.md)
  [`hypergraph_sample_regular()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_uniform.md)
  : Sample Uniform or Regular Random Hypergraphs

### Measures and centralities

- [`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md)
  [`print(`*`<net_hypergraph_measures>`*`)`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md)
  : Structural measures for a hypergraph
- [`hypergraph_centrality()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_centrality.md)
  : Hypergraph eigenvector centralities
- [`hg_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md)
  [`hypergraph_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md)
  : Hypergraph PageRank with edge-dependent vertex weights
- [`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
  [`hypergraph_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
  : Hyperedge-level measures, as a tidy table
- [`hg_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md)
  [`hypergraph_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md)
  : s-betweenness and s-closeness centrality of hyperedges
- [`hg_growth()`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  [`hypergraph_growth()`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  [`plot(`*`<hypernets_series>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  : Growth of a temporal hypergraph over time
- [`hg_representations()`](https://mohsaqr.github.io/hypernets/reference/hg_representations.md)
  [`hypergraph_representations()`](https://mohsaqr.github.io/hypernets/reference/hg_representations.md)
  : Graph and hypergraph representations of the same data
- [`hg_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  [`hypergraph_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  [`as.data.frame(`*`<hypernets_motifs>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  [`plot(`*`<hypernets_motifs>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  : Y, T and O motif census with a hypergraph configuration null

### Spectral methods

- [`hypergraph_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_laplacian.md)
  : Normalized hypergraph Laplacian
- [`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md)
  : Spectral or symmetric-NMF clustering of hypergraph vertices
- [`hypergraph_joint_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_joint_cluster.md)
  : Joint NMF clustering with an auxiliary vertex relation
- [`hg_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md)
  [`hypergraph_embed()`](https://mohsaqr.github.io/hypernets/reference/hg_embed.md)
  : Low-Dimensional Hypergraph Embedding
- [`hypergraph_transduction()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_transduction.md)
  : Transductive label spreading on a hypergraph

### Projections

- [`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md)
  : Clique expansion of a hypergraph
- [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md)
  [`hypergraph_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md)
  : Project a hypergraph onto a weighted graph
- [`hg_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md)
  [`hypergraph_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md)
  : The s-line graph of a hypergraph

### Communities

- [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`print(`*`<hg_communities>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`as.data.frame(`*`<hg_communities>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`plot(`*`<hg_communities>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`hypergraph_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  : Stable Infomap communities of a hypergraph projection
- [`hg_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md)
  [`hypergraph_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md)
  : Quality of a projected-hypergraph partition
- [`hg_compare_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`print(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`as.data.frame(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`summary(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`plot(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  : Compare community structure across representations

### Null models and neural networks

- [`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)
  [`hypergraph_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)
  : Degree-preserving null-model test of hypergraph structure
- [`hg_hypa()`](https://mohsaqr.github.io/hypernets/reference/hg_hypa.md)
  : Hypergeometric anomaly detection for co-occurring node pairs
- [`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md)
  [`hypergraph_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md)
  : Hypergraph neural network classifier (HGNN)
- [`hg_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md)
  [`text_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md)
  : HyperGAT document classifier
- [`heterogeneous_hgat()`](https://mohsaqr.github.io/hypernets/reference/heterogeneous_hgat.md)
  : Heterogeneous graph attention classifier (HGAT)
- [`hypergraph_hypergcn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hypergcn.md)
  [`hg_hypergcn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hypergcn.md)
  : HyperGCN node classifier
- [`hypergraph_hnhn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hnhn.md)
  [`hg_hnhn()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_hnhn.md)
  : HNHN node classifier
- [`hypergraph_allset()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_allset.md)
  [`hg_allset()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_allset.md)
  [`hypergraph_alldeepsets()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_allset.md)
  [`hypergraph_allset_transformer()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_allset.md)
  : AllSet hypergraph node classifier

## Text hypergraphs

A corpus as a weighted document-word hypergraph, analysed with tidy
verbs.

- [`clean_text()`](https://mohsaqr.github.io/hypernets/reference/clean_text.md)
  : Clean a corpus before building a text hypergraph
- [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  : Build a weighted hypergraph from a text corpus
- [`stop_words_en()`](https://mohsaqr.github.io/hypernets/reference/stop_words_en.md)
  : English function-word stop list
- [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md)
  : Structural measures of a hypergraph, as tidy tables
- [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md)
  : Hypergraph node centralities, as a tidy table
- [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  : Spectral clustering of a hypergraph, as a tidy table
- [`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  [`print(`*`<hypernets_keywords>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  [`plot(`*`<hypernets_keywords>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  [`hypergraph_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  : Characteristic words (keywords) per cluster
- [`hg_relations()`](https://mohsaqr.github.io/hypernets/reference/hg_relations.md)
  [`hypergraph_relations()`](https://mohsaqr.github.io/hypernets/reference/hg_relations.md)
  : Relations between topics: a weighted network of shared vocabulary
- [`hg_topic_sizes()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
  [`plot(`*`<hypernets_topic_sizes>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
  [`hypergraph_topic_sizes()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
  : Topic sizes
- [`hg_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
  [`plot(`*`<hypernets_topic_quality>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
  [`hypergraph_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
  : Topic coherence and exclusivity
- [`hg_membership()`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
  [`plot(`*`<hypernets_membership>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
  [`hypergraph_membership()`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
  : Soft topic membership from the spectral embedding
- [`hg_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_stability.md)
  [`hypergraph_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_stability.md)
  : Seed stability of a hypergraph clustering across resolutions
- [`hg_agreement()`](https://mohsaqr.github.io/hypernets/reference/hg_agreement.md)
  [`hypergraph_agreement()`](https://mohsaqr.github.io/hypernets/reference/hg_agreement.md)
  : Agreement between two labelings of the same nodes
- [`hg_seeds()`](https://mohsaqr.github.io/hypernets/reference/hg_seeds.md)
  [`hypergraph_seeds()`](https://mohsaqr.github.io/hypernets/reference/hg_seeds.md)
  : Seed labels from a clustering, for spreading or training
- [`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md)
  [`hypergraph_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md)
  : Transductive label spreading on a hypergraph, as a tidy table

## Methods and accessors

- [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  [`print(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  [`summary(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  : Higher-order hypergraph from a network's clique structure

- [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`print(`*`<hg_communities>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`as.data.frame(`*`<hg_communities>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`plot(`*`<hg_communities>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  [`hypergraph_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  : Stable Infomap communities of a hypergraph projection

- [`hg_compare_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`print(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`as.data.frame(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`summary(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  [`plot(`*`<hypernets_community_comparison>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md)
  : Compare community structure across representations

- [`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  [`print(`*`<hypernets_keywords>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  [`plot(`*`<hypernets_keywords>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  [`hypergraph_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  : Characteristic words (keywords) per cluster

- [`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md)
  [`print(`*`<net_hypergraph_measures>`*`)`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md)
  : Structural measures for a hypergraph

- [`print(`*`<net_hon>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_hon.md)
  : Print Method for net_hon

- [`print(`*`<net_hon_boot>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_hon_boot.md)
  : Print method for net_hon_boot

- [`print(`*`<net_hon_compare>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_hon_compare.md)
  : Print method for net_hon_compare

- [`print(`*`<net_honem>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_honem.md)
  : Print Method for net_honem

- [`print(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_hypa.md)
  : Print Method for net_hypa

- [`print(`*`<net_hypergraph_cluster>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_hypergraph_cluster.md)
  : Print method for net_hypergraph_cluster

- [`print(`*`<net_hypergraph_transduction>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_hypergraph_transduction.md)
  : Print method for net_hypergraph_transduction

- [`print(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_markov_order.md)
  : Print Method for net_markov_order

- [`print(`*`<net_markov_order_group>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_markov_order_group.md)
  :

  Print method for `net_markov_order_group`

- [`print(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_mogen.md)
  : Print Method for net_mogen

- [`print(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_path_dependence.md)
  :

  Print method for `net_path_dependence`

- [`print(`*`<net_persistence_landscape>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_persistence_landscape.md)
  : Print Persistence Landscape

- [`print(`*`<net_persistent_homology>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_persistent_homology.md)
  : Print persistent homology results

- [`print(`*`<net_q_analysis>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_q_analysis.md)
  : Print Q-analysis results

- [`print(`*`<net_simplicial>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.net_simplicial.md)
  : Print a simplicial complex

- [`print(`*`<summary.net_path_dependence>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.summary.net_path_dependence.md)
  :

  Print method for `summary.net_path_dependence`

- [`print(`*`<text_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/print.text_hypergraph.md)
  : Print a text hypergraph

- [`summary(`*`<net_hon>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_hon.md)
  : Summary Method for net_hon

- [`summary(`*`<net_hon_boot>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_hon_boot.md)
  : Summary method for net_hon_boot

- [`summary(`*`<net_hon_compare>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_hon_compare.md)
  : Summary method for net_hon_compare

- [`summary(`*`<net_honem>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_honem.md)
  : Summary Method for net_honem

- [`summary(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_hypa.md)
  : Summary Method for net_hypa

- [`summary(`*`<net_hypergraph_cluster>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_hypergraph_cluster.md)
  : Summary method for net_hypergraph_cluster

- [`summary(`*`<net_hypergraph_transduction>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_hypergraph_transduction.md)
  : Summary method for net_hypergraph_transduction

- [`summary(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_markov_order.md)
  : Summary Method for net_markov_order

- [`summary(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_mogen.md)
  : Summary Method for net_mogen

- [`summary(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/hypernets/reference/summary.net_path_dependence.md)
  :

  Summary method for `net_path_dependence`

- [`hg_growth()`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  [`hypergraph_growth()`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  [`plot(`*`<hypernets_series>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  : Growth of a temporal hypergraph over time

- [`hg_membership()`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
  [`plot(`*`<hypernets_membership>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
  [`hypergraph_membership()`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
  : Soft topic membership from the spectral embedding

- [`hg_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  [`hypergraph_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  [`as.data.frame(`*`<hypernets_motifs>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  [`plot(`*`<hypernets_motifs>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md)
  : Y, T and O motif census with a hypergraph configuration null

- [`hg_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
  [`plot(`*`<hypernets_topic_quality>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
  [`hypergraph_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
  : Topic coherence and exclusivity

- [`hg_topic_sizes()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
  [`plot(`*`<hypernets_topic_sizes>`*`)`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
  [`hypergraph_topic_sizes()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
  : Topic sizes

- [`plot(`*`<net_hon_boot>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_hon_boot.md)
  : Plot method for net_hon_boot

- [`plot(`*`<net_hon_compare>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_hon_compare.md)
  : Plot method for net_hon_compare

- [`plot(`*`<net_honem>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_honem.md)
  : Plot Method for net_honem

- [`plot(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_hypergraph.md)
  : Plot a hypergraph with hyperedges as blobs

- [`plot(`*`<net_hypergraph_cluster>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_hypergraph_cluster.md)
  : Plot method for net_hypergraph_cluster

- [`plot(`*`<net_hypergraph_transduction>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_hypergraph_transduction.md)
  : Plot method for net_hypergraph_transduction

- [`plot(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_markov_order.md)
  : Plot Method for net_markov_order

- [`plot(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_mogen.md)
  : Plot Method for net_mogen

- [`plot(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_path_dependence.md)
  :

  Plot method for `net_path_dependence`

- [`plot(`*`<net_persistence_landscape>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_persistence_landscape.md)
  : Plot Persistence Landscape

- [`plot(`*`<net_persistent_homology>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_persistent_homology.md)
  : Plot Persistent Homology

- [`plot(`*`<net_q_analysis>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_q_analysis.md)
  : Plot Q-Analysis

- [`plot(`*`<net_simplicial>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_simplicial.md)
  : Plot a Simplicial Complex

- [`plot(`*`<net_temporal_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/plot.net_temporal_hypergraph.md)
  : Plot a temporal-hypergraph snapshot through cograph

- [`as.data.frame(`*`<net_hon>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hon.md)
  : Coerce a net_hon to a tidy table

- [`as.data.frame(`*`<net_hon_boot>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hon_boot.md)
  : Coerce a net_hon_boot to its tidy inference table

- [`as.data.frame(`*`<net_hon_compare>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hon_compare.md)
  : Coerce a net_hon_compare to its tidy edge table

- [`as.data.frame(`*`<net_honem>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_honem.md)
  : Coerce a net_honem to a tidy table

- [`as.data.frame(`*`<net_hypa>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hypa.md)
  : Coerce a net_hypa to a tidy table

- [`as.data.frame(`*`<net_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hypergraph.md)
  : Coerce a net_hypergraph to a tidy hyperedge table

- [`as.data.frame(`*`<net_hypergraph_cluster>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hypergraph_cluster.md)
  : Coerce a net_hypergraph_cluster to a data.frame

- [`as.data.frame(`*`<net_hypergraph_measures>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hypergraph_measures.md)
  : Coerce a net_hypergraph_measures to a tidy table

- [`as.data.frame(`*`<net_hypergraph_transduction>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hypergraph_transduction.md)
  : Coerce a net_hypergraph_transduction to a data.frame

- [`as.data.frame(`*`<net_markov_order>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_markov_order.md)
  : Coerce a net_markov_order to a tidy table

- [`as.data.frame(`*`<net_mogen>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_mogen.md)
  : Coerce a net_mogen to a tidy table

- [`as.data.frame(`*`<net_path_dependence>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_path_dependence.md)
  : Coerce a net_path_dependence to its tidy per-context table

- [`as.data.frame(`*`<net_persistence_landscape>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_persistence_landscape.md)
  : Coerce a net_persistence_landscape to its tidy table

- [`as.data.frame(`*`<net_persistent_homology>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_persistent_homology.md)
  : Coerce a net_persistent_homology to a tidy table

- [`as.data.frame(`*`<net_q_analysis>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_q_analysis.md)
  : Coerce a net_q_analysis to a tidy table

- [`as.data.frame(`*`<net_simplicial>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_simplicial.md)
  : Coerce a net_simplicial to a tidy table

- [`as.data.frame(`*`<net_temporal_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_temporal_hypergraph.md)
  : Tidy tables of a temporal hypergraph

- [`as.data.frame(`*`<text_hypergraph>`*`)`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.text_hypergraph.md)
  : Tidy tables of a text hypergraph

## Data

- [`human_long`](https://mohsaqr.github.io/hypernets/reference/long-data.md)
  [`ai_long`](https://mohsaqr.github.io/hypernets/reference/long-data.md)
  : Human-AI Vibe Coding Interaction Data (Long Format)
- [`covid_abstracts`](https://mohsaqr.github.io/hypernets/reference/covid_abstracts.md)
  : COVID-19 education research abstracts
- [`covid_embeddings`](https://mohsaqr.github.io/hypernets/reference/covid_embeddings.md)
  : Sentence embeddings of the COVID-19 abstracts
- [`icsid_tribunals`](https://mohsaqr.github.io/hypernets/reference/icsid_tribunals.md)
  : ICSID arbitration tribunals, one row per seat
- [`gfcc_decisions`](https://mohsaqr.github.io/hypernets/reference/gfcc_decisions.md)
  : Decisions of the German Federal Constitutional Court
- [`gfcc_citations`](https://mohsaqr.github.io/hypernets/reference/gfcc_citations.md)
  : Citation blocks of the German Federal Constitutional Court
