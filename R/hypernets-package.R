#' hypernets: Higher-Order Network Analysis
#'
#' A higher-order network is one in which a relation reaches beyond a single
#' pair of nodes at a single moment: either it *binds more than two nodes at
#' once*, or it *depends on more than the current node*. The literature
#' (Battiston et al. 2020; Bianconi 2021) organizes that idea into three
#' structure families, and hypernets implements all three behind one taxonomy.
#'
#' @section Ownership:
#'
#' The memory-network and simplicial-complex estimators are imported from a
#' sibling estimation package; hypernets wraps them under its own names
#' ([hon()], [honem()], [mogen()], [hypa()], [markov_order()],
#' [memory()], [hg_markov_stability()], [simplicial()],
#' [hg_homology()], [hg_landscape()], [hg_bottleneck()], [hg_betti()],
#' [hg_euler()], [hg_qanalysis()], [hg_degree()]) and returns their result
#' objects unchanged, so every number is the estimator's own. The wrappers
#' add one sequence-input contract (long, wide, list or model input; see
#' [sequence-input]), the [hg_get()] reader for every result class, and the
#' verbs hypernets computes itself on top of them ([hg_bootstrap()],
#' [hg_compare()], [hg_centrality()], [hg_communities()],
#' [hg_wasserstein()]). The hypergraph and text families are hypernets' own.
#'
#' @section The three families:
#'
#' \describe{
#'   \item{**Memory networks**}{A node is a state *plus the memory of how it
#'     was reached*, so a relation depends on history rather than only on the
#'     present state. Built from categorical sequences.
#'     Constructors [hon()], [honem()], [hypa()],
#'     [mogen()]; diagnostics [markov_order()],
#'     [memory()]; inference [hg_bootstrap()], [hg_compare()];
#'     measures [hg_centrality()].}
#'   \item{**Simplicial complexes**}{A relation is a set of nodes that are
#'     *all* mutually related, together with every one of its subsets, which
#'     gives the object a geometry and hence a topology.
#'     Constructor [simplicial()]; measures [hg_betti()],
#'     [hg_euler()], [hg_degree()], [hg_qanalysis()];
#'     topology [hg_homology()], [hg_landscape()],
#'     [hg_bottleneck()], [hg_wasserstein()].}
#'   \item{**Hypergraphs**}{A relation is an arbitrary set of nodes bound as
#'     a unit, with no requirement that its subsets also be relations.
#'     Constructors [network_hypergraph()], [window_hypergraph()],
#'     [group_hypergraph()], [temporal_hypergraph()]; random models
#'     [hg_sample_gnp()], [hg_sample_sbm()], [hg_sample_uniform()]; measures [hg_measures()],
#'     [hg_centrality()]; spectral methods [hg_laplacian()],
#'     [hg_cluster()], [hg_classify()]; PageRank
#'     [hg_pagerank()]; embeddings [hg_embed()]; projections [hg_clique_expansion()],
#'     [hg_project()], [hg_line_graph()],
#'     [dual_hypergraph()]; temporal views [hg_snapshot()],
#'     [hg_snapshots()]; hyperedge tables [hg_edges()] and
#'     s-centrality [hg_edge_centrality()]; communities
#'     [hg_communities()] and [hg_community_quality()]; motifs
#'     [hg_motifs()]; null models [hg_null_test()]; neural networks
#'     [hg_neural()], [text_hypergat()], [heterogeneous_hgat()],
#'     [hg_hypergcn()], [hg_hnhn()], [hg_allset()]; embedding constructor
#'     [knn_hypergraph()].}
#'   \item{**Text hypergraphs**}{A corpus is a bipartite document-word
#'     structure, which is a hypergraph in either orientation: documents as
#'     nodes bound by shared words, or words as nodes bound by shared
#'     documents. Constructor [text_hypergraph()] (bag of words, token
#'     windows, or embedding nearest neighbours); tidy verbs
#'     [hg_measures()], [hg_centrality()], [hg_cluster()], [hg_keywords()],
#'     [hg_classify()], [hg_stability()], [hg_agreement()], [hg_seeds()].
#'     Every verb also accepts any `net_hg`.}
#' }
#'
#' @section Verb grammar:
#'
#' The same naming rules hold across all families:
#'
#' \describe{
#'   \item{Constructors are nouns}{[hon()], [honem()], [mogen()],
#'     [simplicial()], [markov_order()], [memory()],
#'     [text_hypergraph()],
#'     [window_hypergraph()], [group_hypergraph()], ... Where a family admits
#'     several construction routes, they are selected with `type =`.}
#'   \item{Everything else is `hg_*()`}{Measures return a tidy
#'     `data.frame`, one row per node or per structure
#'     ([hg_centrality()], [hg_measures()], [hg_degree()]); inference
#'     returns a result object carrying estimates, intervals and p-values
#'     ([hg_bootstrap()], [hg_compare()], [hg_null_test()]). One verb
#'     names one idea and dispatches on its input: [hg_centrality()],
#'     [hg_communities()] and [hypa()] take a memory network or a
#'     hypergraph.}
#'   \item{[hg_get()] is the one reader}{Every result object hands over its
#'     tables through `hg_get(x, what = )` -- never through `$`. `what =`
#'     selects a secondary table; filters, `sort_by` and `top` are named
#'     arguments. hypernets' result classes also have `print()`,
#'     `summary()` and `plot()` methods.}
#' }
#'
#' @section Crossing between families:
#'
#' The families are entry points into one another, not islands. The text
#' family is the corpus front end of the hypergraph family: a
#' `text_hypergraph` *is* a `net_hg`, so every hypergraph verb takes
#' it, and every `hg_*()` verb takes any `net_hg` in return.
#' [simplicial()] with `type = "pathway"` turns a memory network into a
#' simplicial complex; [window_hypergraph()] turns the same sequences into a
#' hypergraph; [hg_clique_expansion()] projects a hypergraph back to a pairwise
#' network that any of the first-order tools accept.
#' `hg_get(x, what = "pathways")` on a memory network hands its
#' sequence-derived path strings to the other two.
#'
#' @references
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
#' Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
#'
#' Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University Press.
#'
#' @keywords internal
#' @importFrom cograph splot
#' @importFrom Nestimate build_simplicial
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom ggplot2 .data
"_PACKAGE"
