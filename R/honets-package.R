#' honets: Higher-Order Network Analysis
#'
#' A higher-order network is one in which a relation reaches beyond a single
#' pair of nodes at a single moment: either it *binds more than two nodes at
#' once*, or it *depends on more than the current node*. The literature
#' (Battiston et al. 2020; Bianconi 2021) organizes that idea into three
#' structure families, and honets implements all three behind one taxonomy.
#'
#' @section The three families:
#'
#' \describe{
#'   \item{**Memory networks**}{A node is a state *plus the memory of how it
#'     was reached*, so a relation depends on history rather than only on the
#'     present state. Built from categorical sequences.
#'     Constructors [build_hon()], [build_honem()], [build_hypa()],
#'     [build_mogen()]; diagnostics [markov_order_test()],
#'     [path_dependence()]; inference [bootstrap_hon()], [compare_hon()];
#'     measures [hon_centrality()].}
#'   \item{**Simplicial complexes**}{A relation is a set of nodes that are
#'     *all* mutually related, together with every one of its subsets, which
#'     gives the object a geometry and hence a topology.
#'     Constructor [build_simplicial()]; measures [betti_numbers()],
#'     [euler_characteristic()], [simplicial_degree()], [q_analysis()];
#'     topology [persistent_homology()], [persistence_landscape()],
#'     [bottleneck_distance()].}
#'   \item{**Hypergraphs**}{A relation is an arbitrary set of nodes bound as
#'     a unit, with no requirement that its subsets also be relations.
#'     Constructors [build_hypergraph()], [window_hypergraph()],
#'     [group_hypergraph()], [temporal_hypergraph()]; measures [hypergraph_measures()],
#'     [hypergraph_centrality()]; spectral methods [hypergraph_laplacian()],
#'     [hypergraph_cluster()], [hypergraph_transduction()]; PageRank
#'     [hypergraph_pagerank()]; projections [clique_expansion()],
#'     [hypergraph_project()], [hypergraph_line_graph()],
#'     [dual_hypergraph()]; temporal views [hypergraph_snapshot()],
#'     [hypergraph_snapshots()]; hyperedge tables [hypergraph_edges()] and
#'     s-centrality [hypergraph_edge_centrality()]; communities
#'     [hypergraph_communities()] and [hypergraph_community_quality()]; motifs
#'     [hypergraph_motifs()]; null models [hypergraph_null_test()]; neural networks
#'     [hypergraph_neural()], [text_hypergat()]; embedding constructor
#'     [knn_hypergraph()].}
#'   \item{**Text hypergraphs**}{A corpus is a bipartite document-word
#'     structure, which is a hypergraph in either orientation: documents as
#'     nodes bound by shared words, or words as nodes bound by shared
#'     documents. Constructor [text_hypergraph()] (bag of words, token
#'     windows, or embedding nearest neighbours); tidy verbs
#'     [hg_measures()], [hg_centrality()], [hg_cluster()], [hg_keywords()],
#'     [hg_classify()], [hg_stability()], [hg_agreement()], [hg_seeds()].
#'     Every verb also accepts any `net_hypergraph`.}
#' }
#'
#' Long-form hypergraph verbs also have direct compact `hg_*` aliases; for
#' example, [hypergraph_edges()] and [hg_edges()] are the same function.
#' Existing engine names ([hypergraph_measures()], [hypergraph_centrality()],
#' [hypergraph_cluster()]) remain distinct from their tidy `hg_*` views.
#'
#' @section Verb grammar:
#'
#' The same naming rules hold across all three families:
#'
#' \describe{
#'   \item{`build_*()`}{Constructors. Always return an object of class
#'     `net_*`. Where a family admits several construction routes, they are
#'     selected with `type =`, never with a differently named argument.}
#'   \item{`*_centrality()`, `*_measures()`, `*_degree()`}{Measures. Always
#'     return a tidy `data.frame`, one row per node or per structure.}
#'   \item{`bootstrap_*()`, `compare_*()`, `*_test()`}{Inference. Return a
#'     `net_*` result object carrying estimates, intervals and p-values.}
#'   \item{`print()`, `summary()`, `plot()`, `as.data.frame()`}{Every `net_*`
#'     result object supports all four. `as.data.frame()` is the tidy
#'     accessor -- reaching into a result with `$` is never required, and a
#'     secondary table is selected with `what =`.}
#' }
#'
#' @section Crossing between families:
#'
#' The families are entry points into one another, not islands. The text
#' family is the corpus front end of the hypergraph family: a
#' `text_hypergraph` *is* a `net_hypergraph`, so every hypergraph verb takes
#' it, and every `hg_*()` verb takes any `net_hypergraph` in return.
#' [build_simplicial()] with `type = "pathway"` turns a memory network into a
#' simplicial complex; [window_hypergraph()] turns the same sequences into a
#' hypergraph; [clique_expansion()] projects a hypergraph back to a pairwise
#' network that any of the first-order tools accept. [pathways()] is the
#' shared accessor that hands sequence-derived path strings from the memory
#' family to the other two.
#'
#' @references
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
#' Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
#'
#' Bianconi, G. (2021). *Higher-Order Networks*. Cambridge University Press.
#'
#' @keywords internal
#' @importFrom cograph plot_simplicial splot
#' @importFrom stats setNames
#' @importFrom utils head
#' @importFrom ggplot2 .data
"_PACKAGE"

## ggplot2 non-standard-evaluation column names used in plot methods.
## New plot code uses the .data pronoun instead; these cover the bare-name
## aes() calls inherited with the simplicial family.
utils::globalVariables(c(
  # memory family
  "KL", "context", "label", "x", "xend", "y", "yend",
  "top_o1", "top_ok",
  # simplicial family
  "betti", "birth", "components", "count", "death", "dim", "k", "node",
  "dim_label", "persistence", "q", "t", "threshold", "total", "value"
))
