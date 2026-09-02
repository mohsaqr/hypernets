test_that("long and short hypergraph API names are direct aliases", {
  aliases <- list(
    hypergraph_agreement = hg_agreement,
    hypergraph_classify = hg_classify,
    hypergraph_communities = hg_communities,
    hypergraph_community_quality = hg_community_quality,
    hypergraph_edge_centrality = hg_edge_centrality,
    hypergraph_edges = hg_edges,
    hypergraph_embed = hg_embed,
    hypergraph_keywords = hg_keywords,
    hypergraph_line_graph = hg_line_graph,
    hypergraph_motifs = hg_motifs,
    hypergraph_neural = hg_neural,
    hypergraph_hypergcn = hg_hypergcn,
    hypergraph_hnhn = hg_hnhn,
    hypergraph_allset = hg_allset,
    hypergraph_null_test = hg_null_test,
    hypergraph_pagerank = hg_pagerank,
    hypergraph_project = hg_project,
    hypergraph_snapshot = hg_snapshot,
    hypergraph_snapshots = hg_snapshots,
    hypergraph_seeds = hg_seeds,
    hypergraph_stability = hg_stability,
    text_hypergat = hg_hypergat
  )
  short <- list(
    hg_agreement, hg_classify, hg_communities, hg_community_quality,
    hg_edge_centrality, hg_edges, hg_embed, hg_keywords, hg_line_graph, hg_motifs,
    hg_neural, hg_hypergcn, hg_hnhn, hg_allset, hg_null_test, hg_pagerank, hg_project, hg_snapshot, hg_snapshots, hg_seeds,
    hg_stability, hg_hypergat
  )

  identical_alias <- Map(identical, aliases, short)
  expect_true(all(unlist(identical_alias)),
              info = paste(names(aliases)[!unlist(identical_alias)],
                           collapse = ", "))
})
