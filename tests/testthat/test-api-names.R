# hypernets imports Nestimate (ROADMAP Phase 0b), so both namespaces load in
# every session. Nestimate keeps its hypergraph verbs and the net_hypergraph
# class; hypernets' hypergraph verbs are hg_* and its class is net_hg. These
# guards fail if an export or an S3 registration reintroduces a collision.

nestimate_hypergraph_names <- c(
  "build_hypergraph", "bipartite_groups", "hypergraph_centrality",
  "hypergraph_measures", "hypergraph_laplacian", "hypergraph_cluster",
  "hypergraph_transduction", "clique_expansion"
)

test_that("no export reuses a Nestimate hypergraph verb name", {
  exports <- getNamespaceExports("hypernets")
  expect_length(intersect(exports, nestimate_hypergraph_names), 0L)
  expect_length(grep("^hypergraph_", exports, value = TRUE), 0L)
})

test_that("no S3 method is registered on Nestimate's hypergraph classes", {
  methods <- getNamespaceInfo("hypernets", "S3methods")[, 3L]
  expect_length(grep("\\.net_hypergraph", methods, value = TRUE), 0L)
  expect_length(grep("\\.hypergraph_measures$", methods, value = TRUE), 0L)
})

test_that("hypergraph constructors return net_hg", {
  hg <- group_hypergraph(data.frame(p = c("a", "b", "c", "a"),
                                    g = c("x", "x", "y", "y")),
                         actor = "p", group = "g")
  expect_s3_class(hg, "net_hg")
  expect_false(inherits(hg, "net_hypergraph"))
})

test_that("text_hypergat is the text-family name of hg_hypergat", {
  expect_identical(text_hypergat, hg_hypergat)
})

test_that("no hypernets export shares a name with a Nestimate export", {
  ours <- getNamespaceExports("hypernets")
  theirs <- getNamespaceExports("Nestimate")
  expect_length(intersect(ours, theirs), 0L)
})

test_that("constructors are nouns and no build_* or hon_* export remains", {
  exports <- getNamespaceExports("hypernets")
  expect_length(grep("^build_", exports, value = TRUE), 0L)
  expect_length(grep("^hon_", exports, value = TRUE), 0L)
  expect_true(all(c("hon", "honem", "mogen", "simplicial", "markov_order",
                    "memory", "hg_get") %in%
                    exports))
  gone <- c("bootstrap_hon", "compare_hon", "hon_centrality",
            "hon_communities", "wasserstein_distance", "markov_order_test",
            "markov_stability", "path_dependence", "persistent_homology",
            "persistence_landscape", "bottleneck_distance", "betti_numbers",
            "euler_characteristic", "q_analysis", "simplicial_degree",
            "verify_simplicial", "mogen_transitions", "path_counts",
            "pathways")
  expect_length(intersect(exports, gone), 0L)
})
