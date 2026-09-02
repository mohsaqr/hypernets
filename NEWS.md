# honets 0.3.0

* **The complete Hayashi clustering family is implemented.**
  `hypergraph_cluster(algorithm = "symnmf")` adds Algorithm 2, RDC-Sym,
  alongside the existing RDC-Spec implementation. New
  `hypergraph_joint_cluster()` implements the patent experiment's J-NMF
  (Eq. 18) and JS-NMF (Eq. 19) objectives when an auxiliary node-relation
  matrix is available. All NMF paths expose their objective trace,
  convergence state, iteration count, chosen restart and fitted factors.
* **Full HyperGAT semantic hyperedges.** `hg_hypergat(semantic = "lda")`
  now implements the Ding et al. (2020) LDA path in native R: online
  variational Bayes on labeled training documents, class-count topics, and
  per-document topic edges from the top words. Precomputed keyword lists are
  accepted for exact replay of official preprocessing artifacts; the default
  remains the backward-compatible sentence-only ablation.
* **The neural paper set is complete.** New `hypergraph_hypergcn()` provides
  dynamic, fast and one-edge HyperGCN; `hypergraph_hnhn()` exposes the HNHN
  alpha/beta normalization exponents; and `hypergraph_allset()` implements
  both AllDeepSets and AllSetTransformer. `heterogeneous_hgat()` implements
  Linmei et al.'s distinct heterogeneous document/topic/entity graph model
  with node- and type-level attention. Each method has direct equation or
  construction invariants plus end-to-end torch tests.

* **cograph is the shared graph and plotting engine.** It is promoted from
  Suggests to Imports: honets owns hypergraph construction, incidence algebra
  and hypergraph-specific transformations, then uses cograph for ordinary
  graph algorithms and rendering. Dynet remains a separate peer and is not a
  dependency.
* **Dual hypergraph verb names.** Descriptive names such as
  `hypergraph_edges()`, `hypergraph_project()`, `hypergraph_pagerank()` and
  `hypergraph_classify()` are exported as direct bindings to their compact
  `hg_*()` forms; no implementation is duplicated. The raw-text attention
  model is also available as `text_hypergat()`. Existing engine names
  `hypergraph_measures()`, `hypergraph_centrality()` and
  `hypergraph_cluster()` keep their established meanings.
* **The full *Legal hypergraphs* method layer is implemented.** New temporal
  hypergraphs and snapshots; binary/multi and self-association projections;
  s-betweenness and s-closeness of hyperedges; log-subhypergraph centrality;
  induced Y/T/O motif censuses with the HypergraphX configuration MCMC;
  repeated Infomap with AMI-medoid selection; ARI/AMI/NMI agreement; and
  coverage, weighted coverage, performance and modularity. The authors'
  released ICSID data reproduce the published aggregate motif counts exactly:
  Y = 478, T = 7, O = 0.
* **The full Legal Hypergraphs archive reproduction now passes.** A standalone
  runner rebuilds the sparse 3,618-node GFCC citation-block hypergraph,
  verifies all four association projections and all eight archived Figure 8
  medoids, and reproduces the active/aggregate ICSID Figure 6 centrality
  rankings and values. `group_hypergraph()` now accepts an explicit node
  universe and sparse incidence, and normalized hyperedge centrality follows
  the exact NetworkX/HypergraphX convention on disconnected line graphs.

## The text family: texthypergraph folds into honets

The `texthypergraph` package is retired and its whole surface now lives here,
with its history of tests intact (491 shipped expectations came across; the
merged suite passes in full). honets gains a fourth family, **text
hypergraphs**, and the hypergraph family gains every generic method that
texthypergraph had built under its frozen-Nestimate contract.

* **New family — text hypergraphs.** `text_hypergraph()` (bag of words with
  smoothed tf-idf, token windows, or embedding kNN; dense or sparse),
  `stop_words_en()`, and the tidy verbs `hg_measures()`, `hg_centrality()`,
  `hg_cluster()`, `hg_keywords()`, `hg_classify()`, `hg_stability()`,
  `hg_agreement()`, `hg_seeds()`. Data: `covid_abstracts`,
  `covid_embeddings`. Vignettes `text-hypergraphs` and `text-constructions`.
* **Hypergraph family additions.** `knn_hypergraph()`, `dual_hypergraph()`,
  `hg_pagerank()` (personalized EDVW PageRank, sparse-capable),
  `hg_project()` (clique and association weightings), `hg_line_graph()`,
  `hg_edges()`, `hg_null_test()` (swap and configuration nulls), and the
  neural tier `hg_neural()` (HGNN) and `hg_hypergat()` (HyperGAT) in native
  torch. Sparse `Matrix::dgCMatrix` incidences are supported end to end by
  `hg_cluster()`, `hg_classify()`, `hg_pagerank()` and `hg_measures()`.
* **`hypergraph_transduction(normalization = )`.** New argument, default
  `"none"` (the raw Zhou 2006 argmax, unchanged behaviour). `"class_mass"`
  divides each class column by its total spread mass before the argmax (Zhu,
  Ghahramani & Lafferty 2003); without it, class-imbalanced seeds collapse
  every prediction onto the majority class (on R8, every test document is
  predicted "earn"). The result object records `$normalization`.
* **Nestimate dependency removed.** texthypergraph delegated its incidence
  construction and measures to Nestimate; those engines already lived here,
  and `group_hypergraph()` is `identical()` to `Nestimate::bipartite_groups()`
  up to the `member`/`player` argument name (tested). The package now
  imports only cograph, ggplot2, graphics, grid, Matrix, methods, parallel,
  RSpectra, stats and utils.
* **Collisions resolved without a value change.** texthypergraph carried a
  verbatim copy of the spectral trio; honets' copies were kept (they carry
  the plot methods, `top =`, and scalar `edge_weights`), and only the
  normalization argument was ported. `hg_pagerank()` agrees with
  `hypergraph_centrality(type = "pagerank")` to `1e-10` (tested);
  `hg_project(method = "clique")` equals `clique_expansion()` (tested);
  `text_hypergraph(construction = "window")` and `window_hypergraph()`
  produce the same incidence matrix on their shared domain (tested).
* **Condition classes** of the incoming code are `honets_*`
  (`honets_bad_input`, `honets_no_converge`,
  `honets_hypergraph_disconnected`, `honets_empty_corpus`,
  `honets_dropped_documents`, `honets_missing_embeddings`,
  `honets_missing_torch`, `honets_nonpositive_similarity`,
  `honets_sparse_unsupported`, `honets_sparse_too_large`,
  `honets_configuration_collapse`).
* **Infrastructure.** GitHub Actions (`R-CMD-check`, `pkgdown`), a pkgdown
  reference index covering every topic, `VignetteBuilder: knitr`, the
  benchmark harness under `benchmarks/` (build-ignored), and the
  text-family equivalence suites (HyperNetX, HyperG, DHG, the official
  HyperGAT code) under `local_testing_and_equivalence/`.

# honets 0.2.1

## Every accessor takes `top =`

Every verb that can return many rows now takes a `top =` argument, so a
caller never has to write `head()` around a result:

```r
as.data.frame(mo, what = "transitions", order = 2, top = 4)   # not head(..., 4)
```

`top` is applied **last** - after `what`, after every filter (`order_min`,
`min_count`, `dim`, `k`, `dimension`, `significant`), and after `sort_by` -
so `sort_by` and `top` compose: `top = n` is the first `n` rows of the table
as ordered. `top = NULL` (the default) returns everything, and the default
return of every accessor is unchanged. Semantics match the `top` that already
shipped on `path_counts()` and `pathways()`.

Added to: `as.data.frame()` for `net_hon`, `net_honem`, `net_hypa`,
`net_mogen`, `net_markov_order`, `net_path_dependence`, `net_hon_boot`,
`net_hon_compare`, `net_simplicial`, `net_q_analysis`,
`net_persistent_homology`, `net_persistence_landscape`, `net_hypergraph`,
`net_hypergraph_measures`, `net_hypergraph_cluster`,
`net_hypergraph_transduction`; and to `mogen_transitions()`,
`hon_centrality()`, `simplicial_degree()`, `hypergraph_centrality()`.

`path_counts(top =)` now validates its argument like the rest of the family:
a non-whole value such as `top = 2.5` is an error rather than a silent
truncation to 2, matching how `k` already behaved.

## Bug fix

* `plot.net_path_dependence()` clipped the label of its highest-KL context -
  the row the plot exists to show. The modal-flip labels are drawn to the
  right of each point, and the panel did not extend past the largest value,
  so ggplot cut the label off. The x scale now leaves room for it.

## Documentation

`docs/` (build-ignored) is reorganised from 23 per-verb vignettes into **four
documents**, one per structure family plus an overview:
`overview`, `memory-networks`, `simplicial-complexes`, `hypergraphs`. Each
section is one verb, worked end to end on the same data; the per-verb sources
are archived under `docs/_sections/`.

They gain **37 figures**, drawn with cograph. honets results are dual-classed
`cograph_network`, so `cograph::splot()` and `cograph::plot_simplicial()`
take them with no conversion step. `cograph` is added to `Suggests` and every
plot chunk is guarded on it.

The prose was also rewritten to remove 37 `head()` and 12 `subset()` calls
that subset a returned table on the public surface - the idiom `top =` now
replaces.

# honets 0.2.0

honets becomes **the** higher-order networks package: one package covering
all three structure families of the higher-order literature (Battiston et al.
2020) — memory networks, simplicial complexes, and hypergraphs — under a
single taxonomy.

## Absorbed families

* **Simplicial complexes and topological data analysis**, moved verbatim
  from Nestimate 0.9.0: `build_simplicial()` (clique, Vietoris-Rips and
  pathway complexes), `betti_numbers()`, `euler_characteristic()`,
  `persistent_homology()`, `persistence_landscape()`,
  `bottleneck_distance()`, `simplicial_degree()`, `q_analysis()`,
  `verify_simplicial()`.
* **Hypergraphs**, moved verbatim from Nestimate 0.9.0 by way of the
  short-lived `hypernets` package (0.1.2, never released), which is folded in
  and retired: `build_hypergraph()`, `window_hypergraph()`,
  `group_hypergraph()`, `hypergraph_measures()`, `hypergraph_centrality()`,
  `hypergraph_laplacian()`, `hypergraph_cluster()`,
  `hypergraph_transduction()`, `clique_expansion()`.
* The consolidation **deletes 223 lines of duplication** from hypernets'
  323-line `utils.R` -- 171 of them the clique-enumeration closure copied out
  of Nestimate's `simplicial.R`, the rest a `build_simplicial()` shim and a
  second copy of honets' own `.extract_edges_from_matrix()`. hypernets had to
  carry that closure because `build_hypergraph()` needs clique enumeration;
  with both families in one package, `build_hypergraph()` calls the real
  `build_simplicial()` again. The seven copied helpers were verified
  byte-identical to Nestimate's originals before the copy was removed.

## Taxonomy: renames

All renames are to the **surface only**. No computed value changed anywhere —
the `identical()` contracts against Nestimate 0.9.0 still hold for both
absorbed families, with only these names normalised away
(`local_testing_and_equivalence/test-identity-nestimate-{hypergraph,simplicial}.R`).

Every result class now carries a `net_*` class:

| Was | Is |
|---|---|
| `simplicial_complex` | `net_simplicial` |
| `persistent_homology` (class) | `net_persistent_homology` |
| `q_analysis` (class) | `net_q_analysis` |
| `persistence_landscape` (class) | `net_persistence_landscape` |
| `hypergraph_measures` (class) | `net_hypergraph_measures` |

Constructors and arguments:

| Was | Is | Why |
|---|---|---|
| `bipartite_groups()` | `group_hypergraph()` | it returns a hypergraph, not bipartite groups; now matches `window_hypergraph()` |
| `bipartite_groups(player =)` | `group_hypergraph(member =)` | the column names hypergraph nodes, which need not be people |
| `build_hypergraph(method =)` | `build_hypergraph(type =)` | same construction axis as `build_simplicial(type =)` |
| `params$method` | `params$type`, plus `params$source` | every hypergraph constructor now records `source`, so a result says how it was built |

Classed conditions are now uniformly `honets_*`:
`hypernets_no_converge` and `nestimate_hypergraph_disconnected` became
`honets_no_converge` (shared with the memory family's power iteration) and
`honets_hypergraph_disconnected`.

## Taxonomy: complete tidy-accessor coverage

Every `net_*` result class now has an `as.data.frame()` method, so no result
object requires reaching in with `$`. Eleven are new, each with a `what =`
argument for its secondary table:

* `net_hon` (`"rules"` / `"nodes"`, plus `order_min` and `sort_by`)
* `net_honem` (`"embeddings"` / `"variance"`)
* `net_hypa` (`"scores"` / `"over"` / `"under"`, plus `sort_by`)
* `net_mogen` (`"orders"` / `"transitions"`)
* `net_markov_order` (`"orders"` / `"null"`)
* `net_path_dependence` (plus `min_count`, `sort_by`)
* `net_simplicial` (`"simplices"` / `"f_vector"`, plus `dim`)
* `net_q_analysis` (`"q_levels"` / `"nodes"`)
* `net_persistent_homology` (`"persistence"` / `"betti"`, plus `dimension`,
  `sort_by`)
* `net_persistence_landscape` (plus `k`)
* `net_hypergraph_measures` (`"nodes"` / `"edges"` / `"global"`, plus
  `sort_by`)

A regression test asserts the coverage, so a future `net_*` class without an
accessor fails the suite.

## Structure

* `R/` is organised by family: `memory_*.R` (8 files), `simplicial_*.R` (4),
  `hypergraph_*.R` (7), plus shared `utils.R`, `pathways.R`, `data.R` and the
  package doc. Test files follow the same names.
* Nestimate's 1,561-line `simplicial.R` was split along its real dependency
  seams into `simplicial.R` (construction and structural measures),
  `simplicial_filtration.R` (the filtration and Z/2 boundary-reduction layer
  that both `build_simplicial(type = "vr")` and `persistent_homology()` sit
  on) and `simplicial_homology.R`. Code unchanged; the split was verified
  line-for-line content-preserving.

## Other

* `group_hypergraph()`'s weighted branch and `build_hypergraph()`'s incidence
  fill are vectorised (they were `for` loops accumulating into a matrix).
  Duplicate `(member, group)` cells are summed before assignment, which
  index assignment alone would not do.
* Package-level documentation (`?honets`) now states the three-family
  taxonomy, the verb grammar, and how the families cross into one another.

# honets 0.1.5

* New verb `hon_centrality()` (roadmap item A2): PageRank, betweenness
  and closeness computed on the higher-order topology and projected back
  onto first-order states (Scholtes, Wider & Garas 2016). Semantics
  follow pathpy 2.2.0 - that paper's reference implementation -
  generalized from fixed-order to the variable-order networks
  `build_hon()` produces, and verified against it: betweenness and
  closeness match exactly (< 1e-10) on second- and third-order
  topologies, PageRank to pathpy's own `tol = 1e-6`. The underlying
  kernels additionally match `igraph` on the same topologies.
  `project = FALSE` reports the centralities of the memory contexts
  themselves; `projection =` chooses how a higher-order node's PageRank
  is distributed ("scaled", "last", "first", "all"); `sort_by =` returns
  the table ranked.
* `build_hon()` gained the long-format interface already used by the
  inference verbs: `action`, `actor` and `time` column names, so an
  event table needs no manual splitting. Existing calls are unaffected
  (the arguments default to `NULL`) and produce byte-identical networks.
* New tutorial `Tutorial_docs/hon_centrality.html`: why the first-order
  network of a dense corpus cannot rank its states at all (complete
  digraph, uniform PageRank), what the higher-order ranking recovers,
  which contexts carry the flow, and when betweenness and closeness are
  saturated.

# honets 0.1.4

* New inference verbs for higher-order rules (roadmap item A1):
  `bootstrap_hon()` — sequence bootstrap with percentile CIs for rule
  probabilities and per-rule extraction *support*; `compare_hon()` —
  two-sample permutation comparison with per-edge BH adjustment and a
  pooled-count-weighted global test. Both precompute per-sequence counts
  once and rebuild replicates from reweighted counts (proven identical
  to re-counting the resampled multiset), draw all randomness serially
  (parallel runs reproduce serial results under a seed), accept long
  format (`action`/`actor`/`time`), and ship with `print`/`summary`/
  `plot`/`as.data.frame` methods (accessor filters `min_support`,
  `order_min`, `significant`, `sort_by`).
* Bundled example data `human_long` and `ai_long` (coded human-AI pair
  programming sessions, long format).
* New tutorial `Tutorial_docs/hon_inference.html`: rule stability by
  order, support-filtered reporting, and a diffuse early-vs-late shift
  (significant global test, no significant single edge).
* Internal: `.hon_extract_rules()` split into a counts-based core
  (`.hon_extract_rules_count()`); behavior unchanged (full equivalence
  suite re-verified).

# honets 0.1.3

* Roadmap: hypernets B2 (EDVW hypergraph PageRank) marked done in
  `EXPANSION-PLAN.md`. No package code changed.

# honets 0.1.2

* Roadmap: hypernets B1 (windowed sequence hyperedges) marked done in
  `EXPANSION-PLAN.md` with the shipped design recorded. No package code
  changed.

# honets 0.1.1

* Added the consolidated family expansion roadmap (`EXPANSION-PLAN.md`,
  build-ignored): honets higher-order features A1–A4 and the hypernets
  hypergraph sibling (scaffolded 2026-08-25). No package code changed.

# honets 0.1.0

* Initial release. Code moved from Nestimate 0.9.0 (delegation T0): `build_hon()`,
  `build_honem()`, `build_hypa()`, `build_mogen()`, `mogen_transitions()`,
  `path_counts()`, `markov_order_test()`, `path_dependence()`, and the
  `pathways()` generic with methods for `net_hon`, `net_hypa`, and `net_mogen`.
  Numbers are identical to the Nestimate implementations (same code, same RNG
  streams).
* Corrected the HONEM reference (Saebi, Ciampaglia, Kaplan & Chawla 2020,
  \doi{10.1089/big.2019.0169}); the author list previously cited was wrong.
