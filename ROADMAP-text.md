# hypernets text and hypergraph roadmap

This file preserves the release history of the retired `texthypergraph`
package and records its continuation inside hypernets. Sections v0.1-v0.6 are
historical implementation notes; v0.7 and later describe the maintained
hypernets package.

**What the package is:** hypergraph text analysis in R — corpus in, weighted
text hypergraph out, spectral/statistical analysis through tidy verbs. First
of its kind on CRAN (gap verified 2026-08-24: no hypergraph text
classification exists in R at all, not even non-neural Zhou 2006).

**Current ownership contract (no duplication):**

- hypernets owns all higher-order and text-hypergraph construction, analysis,
  inference and result classes. The retired package is an archive only.
- Nestimate is neither an Import nor a Suggest. Historical identity checks
  remain local equivalence tests; no hypernets runtime path calls Nestimate.
- cograph is the imported graph and plotting engine at the graph-shaped
  boundary. Dynet remains a separate temporal-network peer.
- `sbert` (Suggests) is the native embedding front-end; every verb also
  accepts a precomputed `embeddings` matrix so the package runs offline.
- Oracles (local-only, never dependencies): HyperNetX (EDVW Laplacian,
  verified author-adjacent), XGI (tensor centralities), HyperG (unweighted
  constructions/samplers), gudhi (wasserstein). Full map: `repos/README.md`.

**Filter (house rules):** deterministic, base R + the two family Imports,
tidy one-verb APIs, every method oracle-tested before release, effect sizes +
CIs in vignettes, multi-seed sensitivity for every arbitrary choice (k,
window, threshold).

## v0.1 — Scaffold + the bridge (no new math)

**Status 2026-08-25: v0.1 COMPLETE.** Skeleton, `text_hypergraph()`
(+ curly-apostrophe normalization), `stop_words_en()`, four delegating verbs
(+ `sort_by`/`n` on centrality), `covid_abstracts` dataset (165 abstracts,
2020-2024), full worked vignette (clustering + k-sensitivity + few-label
classification with verified numbers), package README. 74 tests pass, engine
parity `identical()`, R CMD check 0 errors / 0 warnings / 1 environmental
NOTE.

- Package skeleton: DESCRIPTION, testthat 3e, pkgdown, CI matrix;
  papers/repos/TODO/ROADMAP `.Rbuildignore`d.
- `text_hypergraph(x, column, weight = c("n", "tfidf"))` — document–word
  incidence from raw text or any long table (base-R tokenization; quanteda
  dfm / tidytext tables accepted as inputs, never dependencies).
- Thin tidy verbs delegating to Nestimate: `hg_measures()`, `hg_centrality()`,
  `hg_cluster()`, `hg_classify()` (transduction), `hg_embed()` — one call,
  named arguments, tidy data.frame returns, print/summary/plot +
  `as.data.frame()` per house rules.
- Bundled corpus + the bridge vignette (the Phase-1 vignette, now living
  here instead of Nestimate).
- Gate: `--as-cran` clean; parity of every delegated verb against direct
  Nestimate calls (`identical()`).

## v0.2 — Text-native constructions

**Status 2026-08-25: v0.2 COMPLETE.** Both constructions shipped:
`text_hypergraph(construction = "window")` (sliding/tumbling, set-valued
window hyperedges; w = 2 off-diagonal parity with `Nestimate::wtna()`
shipped as a package test — diagonals differ by design, wtna counts
within-window repeats) and `construction = "knn"` + `knn_hypergraph()`
(cosine-weighted, deterministic tie-breaks, classed refusal of non-positive
similarities; binary support verified identical to `HyperG::knn_hypergraph`
in local_testing_and_equivalence/). `covid_embeddings` bundled (sbert
all-MiniLM-L6-v2, 165x384). `vignette("constructions")` demonstrates both
on the real corpus with w and k sensitivity checks. 136 tests; check
0 / 0 / 1 environmental NOTE. Version 0.2.0.

- **Windowed sequence hyperedges** — `text_hypergraph(construction =
  "window", w = )`: the HyperGAT construction, weighted/windowed extension.
  No upstream oracle; invariant-gated (window-count conservation; `w = 2`
  reduces exactly to `wtna_cooccurrence`).
- **kNN embedding hypergraph** — `knn_hypergraph(embeddings, k, weight =
  "cosine")` + `text_hypergraph(construction = "knn", model = )` calling
  sbert when installed. Construction oracle: `HyperG::knn_hypergraph`
  (unweighted). Vignette-first before the API freezes; k as a sensitivity
  axis, encoder-dependence reported per model.
- Gate: construction oracles/invariants + a real-corpus vignette for each.

## v0.3 — Theory completions (implemented HERE; Nestimate frozen)

**Status 2026-08-25: v0.3 COMPLETE** (hg_pagerank + dual_hypergraph +
hg_null_test; 180 tests; check 0/0/1 environmental NOTE; version 0.3.0).

- [x] DONE 2026-08-25: **`hg_pagerank()`** (Chitra & Raphael 2019 EDVW walk
  + damping/personalization). Four verification layers: direct linear-solve
  reference; collapse-theorem closed form (1e-12); Nestimate `$pi` parity
  (1e-12); HyperNetX `get_pi(prob_trans(weights=TRUE))` parity (1.4e-16).
- [x] DONE 2026-08-25: **`dual_hypergraph()`** — transpose identity,
  involution, and dual == opposite-orientation-construction tests shipped;
  HyperG::dual_hypergraph support parity PASS (local).
- [x] DONE 2026-08-25: **`hg_null_test()`** — degree-preserving
  checkerboard null (Gotelli 2000), permutation p (Phipson & Smyth 2010),
  null quantiles + z; margin-conservation and blocky-structure detection
  tests shipped; avg_edge_size invariance as built-in falsification.
- Gate: same discipline — hand-computed fixtures, invariants, oracle
  scripts, real-corpus vignette.

## v0.4 — Scale: sparse core (MANDATE 2026-08-25)

**Status 2026-08-25: sparse core SHIPPED** (bag construction, PageRank,
transduction via CG, RSpectra spectral clustering, measures, dual; window/
knn/tensor-centrality/null-test sparse paths still dense-only, guarded by
classed errors). Scale demo: 20,000 docs x 17,576 words — build 12.5s,
classify 0.6s (accuracy 1.000 on planted blocks), cluster 0.3s (perfect
partition), 14.4 MB sparse vs 2.8 GB dense equivalent. 291 tests; check
0/0/0. Version 0.4.0. The 20NG-class benchmarks (v0.5) are now feasible.

User verdict on v0.1-0.3: correct foundation, toy scale. The dense incidence
caps the package at hundreds of documents while the literature base
benchmarks on 7k-18k (Ohsumed, 20NG). v0.4 removes the ceiling:

- Sparse incidence (`Matrix::dgCMatrix`) as a first-class representation
  (`text_hypergraph(sparse = TRUE)`), same construction semantics.
- Sparse engines: EDVW transition + PageRank (sparse mat-vec), transduction
  via conjugate gradient on the operator (never materializing `solve`),
  spectral clustering via `RSpectra` partial eigenpairs on the similarity
  operator, structural measures on sparse cross-products.
- **Oracle: the shipped dense engines themselves** (Nestimate/HyperNetX-
  parity-verified) — sparse and dense must agree to 1e-8 on every small
  corpus; scale gate: 20k x 50k corpora construct + classify in minutes.

## v0.5 — Benchmarks: the papers' own yardstick

**Status 2026-08-25: SHIPPED.** Harness in `benchmarks/` (TextGCN
corpora/splits — the exact files HyperGAT preprocesses; split sizes match
Ding Table 1 on all five datasets, vocab exactly on R8/R52). Results
(`benchmarks/RESULTS.md`, pkgdown article `vignettes/articles/`):
20NG 0.8477 and MR 0.7684 beat the tf-idf centroid (0.7796 / 0.6851) and
fastText (published 0.7938); MR edges published TextGCN-transductive
(0.7674); skewed corpora (R8/R52/Ohsumed) go to the centroid — honest
both ways. Found + fixed a real method defect: raw Zhou argmax collapses
onto the majority class under imbalanced seeds (R8 0.4947 = the
majority-class rate) → new `normalization = "class_mass"` argument
(Zhu et al. 2003 CMN) on `hg_classify()`/`hypergraph_transduction()`,
fixture + invariance + parity + mutation tested. Low-label study
(stratified 1-20%, 5 draws): transduction leads the centroid at every
fraction on MR only. Version 0.5.0.

- R8 / R52 / MR / Ohsumed / 20NG harness (HyperGAT's Table 2 datasets);
  train/test splits as published.
- Zhou/EDVW transduction + tf-idf baselines vs the published accuracy
  tables; accuracy/F1 with bootstrap CIs; honest reporting either way.
- Results as a pkgdown article; harness reused by v0.6.

## v0.6 — Neural tier: hypergraph GNNs natively in R ({torch}, Suggests)

**Status 2026-08-26 (2): HyperGAT SHIPPED — v0.6 complete** (`hg_hypergat()`,
version 0.6.4). Document-level hypergraphs (unique words as vertices,
sentences as hyperedges), dual attention mirroring the official code
(kaize0409): forward parity 2.4e-7 with fixed weights; independent
plain-R layer-math oracle in the tests; word-order invariance;
mutation-checked. Sentence-hyperedge and full LDA semantic-hyperedge paths
are implemented; the original benchmark rows remain sentence-only and their
comparator is the paper's "w/o semantic" ablation. No WordNet lemmatizer
(documented). Benchmarks (official
corpus files, official protocol, 3 seeds): **R8 0.9665 +/- 0.003** (1.3
pts from full published HyperGAT, at published TextGCN level), **R52
0.9433 +/- 0.002 — above published transductive TextGCN (0.9356)**;
~2-4 min per training run on CPU. Package-best on both.

**Status 2026-08-26: HGNN SHIPPED** (`hg_neural()`, torch in Suggests,
version 0.6.0). Verified four ways: propagation == Zhou operator from the
in-package spectral core (1e-12); hand-computed weighted factorization;
forward parity vs the official DHG `HGNNConv` at 2.4e-7
(local_testing_and_equivalence/test-equiv-hgnn-dhg.R — note DHG
deduplicates identical hyperedges); mutation check. Benchmarks (3 seeds):
R8 0.9539, R52 0.8440, MR 0.7692 — package-best on all three; Ohsumed and
20NG stay below closed-form transduction at every configuration tested
(corpus-level doc-node design oversmooths; see the pkgdown article).
HyperGAT (document-level hypergraphs + dual attention) is the open half
of this stage and the expected answer to the 20NG/Ohsumed gap.

The no-neural rule is LIFTED (user decision 2026-08-25). First hypergraph
GNNs in R, no Python:

- HGNN (Feng 2019): spectral convolution on the Zhou Laplacian — the
  simplest layer, first.
- HyperGAT (Ding 2020): dual node-level/edge-level attention.
- Equivalence discipline for neural code: forward-pass parity against the
  official PyTorch implementations with fixed weights (reticulate,
  local_testing_and_equivalence/); training curves + benchmark accuracy
  vs the papers' tables via the v0.5 harness.

## v0.7 — Analytics depth: the projection tier

**Status 2026-09-02: IMPLEMENTED.** Scope set after reading [Coupette, Hartung &
Katz (2024), "Legal hypergraphs"](https://doi.org/10.1098/rsta.2023.0141),
Phil. Trans. R. Soc. A 382(2270):20230141. Hyperedge measures, projections,
s-line graphs, temporal snapshots, centralities, motifs, communities and
configuration nulls are shipped and tested against the paper's oracle code.

The open-access paper is archived in the source-tree
[research library](papers/2024-PhilTrans-LegalHypergraphs-Coupette.pdf), with
a [method and reproducibility note](repos/legal-hypergraphs.md).

**The gap this closes.** Every construction shipped through v0.6 derives its
hyperedges from *within-document word co-occurrence* (`bag` / `window` /
`knn`), and every verb reaches its answer through the incidence matrix and a
Laplacian eigen-decomposition. The package therefore has no combinatorial
tier: no hyperedge-level measures, no principled projection, no paths, no
communities beyond spectral. That followed from the v0.5/v0.6 benchmark
trajectory — the TextGCN corpora ship as flat (label, text) pairs carrying no
metadata, so nothing motivated a construction that uses document relations.
It was not a scoping decision: "Out of scope" below excludes none of it.

**Package-boundary contract (settled 2026-09-01; clarified 2026-09-02).**

```
hypernets           owns the hypergraph layer: construction, incidence and
                 spectral algebra, projections, hypergraph nulls and motifs,
                 plus any temporal hypergraph representation

cograph          imported graph and plotting engine; once hypernets produces a
                 projection or s-line graph, cograph owns ordinary paths,
                 centralities, communities, comparisons and rendering

Dynet            separate temporal-network peer; never an hypernets dependency
                 or computational layer. Similar start/end/step/window names
                 are interoperability, not delegation
```

**Dependency decision.** Dynet is absent from both Imports and Suggests.
`cograph` is an explicit Import: it is already the rendering contract for
hypernets objects, so making the graph algorithms optional would add guards
without delivering meaningful independence. Hypergraph-specific mathematics
stays in hypernets; cograph begins at the graph-shaped boundary.

### Items — each names its equivalence oracle before implementation

- [x] **DONE 2026-09-01 `hg_project()`** — hypergraph to weighted graph. `method = "clique"`
  (unnormalised co-occurrence, `identical()` parity with
  `Nestimate::clique_expansion()`) and `method = "association"`, the
  random-walk-correct `w({u,v}) = sum_e 1/(|e|-1)` of Coupette et al. (2024).
  Verified 2026-09-01: `clique_expansion()` is `tcrossprod(incidence)` with a
  zeroed diagonal and applies no `|e|` normalisation, so the association
  weighting exists nowhere in the stack. *Oracle*: the paper's Zenodo
  reproducibility package (doi:10.5281/zenodo.8081507). *Invariant*: each
  hyperedge contributes total weight exactly 1 to the edges incident on each
  of its own members.
- [x] **DONE 2026-09-01 `hg_line_graph(s = 1)`** — the s-line graph (hyperedges as vertices,
  adjacent when `|e_i & e_j| >= s`). `dual_hypergraph()` is the s = 1 case;
  this generalises it. *Oracle*: HyperNetX `s_line_graph`.
- [x] **DONE 2026-09-01 `hg_edges()`** — the missing hyperedge-level tier: cardinality,
  degree, node-neighbourhood size, one row per hyperedge, with
  `what = "distribution"` for degree/cardinality CCDFs. Everything the
  package returns today is keyed on vertices. *Oracle*: XGI edge-stat suite.
- [x] **DONE 2026-09-02 `hg_edge_centrality()`** — s-betweenness and s-closeness via
  `hg_line_graph()` handed to cograph's statnet-verified graph kernels. hypernets
  owns the s-line transformation; cograph owns the ordinary shortest-path
  centrality after that boundary.
- [x] **DONE 2026-09-01 `hg_null_test(method = "configuration")` plus
  DONE 2026-09-02 `hg_motifs()`** — the hypergraph
  configuration model beside the shipped degree-preserving checkerboard.
  configuration null plus induced Y/T/O counts, relative abundance and the
  normalized motif profile.
  *Oracle*: XGI Chung-Lu / configuration samplers. *Invariant*: expected
  degree sequence preserved.
- [x] **DONE 2026-09-02 `hg_communities()` and
  `hg_community_quality()`** — Infomap over
  `hg_project(method = "association")`, including the paper's optional
  self-association weights, repeated seeded fits and medoid selection. hypernets
  builds the paper-specific association graph; cograph performs the ordinary
  graph clustering and renders the result.
- [x] **DONE 2026-09-02 `hg_agreement(method = "ami")`** — adjusted mutual
  information and NMI beside the shipped ARI. *Oracle*: scikit-learn; exact
  fixture parity to 1e-12.
- [x] **DONE 2026-09-02 `temporal_hypergraph()` plus snapshots** — long or
  wide relational memberships, growing or interval evolution, active,
  cumulative and aggregate views, with optional multi-edge collapse. This is
  native hypernets code and introduces no Dynet dependency. *Oracle*: the paper's
  temporal GFCC/ICSID construction.
- [x] **DONE 2026-09-02:** worst-community weighted conductance added to
  `hg_community_quality()`. Larger bundled corpora remain carried forward.
  **Resolved 2026-09-02:** `hg_embed()` and `hypergraph_embed()` now expose
  the same documented wrapper over the shipped spectral/NMF coordinates.
- [x] **DONE 2026-09-02 paper-oracle validation** — reconstructed the authors'
  ICSID aggregate (441 nodes, 742 cases) from Zenodo 8081507 and reproduced
  the published observed motif counts exactly: Y = 478, T = 7, O = 0.
  Formula fixtures also cover HypergraphX 1.5 s-centrality, subhypergraph
  centrality and configuration-MCMC conventions.
  Not a vignette: legal citation networks sit outside the ownership contract.
  It is the equivalence oracle for the entire projection tier, which has none
  today — HyperNetX/XGI/HyperG cover construction and spectral only.


## v1.0 — CRAN + the paper

- CRAN submission (the first-in-R claims re-verified against CRAN at
  submission).
- Flagship paper + vignette: **first fully-native R hypergraph text
  analysis** — sbert embeddings → kNN/windowed hyperedges → EDVW spectral
  clustering + transductive classification → tensor centralities; zero
  Python at runtime, oracle-verified against HyperNetX/XGI.
- [x] DONE 2026-09-05: BERTopic as the reported baseline (actual package
  0.17.4, ten seeds, paired effects + bootstrap CIs; `benchmarks/RESULTS.md`).

## Scope boundary

Native hypergraph neural models are in scope and already include HGNN and
HyperGAT. PLM training/wrappers, hMETIS/KaHyPar partitioning and retrieval
stacks remain outside hypernets; pretrained embeddings enter through sbert or a
caller-supplied matrix. Nestimate is retained only as a historical oracle,
never as a runtime dependency.
