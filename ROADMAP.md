# hypernets roadmap and timeline

Single source of truth, written 2026-09-29 from three audits (papers,
reference repositories, every planning file against NAMESPACE and git).
It supersedes `EXPANSION-PLAN.md`, `ROADMAP-analytics.md`,
`ROADMAP-text.md`, `TODO-text.md`, `MIGRATION-PLAN.md` and
`AGENT-NOTE-STBM-HYPERTEXT-TOPIC-MODELS.md`, which stay as history.
Sidelined and rejected work is registered in `workinprogress/README.md`.

## Rules every item obeys

1. **No method ships without a published reference.** A construction the
   literature does not define goes to `workinprogress/`, not the package.
2. **Every numerical verb gets an independent oracle** (another package,
   the authors' code, or a hand-computed published example). Identity with
   Nestimate is provenance, not an oracle.
3. Tidy verbs, named arguments, `net_*` + four methods (CLAUDE.md taxonomy).
   Classes of verbs grounded in Nestimate are Nestimate's classes.
4. Plots delegate to cograph or the package's own hypergraph plot; no
   hand-built layers.
5. Dates below are targets proposed on 2026-09-29, not commitments.

## Where we are (2026-09-29)

hypernets 0.5.1, committed as 3f1cca2 (not pushed). 115 exports.
NOT_CRAN suite 4054 pass / 0 fail; `R CMD check --as-cran` 0/0/0;
equivalence suite green; Eventdata26 `plot_blobs()` parity 371/0.

| Version | Date | Added |
|---|---|---|
| 0.1.0 | 2026-08-24 | Memory family moved verbatim from Nestimate 0.9.0 |
| 0.1.4–0.1.5 | 2026-08-25 | `bootstrap_hon`, `compare_hon`, `hon_centrality`, long-format input |
| 0.2.0–0.2.1 | 2026-08-26 | Simplicial and hypergraph families absorbed; taxonomy; `top =` |
| 0.3.0–0.3.12 | 2026-09-02 → 09-07 | Text family folded in; Legal Hypergraphs reproduction; NMF clustering; neural set (HGNN, HyperGAT, HyperGCN, HNHN, AllSet, HGAT); Wasserstein; random generators; topic verbs |
| 0.4.0 | 2026-09-07 | Renamed honets → hypernets |
| 0.4.1–0.4.8 | 2026-09-07 → 09-20 | `hg_hypa`, pebble plots, sparse overflow fix, component packing |
| 0.5.0 | 2026-09-20 | `markov_stability`, `hg_sequences`, long-format guard (`hon_outcome`, later withdrawn) |
| 0.5.1 | 2026-09-29 | `hon_communities` (map equation, Infomap oracle); blob/overlay plots with Eventdata26 parity; `ring_sequences`; `hon_outcome` withdrawn; Tier C rejected |

## Phase 0 — close 0.5.1 honestly (by 2026-10-03)

| Item | Why | Evidence / reference |
|---|---|---|
| Add a primary reference to `markov_order_test()` | It has none; breaks rule 1 | `R/memory_markov_order.R` has 0 `@references`; Anderson & Goodman (1957), Scholtes (2017) |
| Document the field of `betti_numbers()` (rationals, `qr()` rank) vs `persistent_homology()` (Z/2) | They disagree only with torsion; undocumented | `R/simplicial.R:435`; `R/simplicial_filtration.R:178` |
| HyperGAT benchmark: compare like with like | R8/R52 run with `semantic = "none"` (the paper's "w/o semantic" ablation) but compared against full HyperGAT | `benchmarks/results/hypergat.csv`; Ding et al. (2020) Table 4 |
| Put rule 1 in CLAUDE.md | Recorded only in LEARNINGS/NEWS | audit contradiction 18 |
| Refresh stale notes: `cran-comments.md` counts, HANDOFF, NEWS bullet saying `hon_communities` plot is cograph-based, NEWS missing 0.4.6 | Files contradict code | audit contradictions 8–12 |
| Sync `repos/` notes | **Done 2026-09-29**: all oracle claims re-verified by grep; only Dynet was stale (source now `../Dyna` 0.5.1; `../temporal` no longer exists) | repos audit §5 |
| Complete `papers/README.md` bibliography | **Done 2026-09-29**: every roxygen `@references` work has a Crossref-verified entry; LDA, Deep Sets, Set Transformer, GCN, TextGCN, SBERT marked used-not-cited; XLNet unused. Still unverified: Page 1999 report no., Lee 2025 article no., some proceedings pages | papers audit §4 |
| Decide version (0.6.0 recommended) and commit | New export, withdrawn export, plot API | — |
| Pipeline: render `with_blob_plot.Rmd`; rebuild the 3 `mmm_steps` figures | Parity covered 19 of 22 calls | Step-Groups-100k l.275–277 |

## Phase 0b — hypernets imports Nestimate; hypergraphs are hypernets' own (before CRAN)

Decided 2026-09-29 (author). hypernets cannot be independent of the code it
was moved from, so it builds on Nestimate instead of carrying copies.

- **Direction:** hypernets `Imports: Nestimate (>= 0.8.5)`. Nestimate never
  imports or suggests hypernets, and **Nestimate does not change**.
- **Grounded in Nestimate** (hypernets deletes its copies and calls or
  re-exports Nestimate's exported verbs): memory builders (`build_hon`,
  `build_honem`, `build_hypa`, `build_mogen`, `markov_order_test`,
  `path_dependence`, `mogen_transitions`, `path_counts`, `pathways`), the
  simplicial family (`build_simplicial`, `persistent_homology`, `q_analysis`,
  `betti_numbers`, `euler_characteristic`, `simplicial_degree`,
  `verify_simplicial`, `bottleneck_distance`, `persistence_landscape`),
  `markov_stability`, `clique_expansion`. All are exported by CRAN 0.8.5
  (checked from the tarball).
- **hypernets' speciality:** the hypergraph family (every hypergraph
  constructor, the Laplacian / clustering / transduction / centrality /
  measures engines, `hg_*`, neural set, plots) and the text family; plus
  the verbs it adds on top of Nestimate's objects (`bootstrap_hon`,
  `compare_hon`, `hon_centrality`, `hon_communities`, `wasserstein_distance`).

Measured 2026-09-29 with both packages loaded: **26 functions masked and 32
S3 methods overwritten**; a Nestimate `bipartite_groups()` hypergraph printed
through hypernets' method loses its "Source" line. With the import, both load
every time, so the collision must be removed, not tolerated.

| # | Item | Why |
|---|---|---|
| G1 | **Done 2026-09-29.** `hg_*` is the only public hypergraph naming: 27 `hypergraph_*` aliases dropped; `hg_allset` / `hg_hypergcn` / `hg_hnhn` / `hg_snapshot(s)` / `hg_laplacian` / `hg_joint_cluster` primary; engines `hypergraph_centrality` / `_cluster` / `_measures` / `_transduction` internal (`.hg_*_fit`) behind `hg_centrality` / `hg_cluster` / `hg_measures` / `hg_classify`, which gained the engines' arguments and docs; `build_hypergraph()` -> `network_hypergraph()`; `clique_expansion()` -> `hg_clique_expansion()`; class `net_hypergraph*` -> `net_hg*`. Exports 115 -> 77. Guard: `tests/testthat/test-api-names.R`. Tests 4064/0; `--as-cran` 0/0/1 (network-time note only); R8 guard 0.8451; equivalence suite green after two pre-existing test-harness fixes (source-Nestimate `load_all()` left attached; blob reference pinned to pipeline commit cc58acf) | The inherited `hypergraph_*` names are Nestimate's; the same names overwrite each other's S3 methods and mask verbs in every session. hypernets is unreleased, so the rename breaks no one |
| G2 | **Done 2026-09-29.** Delete hypernets' copies of the grounded verbs; `importFrom` + re-export them; results and classes become Nestimate's (`simplicial_complex`, `persistent_homology`, ...) | One implementation of each |
| G3 | **Done 2026-09-29.** Keep hypernets' additions to grounded objects only where Nestimate has no method: `as.data.frame()` for `net_hon`, `net_honem`, `net_hypa`, `net_mogen`, `net_markov_order`, `net_path_dependence`, `simplicial_complex`, `persistent_homology`, `q_analysis`, `persistence_landscape`, `net_markov_stability` | Registering a method Nestimate lacks overwrites nothing; Rule 0 needs the accessors |
| G4 | **Done 2026-09-29.** hypernets-only verbs that used Nestimate *internals* (`bootstrap_hon`, `compare_hon`: `.hon_build_*`, `.hon_extract_rules_count`; `temporal_hypergraph`, `window_hypergraph`: `.coerce_sequence_input`; `wasserstein_distance`: `.ph_as_diagram`; `hypergraph_joint_cluster`) keep private copies of those internals, with an identity test against Nestimate's | `:::` is not allowed on CRAN and Nestimate does not change |
| G5 | **Done 2026-09-29.** Losses to record in NEWS (Nestimate lacks them): `build_hon(action/actor/time)` long input, the long-format guard, `top =` on `mogen_transitions` / `path_counts` / `simplicial_degree`, hypernets' `markov_stability` rewrite (reducibility check, landscape plot), `net_*` class names. Upstreaming them to Nestimate later is optional and additive | Nestimate keeps what it has |
| G6 | **Done 2026-09-29.** Against CRAN 0.8.5 (tolerance 0, 23 calls over memory, simplicial, `markov_stability`, `clique_expansion`): 22 identical; `markov_stability` differs only because 0.8.5 rounds `$stability` (4 / 2 decimals; passage times identical). Floor: `Nestimate (>= 0.8.5)`. Script: `local_testing_and_equivalence/g6-nestimate-085.R` | hypernets was proven identical to 0.9.0, not to the CRAN version it would import |
| G7 | **Done 2026-09-29.** Rewire hypernets internals that assumed its own classes: `build_hypergraph` calls `build_simplicial(type = "clique")` (now returns `simplicial_complex`), `.hg_extract_adj`, `build_simplicial(type = "pathway")` readers | Class names change back to Nestimate's |
| G8 | **Done 2026-09-29.** CLAUDE.md: replace "Nestimate is not a dependency" and the provenance section with this ownership | Contradicts the decision |

Superseded the same day: "Nestimate imports hypernets" (would have made
Nestimate change and depend on hypernets).

## Phase 1 — oracles for what already ships, then CRAN (October 2026)

Rule 2 applied to the shipped surface. No new methods in this phase.

| Verb(s) | Oracle to add |
|---|---|
| `betti_numbers`, `euler_characteristic`, `bottleneck_distance`, `persistence_landscape` | GUDHI (SimplexTree Betti, `bottleneck_distance`, `representations.Landscape`) |
| `hg_hypergcn`, `hg_hnhn` | **Done 2026-09-29**: DHG 0.9.7 + official HyperGCN code (`de04938`), 6/6 after fixing the symmetric-matrix bug it found; HNHN covered for alpha = beta = 0 only (official code unusable as an oracle beyond that) |
| `hg_allset`, `heterogeneous_hgat` | open: official AllSet needs torch_geometric/torch_scatter and differs architecturally (unscaled leaky-ReLU attention, biased K/V, ReLU in the block); official HGAT code departs from the paper equations -- **decision needed**: follow the code or the paper |
| `markov_stability` | `markovchain::steadyStates`, `meanFirstPassageTime` |
| `build_hypa` | the HYPA reference code (not only BiasedUrn / Monte Carlo) |
| `hg_topic_quality` | **Done 2026-09-29**: stm 1.3.7 UMass/FREX (max diff 0), text2vec 0.6.6 NPMI (1e-8); `hg_stability` vs fpc `clusterboot` (1.1e-16); `hg_cocluster` vs sklearn (1e-8) |
| `text_hypergraph` incidence | **Done 2026-09-29**: official `Data.get_slice()`, exact on 40 cells |
| Zhou Laplacian | **Done 2026-09-29**: XGI 0.10.2 and HNX 2.4.3 to 2.2e-16 (XGI's weighted variant uses unweighted degree; exact mapping tested) |
| `hg_motifs` Y/T/O convention | **Done 2026-09-29**: hypergraphx 1.8.0, 72/72 incl. the configuration null replayed draw for draw |

Then: win-builder (devel + release), GitHub Actions on Linux/macOS/Windows,
fill `cran-comments.md`, submit. (Needs the author's go: outward-facing.)

## Phase 2 — referenced structural gaps (November 2026)

Ranked by value × cost; each has its reference and oracle before any code.

| # | Item | Reference | Oracle |
|---|---|---|---|
| 0 | **Done 2026-09-29** (`hg_mmsbm()`; authors' code to 2.1e-12). **Probabilistic (mixed) membership: Hy-MMSBM** -- moved forward 2026-09-29 at the author's request. Each node gets a membership vector over K communities from a likelihood model of the hypergraph (EM with random starts); replaces the fuzzy-c-means reading in `hg_membership()` as the principled soft membership, and gives topic membership for documents | Ruggeri, Contisciani, Battiston & De Bacco (2023), *Sci. Adv.* 9(28):eadg9159 | authors' code github.com/nickruggeri/Hy-MMSBM: one EM step from fixed u, w compared exactly; log-likelihood at convergence |
| 1 | Hodge / combinatorial Laplacians; export boundary matrices (already computed internally) | Horak & Jost (2013); Lim (2020); Schaub et al. (2020) | XGI `hodge_laplacian`, TopoNetX |
| 2 | Hypergraph → simplicial closure (`build_simplicial(type = "hypergraph")`), simpliciality | Battiston et al. (2020) §II.2; Landry et al. (2023) | XGI simpliciality, HNX `homology_mod2`, GUDHI |
| 3 | **Done 2026-09-29** except E&RV global C2(H) (no oracle; `workinprogress/`): `hg_transitivity()`, `hg_assortativity()`, `hg_degree_correlation()`. Clustering coefficients and degree assortativity | Estrada & Rodríguez-Velázquez (2006); Battiston §III.C | XGI clustering, HGX `degree_correlation` |
| 4 | Node-level s-centralities (betweenness, closeness, harmonic, eccentricity) | Aksoy et al. (2020) | HyperNetX `s_centrality_measures` |
| 5 | **Done 2026-09-29** (`hg_centrality(type = "katz")`). Katz / alpha centrality | Battiston §III.B.3 | XGI, igraph `alpha_centrality` |
| 6 | General order-3/4 and directed hypergraph motifs | Lotito et al. (2022) | hypergraphx `compute_motifs` |
| 7 | **Done 2026-09-29** (`hg_read_hif()` / `hg_write_hif()`, schema-valid, XGI + HNX round-trips). HIF read/write (makes every later oracle test cheaper) | HIF standard (Coll et al.) | round-trip via XGI / HNX |
| 8 | Export the configuration-model sampler already inside `hg_null_test`; statistically validated hyperedges | Chodrow (2020); Musciotto et al. (2021) | HGX `configuration_model`, `get_svh` |
| 8b | **Done 2026-09-29** (`hg_modularity()`, `hg_communities(type = "irmm")`; HNX 2.4.3 to 5.6e-16). Topic follow-ups: Kumar IRMM + Kaminski hypergraph modularity score for topic communities (proposal rank 3; oracle HyperNetX 2.4.3 `kumar`, `modularity`); topic shares by period via `hg_snapshots` (rank 5) | Kumar et al. (2020); Kaminski et al. (2019) | HyperNetX 2.4.3 |
| 9 | Evaluation verbs: accuracy, macro/micro-F1, confusion table; Jaccard/F1/ANC/ACo for clusterings | Li et al. (2022); Hayashi et al. (2020) §6 | scikit-learn |

## Phase 3 — memory family (December 2026)

| Item | Status | Reference | Oracle |
|---|---|---|---|
| Tier B: `compare_hon` for k groups, covariates, permutation within clusters | open | Anderson & ter Braak (2003); Winkler et al. (2015) | frozen two-group fixture; KS uniformity under the null |
| Group model: one HON per group (condition / cohort / human vs AI) | idea; **reference needed first** | — | — |
| A4: `simulate()` / `predict()` for `net_mogen` | open | Gote & Scholtes (2023); Scholtes (2017) | pathpy `MultiOrderModel` |
| A1b: MOGen inference with order re-selection per replicate | open | Efron & Tibshirani (1993) | coverage simulation |
| Order selection vs pathpy `estimate_order` | open | Scholtes (2017) | pathpy 2.2.0 |
| HON eigenvector / visitation centralities, algebraic connectivity | open | Scholtes et al. (2014) | pathpy 2.2.0 |
| E2: R↔JS equivalence of `markov_order_test` vs `tnaj` | open | — | `validation/` runner |
| `hon_communities` coarse-tuning (Alg. 6) | partial | Edler et al. (2017) | infomap |

## Phase 4 — decisions for 2027 Q1

| Item | Decision needed |
|---|---|
| Benchmark reproductions (Hayashi 20NG G1–G4; HyperGAT 20NG/MR/Ohsumed incl. `semantic = "lda"`; neural set on Cora/Citeseer/Pubmed) | worth the compute? |
| Hypergraph-native communities (Hypergraph-MT (Kumar modularity is Phase 2 8b; Hy-MMSBM moved to Phase 2 item 0)) | Kumar et al. (2020); Ruggeri et al. (2023, Hy-MMSBM, doi:10.1126/sciadv.adg9159); Contisciani et al. (2022, Hypergraph-MT); oracles HNX 2.4.3 (`kumar`, `modularity`), authors' code |
| Random and growth models for simplicial complexes and hypergraphs | Iacopini et al. (2019); Courtney & Bianconi (2016); oracle XGI |
| **Dynamics family** (contagion, Kuramoto, opinion, games; Battiston §V–VIII) | a new family: in scope or not? |
| C1: metric Vietoris–Rips for `build_hypergraph(type = "vr")` | Zomorodian (2010); GUDHI |
| Flagship paper | after CRAN |

## Register: sidelined, rejected, deferred

| Item | Status | Where / why |
|---|---|---|
| `hon_outcome()` | sidelined 2026-09-29 | `workinprogress/hon_outcome/`; "exposure" feature has no reference |
| Per-actor / per-session HON (`build_hon(by =)`) | rejected 2026-09-29 | `workinprogress/hon_per_actor/`; overfits (median session 17 transitions) |
| STBM / hypertext topic models | deferred 2026-08-29 | data-model question unresolved |
| Nestimate long-format guard port | optional | Nestimate keeps what it has (Phase 0b, G5) |
