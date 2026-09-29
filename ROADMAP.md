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
4. Plots delegate to cograph or the package's own hypergraph plot; no
   hand-built layers.
5. Dates below are targets proposed on 2026-09-29, not commitments.

## Where we are (2026-09-29)

hypernets 0.5.1, uncommitted on top of 664e379. 115 exports.
NOT_CRAN suite 3984 pass / 0 fail; `R CMD check --as-cran` 0/0/0;
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
| Sync `repos/` notes | 7 stale: Wasserstein "open", BERTopic "pending", HGX 1.5→1.8, Dynet 0.3.53→0.4.10, pathpyG unrecorded, GUDHI/SimplicialComplex claimed as oracles but never used | repos audit §5 |
| Complete `papers/README.md` bibliography | ~40 cited works missing; LDA, Deep Sets, GCN, TextGCN, SBERT used but uncited; XLNet unused | papers audit §4 |
| Decide version (0.6.0 recommended) and commit | New export, withdrawn export, plot API | — |
| Pipeline: render `with_blob_plot.Rmd`; rebuild the 3 `mmm_steps` figures | Parity covered 19 of 22 calls | Step-Groups-100k l.275–277 |

## Phase 1 — oracles for what already ships, then CRAN (October 2026)

Rule 2 applied to the shipped surface. No new methods in this phase.

| Verb(s) | Oracle to add |
|---|---|
| `betti_numbers`, `euler_characteristic`, `bottleneck_distance`, `persistence_landscape` | GUDHI (SimplexTree Betti, `bottleneck_distance`, `representations.Landscape`) |
| `hypergraph_hypergcn`, `hypergraph_hnhn` | DHG 0.9.7 `HyperGCNConv`, `HNHNConv` (same pattern as the HGNN test) |
| `hypergraph_allset`, `heterogeneous_hgat` | official AllSet and HGAT code |
| `markov_stability` | `markovchain::steadyStates`, `meanFirstPassageTime` |
| `build_hypa` | the HYPA reference code (not only BiasedUrn / Monte Carlo) |
| `hg_topic_quality` | `stm::semanticCoherence`, `exclusivity` |
| `text_hypergraph` incidence | official HyperGAT `utils.py::get_slice()` |
| Zhou Laplacian | HyperNetX / XGI normalized Laplacian |
| `hg_motifs` Y/T/O convention | hypergraphx 1.8 `compute_motifs` |

Then: win-builder (devel + release), GitHub Actions on Linux/macOS/Windows,
fill `cran-comments.md`, submit. (Needs the author's go: outward-facing.)

## Phase 2 — referenced structural gaps (November 2026)

Ranked by value × cost; each has its reference and oracle before any code.

| # | Item | Reference | Oracle |
|---|---|---|---|
| 1 | Hodge / combinatorial Laplacians; export boundary matrices (already computed internally) | Horak & Jost (2013); Lim (2020); Schaub et al. (2020) | XGI `hodge_laplacian`, TopoNetX |
| 2 | Hypergraph → simplicial closure (`build_simplicial(type = "hypergraph")`), simpliciality | Battiston et al. (2020) §II.2; Landry et al. (2023) | XGI simpliciality, HNX `homology_mod2`, GUDHI |
| 3 | Clustering coefficients and degree assortativity | Estrada & Rodríguez-Velázquez (2006); Battiston §III.3 | XGI clustering, HGX `degree_correlation` |
| 4 | Node-level s-centralities (betweenness, closeness, harmonic, eccentricity) | Aksoy et al. (2020) | HyperNetX `s_centrality_measures` |
| 5 | Katz / alpha centrality | Battiston §III.2.3 | XGI, igraph `alpha_centrality` |
| 6 | General order-3/4 and directed hypergraph motifs | Lotito et al. (2022) | hypergraphx `compute_motifs` |
| 7 | HIF read/write (makes every later oracle test cheaper) | HIF standard (Coll et al.) | round-trip via XGI / HNX |
| 8 | Export the configuration-model sampler already inside `hg_null_test`; statistically validated hyperedges | Chodrow (2020); Musciotto et al. (2021) | HGX `configuration_model`, `get_svh` |
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
| Hypergraph-native communities (Kumar modularity, Hy-MMSBM) | Kumar et al. (2020); Contisciani et al. (2022); oracles HNX, HGX |
| Random and growth models for simplicial complexes and hypergraphs | Iacopini et al. (2019); Courtney & Bianconi (2016); oracle XGI |
| **Dynamics family** (contagion, Kuramoto, opinion, games; Battiston §V–VIII) | a new family: in scope or not? |
| C1: metric Vietoris–Rips for `build_hypergraph(type = "vr")` | Zomorodian (2010); GUDHI |
| Flagship paper | after CRAN |
| Nestimate delegation (T2) | conflicts with "hypernets is independent" in CLAUDE.md |

## Register: sidelined, rejected, deferred

| Item | Status | Where / why |
|---|---|---|
| `hon_outcome()` | sidelined 2026-09-29 | `workinprogress/hon_outcome/`; "exposure" feature has no reference |
| Per-actor / per-session HON (`build_hon(by =)`) | rejected 2026-09-29 | `workinprogress/hon_per_actor/`; overfits (median session 17 transitions) |
| STBM / hypertext topic models | deferred 2026-08-29 | data-model question unresolved |
| Nestimate long-format guard port | blocked | must be done from a Nestimate session |
