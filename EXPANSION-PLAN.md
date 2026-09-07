# Expansion Plan — higher-order networks

Owner decisions taken 2026-08-25, **substantially revised 2026-08-26**.
This is the consolidated roadmap reconciling
`Nestimate/todo/COVERAGE-CATCHUP.md`, `../texthypergraph/TODO.md`, and
the HON-1..15 / EG-1..13 sections of `Nestimate/TODO.md`.

## Decision reversed 2026-08-26: one package, not two

Decisions 1 and 2 below were made on 2026-08-25 and **overturned the
next day**. They are kept here with the evidence that overturned them,
because the reasoning is the useful part.

> ~~1. **hypernets is the fourth delegation sibling**, scaffolded now.~~
> ~~2. **hypernets stays sequence/path only.** No hypergraph code here —
> the two paradigms share no infrastructure (ordered k-grams vs
> unordered multi-way co-membership).~~

**What disproved decision 2, within hours of taking it:**

- [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
  needs clique enumeration, which lives in Nestimate’s `simplicial.R`.
  Scaffolding hypernets therefore required **copying a 171-line
  closure** into `hypernets/R/utils.R` behind an internal
  [`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md)
  shim. “Zero shared infrastructure” was false on day one.
- B1
  ([`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)),
  the first feature built after the split, takes **the same
  `action`/`actor`/`time` sequence input as
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)**.
  It is a sequence verb that emits hyperedges, and it had no natural
  home: it went into hypernets, but hypernets had an equal claim.
- Nestimate has **zero internal consumers** of its own 1,956-line
  simplicial /TDA family once `hypergraph.R` departs — while
  `build_simplicial(type = "pathway")` dispatches on hypernets’
  `net_hon`, `net_hypa` and `net_mogen`. The TDA layer already read
  hypernets classes.
- The literature the plan itself cites (Battiston et al. 2020; Bianconi
  2021; Tian & Zafarani 2024) treats memory networks, simplicial
  complexes and hypergraphs as **one field**, called higher-order
  networks.
- `hypernets` and `hypernets` are three characters apart, both ending in
  `nets`, both higher-order packages by the same author. That is a
  usability defect independent of the architecture.

The split had put a boundary through the middle of a field rather than
between fields.

**Decision now in force:** hypernets is *the* higher-order networks
package, covering all three structure families under one taxonomy (see
[`?hypernets`](https://mohsaqr.github.io/hypernets/reference/hypernets-package.md)
and `CLAUDE.md`). hypernets 0.1.2 is folded in and retired — nothing was
thrown away: its ten R files, tests, equivalence suite and both
tutorials moved across, and 223 lines of duplication were **deleted**
rather than maintained. Nestimate imports one sibling instead of two.

## Family map after the consolidation

| Package | Scope | Status |
|----|----|----|
| Nestimate | estimation hub; re-exports delegated verbs | 0.9.0, untouched |
| psychnets | cross-sectional psychometric networks | delegated |
| idiographic | person-specific temporal models | delegated |
| **hypernets** | **higher-order networks: memory + simplicial + hypergraph** | **0.2.0 (T0 done for all three families)** |
| ~~hypernets~~ | ~~hypergraphs~~ | **retired 2026-08-26, folded into hypernets** |

## Track A — hypernets (higher-order / sequence paradigm)

New verbs stay in the existing one-file-per-method layout, each with
shipped tests + an entry in `local_testing_and_equivalence/` where an
oracle exists.

| \# | Feature | Verb sketch | Oracle / validation | Reference |
|----|----|----|----|----|
| A1 | **DONE 2026-08-25 (HON part).** Bootstrap CIs + rule support (`bootstrap_hon(data, n_boot, level, ...)`) and two-sample permutation comparison (`compare_hon(x, y, n_perm, ...)` — two-cohort interface instead of the sketched group vector; per-edge BH + global test). Per-sequence counts precomputed once, replicates reweight counts (never re-scan); RNG drawn serially so parallel == serial under a seed; long-format input (`action`/`actor`/`time`); `human_long`/`ai_long` bundled; tutorial `hon_inference.html`. **MOGen-transition inference deferred to a new item A1b** (needs per-replicate order re-selection — different machinery). | `bootstrap_hon(data, n_boot)`, `compare_hon(x, y)` | exact multiset-recount equivalence of the weighted aggregation; coverage of a known conditional probability (8+/10 seeded runs); null calibration of the global permutation test; planted-difference power — all green | Efron & Tibshirani (1993); Good (2005); Benjamini & Hochberg (1995) |
| A2 | **DONE 2026-08-25.** Higher-order centralities projected to first-order states, plus `project = FALSE` for context-level values, `projection =` (scaled/last/first/all) and `sort_by =`. [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md) also gained the long-format `action`/`actor`/`time` interface. | `hon_centrality(hon, type = c("pagerank", "betweenness", "closeness"))` | pathpy 2.2.0 on identical topologies (k = 1,2,3 + 8 random corpora): betweenness/closeness exact \< 1e-10, PageRank to pathpy’s own tol = 1e-6; igraph for the unprojected kernels. Found and fixed a duplicate-cell bug in the distance projection. | Scholtes, Wider & Garas (2016) EPJ B 89:61 |
| A3 | Memory-network community detection (map equation on the HON state graph) | `hon_communities(hon, method = "map_equation")` | pathpy/Infomap on second-order state networks | Rosvall et al. (2014) Nat. Commun. 5:4630 |
| A4 | Simulation + prediction from fitted models | `simulate(net_mogen, n)`, `predict(net_mogen, newdata)` S3 methods | pathpy `MultiOrderModel` likelihood/prediction; RNG-seeded reproducibility tests | Gote & Scholtes (2023), Scholtes (2017) |

Suggested order: A1 first (pure house paradigm, no new theory, highest
research value), then A2 (small, oracle exists), A3, A4.

**CRAN sequencing:** hold the T1 CRAN submission until A1–A2 land, then
submit as hypernets 0.2.0 — one submission instead of two.

## Track B — hypergraph family (was: hypernets)

All Track B work now lands in hypernets under the `hypergraph_*.R`
files. The “hypernets” labels below are historical.

| \# | Feature | Verb sketch | Oracle / validation | Reference |
|----|----|----|----|----|
| B0 | **T0 scaffold — DONE this session:** moved `hypergraph.R`, `hypergraph_centrality.R`, `hypergraph_laplacian.R`, `hypergraph_measures.R`, `bipartite_groups.R`, `clique_expansion.R` (+ helpers `.wrap_netobject`, `.validate_mcml_matrix`, `.extract_edges_from_matrix`, and the simplicial-clique closure behind an internal [`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md) shim) verbatim from Nestimate 0.9.0 | 8 exports, 11 S3 methods | exact [`identical()`](https://rdrr.io/r/base/identical.html) identity vs Nestimate proven (61 assertions, incl. RNG paths); 281 shipped tests + 682 equiv assertions green; quick check 0 errors / 0 warnings | — |
| B1 | **DONE 2026-08-25.** Windowed sequence hyperedges (tumbling/sliding window w, hyperedge = distinct states in window, weight = window count; incidence cells = occurrence totals, i.e. EDVW). Shipped as a dedicated constructor `window_hypergraph(data, window, step, action, actor, time, min_size)` rather than an input mode of [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md) (whose adjacency-typed first argument and clique args don’t transfer); Laplacian family defaults `edge_weights` to the window counts; new tidy [`as.data.frame.net_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_hypergraph.md) accessor. | `window_hypergraph(data, window = w)` | wtna co-occurrence oracle (w = 2 tumbling, off-diagonals), bipartite_groups whole-sequence reduction, conservation invariants — all green; identity vs Nestimate re-verified after the edge-weight default change | Ding et al. (2020) EMNLP (HyperGAT construction) |
| B2 | **DONE 2026-08-25.** EDVW random-walk centrality / hypergraph PageRank. Shipped as `hypergraph_centrality(hg, type = "pagerank", damping =, edge_weights =)` (opt-in; vertex weights come from the incidence itself, so no `vertex_weights` arg needed); EDVW transition matrix shared with the random-walk Laplacian; scalar `edge_weights` recycling; tutorial `hypergraph_pagerank.html`. | `hypergraph_centrality(hg, type = "pagerank")` | igraph collapse-theorem oracle (25 configs), dense fixed-point solve, damping→1 vs Laplacian stationary pi — all green; identity vs Nestimate re-verified after the shared-helper refactor | Chitra & Raphael (2019) ICML |
| B3 | Bayesian hypergraph reconstruction suite (Nestimate TODO HON-11..15): `reconstruct_hypergraph`, `hyperedge_significance`, `hypergraph_order_select`, `compare_hypergraphs`, soft clustering | as named in `Nestimate/TODO.md` | per-item oracles listed there | Young, Petri & Peixoto (2021) lineage per TODO |
| B4 | kNN embedding hypergraph + NLP bridge vignette (quanteda/tidytext → [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md) → Laplacian machinery) | `knn_hypergraph(embeddings, k, weight = "cosine")` | `HyperG::knn_hypergraph` for unweighted construction; vignette-first before freezing the API | texthypergraph/TODO.md |

Small carried-over items:
[`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md)
(hypernets, when a use-case appears), XGI second oracle for shipped
centralities (hypernets equiv suite),
[`wasserstein_distance()`](https://mohsaqr.github.io/hypernets/reference/wasserstein_distance.md)
(**done 2026-09-02** in hypernets’ simplicial family), random hypergraph
samplers (**done 2026-09-02** with HyperG parity; hypernets owns
simulation as part of its hypergraph family), per the
simulation/computation split).

Suggested order: B0 → B1 → B2 (all done) → B3 → B4.

## Track C — simplicial / TDA family (absorbed 2026-08-26)

Moved verbatim from Nestimate 0.9.0 in the consolidation. Nine exports,
identity-tested. Carried-over items from `Nestimate/TODO.md`:

| \# | Feature | Notes |
|----|----|----|
| C1 | A genuine metric Vietoris-Rips construction for `build_hypergraph(type = "vr")` | currently raises rather than silently aliasing `"clique"` |
| C2 | **DONE 2026-09-02.** [`wasserstein_distance()`](https://mohsaqr.github.io/hypernets/reference/wasserstein_distance.md) for persistence diagrams | finite Wasserstein order and configurable ground metric; native exact Hungarian assignment |
| C3 | Convert the inherited bare-name `aes()` plot code to the `.data` pronoun | the family is the only remaining user of [`globalVariables()`](https://rdrr.io/r/utils/globalVariables.html) here |

## Delegation tiers

1.  **T0 — DONE.** Memory family 2026-08-24 (hypernets 0.1.0);
    simplicial and hypergraph families 2026-08-26 (hypernets 0.2.0, the
    latter via hypernets 0.1.2). All three identity-tested against
    Nestimate 0.9.0. Nestimate untouched.
2.  **T1 — CRAN submission of hypernets 0.2.0.** `--as-cran` is already
    clean (0 errors, 0 warnings). Remaining: `cran-comments.md`, a
    decision on whether the pre-existing `test-hypa.R` warnings are
    worth clearing first.
3.  **T2 (Nestimate delegation release):** delete the memory, simplicial
    and hypergraph R files + moved helpers where unused elsewhere;
    `Imports: hypernets`; thin forwarders (gimme pattern), mapping the
    renamed surface — `bipartite_groups()` →
    [`hypernets::group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
    `build_hypergraph(method =)` → `type =`, and the `net_*` result
    classes; move equivalence suites; htna gate + `--as-cran`.

## Bookkeeping done / to do

`Nestimate/HYPERNETS-DELEGATION-PLAN.md` §6b — annotate: revisited
2026-08-25 (hypernets scaffolded), **reversed 2026-08-26** (folded into
hypernets); pointer here.

`Nestimate/todo/COVERAGE-CATCHUP.md` + `../texthypergraph/TODO.md` —
retarget the hypergraph and TDA items’ repo labels Nestimate → hypernets
at T2 (not before; until T2 the shipping copies are Nestimate’s).

Add hypernets to the `saqr_AR` skill family list
(`../Writing/saqr_Coding conventions.md`); remove hypernets.

Archive the `hypernets` GitHub repo with a README pointing at hypernets.
