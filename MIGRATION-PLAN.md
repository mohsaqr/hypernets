# Migration plan — texthypergraph folds into hypernets

**Status (2026-09-01): EXECUTED in the working tree, uncommitted.** See
`HANDOFF.md` / `CHANGES.md` for what was done and the measured results.
Deviations from the plan below, found while executing it:

- `Nestimate::bipartite_groups()` does not exist in hypernets under that name;
  it is `group_hypergraph(member =)` — a pure rename with `identical()`
  output (verified on binary and weighted builds). Four call sites and all
  test fixtures were renamed.
- The spectral S3 methods: **hypernets' were kept, not ours** — they are the
  superset (plot methods, `top =`, scalar `edge_weights`, window-count
  defaults). Only the `normalization` argument and its helper
  (`.hl_score_predictions()`) were ported.
- `hg_pagerank()` stays as the tidy verb (personalization, sparse path,
  convergence warning) beside `hypergraph_centrality(type = "pagerank")`;
  parity to `1e-10` is tested rather than merging the two engines.
- Brought in by file copy in the working tree, not `git subtree add` —
  both trees were uncommitted and the user asked for a clean merge; the
  texthypergraph history stays in its own (retired) repo.
- CI moved *to* hypernets (it had none): `R-CMD-check.yaml`, `pkgdown.yaml`,
  `_pkgdown.yml` (destination `pkgdown/`, since `docs/` holds the hand-knit
  family documents).
- `Tutorial_docs/hon_inference.{Rmd,html}` (deleted on disk, tracked at
  HEAD) was left for the commit step.

---


**Decision (2026-09-01):** one package. `hypernets` is the home and keeps its
name; `texthypergraph` becomes its text family. Nothing moves until this plan
is approved.

## Why one package

- `texthypergraph` is 52% generic hypergraph math, and ~1,100 of those lines
  duplicate hypernets exactly (verified numerically identical this session, not
  assumed).
- **The merge drops the `Nestimate` dependency.** texthypergraph imports three
  Nestimate functions in code — `bipartite_groups`, `hypergraph_measures`,
  `hypergraph_centrality` (`wtna` appears only in comments). hypernets already
  has all three. CRAN's Nestimate is 0.8.5 while texthypergraph requires
  >= 0.9.0, so today texthypergraph **cannot** be submitted; merged, the
  dependency disappears and the blocker with it.
- The families are loosely coupled: memory and simplicial call nothing in the
  hypergraph family, and the only cross-family call is
  `build_hypergraph()` -> `build_simplicial()`. Four families coexist with
  file prefixes, no tangle.

## Prerequisite — both trees are dirty. Commit before moving anything.

| Repo | State |
|---|---|
| hypernets | hypergraph fold-in **untracked** (13 R, 12 tests, 47 man); rename in flight (`hon.R` -> `memory_hon.R` x7); 48 modified; `Tutorial_docs/hon_inference.{Rmd,html}` deleted on disk but tracked at HEAD. HEAD has **zero** hypergraph exports. |
| texthypergraph | 44 files uncommitted — everything since 0.6.2 (HyperGAT, agreement API, projection tier). Version 0.6.4. |

Also note: until hypernets commits the fold-in, the only *committed* copy of the
hypergraph engine is `mohsaqr/hypernet_retired` (private). Do not delete that
repo until this migration is committed.

## Name collisions — 3 exports, 6 S3 methods

All three are the spectral trio texthypergraph migrated on 2026-08-25.

| Export | Resolution |
|---|---|
| `hypergraph_laplacian` | identical signature and values -> **keep hypernets'**, delete ours |
| `hypergraph_cluster` | identical signature and values -> **keep hypernets'**, delete ours |
| `hypergraph_transduction` | **KEEP OURS.** Ours has `normalization = c("none", "class_mass")`; hypernets' does not. That is the Zhu et al. (2003) class-mass fix, without which Zhou transduction predicts the majority class for every node under imbalanced seeds (R8: all 2,189 test docs -> "earn"). Dropping it silently regresses every benchmark. |

S3 methods colliding (`print`/`summary`/`as.data.frame` for
`net_hypergraph_cluster` and `net_hypergraph_transduction`): keep one set,
ours, since they are the tidier `as.data.frame(x, what =)` versions written
under Rule 0.

## Duplicate collapses — ~1,100 lines

Each pair is already verified numerically identical, so collapsing is a
deletion plus a parity test, never a rewrite.

| texthypergraph | hypernets | Keep |
|---|---|---|
| `R/spectral.R` (592) | `R/hypergraph_laplacian.R` (746) | hypernets', plus our `normalization` arg ported onto `hypergraph_transduction` |
| `R/pagerank.R` (196) | hypernets' B2 hypergraph PageRank | compare first — ours has `personalized =` and `sort_by =`; keep the superset |
| `hg_project(method = "clique")` | `R/hypergraph_expansion.R` (93) | both — `clique_expansion()` stays as the engine, `hg_project()` as the tidy verb, with an `identical()` parity test |
| `text_hypergraph(construction = "window")` | `R/hypergraph_window.R` (312) | both — ours tokenises text, hypernets' takes sequences; assert the reduction identity |

## File moves

Straight copies, prefixed `text_` where they are text-specific:

    R/text_hypergraph.R    -> R/text_hypergraph.R      (540, constructor)
    R/verbs.R              -> R/text_verbs.R           (340)
    R/hypergat.R           -> R/text_hypergat.R        (449)
    R/neural.R             -> R/text_neural.R          (265, HGNN — generic, but torch-gated)
    R/knn.R                -> R/hypergraph_knn.R       (105)
    R/stop_words.R         -> R/text_stopwords.R       (29)
    R/project.R            -> R/hypergraph_project.R   (186)
    R/edges.R              -> R/hypergraph_edges.R     (89)
    R/dual.R               -> R/hypergraph_dual.R      (53)
    R/null_test.R          -> R/hypergraph_null.R      (242)
    R/agreement.R          -> R/agreement.R            (207, generic — used by all families)
    R/sparse.R             -> R/hypergraph_sparse.R    (364)
    R/spectral.R           -> DELETE (except the normalization arg)
    R/pagerank.R           -> merge into hypernets' PageRank
    R/data.R               -> merged into R/data.R

Tests move alongside, one file per source file, names matching.

## DESCRIPTION union

    Imports:  ggplot2, grid, graphics, parallel, stats, utils,
              Matrix, methods, RSpectra          # Nestimate REMOVED
    Suggests: cograph, igraph, gridExtra, knitr, rmarkdown, testthat,
              sbert, torch
    Additional_repositories: https://mohsaqr.r-universe.dev   # for sbert

Title/Description need rewriting to cover four families including text.

## Non-R assets

- `data/`: no collision — `covid_abstracts`, `covid_embeddings` join
  `ai_long`, `human_long`. `data-raw/` moves too.
- `Tutorial_docs/`: 3 texthypergraph tutorials join the 2 ported from
  hypernets. Resolve the deleted `hon_inference.{Rmd,html}` first.
- `benchmarks/` (R8/R52/MR/Ohsumed/20NG harness) and
  `vignettes/articles/benchmarks.Rmd` move as-is; add `benchmarks/` to
  `.Rbuildignore`.
- `local_testing_and_equivalence/`: merge both, keeping every oracle. The
  Nestimate identity tests stay meaningful as a regression guard.
- `_pkgdown.yml`: reference index gains the text and projection sections;
  pkgdown CI fails on unindexed topics, so do this in the same commit.

## Order

1. Commit hypernets' fold-in + rename (resolve the `hon_inference` deletion).
2. Commit texthypergraph 0.6.5.
3. Bring texthypergraph in with history preserved (`git subtree add`), not a
   file copy.
4. Apply the collision and duplicate resolutions above.
5. `devtools::document()`, full suite, `R CMD check --as-cran`.
   Baselines to meet or beat: **hypernets 1993 pass / 0 fail**,
   **texthypergraph 491 pass / 0 fail**, both checks 0/0/0.
6. Update `CLAUDE.md`, `ROADMAP.md`, `TODO.md`; retire texthypergraph's repo
   the same way hypernets was (private + `_retired` suffix), and only then
   consider releasing `hypernet_retired`.

## Risks

- **Silent regression via the transduction arg** — the single most likely way
  to lose real capability. Guard it with the R8 benchmark before and after.
- **Check time**: ~14,300 R lines, ~2,484 tests, torch and sbert in Suggests.
  `--as-cran` will be slow; keep heavy paths behind `skip_on_cran()`.
- **One release cadence** for four unrelated families is the accepted cost of
  this decision.

## Not in scope

Moving `build_hypergraph()` out of hypernets (it needs `build_simplicial()` —
this is what killed hypernets), any new methods, and the remaining v0.7
roadmap items.
