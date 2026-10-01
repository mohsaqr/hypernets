# workinprogress/

Code taken out of the shipped package but kept for later. Build-ignored
(`^workinprogress$` in `.Rbuildignore`); nothing here is loaded, exported,
tested or documented by the package.

## hon_outcome/ (sidelined 2026-09-29, from 0.5.1)

`hon_outcome()` and its `net_outcome` class: regress an actor-level outcome on
per-actor features, with CR1/CR3 cluster-robust SEs and BH correction.

**Why sidelined:** its only feature construction -- the visit-weighted mean of
a higher-order node value (centrality or HONEM coordinate) over an actor's
visits, called "exposure" -- has no published basis. No reference defines
it, and "exposure" already means something else in network science
(Valente's network exposure). The author's rule: no method ships without a
reference.

**State when sidelined:** working and tested (168 tests pass). CR3 for GLMs
fixed to match `sandwich::vcovCL(type = "HC3")`; per-SD effects and
25th/75th-percentile contrasts added; review page in `review/`.

**To revive:** replace the exposure features with a referenced construction
(e.g. per-actor state distribution, Gabadinho et al. 2011), move `R/` back to
`R/memory_outcome.R` and `tests/` to `tests/testthat/`, re-add `sandwich` to
Suggests and `hon_outcome` to `_pkgdown.yml`, run `devtools::document()`.

## hon_per_actor/ (removed 2026-09-29, never released)

`build_hon(by =)` returning `net_hon_group`: an independent BuildHON model per
unit (actor, project, session). **Rejected, not just sidelined:** on
`human_long` a median session (17 transitions) got an 8-node / 14-edge model
and 234 of 425 sessions were given second- to fourth-order memory from that
little data -- noise, not memory. Per-unit models also have different node
sets, so units cannot be compared. Nestimate never did this: its sessions are
the *trajectories* of one pooled model (transitions never cross a session
boundary), which is what `build_hon(actor = "session_id")` already does.

Kept: `R/memory_hon_group.R`, its tests, the frozen `build_hon()` fixture,
`build_hon_by.patch` (the `build_hon()` / `.hi_parse()` edits), docs, review page.

**What may be wanted instead** (author, 2026-09-29): a *group* model -- one
HON per group (condition, cohort, human vs AI), where each group has enough
data to estimate memory; the higher-order analogue of TNA's group models.
Reuse the per-unit fitting loop here, but with groups, not individuals, and
find a reference before building.

## hypergraph global clustering C2(H) (Estrada & Rodriguez-Velazquez 2006) -- sidelined 2026-09-29

**What it is.** A global clustering coefficient for hypergraphs: triangles of
the clique expansion not contained in a single hyperedge, over 2-paths minus
three per triple inside one hyperedge; and its dual C2(H*).

**Why it is out.** The formula is well defined and published (Physica A 364,
581-594, doi:10.1016/j.physa.2005.12.002), but (1) no library implements it,
so there is no independent oracle (ROADMAP rule 2), and (2) the paper's own
worked example (Fig. 1, "C2 = 0.25 of 18 triples") could not be rebuilt from
the text, so it cannot serve as the oracle either.

**How to revive.** Find the authors' code, a later paper that reproduces a
value, or reconstruct Fig. 1's hypergraph exactly; then implement it as
`hg_transitivity(type = "estrada")` next to the local variants in
`R/hypergraph_transitivity.R`.
