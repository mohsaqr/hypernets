# ROADMAP — analytical capability

Status 2026-09-20. Companion to `ROADMAP-text.md` (constructions) and
`EXPANSION-PLAN.md` (family structure). Those two describe what hypernets can
*represent*. This one describes what it can *conclude*, which is currently
almost nothing.

## The gap, with evidence

hypernets 0.4.8 exports 112 verbs. Grepping `R/` for `lm(`, `glm`, `lmer`,
mixed effects, cluster-robust inference or `predict` returns **nothing** — the
apparent hits are prose ("outcome states", "valid outcomes"). No verb anywhere
accepts an outcome variable or a covariate.

The whole inferential surface is four verbs, and none relates structure to
anything outside itself:

| Verb | Question it answers |
|---|---|
| `bootstrap_hon()` | are these edges stable under resampling? |
| `compare_hon(x, y)` | do **exactly two** groups differ? |
| `hg_null_test()` | does this structure differ from a random null? |
| `markov_order_test()` | how far back does memory reach? |

So the package characterises higher-order structure in great detail and cannot
answer the questions its users actually have: *does this pattern predict the
outcome? does it differ across conditions, controlling for prior attainment,
with students nested in classes? which people show which structure?*

It is descriptive and unsupervised end to end. Everything upstream already
produces the inputs an inferential layer would need — `build_honem()` returns
per-node vector features, `hg_embed()` spectral embeddings, `hg_cluster()`
per-document assignments — and all of it dead-ends.

## Design constraints (non-negotiable, from CLAUDE.md)

1. **Taxonomy.** `build_*()` constructs and returns `net_*`; measures return a
   tidy one-row-per-observation `data.frame`; inference returns a `net_*`
   carrying estimates, intervals and p-values. Alternative routes are chosen
   with `type =`, never `method =`.
2. **Rule 0.** One call, named arguments, tidy return. If a workflow forces the
   caller to `merge()`/`order()`/`split()` by hand, the API is wrong — add the
   verb. Every `net_*` gets `print`, `summary`, `plot`, `as.data.frame`.
3. **Statistics.** Effect sizes with confidence intervals, never bare
   p-values. `p.adjust(method = "BH")` whenever many tests are computed, and
   say so in the result. Nested data gets a model that respects the nesting.
   Every stochastic result reported from multiple seeds.
4. **Conditions** are `hypernets_*`; reuse an existing class over minting a
   near-synonym.
5. **Dependencies.** Base R first. A mixed-effects backend is the one place a
   new Suggests may be justified; it must be guarded at every use site.

## Tier A — structure to outcome

The multiplier. Converts the package from "describes structure" to "explains
learning".

```r
hon_outcome(fit, outcome = scores, by = "student",
            features = c("centrality", "honem"), nested_in = "class")
```

- **Input**: a `net_hon` / `net_mogen` / `net_hypa`, or a hypergraph, plus an
  outcome keyed by actor. Features come from the verbs that already exist —
  `hon_centrality()`, `build_honem()`, `hg_embed()`, `hg_centrality()` — so
  this is an assembly and inference layer, not new mathematics.
- **Return**: `net_outcome`, whose `as.data.frame()` is one row per feature
  with `estimate`, `conf_low`, `conf_high`, `p`, `p_adj`, `n`. `summary()`
  reports model fit and what was corrected for; `plot()` is a coefficient
  plot with intervals.
- **Inference**: OLS/GLM by default; `nested_in =` switches to cluster-robust
  standard errors (base R, no dependency) or a mixed model where a backend is
  available and permitted. Multiplicity across features is BH-corrected and
  the correction is named in the result.
- **Naming**: `hon_outcome()` / `hg_outcome()` pairs with the existing
  `hon_centrality()` / `hg_centrality()`, and keeps the noun "outcome" that
  `codyna` already uses (`analyzeOutcome`). **Open decision** — the taxonomy
  says inference verbs are `bootstrap_*`, `compare_*`, `*_test`, and an
  outcome model is none of those.
- **Verification**: coefficients and standard errors must match `stats::lm` /
  `stats::glm` on the same design matrix to `< 1e-10`; cluster-robust SEs
  against a hand-computed sandwich estimator and against `sandwich::vcovCL`
  where installed; a simulation showing nominal coverage of the intervals at
  a known effect size over multiple seeds.

**Risk to name honestly.** The features are high-dimensional and collinear,
the sample is actors rather than transitions, and it is easy to ship a verb
that produces confident nonsense. The coverage simulation is not optional
garnish — it is the thing that makes this verb publishable rather than
dangerous.

## Tier B — comparison beyond two pooled groups

`compare_hon(x, y, ...)` is hard-wired to two groups, pooled, no covariates.
Learning data is many groups, nested, with covariates. Pooled two-group
permutation is the wrong model for students in classes, and CLAUDE.md's own
statistics rules say so.

```r
compare_hon(fit, by = "condition", covariates = ~ prior_score,
            nested_in = "class")
```

- Generalise `by =` to k groups; keep the present two-argument form working.
- Permutation respects the nesting: shuffle labels **within** clusters, not
  across them, or the null is wrong.
- Return keeps the existing `net_*` shape with an added contrast column;
  BH across contrasts.
- **Verification**: with two groups and no covariates, results must be
  `identical()` to today's `compare_hon()` — a frozen regression fixture. The
  within-cluster permutation gets the same uniformity check that
  `markov_order_test()` already has (KS against U(0,1) under a true null).

## Tier C — idiographic higher-order models

Everything currently pools across people, discarding the variation learning
research exists to study. One model per actor, then compare, cluster and
relate those models to outcomes. The sibling `idiographic` package sets the
precedent.

```r
fits <- build_hon(sequences, by = "student")   # net_hon_group
as.data.frame(fits)                            # one row per student per edge
hg_cluster(fits)                               # types of learner
hon_outcome(fits, outcome = scores)            # Tier A consumes it directly
```

- A `net_hon_group` container with the usual four methods, so per-actor models
  are first-class rather than a list the user loops over.
- Sparse-data guards: most actors have short sequences. Report per-actor
  support and raise `hypernets_insufficient_support` rather than returning a
  confidently wrong per-person model.

## Tier D — the text-to-sequence bridge (prerequisite for the workflow)

Demonstrated end to end on 2026-09-20 and found to work **except for the
join**. `text_hypergraph()` already carries document covariates through
(`student`, `turn` survive into `as.data.frame(hg, what = "documents")`), and
`hg_cluster()` returns one row per document with its topic. Getting from there
to per-actor topic sequences currently requires the caller to write the join
by hand: merge the documents table onto the cluster assignments, sort the
result by actor and then by turn, and split the topic column by actor. Three
lines, every one of them the bracket-and-`order()` ritual Rule 0 forbids on a
public surface — which is the argument for the verb below. (Reproduced in
`tmp/la_workflow.R`, outside the public surface, for the round-trip fixture.)

```r
hg_sequences(hg, topics, actor = "student", order_by = "turn")
```

- **Return**: the long sequence table `actor` / `time` / `action`. CORRECTED
  2026-09-20 after implementation — the first draft of this section was wrong
  twice, and both errors were caught by running the code:
  - Only `build_hon()`, `bootstrap_hon()`, `compare_hon()` and
    `window_hypergraph()` take `action`/`actor`/`time`. `build_mogen()`,
    `build_hypa()`, `markov_order_test()` and `path_dependence()` do **not**.
    Giving those four the long-form arguments is a real gap, and memory-family
    surgery beyond this tier.
  - A bare `build_hon(seqs)` does **not** work and does **not** error. It
    silently reads the long table as a wide one, one trajectory per row, so
    `actor` ids and timestamps become states: a 2-actor x 4-turn table yields
    states `1, 2, 3, 4, A, B, C, s1, s2` and 8 "trajectories" instead of 3
    states and 2. See the defect below.
- This closes the one missing edge in the cross-family design: `pathways()`
  bridges memory to the other families and `window_hypergraph()` bridges
  sequences to hypergraphs, but the text family dead-ends at clustering.
- **Verification**: round-trip against the hand-written join above on a fixed
  corpus (`identical()`); sequence construction cross-checked against
  `TraMineR` / `tna` / `codyna` rather than a Python oracle.

## Defect found while building Tier D — silent misread of a long table

`build_hon()` cannot tell a long sequence table from a wide one and assumes
wide, without an error or a warning. Verified 2026-09-20:

| call | nodes | trajectories |
|---|---|---|
| `build_hon(long, action =, actor =, time =)` | 3 | 2 |
| `build_hon(long)` | 9 | 8 |

The bare call promotes timestamps and actor ids into behavioural states. This
violates the package's own "no silent failure" rule and it is reachable by an
ordinary user mistake — more reachable once `hg_sequences()` hands people a
long table.

The fix belongs in `.coerce_sequence_input()` (`R/utils.R`), which is **shared
with Nestimate under the provenance contract** ("keep them in sync"), so it is
an author decision, not a drive-by change. Options: detect an `actor`/`time`/
`action`-shaped frame and route it long; or refuse a frame whose columns match
those names unless the arguments are given. Either way both copies move
together, and the memory identity test now in
`local_testing_and_equivalence/test-identity-nestimate-memory.R` will catch a
one-sided change.

## Tier E — Markov order and the memory family's loose ends

Small next to A–D, but they sit on the same code and are cheap.

1. **`markov_stability()` is missing.** `Nestimate` 0.9.12 exports it and the
   JS `tnaj` has `markovStability`; hypernets is the only one of the three
   without it. Decide: port it, or state that hypernets deliberately does not
   carry it. `extract_pathways` is in the same position.
2. **R↔JS equivalence for the order test.** `hypernets::markov_order_test()`
   and `tnaj::markovOrderTest()` are the same method in two languages on one
   machine, never compared. CLAUDE.md mandates cross-language numerical
   equivalence, and the workspace already has a `validation/` runner
   (R→JSON→TS) built for exactly this.
3. **The cross-package identity test does not run.**
   `Nestimate/local_testing_and_equivalence/test-equiv-honets.R:12` is
   `skip_if_not_installed("honets")` — a package renamed to `hypernets` at
   0.4.0 — so it has skipped silently since the rename. Meanwhile Nestimate
   touched the memory family on 2026-09-13 (0.9.5) and hypernets' copy last
   moved 2026-09-06. Re-home the test **into hypernets** comparing against
   installed `Nestimate`, rather than editing Nestimate (forbidden from here).
   A read-only diff on 2026-09-20 found `markov_order.R` differing only in
   `invisible(out)` and hypernets' added `as.data.frame.net_markov_order()` —
   no mathematics — but the other five memory files were not audited.
4. **The equivalence gate is misdocumented.** `CLAUDE.md:114` says
   `HYPERNETS_EQUIV_TESTS`; `helper-equiv-gate.R` reads `HONETS_EQUIV_TESTS`.
   Following the documentation skips all 31 oracle files and reports
   `[ FAIL 0 | PASS 0 ]`. One-line fix, and until it is fixed nobody knows
   what that suite says.

## Sequencing

| Order | Item | Why here |
|---|---|---|
| 1 | E4, then run the suite | One line, and it tells us whether anything is already broken before new work lands on top |
| 2 | **Tier A** | The multiplier; its inputs already exist and go nowhere |
| 3 | **Tier D** | Small, and it makes Tier A reachable from a transcript |
| 4 | **Tier B** | Generalises an existing verb; fixture-protected |
| 5 | E3 | Re-home the identity test once the memory family is being touched again |
| 6 | **Tier C** | Largest surface; benefits from A and B being settled |
| 7 | E1, E2 | Decide `markov_stability`; R↔JS equivalence |

## Open questions for the author

- Naming: `hon_outcome()` / `hg_outcome()`, or fold into the existing
  `compare_*` grammar? The taxonomy does not currently have a slot for a
  supervised verb.
- Mixed-effects backend: accept a Suggests dependency, or ship cluster-robust
  standard errors only and stay base-R? Cluster-robust covers most learning
  designs and costs nothing.
- Does Tier C belong here or in `idiographic`? The sibling exists; duplicating
  it would repeat the Nestimate divergence this roadmap is trying to close.
- Scope: are hypergraph-side outcomes (`hg_outcome()`) wanted in the same pass,
  or is the memory family enough to prove the design?
