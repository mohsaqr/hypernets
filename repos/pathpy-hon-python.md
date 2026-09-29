# HON-family Python ecosystem — pathpy / pyHON / HYPA / HONEM

One note for the *sequence* higher-order side (the hypernets paradigm, distinct
from the hypergraph packages in this folder):

- **pathpy** (Scholtes group) — the only real package (`pip install pathpy`):
  path data, `HigherOrderNetwork`, `MultiOrderModel` (MOGen's basis), model
  selection; 2.2.0 is still the PyPI release and the oracle version.
  **pathpyG** is the newer PyTorch successor: latest release v0.2.0-alpha
  (GitHub pre-release, 2025-10-02), **not on PyPI** (checked 2026-09-29);
  not used by any hypernets test.
- **pyHON** (github.com/xyjprc/hon) — Xu et al. reference BuildHON/BuildHON+;
  research code, not a package.
- **HYPA** (github.com/tlarock/hypa) — LaRock et al. reference; research code.
  **Not used as an oracle yet**: `build_hypa()` is checked against a
  Monte-Carlo configuration-model null and `BiasedUrn::pMWNCHypergeo`
  (Wallenius), in `local_testing_and_equivalence/test-equiv-hypa.R`.
- **HONEM** — Saebi et al. research code.
- No Python package implements a permutation-based Markov-order test
  (`markov_order_test()` has no upstream equivalent).

**Status for us**: already integrated — hypernets' local equivalence
suite uses pyHON, pyMOGen, and pathpy as oracles via reticulate
(`helper-python-equiv.R`; MOGen matches pathpy.MultiOrderModel at machine
precision). `hon_communities()` is checked against the Python `infomap`
package (tested with 2.15.1; see [`infomap.md`](infomap.md)). The parity
proofs live in `hypernets/local_testing_and_equivalence/`; the HYPA reference
code is the one planned addition (ROADMAP.md Phase 1).
