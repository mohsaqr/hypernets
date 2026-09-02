# HON-family Python ecosystem — pathpy / pyHON / HYPA / HONEM

One note for the *sequence* higher-order side (the honets paradigm, distinct
from the hypergraph packages in this folder):

- **pathpy** (Scholtes group) — the only real package (`pip install pathpy`):
  path data, `HigherOrderNetwork`, `MultiOrderModel` (MOGen's basis), model
  selection. **pathpyG** is the newer PyTorch successor.
- **pyHON** (github.com/xyjprc/hon) — Xu et al. reference BuildHON/BuildHON+;
  research code, not a package.
- **HYPA** (github.com/tlarock/hypa) — LaRock et al. reference; research code.
- **HONEM** — Saebi et al. research code.
- No Python package implements a permutation-based Markov-order test
  (`markov_order_test()` has no upstream equivalent).

**Status for us**: already integrated — honets' local equivalence
suite uses pyHON, pyMOGen, and pathpy as oracles via reticulate
(`helper-python-equiv.R`; MOGen matches pathpy.MultiOrderModel at machine
precision). Nothing new to adopt; the parity proofs live in
`honets/local_testing_and_equivalence/`.
