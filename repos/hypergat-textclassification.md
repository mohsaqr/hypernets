# HyperGAT_TextClassification — official code (Python, PyTorch)

- **What**: official implementation of Ding et al. (2020), "Be More with
  Less: Hypergraph Attention Networks for Inductive Text Classification"
  (`../papers/2020-EMNLP-HyperGAT-Ding.pdf`). Research code, not a package.
- **Verified 2026-08-24** (reading the source): hypergraph construction lives
  in `utils.py::get_slice()` — incidence matrix with each SENTENCE as one
  hyperedge (word-row, sentence-column), LDA topic hyperedges appended as
  extra incidence columns when enabled (`generate_lda.py` precomputes
  keywords). **The paper's sliding-window alternative is NOT implemented** —
  the honets sliding-window construction is our extension with no
  upstream oracle.
- **Stack**: PyTorch 1.4, Python 3.6 era; datasets 20NG/R8/R52/Ohsumed/MR.
  Benchmark numbers we quote are from its paper's Table 2 (read 2026-08-25).
- **Role for us**: the sentence-only dual-attention model is shipped as
  `hg_hypergat()`/`text_hypergat()` with forward parity. The remaining work
  is its `--use_LDA` semantic-hyperedge path and the 20NG/MR/Ohsumed runs.
- **Links**: https://github.com/kaize0409/HyperGAT_TextClassification
