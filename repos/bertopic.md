# BERTopic — Python native + CRAN wrapper

- **Python (canonical)**: Grootendorst's transformer topic modeling —
  SBERT-class embeddings -> UMAP -> HDBSCAN -> class-based TF-IDF topic
  descriptions; modular (each stage swappable), widely replicated as
  better-coherence-than-LDA. The de-facto LDA replacement (tier 1).
- **R wrapper — verified 2026-08-24 via CRAN_package_db**: CRAN package
  `BERTopic` 0.1.0 — reticulate interface to the Python package
  (fit/transform, update/reduce topics, topic- and document-level info,
  interactive visualizations). NOT a native implementation; Python managed
  through reticulate.
- **Role for us**: the practical baseline any hypergraph-clustering claim
  must beat or complement. The actual-package benchmark remains open because
  the current environment has no Python `bertopic` installation. A ten-seed
  TF-IDF/SVD -> UMAP -> HDBSCAN comparison is retained under its accurate
  name as a separate baseline, not as a BERTopic result. UMAP/HDBSCAN are
  stochastic and BERTopic has no generative likelihood.
- **Links**: https://maartengr.github.io/BERTopic ·
  https://cran.r-project.org/package=BERTopic
