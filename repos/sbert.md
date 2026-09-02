# sbert — native embedding frontend

- **Relationship**: suggested dependency used only when
  `text_hypergraph(construction = "knn")` must compute embeddings. Callers can
  always pass a numeric embedding matrix and run offline without sbert.
- **Boundary**: sbert owns encoders and embedding inference; honets owns the
  kNN hypergraph and all downstream analysis.
- **Collision audit 2026-09-02**: no exported name overlaps with honets.
- **Version checked locally**: 0.5.2.

