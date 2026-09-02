# DHG — DeepHypergraph (Python, PyTorch)

- **What**: the packaged home of hypergraph neural networks (iMoonLab):
  HGNN, HGNN+, HyperGCN, UniGNN-family and more as ready PyTorch models,
  plus hypergraph data structures and spectral utilities.
- **Verified 2026-09-02**: official repository and documentation list HGNN,
  HGNN+, HyperGCN, HNHN and the UniGNN family among the high-order models.
- **Relatives**: PyTorch Geometric ships a `HypergraphConv` layer;
  TopoNetX / TopoModelX (pyt-team) cover topological deep learning
  (simplicial/cell/hypergraph message passing).
- **Role for us**: HGNN is checked against DHG's `HGNNConv`. honets now also
  ships HyperGCN (dynamic, fast and one-edge mediator variants) and HNHN
  (paper alpha/beta normalization). Direct paper-equation fixtures are the
  durable oracle; isolated DHG parity fixtures remain useful follow-up work.
- **Links**: https://github.com/iMoonLab/DeepHypergraph ·
  https://deephypergraph.readthedocs.io · `pip install dhg`
