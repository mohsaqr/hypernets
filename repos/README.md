# Repo index — oracles vs guiding references

One line per surveyed repo/package (each has a full note in this folder).
Compiled 2026-08-25. "Verified" = source/docs checked on the stated date;
unverified rows are background knowledge — re-check before relying on a
specific function.

## Ours (anchors)

| Repo | What | Note |
|---|---|---|
| Nestimate hypergraph module | The anchor all comparisons target: `build_hypergraph`, `bipartite_groups`, `hypergraph_measures`, `hypergraph_centrality` (CEC + tensor Z/H), `clique_expansion`, + the Zhou/Hayashi Laplacian trio (2026-08-25, HyperNetX parity < 1e-12) | `nestimate-hypergraph-module.md` |
| sbert (R, ../SBERT) | Native Python-free embedding front-end for the kNN-hypergraph pipeline (replaces the `text`/reticulate route); 13 ONNX + 1 static pinned models, SentenceTransformers parity ~1e-7 | — (own repo docs) |

## Equivalence oracles (numerical parity targets; reticulate/local-only, never dependencies)

| Repo | Lang | Oracle for | Verified |
|---|---|---|---|
| HyperNetX (PNNL) | Py | **THE** oracle for the Zhou/Hayashi weighted Laplacian trio (`prob_trans`/`get_pi`/`norm_lap`/`spec_clus` implement Hayashi EDVW exactly, author-adjacent); already used, parity < 1e-12. Homology module = possible extra TDA oracle | 2026-08-25 |
| XGI | Py | Second independent oracle for shipped `hypergraph_centrality()` (`h_/z_/clique_eigenvector_centrality`, Katz) — cross-check parked in TODO. Also the measure-naming reference | 2026-08-25 |
| HyperG (CRAN) | R | Unweighted cases: `hypergraph_laplacian_matrix`/`cluster_spectral` (Laplacian item), `knn_hypergraph` (kNN construction), `sample_*` (Saqrlab samplers), `dual_hypergraph` (dual accessor) | 2026-08-24 |
| SimplicialComplex (CRAN) | R | Parked `wasserstein_distance()`; second oracle for shipped `bottleneck_distance`. Name collisions with Nestimate are benign (last-loaded wins) | 2026-08-24 |
| gudhi / ripser.py / giotto-tda | Py | Candidate oracles for `wasserstein_distance()` (gudhi = field reference); Nestimate already uses `TDAstats` | background |
| pathpy / pyHON / HYPA / HONEM | Py | **Already integrated** in honets/Nestimate `local_testing_and_equivalence/` (HON sequence side; MOGen matches pathpy at machine precision). Nothing new to adopt | integrated |

## Guiding / reference repos (design, construction, or baseline — not parity targets)

| Repo | Lang | Guides | Verified |
|---|---|---|---|
| HyperGAT_TextClassification | Py | Construction reference for windowed sequence hyperedges (`utils.py::get_slice()`: sentence = hyperedge + LDA-topic edges). **Sliding-window variant NOT in the code** — our windowed item has no upstream oracle, invariant-gated instead | 2026-08-24 |
| BERTopic (Py + CRAN wrapper) | Py/R | The practical baseline any hypergraph-clustering claim must beat or complement (stochastic — multi-seed reporting required) | 2026-08-24 (wrapper) |
| text (CRAN) | R | Was the assumed embedding step (reticulate → HF); **superseded by sbert** for our pipeline. Remains the only R route to XLNet-class PLMs | 2026-08-24 (presence) |
| AllSet (official) | Py | Neural tier, reference-only; strongest general hypergraph-NN framework if carm-ml ever grows one; its 10-dataset benchmark harness is a design reference | paper 2026-08-24; repo unaudited |
| DHG / DeepHypergraph (+ PyG `HypergraphConv`, TopoNetX) | Py | Neural tier, reference-only: the packaged HGNN/HyperGCN/UniGNN zoo to study if a GNN story lands | background |
| HypergraphX (Battiston) | Py | Tertiary reference; candidate extra oracle for communities/motifs if those ever land — XGI + HyperNetX cover current needs | background |
| [Legal Hypergraphs](legal-hypergraphs.md) | R/data | Edge measures, association projection, s-line centralities, configuration nulls and temporal analysis; published reproducibility archive is the planned end-to-end oracle | paper + archive 2026-09-02 |
