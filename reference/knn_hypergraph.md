# Build a k-nearest-neighbor hypergraph from embeddings

Every row of `embeddings` becomes one hyperedge containing the row
itself and its `k` nearest neighbors by cosine similarity (ties broken
deterministically by row name). With `weight = "cosine"` the incidence
entries are the similarities to the hyperedge's center (1 for the center
itself); `weight = "binary"` gives the unweighted membership hypergraph
(the `HyperG::knn_hypergraph()` construction, which this generalizes).

## Usage

``` r
knn_hypergraph(embeddings, k, weight = c("cosine", "binary"))
```

## Arguments

- embeddings:

  Numeric matrix, one row per item, with unique non-empty rownames (the
  item IDs). No missing values; no all-zero rows.

- k:

  Number of neighbors per hyperedge (between 1 and
  `nrow(embeddings) - 1`).

- weight:

  `"cosine"` (default) or `"binary"`.

## Value

A `net_hypergraph` (from
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md))
with one hyperedge per item, each of size `k + 1`, plus a `knn` field
recording `k` and the weighting. Accepted by
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
[`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
and
[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md).

## Details

Cosine similarity and Euclidean distance order neighbors identically on
L2-normalized embeddings, so normalized encoder output (the default of
`sbert::encode()`) gives the conventional kNN structure.

## Conditions

Raises `hypernets_bad_input` for broken contracts, and
`hypernets_nonpositive_similarity` when `weight = "cosine"` selects a
neighbor with similarity \<= 0 – a non-positive incidence weight would
invalidate the random-walk machinery downstream; use `weight = "binary"`
or a smaller `k` instead.

## Examples

``` r
set.seed(1)
emb <- matrix(rnorm(20, mean = 3), nrow = 5,
              dimnames = list(paste0("doc_", 1:5), NULL))
hg <- knn_hypergraph(emb, k = 2)
hg_measures(hg, what = "edges")
#>    edge size
#> 1 doc_1    3
#> 2 doc_2    3
#> 3 doc_3    3
#> 4 doc_4    3
#> 5 doc_5    3
```
