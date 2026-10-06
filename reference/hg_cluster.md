# Spectral or symmetric-NMF clustering of hypergraph vertices

Partitions the nodes of a hypergraph into `k` clusters with either of
Hayashi et al.'s (2020) representative-digraph algorithms.
`algorithm = "spectral"` (RDC-Spec) row-normalizes the `k` smallest
Laplacian eigenvectors and applies k-means. `algorithm = "symnmf"`
(RDC-Sym) computes a rank-`k` non-negative factorization `T ~= U U'` of
the normalized similarity `T = I - L`, then assigns each vertex to the
largest entry in its row of `U`, exactly as Algorithm 2 specifies. With
`type = "random_walk"` and a weighted incidence (e.g. from
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
with `weight =`), the edge-dependent vertex weights genuinely change the
partition - with edge-independent weights the walk collapses to a graph
random walk (Chitra & Raphael 2019).

## Usage

``` r
hg_cluster(
  hg,
  k,
  type = c("zhou", "random_walk"),
  edge_weights = NULL,
  seed = NULL,
  nstart = 25L,
  what = c("clusters", "embedding", "eigenvalues", "membership"),
  n = Inf,
  algorithm = c("spectral", "symnmf"),
  max_iter = 500L,
  tol = 1e-06
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hg`).

- k:

  Number of clusters (explicit by design; there is no correct default).

- type:

  The Laplacian, as in
  [`hg_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hg_laplacian.md):
  `"zhou"` (Zhou et al. 2006) or `"random_walk"` (Hayashi et al. 2020).

- edge_weights:

  Hyperedge weights for the Laplacian: `NULL` (default, every hyperedge
  counts once), a positive numeric vector with one entry per hyperedge,
  or `"idf"` to weight each word hyperedge by its smoothed inverse
  document frequency (needs `weight = "tfidf"`), so a word used in most
  documents binds them less than a rare shared word. Note that the
  default `type = "zhou"` reads only which documents a hyperedge binds;
  the incidence cells (`weight =`) enter the `"random_walk"` type only.

- seed:

  Random seed passed to the selected solver; set it for a reproducible
  partition.

- nstart:

  Number of solver starts (default `25L`).

- what:

  What to return: `"clusters"` (default) for the partition,
  `"embedding"` for the partition plus the row-normalized spectral
  embedding used by k-means (`dim1..dimk` – plot these to map the
  corpus) and the stationary weight `pi`, or `"eigenvalues"` for the
  Laplacian spectrum (dense engines return all `n` values; the sparse
  engine returns the `k + 1` it computed – raise `k` for an eigengap
  scan), or, with `algorithm = "symnmf"`, `"membership"` for the graded
  membership of every node in every cluster: the node's row of the
  non-negative factor, normalised to sum to one (Kuang, Ding & Park
  2012).

- n:

  For `what = "eigenvalues"`, how many rows to keep, leading first
  (default `Inf`, all of them). The spectrum carries one eigenvalue per
  node, so a corpus of a few thousand documents returns a few thousand
  rows, and only the leading ones carry the gap that decides how many
  groups the structure supports. Ignored for the other values of `what`.

- algorithm:

  Character. `"spectral"` (default; RDC-Spec) or `"symnmf"` (RDC-Sym).

- max_iter:

  Maximum multiplicative-update iterations for `algorithm = "symnmf"`.

- tol:

  Relative objective tolerance for `algorithm = "symnmf"`.

## Value

A base `data.frame`. For `what = "clusters"`: one row per node, columns
`node` and `cluster`. For `what = "embedding"`: `node`, `cluster`, `pi`,
`dim1..dimk`. For `what = "eigenvalues"`: one row per eigenvalue,
columns `index`, `value` (ascending) and `gap` (the distance to the next
eigenvalue – large gaps indicate supported cluster counts; `NA` on the
last row). For `what = "membership"`: one row per node and cluster,
columns `node`, `cluster` and `membership` (summing to one over a node's
clusters); a node's largest membership is its cluster in
`what = "clusters"`.

## Details

Both solvers are stochastic: `nstart` initializations are used and a
`seed` fixes the result. Report stability across seeds for consequential
results.

`type = "random_walk"` (Hayashi et al. 2020) is the natural choice for
tf-idf-weighted text hypergraphs. `what =` returns the cluster table,
the embedding or the leading eigenvalues as a tidy data.frame.

## References

Kuang, D., Ding, C., & Park, H. (2012). Symmetric nonnegative matrix
factorization for graph clustering. *Proceedings of the 2012 SIAM
International Conference on Data Mining*, 106-117.
[doi:10.1137/1.9781611972825.10](https://doi.org/10.1137/1.9781611972825.10)

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
edge-dependent vertex weights. *ICML 2019*.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
hg_cluster(hg, k = 2, type = "random_walk", seed = 1)
#>        node   cluster
#> 1 cooking_1 Cluster 1
#> 2 cooking_2 Cluster 1
#> 3   space_1 Cluster 2
#> 4   space_2 Cluster 2
hg_cluster(hg, k = 2, seed = 1, what = "embedding")
#>        node   cluster   pi      dim1       dim2
#> 1 cooking_1 Cluster 1 0.20 0.5695852  0.8219323
#> 2 cooking_2 Cluster 1 0.30 0.8513708  0.5245643
#> 3   space_1 Cluster 2 0.25 0.6550598 -0.7555770
#> 4   space_2 Cluster 2 0.25 0.8037070 -0.5950253
hg_cluster(hg, k = 2, seed = 1, what = "eigenvalues", n = 5)
#>   index         value        gap
#> 1     1 -2.081668e-17 0.07162811
#> 2     2  7.162811e-02 0.17318608
#> 3     3  2.448142e-01 0.23041019
#> 4     4  4.752244e-01         NA
```
