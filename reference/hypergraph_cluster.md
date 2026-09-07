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
hypergraph_cluster(
  hg,
  k,
  type = c("zhou", "random_walk"),
  edge_weights = NULL,
  nstart = 25L,
  seed = NULL,
  algorithm = c("spectral", "symnmf"),
  max_iter = 500L,
  tol = 1e-06
)
```

## Arguments

- hg:

  A connected `net_hypergraph`.

- k:

  Integer number of clusters, `2 <= k <= n_nodes - 1`.

- type, edge_weights:

  Passed to
  [`hypergraph_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_laplacian.md).

- nstart:

  Integer. Random restarts (default 25).

- seed:

  Optional integer seed for initialization.

- algorithm:

  Character. `"spectral"` (default; RDC-Spec) or `"symnmf"` (RDC-Sym).

- max_iter:

  Maximum multiplicative-update iterations for `algorithm = "symnmf"`.

- tol:

  Relative objective tolerance for `algorithm = "symnmf"`.

## Value

An object of class `net_hypergraph_cluster`: a list with `$clusters`
(data.frame, one row per node: `node`, `cluster` - labels `"Cluster 1"`,
`"Cluster 2"`, ... ordered by first appearance), `$embedding` (node x k
row-normalized spectral embedding used by k-means, dims `dim1..dimk`),
`$k`, `$type`, `$eigenvalues` (full Laplacian spectrum, increasing),
`$eigengap` (gap after the k-th eigenvalue), `$sizes` (data.frame
`cluster`/`size`), `$pi` (stationary distribution) and `$params`. Has
`print`, `summary`, `plot` and `as.data.frame` methods;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
one row per node with `node`, `cluster`, the stationary probability
`pi`, and the embedding coordinates.

## Details

Both solvers are stochastic: `nstart` initializations are used and a
`seed` fixes the result. Report stability across seeds for consequential
results.

## References

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
edge-dependent vertex weights. *ICML 2019*.

## Examples

``` r
events <- data.frame(
  person = c("a", "b", "c", "a", "b", "c", "d", "e", "f",
             "d", "e", "f", "c", "d"),
  meeting = c("m1", "m1", "m1", "m2", "m2", "m2", "m3", "m3", "m3",
              "m4", "m4", "m4", "m5", "m5")
)
hg <- group_hypergraph(events, actor = "person", group = "meeting")
cl <- hypergraph_cluster(hg, k = 2, seed = 1)
cl
#> Hypergraph RDC-Spec spectral clustering (zhou Laplacian)
#>   Nodes: 6 | Hyperedges: 5 | k: 2
#>   Cluster sizes: Cluster 1 = 3, Cluster 2 = 3
#>   Eigengap after k: 0.6667
as.data.frame(cl)
#>   node   cluster        pi       dim1       dim2
#> 1    a Cluster 1 0.1428571 -0.6575959  0.7533708
#> 2    b Cluster 1 0.1428571 -0.6575959  0.7533708
#> 3    c Cluster 1 0.2142857 -0.7947194  0.6069770
#> 4    d Cluster 2 0.2142857 -0.7947194 -0.6069770
#> 5    e Cluster 2 0.1428571 -0.6575959 -0.7533708
#> 6    f Cluster 2 0.1428571 -0.6575959 -0.7533708
```
