# Spectral clustering of a hypergraph, as a tidy table

Calls the in-package
[`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md)
engine (Zhou et al. 2006 normalized Laplacian, or the Hayashi et al.
2020 random-walk Laplacian with edge-dependent vertex weights – the
natural choice for tf-idf-weighted text hypergraphs).

## Usage

``` r
hg_cluster(
  hg,
  k,
  type = c("zhou", "random_walk"),
  edge_weights = NULL,
  seed = NULL,
  nstart = 25L,
  what = c("clusters", "embedding", "eigenvalues"),
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
  (or any hypernets `net_hypergraph`).

- k:

  Number of clusters (explicit by design; there is no correct default).

- type:

  `"zhou"` or `"random_walk"`, as in
  [`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md).

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
  scan).

- n:

  For `what = "eigenvalues"`, how many rows to keep, leading first
  (default `Inf`, all of them). The spectrum carries one eigenvalue per
  node, so a corpus of a few thousand documents returns a few thousand
  rows, and only the leading ones carry the gap that decides how many
  groups the structure supports. Ignored for the other values of `what`.

- algorithm:

  `"spectral"` (RDC-Spec) or `"symnmf"` (RDC-Sym), as in
  [`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md).
  SymNMF currently requires a dense incidence matrix because its paper
  objective factorizes a dense node similarity.

- max_iter, tol:

  SymNMF convergence controls passed to
  [`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md).

## Value

A base `data.frame`. For `what = "clusters"`: one row per node, columns
`node` and `cluster`. For `what = "embedding"`: `node`, `cluster`, `pi`,
`dim1..dimk`. For `what = "eigenvalues"`: one row per eigenvalue,
columns `index`, `value` (ascending) and `gap` (the distance to the next
eigenvalue – large gaps indicate supported cluster counts; `NA` on the
last row).

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
#> 1     1 -5.551115e-17 0.07162811
#> 2     2  7.162811e-02 0.17318608
#> 3     3  2.448142e-01 0.23041019
#> 4     4  4.752244e-01         NA
```
