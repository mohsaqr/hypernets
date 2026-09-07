# Low-Dimensional Hypergraph Embedding

Returns the node coordinates computed by honets' existing hypergraph
spectral or symmetric-NMF engine without exposing the incidental k-means
assignments produced by
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md).
This is the direct embedding verb; `hg_embed()` and `hypergraph_embed()`
are identical names for it.

## Usage

``` r
hg_embed(
  hg,
  dimensions = 2L,
  type = c("zhou", "random_walk"),
  method = c("spectral", "symnmf"),
  seed = NULL,
  nstart = 25L,
  max_iter = 500L,
  tol = 1e-06
)

hypergraph_embed(
  hg,
  dimensions = 2L,
  type = c("zhou", "random_walk"),
  method = c("spectral", "symnmf"),
  seed = NULL,
  nstart = 25L,
  max_iter = 500L,
  tol = 1e-06
)
```

## Arguments

- hg:

  Any honets `net_hypergraph`.

- dimensions:

  Number of embedding coordinates, between 2 and `n_nodes - 1`.

- type:

  Laplacian type: `"zhou"` or `"random_walk"`.

- method:

  `"spectral"` (default) or Hayashi et al.'s `"symnmf"`. SymNMF requires
  a dense incidence matrix.

- seed:

  Optional initialization seed.

- nstart:

  Number of k-means starts for the shared engine; for spectral
  embeddings this affects only the discarded cluster assignment, not the
  returned coordinates. For SymNMF it is the number of factor restarts.

- max_iter, tol:

  SymNMF iteration controls.

## Value

A data frame with `node`, stationary mass `pi`, and `dim1` through
`dim<dimensions>`.

## Examples

``` r
h <- group_hypergraph(
  data.frame(node = c("a", "b", "c", "b", "c", "d"),
             edge = rep(c("e1", "e2"), each = 3)),
  "node", "edge"
)
hg_embed(h, dimensions = 2)
#>   node        pi dim1          dim2
#> 1    a 0.1666667 -0.5  8.660254e-01
#> 2    b 0.3333333 -1.0 -1.145529e-15
#> 3    c 0.3333333 -1.0 -8.735814e-16
#> 4    d 0.1666667 -0.5 -8.660254e-01
```
