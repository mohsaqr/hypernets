# Coerce a net_hypergraph_cluster to a data.frame

Coerce a net_hypergraph_cluster to a data.frame

## Usage

``` r
# S3 method for class 'net_hypergraph_cluster'
as.data.frame(x, ..., top = NULL)
```

## Arguments

- x:

  A `net_hypergraph_cluster` object.

- ...:

  Additional arguments (ignored).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

The tidy assignment table: one row per node, columns `node`, `cluster`,
`pi` (stationary probability of the node under the Laplacian's random
walk) and the spectral-embedding or NMF-factor coordinates `dim1..dimk`.
