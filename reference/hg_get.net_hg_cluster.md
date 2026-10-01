# Assignment table of a hypergraph clustering fit

Assignment table of a hypergraph clustering fit

## Usage

``` r
# S3 method for class 'net_hg_cluster'
hg_get(x, what = "assignments", ..., top = NULL)
```

## Arguments

- x:

  A `net_hg_cluster` object.

- what:

  `"assignments"`, the only table.

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
