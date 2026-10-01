# Tables of a hypergraph measures fit

Tables of a hypergraph measures fit

## Usage

``` r
# S3 method for class 'net_hg_measures'
hg_get(
  x,
  what = c("nodes", "edges", "global"),
  ...,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hg_measures` object.

- what:

  `"nodes"` (default) for the per-node measures, `"edges"` for the
  per-hyperedge measures, or `"global"` for the whole-hypergraph
  summary.

- ...:

  Additional arguments (ignored).

- sort_by:

  `NULL` (node/edge order, default), or a column name to sort by,
  largest first. Ignored when `what = "global"`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "nodes"`, one row per node: `node`,
`hyperdegree`, `node_strength`, `max_edge_size`. For `what = "edges"`,
one row per hyperedge: `hyperedge`, `size`. For `what = "global"`, one
row per statistic: `measure`, `value`.
