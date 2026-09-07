# Coerce a net_hypergraph_measures to a tidy table

Coerce a net_hypergraph_measures to a tidy table

## Usage

``` r
# S3 method for class 'net_hypergraph_measures'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("nodes", "edges", "global"),
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hypergraph_measures` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"nodes"` (default) for the per-node measures, `"edges"` for the
  per-hyperedge measures, or `"global"` for the whole-hypergraph
  summary.

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
