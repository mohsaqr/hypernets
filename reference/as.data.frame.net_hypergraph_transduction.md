# Coerce a net_hypergraph_transduction to a data.frame

Coerce a net_hypergraph_transduction to a data.frame

## Usage

``` r
# S3 method for class 'net_hypergraph_transduction'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("predictions", "scores"),
  top = NULL
)
```

## Arguments

- x:

  A `net_hypergraph_transduction` object.

- row.names:

  Ignored (present for S3 consistency with the generic).

- optional:

  Ignored (present for S3 consistency with the generic).

- ...:

  Additional arguments (ignored).

- what:

  Character. `"predictions"` (default) for the one-row-per-node table,
  `"scores"` for the tidy long score table (one row per node x class:
  `node`, `class`, `score`).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame as selected by `what`.
