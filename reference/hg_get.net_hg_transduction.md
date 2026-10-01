# Tables of a hypergraph transduction fit

Tables of a hypergraph transduction fit

## Usage

``` r
# S3 method for class 'net_hg_transduction'
hg_get(x, what = c("predictions", "scores"), ..., top = NULL)
```

## Arguments

- x:

  A `net_hg_transduction` object.

- what:

  Character. `"predictions"` (default) for the one-row-per-node table,
  `"scores"` for the tidy long score table (one row per node x class:
  `node`, `class`, `score`).

- ...:

  Additional arguments (ignored).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame as selected by `what`.
