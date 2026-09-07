# Coerce a net_markov_order to a tidy table

Coerce a net_markov_order to a tidy table

## Usage

``` r
# S3 method for class 'net_markov_order'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("orders", "null"),
  top = NULL
)
```

## Arguments

- x:

  A `net_markov_order` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"orders"` (default) for the per-order test table, or `"null"` for the
  permutation null distribution behind it.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "orders"`, one row per candidate order with
the log-likelihood, information criteria, likelihood-ratio statistic and
permutation p-value. For `what = "null"`, one row per permutation
replicate: `order`, `replicate`, `g2` (the replicate's likelihood-ratio
statistic under the null).
