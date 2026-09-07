# Coerce a net_path_dependence to its tidy per-context table

Coerce a net_path_dependence to its tidy per-context table

## Usage

``` r
# S3 method for class 'net_path_dependence'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  min_count = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_path_dependence` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- min_count:

  Integer or `NULL`. Keep only contexts observed at least this many
  times.

- sort_by:

  `NULL` (construction order, default) or `"KL"` - sort by the order-k
  vs order-1 divergence, largest first.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per context, carrying the context label, its
observed count, and the Kullback-Leibler divergence between its order-k
and order-1 continuation distributions.
