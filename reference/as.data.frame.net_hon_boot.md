# Coerce a net_hon_boot to its tidy inference table

Coerce a net_hon_boot to its tidy inference table

## Usage

``` r
# S3 method for class 'net_hon_boot'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  min_support = NULL,
  order_min = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hon_boot` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- min_support:

  Numeric in `[0, 1]` or NULL. Keep only rule edges with at least this
  bootstrap support.

- order_min:

  Integer or NULL. Keep only rule edges of at least this order (e.g. `2`
  for the genuinely higher-order rules).

- sort_by:

  `NULL` (order/from/to, default), `"count"`, `"probability"`, or
  `"support"` - sort the table by that column, largest first (ties
  broken by from/to so the order is deterministic).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per rule edge: `from`, `to`, `order`, `count`,
`probability`, `ci_lower`, `ci_upper`, `support`, `n_boot_used`.
