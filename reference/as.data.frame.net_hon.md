# Coerce a net_hon to a tidy table

Coerce a net_hon to a tidy table

## Usage

``` r
# S3 method for class 'net_hon'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("rules", "nodes"),
  order_min = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hon` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"rules"` (default) for the extracted higher-order rules, or `"nodes"`
  for the node table of the higher-order topology.

- order_min:

  Integer or `NULL`. Keep only rules whose source context is at least
  this order (e.g. `2` for the genuinely higher-order rules). Ignored
  when `what = "nodes"`.

- sort_by:

  `NULL` (construction order, default), `"count"` or `"probability"` -
  sort largest first, ties broken by `from`/`to` so the order is
  deterministic. Ignored when `what = "nodes"`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "rules"`, one row per rule edge with columns
`path`, `from`, `to`, `count`, `probability`, `from_order`, `to_order`.
For `what = "nodes"`, one row per higher-order node with columns `id`,
`label`, `name`.
