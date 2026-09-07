# Coerce a net_mogen to a tidy table

Coerce a net_mogen to a tidy table

## Usage

``` r
# S3 method for class 'net_mogen'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("orders", "transitions"),
  order = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_mogen` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"orders"` (default) for the per-order model-selection table, or
  `"transitions"` for the fitted transition probabilities (the same
  table
  [`mogen_transitions()`](https://mohsaqr.github.io/hypernets/reference/mogen_transitions.md)
  returns).

- order:

  Integer or `NULL`. Only for `what = "transitions"`: which layer's
  transitions to return. Default `NULL` uses the selected order.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "orders"`, one row per candidate order:
`order`, `log_likelihood`, `aic`, `bic`, `dof`, `layer_dof`, and
`optimal` (logical, `TRUE` for the selected order). For
`what = "transitions"`, one row per transition.
