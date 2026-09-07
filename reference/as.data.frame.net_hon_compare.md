# Coerce a net_hon_compare to its tidy edge table

Coerce a net_hon_compare to its tidy edge table

## Usage

``` r
# S3 method for class 'net_hon_compare'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  significant = FALSE,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hon_compare` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- significant:

  Logical. `TRUE` restricts to edges whose adjusted p-value falls below
  the object's `alpha`. Default `FALSE` (all edges).

- sort_by:

  `NULL` (order/from/to, default), `"abs_diff"`, `"count"`, or
  `"p_adj"` - sort by absolute difference or count (largest first) or
  adjusted p-value (smallest first), ties broken by from/to.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per pooled rule edge (see
[`compare_hon()`](https://mohsaqr.github.io/hypernets/reference/compare_hon.md)
for the columns).
