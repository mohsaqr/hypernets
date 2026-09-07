# Coerce a net_persistence_landscape to its tidy table

Coerce a net_persistence_landscape to its tidy table

## Usage

``` r
# S3 method for class 'net_persistence_landscape'
as.data.frame(x, row.names = NULL, optional = FALSE, ..., k = NULL, top = NULL)
```

## Arguments

- x:

  A `net_persistence_landscape` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- k:

  Integer or `NULL`. Keep only this landscape level.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per (level, grid point): `k` (landscape level),
`t` (filtration value), `value` (landscape height).
