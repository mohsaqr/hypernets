# Coerce a net_persistent_homology to a tidy table

Coerce a net_persistent_homology to a tidy table

## Usage

``` r
# S3 method for class 'net_persistent_homology'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("persistence", "betti"),
  dimension = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_persistent_homology` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"persistence"` (default) for the persistence diagram, or `"betti"`
  for the Betti curves across the filtration.

- dimension:

  Integer or `NULL`. Keep only this homological dimension.

- sort_by:

  `NULL` (construction order, default) or `"persistence"` -
  longest-lived features first. Only for `what = "persistence"`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "persistence"`, one row per feature:
`dimension`, `birth`, `death`, `persistence`. For `what = "betti"`, one
row per (threshold, dimension): `threshold`, `dimension`, `betti`.
