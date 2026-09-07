# Coerce a net_hypa to a tidy table

Coerce a net_hypa to a tidy table

## Usage

``` r
# S3 method for class 'net_hypa'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("scores", "over", "under"),
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hypa` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"scores"` (default) for every scored path, `"over"` for the
  significantly over-represented paths only, or `"under"` for the
  significantly under-represented ones.

- sort_by:

  `NULL` (construction order, default), `"p"` (sorts by `p_value`,
  smallest first), `"ratio"` (largest observed/expected deviation first)
  or `"observed"` (heaviest paths first). `p_value` is the two-sided
  value; the direction-specific `p_adjusted_under` / `p_adjusted_over`
  columns are left for the caller to read, since a single ordering
  cannot serve both directions.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per path: `path`, `from`, `to`, `observed`,
`expected`, `ratio`, the p-values (`p_value`, `p_under`, `p_over`,
`p_adjusted_under`, `p_adjusted_over`), the `anomaly` label and the path
`order`.
