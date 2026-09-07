# Coerce a net_q_analysis to a tidy table

Coerce a net_q_analysis to a tidy table

## Usage

``` r
# S3 method for class 'net_q_analysis'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("q_levels", "nodes"),
  top = NULL
)
```

## Arguments

- x:

  A `net_q_analysis` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"q_levels"` (default) for the number of q-connected components at
  each level, or `"nodes"` for each node's highest q-level.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "q_levels"`, one row per level: `q`,
`components`. For `what = "nodes"`, one row per node: `node`, `max_q`.
