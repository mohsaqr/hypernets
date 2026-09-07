# Coerce a net_honem to a tidy table

Coerce a net_honem to a tidy table

## Usage

``` r
# S3 method for class 'net_honem'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("embeddings", "variance"),
  top = NULL
)
```

## Arguments

- x:

  A `net_honem` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"embeddings"` (default) for the node coordinates, or `"variance"` for
  the per-dimension singular values.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "embeddings"`, one row per higher-order node:
`node` followed by one column per embedding dimension (`dim1`, `dim2`,
...). For `what = "variance"`, one row per dimension: `dim`,
`singular_value`, `proportion` (that dimension's share of the total
squared singular value).
