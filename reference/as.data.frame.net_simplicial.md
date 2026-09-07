# Coerce a net_simplicial to a tidy table

Coerce a net_simplicial to a tidy table

## Usage

``` r
# S3 method for class 'net_simplicial'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("simplices", "f_vector"),
  dim = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_simplicial` object.

- row.names:

  Ignored (S3 consistency).

- optional:

  Ignored (S3 consistency).

- ...:

  Additional arguments (ignored).

- what:

  `"simplices"` (default) for one row per simplex, or `"f_vector"` for
  the face counts by dimension.

- dim:

  Integer or `NULL`. Only for `what = "simplices"`: keep only simplices
  of this dimension (`0` vertices, `1` edges, `2` triangles).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "simplices"`, one row per simplex: `id`,
`dim`, `size`, `members` (the node names, comma separated). For
`what = "f_vector"`, one row per dimension: `dim`, `count`.
