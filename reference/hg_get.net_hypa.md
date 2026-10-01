# Tables of a path-anomaly (HYPA) fit

Tables of a path-anomaly (HYPA) fit

## Usage

``` r
# S3 method for class 'hypernets_hypa'
hg_get(x, ...)

# S3 method for class 'net_hypa'
hg_get(
  x,
  what = c("anomalies", "scores", "over", "under", "pathways"),
  ...,
  type = NULL,
  order_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hypa` object from
  [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md) on
  sequences.

- ...:

  Unused.

- what:

  `"anomalies"` (default) for the anomalous paths selected by `type`,
  `"scores"` for every scored path, `"over"` or `"under"` for the
  significantly over- or under-represented paths with all their columns,
  or `"pathways"` for the anomalous paths written as pathways
  (`"a b -> c"`: the context states, then the next state), the form
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)
  turns into simplices.

- type:

  For `"anomalies"` and `"pathways"`: `"all"`, `"over"` or `"under"`.
  Default: the `type` given to
  [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md),
  else `"all"`.

- order_by:

  `"sig"` (smallest p-value of the path's own direction first),
  `"ratio"` (largest deviation first: highest ratio for over-represented
  paths, lowest for under-represented ones, largest absolute log ratio
  when both directions are listed), `"freq"` (most observed first) or
  `"path"` (alphabetical). Default: the `order_by` given to
  [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md),
  else `"sig"`. For `"scores"`, `"sig"` sorts by `p_value`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `order_by`. Default `NULL` returns every row.

## Value

A data.frame. For `"anomalies"`, one row per anomalous path: `path`,
`count`, `expected`, `ratio`, `p_adj` (the adjusted p-value of the
path's own direction), `direction`, and `order` when several orders were
scored. For `"scores"`, `"over"` and `"under"`, one row per path with
every column of the fit: `path`, `from`, `to`, `count`, `expected`,
`ratio`, the p-values (`p_value`, `p_under`, `p_over`, `p_adj_under`,
`p_adj_over`), the `direction` (`"over"`, `"under"` or `"normal"`) and
the path `order`. For `"pathways"`, one row per anomalous path with the
column `pathway`.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
fit <- hypa(seqs, order = 2, min_count = 1)
hg_get(fit, type = "over", order_by = "ratio", top = 5)
#> [1] path      count     expected  ratio     p_adj     direction
#> <0 rows> (or 0-length row.names)
```
