# Tables of a higher-order network embedding

Tables of a higher-order network embedding

## Usage

``` r
# S3 method for class 'net_honem'
hg_get(x, what = c("embeddings", "variance"), ..., top = NULL)
```

## Arguments

- x:

  A `net_honem` object from
  [`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md).

- what:

  `"embeddings"` (default) for the node coordinates, or `"variance"` for
  the per-dimension singular values.

- ...:

  Unused.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "embeddings"`, one row per higher-order node:
`node` followed by one column per embedding dimension (`dim1`, `dim2`,
...). For `what = "variance"`, one row per dimension: `dimension`,
`singular_value`, `proportion` (that dimension's share of the total
squared singular value).

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
emb <- honem(hon(seqs, max_order = 2), dim = 2)
hg_get(emb)
#>     node      dim1          dim2
#> 1      a -0.366535 -4.082483e-01
#> 2 a -> b -0.428539 -4.082483e-01
#> 3      b -0.428539  6.106227e-16
#> 4      c -0.366535 -4.082483e-01
#> 5      d -0.366535  4.082483e-01
#> 6      x -0.366535  4.082483e-01
#> 7 x -> b -0.428539  4.082483e-01
hg_get(emb, what = "variance")
#>   dimension singular_value proportion
#> 1         1       1.088329  0.5422208
#> 2         2       1.000000  0.4577792
```
