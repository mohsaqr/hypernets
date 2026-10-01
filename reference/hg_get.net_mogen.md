# Tables of a multi-order generative model

Tables of a multi-order generative model

## Usage

``` r
# S3 method for class 'net_mogen'
hg_get(
  x,
  what = c("orders", "transitions", "paths", "pathways"),
  ...,
  order = NULL,
  k = 2L,
  min_count = 1L,
  top = NULL
)
```

## Arguments

- x:

  A `net_mogen` object from
  [`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md).

- what:

  `"orders"` (default) for the per-order model-selection table,
  `"transitions"` for the fitted transition probabilities of one layer,
  `"paths"` for the frequency of every observed path of `k` states, or
  `"pathways"` for the layer's transitions written as pathways
  (`"a b -> c"`), the form
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)
  turns into simplices.

- ...:

  For `what = "pathways"`, `min_prob` (drop transitions below this
  probability). Otherwise unused.

- order:

  Integer or `NULL`. For `what = "transitions"` and `"pathways"`: which
  layer to read. Default `NULL` uses the selected order.

- k:

  Integer vector, for `what = "paths"`: the path length or lengths in
  states, each between `2` and the model's highest order plus one, e.g.
  `k = 3` or `k = 2:4`. Default `2` (single transitions).

- min_count:

  Integer, for `what = "transitions"` and `"pathways"`: drop transitions
  observed fewer times. Default `1`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "orders"`, one row per candidate order:
`order`, `log_likelihood`, `aic`, `bic`, `df`, `layer_df`, and `optimal`
(logical, `TRUE` for the selected order). For `what = "transitions"`,
one row per transition: `path`, `count`, `probability`, `from`, `to`.
For `what = "paths"`, one row per distinct path of each length in `k`
(`top` applies within each length), grouped by length and most frequent
first (ties alphabetical): `k`, `path`, `count`, `proportion` (share of
all `k`-state paths, rounded to 4 decimals). For `what = "pathways"`,
one row per transition with the column `pathway`.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
mg <- mogen(seqs, max_order = 2)
hg_get(mg)
#>   order log_likelihood      aic      bic df layer_df optimal
#> 1     0     -37.457050 82.91410 87.62631  4        4   FALSE
#> 2     1     -12.712215 35.42443 41.31470  5        1   FALSE
#> 3     2      -7.167038 24.33408 30.22434  5        0    TRUE
hg_get(mg, what = "transitions", top = 5)
#>          path count probability   from to
#> 1 a -> b -> c     4           1 a -> b  c
#> 2 x -> b -> d     4           1 x -> b  d
#> 3 c -> a -> b     2           1 c -> a  b
#> 4 b -> c -> a     2           1 b -> c  a
#> 5 b -> d -> x     2           1 b -> d  x
hg_get(mg, what = "paths", k = 3)
#>   k        path count proportion
#> 1 3 a -> b -> c     4      0.250
#> 2 3 x -> b -> d     4      0.250
#> 3 3 b -> c -> a     2      0.125
#> 4 3 b -> d -> x     2      0.125
#> 5 3 c -> a -> b     2      0.125
#> 6 3 d -> x -> b     2      0.125
```
