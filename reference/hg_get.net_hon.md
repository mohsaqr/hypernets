# Tables of a higher-order network

Tables of a higher-order network

## Usage

``` r
# S3 method for class 'net_hon'
hg_get(
  x,
  what = c("rules", "nodes", "pathways"),
  ...,
  order_min = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hon` object from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md).

- what:

  `"rules"` (default) for the extracted higher-order rules, `"nodes"`
  for the node table of the higher-order topology, or `"pathways"` for
  the genuinely higher-order rules written as pathways (`"a b -> c"`:
  the context states, then the next state), most frequent first – the
  form
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)
  turns into simplices.

- ...:

  For `what = "pathways"`: `min_count`, `min_prob` (drop rules below
  this count or probability) and `order` (keep one context order).
  Otherwise unused.

- order_min:

  Integer or `NULL`. Keep only rules whose source context is at least
  this order (e.g. `2` for the genuinely higher-order rules). Only for
  `what = "rules"`.

- sort_by:

  `NULL` (construction order, default), `"count"` or `"probability"` -
  sort largest first, ties broken by `from`/`to` so the order is
  deterministic. Only for `what = "rules"`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "rules"`, one row per rule edge with columns
`path`, `from`, `to`, `count`, `probability`, `from_order`, `to_order`.
For `what = "nodes"`, one row per higher-order node with columns `id`
and `node`. For `what = "pathways"`, one row per higher-order rule with
the column `pathway`.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
net <- hon(seqs, max_order = 2)
hg_get(net, order_min = 2)
#>          path   from to count probability from_order to_order
#> 1 a -> b -> c a -> b  c     4           1          2        1
#> 2 x -> b -> d x -> b  d     4           1          2        1
hg_get(net, what = "nodes")
#>   id   node
#> 1  1      a
#> 2  2 a -> b
#> 3  3      b
#> 4  4      c
#> 5  5      d
#> 6  6      x
#> 7  7 x -> b
hg_get(net, what = "pathways")
#>    pathway
#> 1 a b -> c
#> 2 x b -> d
```
