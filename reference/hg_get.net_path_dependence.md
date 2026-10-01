# Per-context table of a path-dependence analysis

Per-context table of a path-dependence analysis

## Usage

``` r
# S3 method for class 'net_path_dependence'
hg_get(
  x,
  what = "contexts",
  ...,
  min_count = NULL,
  flips = NULL,
  sort_by = "kl",
  top = NULL
)
```

## Arguments

- x:

  A `net_path_dependence` object from
  [`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md).

- what:

  `"contexts"`, the only table.

- ...:

  Additional arguments (ignored).

- min_count:

  Integer or `NULL`. Keep only contexts observed at least this many
  times.

- flips:

  `NULL` (default, every context), `TRUE` for only the contexts whose
  most probable next state differs from the first-order prediction
  (`flips` is `TRUE`), or `FALSE` for only those where it agrees.

- sort_by:

  `"kl"` (default) sorts by the order-k vs order-1 divergence, largest
  first; `NULL` keeps the construction order.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per context: `context`, `count` (times observed),
`entropy_first_order` and `entropy_order_k` (entropy of the next state
given the last state and given the whole context), `entropy_drop` (their
difference), `kl` (Kullback-Leibler divergence of the order-k from the
first-order distribution), `top_first_order` and `top_order_k` (the most
probable next state under each), and `flips` (`TRUE` when the two
differ).

## Examples

``` r
# one sequence per row
wide <- data.frame(t1 = c("a", "x", "a", "x"), t2 = c("b", "b", "b", "b"),
                   t3 = c("c", "d", "c", "d"), t4 = c("a", "x", "a", "x"),
                   t5 = c("b", "b", "b", "b"), t6 = c("c", "d", "c", "d"))
pd <- memory(wide, order = 2, min_count = 1)
hg_get(pd)
#>   context count entropy_first_order entropy_order_k entropy_drop kl
#> 1  a -> b     4                   1               0            1  1
#> 2  x -> b     4                   1               0            1  1
#> 3  b -> c     2                   0               0            0  0
#> 4  b -> d     2                   0               0            0  0
#> 5  c -> a     2                   0               0            0  0
#> 6  d -> x     2                   0               0            0  0
#>   top_first_order top_order_k flips
#> 1               c           c FALSE
#> 2               c           d  TRUE
#> 3               a           a FALSE
#> 4               x           x FALSE
#> 5               b           b FALSE
#> 6               b           b FALSE
hg_get(pd, flips = TRUE)
#>   context count entropy_first_order entropy_order_k entropy_drop kl
#> 1  x -> b     4                   1               0            1  1
#>   top_first_order top_order_k flips
#> 1               c           d  TRUE
```
