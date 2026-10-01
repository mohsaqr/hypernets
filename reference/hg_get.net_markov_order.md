# Tables of a Markov order test

Tables of a Markov order test

## Usage

``` r
# S3 method for class 'net_markov_order'
hg_get(x, what = c("orders", "null"), ..., top = NULL)
```

## Arguments

- x:

  A `net_markov_order` object from
  [`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md).

- what:

  `"orders"` (default) for the per-order test table, or `"null"` for the
  permutation null distribution behind it.

- ...:

  Additional arguments (ignored).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "orders"`, one row per candidate order:
`order`, `log_likelihood`, `aic`, `bic`, `df`, `g2` (the
likelihood-ratio statistic against the order below), `p_value` (its
permutation p-value), `p_asymptotic` (its chi-squared p-value) and
`significant`. For `what = "null"`, one row per permutation replicate:
`order`, `replicate`, `g2` (the replicate's likelihood-ratio statistic
under the null).

## Examples

``` r
seqs <- list(c("a", "b", "a", "c", "a", "b"), c("b", "a", "c", "a", "b", "a"))
fit <- markov_order(seqs, max_order = 2L, n_perm = 20L, seed = 1L)
hg_get(fit)
#>   order log_likelihood      aic      bic df        g2    p_value p_asymptotic
#> 1     0     -12.136851 28.27370 29.24352 NA        NA         NA           NA
#> 2     1      -5.156818 16.31364 17.76836  4 13.862944 0.04761905  0.007745578
#> 3     2      -2.302585 10.60517 12.05989  1  5.545177 0.52380952  0.018531678
#>   significant
#> 1          NA
#> 2        TRUE
#> 3       FALSE
hg_get(fit, what = "null", top = 5)
#>   order replicate       g2
#> 1     1         1 3.452185
#> 2     1         2 1.726092
#> 3     1         3 3.452185
#> 4     1         4 7.271270
#> 5     1         5 4.498681
```
