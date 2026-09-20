# Tidy accessor for a grouped Markov order test

The per-group test tables stacked into one data.frame, with a `group`
column naming the group each row came from. This is the supported way to
read a `net_markov_order_group`; the print method shows each group in
turn but a comparison across groups needs one table.

## Usage

``` r
# S3 method for class 'net_markov_order_group'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("orders", "null"),
  top = NULL
)
```

## Arguments

- x:

  A `net_markov_order_group`, as returned by
  [`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md)
  on a `netobject_group`.

- row.names, optional, ...:

  Passed on for S3 consistency; unused.

- what:

  `"orders"` (default) for one row per group and order, or `"null"` for
  the permutation null draws, one row per group, order and replicate.

- top:

  Keep only the first `top` rows, applied last.

## Value

A base `data.frame`. For `what = "orders"`, one row per group and order,
with `group` first and then the columns of
[`as.data.frame.net_markov_order()`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_markov_order.md).
For `what = "null"`, one row per group, order and replicate with columns
`group`, `order`, `replicate`, `g2`. Returns a zero-row frame with those
columns when the group is empty.

## Examples

``` r
seqs <- list(c("a", "b", "a", "c"), c("b", "a", "c", "a"))
fit <- markov_order_test(seqs, max_order = 2L, n_perm = 20L, seed = 1L)
as.data.frame(fit)
#>   order    loglik      AIC      BIC df       g2 p_permutation p_asymptotic
#> 1     0 -8.317766 20.63553 20.79442 NA       NA            NA           NA
#> 2     1 -3.988984 13.97797 14.21629  4 8.317766     0.3333333   0.08060755
#> 3     2 -3.178054 12.35611 12.59443  0 0.000000            NA           NA
#>   significant
#> 1          NA
#> 2       FALSE
#> 3       FALSE
```
