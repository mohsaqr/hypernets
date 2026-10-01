# Test the Markov order of sequences

At each order k the likelihood-ratio statistic G2 tests whether the k-th
state back adds information about the next state given the k - 1 most
recent ones; its null distribution comes from an exact within-context
permutation. The order is raised while the test rejects.

## Usage

``` r
markov_order(
  data,
  max_order = 3L,
  n_perm = 500L,
  alpha = 0.05,
  parallel = FALSE,
  n_cores = 2L,
  seed = NULL,
  action = NULL,
  actor = NULL,
  time = NULL,
  session = NULL,
  time_threshold = 900,
  timezone = "UTC"
)
```

## Arguments

- data:

  Sequences in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md);
  a `netobject_group` gives one test per group.

- max_order:

  Integer. Highest order tested. Default `3`.

- n_perm:

  Integer. Permutations per order. Default `500`.

- alpha:

  Significance level of the sequential selection. Default `0.05`.

- parallel:

  Logical. Run the permutations with
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (not on Windows). Default `FALSE`.

- n_cores:

  Integer. Cores when `parallel = TRUE`. Default `2`.

- seed:

  Integer or `NULL`. Seed for the permutations.

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide, list or model input.

## Value

A `net_markov_order` object (a `net_markov_order_group` for a
`netobject_group`). Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the per-order test table (`what = "orders"`, the default) or the
permutation null draws (`"null"`).

## References

Anderson, T. W., & Goodman, L. A. (1957). Statistical inference about
Markov chains. *The Annals of Mathematical Statistics*, 28(1), 89–110.
[doi:10.1214/aoms/1177707039](https://doi.org/10.1214/aoms/1177707039)

Good, P. (2005). *Permutation, Parametric and Bootstrap Tests of
Hypotheses* (3rd ed.). Springer.

## See also

[`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md),
[`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md)

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
```
