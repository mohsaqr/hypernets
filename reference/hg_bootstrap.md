# Bootstrap inference for higher-order network rules

Nonparametric bootstrap over sequences for the rules of a higher-order
network (see
[`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md)):
sequences are resampled with replacement, and for every rule edge of the
observed network the replicate distribution yields a percentile
confidence interval for its conditional probability and a *support* -
the fraction of replicates in which the rule's context is extracted as a
rule at all. Support close to 1 for a higher-order context means the
order elevation is stable under resampling, not an artifact of the
particular sample; first-order contexts have support 1 by construction.

## Usage

``` r
hg_bootstrap(
  data,
  n_boot = 500L,
  level = 0.95,
  max_order = 5L,
  min_freq = 1L,
  collapse_repeats = FALSE,
  action = NULL,
  actor = NULL,
  time = NULL,
  parallel = FALSE,
  n_cores = 2L,
  seed = NULL,
  session = NULL,
  time_threshold = 900,
  timezone = "UTC"
)
```

## Arguments

- data:

  Sequences in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md):
  a long event table (with `action`), a wide data.frame (one sequence
  per row), a list of vectors, or a model object carrying its sequences.
  A group model from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md) with
  `group` is resampled group by group with its own settings, and the
  result is a `net_hon_boot_group` (one bootstrap per group;
  [`summary()`](https://rdrr.io/r/base/summary.html) and
  [`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
  stack them with a `group` column).

- n_boot:

  Integer \>= 2. Bootstrap replicates. Default `500`.

- level:

  Confidence level in (0, 1). Default `0.95` (percentile interval).

- max_order, min_freq, collapse_repeats:

  As in [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md).

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide or list input.

- parallel:

  Logical. Use
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) for
  the replicates (not on Windows). Results are identical to the serial
  run.

- n_cores:

  Integer. Cores when `parallel = TRUE`.

- seed:

  Optional integer seed.

## Value

An object of class `net_hon_boot`: a list with `edges` (the tidy
inference table, one row per rule edge of the observed network: `from`,
`to`, `order`, `count`, `probability`, `ci_lower`, `ci_upper`,
`support`, `n_boot_used`), `n_boot`, `level`, `max_order`, `min_freq`,
`n_trajectories`, and `seed`. Has `print`, `summary` and `plot` methods;
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
returns the inference table (optionally filtered with `min_support =` or
restricted with `order_min =`).

## Details

Per-sequence observation counts are computed once; every replicate is a
weighted aggregation of those counts followed by counts-based rule
extraction (never a re-scan of the data). All randomness is drawn before
any parallel work, so `parallel = TRUE` reproduces the serial result
under the same `seed`.

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances* 2(5),
e1600028.
[doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028)

Efron, B., & Tibshirani, R. J. (1993). *An Introduction to the
Bootstrap*. Chapman & Hall.

## See also

[`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
[`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md),
[`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md)

## Examples

``` r
hg_seqs <- list(
  c("a", "b", "c", "a", "b", "c"),
  c("x", "b", "d", "x", "b", "d"),
  c("a", "b", "c", "a", "b", "c"),
  c("x", "b", "d", "x", "b", "d")
)
bs <- hg_bootstrap(hg_seqs, n_boot = 50, max_order = 2, seed = 1)
bs
#> HON bootstrap: 8 rule edges (2 higher-order) from 4 sequences
#>   50 replicates, 95% percentile CIs
#>   Higher-order rule support: min 0.54, median 0.58, max 0.62
#>    from to order count probability ci_lower ci_upper support n_boot_used
#>       a  b     1     4         1.0  1.00000  1.00000    0.94          47
#>       b  c     1     4         0.5  0.00000  0.94375    1.00          50
#>       b  d     1     4         0.5  0.05625  1.00000    1.00          50
#>       c  a     1     2         1.0  1.00000  1.00000    0.94          47
#>       d  x     1     2         1.0  1.00000  1.00000    0.96          48
#>       x  b     1     4         1.0  1.00000  1.00000    0.96          48
#>  a -> b  c     2     4         1.0  1.00000  1.00000    0.54          47
#>  x -> b  d     2     4         1.0  1.00000  1.00000    0.62          48
hg_get(bs, sort_by = "support", top = 6)
#>   from to order count probability ci_lower ci_upper support n_boot_used
#> 1    b  c     1     4         0.5  0.00000  0.94375    1.00          50
#> 2    b  d     1     4         0.5  0.05625  1.00000    1.00          50
#> 3    d  x     1     2         1.0  1.00000  1.00000    0.96          48
#> 4    x  b     1     4         1.0  1.00000  1.00000    0.96          48
#> 5    a  b     1     4         1.0  1.00000  1.00000    0.94          47
#> 6    c  a     1     2         1.0  1.00000  1.00000    0.94          47
```
