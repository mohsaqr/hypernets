# Bootstrap inference for higher-order network rules

Nonparametric bootstrap over sequences for the rules of a higher-order
network (see
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)):
sequences are resampled with replacement, and for every rule edge of the
observed network the replicate distribution yields a percentile
confidence interval for its conditional probability and a *support* -
the fraction of replicates in which the rule's context is extracted as a
rule at all. Support close to 1 for a higher-order context means the
order elevation is stable under resampling, not an artifact of the
particular sample; first-order contexts have support 1 by construction.

## Usage

``` r
bootstrap_hon(
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
  seed = NULL
)
```

## Arguments

- data:

  Sequence data: wide data.frame (one sequence per row), list of
  vectors, `tna`/`netobject` model objects, or a long data.frame
  together with `action` (and optionally `actor`, `time`).

- n_boot:

  Integer \>= 2. Bootstrap replicates. Default `500`.

- level:

  Confidence level in (0, 1). Default `0.95` (percentile interval).

- max_order, min_freq, collapse_repeats:

  As in
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md).

- action, actor, time:

  Long-format column names: `action` holds the categorical state,
  `actor` groups events into sequences, `time` orders them within an
  actor (row order when `NULL`). Leave `NULL` for wide/list input.

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
`n_trajectories`, and `seed`. Has `print`, `summary`, `plot`, and
`as.data.frame` methods;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the inference table (optionally filtered with `min_support =` or
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

[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md),
[`compare_hon()`](https://mohsaqr.github.io/hypernets/reference/compare_hon.md),
[`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md)

## Examples

``` r
hg_seqs <- list(
  c("a", "b", "c", "a", "b", "c"),
  c("x", "b", "d", "x", "b", "d"),
  c("a", "b", "c", "a", "b", "c"),
  c("x", "b", "d", "x", "b", "d")
)
bs <- bootstrap_hon(hg_seqs, n_boot = 50, max_order = 2, seed = 1)
bs
#> HON bootstrap: 8 rule edges (2 higher-order) from 4 sequences
#>   50 replicates, 95% percentile CIs
#>   Higher-order rule support: min 0.54, median 0.58, max 0.62
#>   Tidy table: as.data.frame(x); higher-order only: as.data.frame(x, order_min = 2)
rules <- as.data.frame(bs)
head(rules)
#>   from to order count probability ci_lower ci_upper support n_boot_used
#> 1    a  b     1     4         1.0  1.00000  1.00000    0.94          47
#> 2    b  c     1     4         0.5  0.00000  0.94375    1.00          50
#> 3    b  d     1     4         0.5  0.05625  1.00000    1.00          50
#> 4    c  a     1     2         1.0  1.00000  1.00000    0.94          47
#> 5    d  x     1     2         1.0  1.00000  1.00000    0.96          48
#> 6    x  b     1     4         1.0  1.00000  1.00000    0.96          48
```
