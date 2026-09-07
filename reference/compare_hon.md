# Two-sample permutation comparison of higher-order network rules

Tests whether two cohorts of sequences differ in their higher-order rule
probabilities. The rule set is extracted from the pooled data (see
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md));
for every pooled rule edge the statistic is the absolute difference of
the two cohorts' conditional probabilities, and its null distribution
comes from permuting cohort labels over sequences. Per-edge p-values are
Benjamini-Hochberg adjusted; a global test aggregates the edge
differences weighted by pooled counts.

## Usage

``` r
compare_hon(
  x,
  y,
  n_perm = 1000L,
  alpha = 0.05,
  max_order = 5L,
  min_freq = 1L,
  collapse_repeats = FALSE,
  action = NULL,
  actor = NULL,
  time = NULL,
  names = c("x", "y"),
  parallel = FALSE,
  n_cores = 2L,
  seed = NULL
)
```

## Arguments

- x, y:

  The two cohorts of sequence data, each in any input format accepted by
  [`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md)
  (wide data.frame, list, `tna`/`netobject`, or long data.frame with
  `action`/`actor`/`time`).

- n_perm:

  Integer \>= 2. Label permutations. Default `1000`.

- alpha:

  Significance level for the `significant` flag on the adjusted
  p-values. Default `0.05`.

- max_order, min_freq, collapse_repeats:

  As in
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md).

- action, actor, time:

  Long-format column names applied to both `x` and `y`; `NULL` for
  wide/list input.

- names:

  Character vector of length 2 naming the cohorts in the output (default
  `c("x", "y")`).

- parallel, n_cores, seed:

  As in
  [`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md).

## Value

An object of class `net_hon_compare`: a list with `edges` (one row per
pooled rule edge: `from`, `to`, `order`, `count`, `count_x`/`count_y`
and `prob_x`/`prob_y` (columns named after `names`), `diff` (prob
difference, first minus second), `p_value`, `p_adj` (BH), `significant`,
`n_perm_used`), `global` (list: `statistic` - the pooled-count-weighted
mean absolute difference - and `p_value`), `names`, `n_perm`, `alpha`,
`max_order`, `min_freq`, `n_trajectories` (per cohort), and `seed`. Has
`print`, `summary`, `plot`, and `as.data.frame` methods;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the edge table (`significant = TRUE` restricts it).

## Details

As in
[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md),
per-sequence counts are precomputed once and every permutation is a
weighted aggregation; permutations are drawn before any parallel work,
so `parallel = TRUE` reproduces the serial result under the same `seed`.
Edges whose context is unobserved in a cohort under some permutation
contribute only their valid permutations (`n_perm_used`).

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances* 2(5),
e1600028.
[doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028)

Good, P. (2005). *Permutation, Parametric and Bootstrap Tests of
Hypotheses* (3rd ed.). Springer.

## See also

[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md),
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md),
[`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md)

## Examples

``` r
first_order  <- replicate(6, sample(c("a", "b", "c"), 12, replace = TRUE),
                          simplify = FALSE)
second_order <- replicate(6, rep(c("a", "b", "c", "b"), 3),
                          simplify = FALSE)
cmp <- compare_hon(first_order, second_order, n_perm = 99,
                   max_order = 2, seed = 1)
cmp
#> HON comparison: x (6 sequences) vs y (6 sequences)
#>   15 pooled rule edges, 99 permutations
#>   Global weighted |diff|: 0.3711, p = 0.01
#>   Significant edges (BH, alpha = 0.05): 7
#>   Tidy table: as.data.frame(x); significant only: as.data.frame(x, significant = TRUE)
head(as.data.frame(cmp))
#>   from to order count count_x count_y    prob_x prob_y        diff p_value
#> 1    a  a     1     8       8       0 0.3809524    0.0  0.38095238    0.02
#> 2    a  b     1    24       6      18 0.2857143    1.0 -0.71428571    0.01
#> 3    a  c     1     7       7       0 0.3333333    0.0  0.33333333    0.01
#> 4    b  a     1    20       8      12 0.3333333    0.4 -0.06666667    0.56
#> 5    b  b     1     9       9       0 0.3750000    0.0  0.37500000    0.01
#> 6    b  c     1    25       7      18 0.2916667    0.6 -0.30833333    0.01
#>        p_adj significant n_perm_used
#> 1 0.02571429        TRUE          99
#> 2 0.01500000        TRUE          99
#> 3 0.01500000        TRUE          99
#> 4 0.56000000       FALSE          99
#> 5 0.01500000        TRUE          99
#> 6 0.01500000        TRUE          99
```
