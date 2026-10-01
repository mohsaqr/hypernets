# Permutation comparison of the rules of two groups

Tests whether two groups of sequences differ in their higher-order rule
probabilities. The rule set is extracted from the pooled data of the two
groups; for every pooled rule the statistic is the absolute difference
of the two groups' conditional probabilities, and its null distribution
comes from permuting the group labels over sequences. Per-rule p-values
are Benjamini-Hochberg adjusted; a global test aggregates the rule
differences weighted by pooled counts.

## Usage

``` r
hg_compare(
  x,
  groups = NULL,
  n_perm = 1000L,
  alpha = 0.05,
  parallel = FALSE,
  n_cores = 2L,
  seed = NULL
)
```

## Arguments

- x:

  A group model from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md) with
  `group`. The estimation settings (`max_order`, `min_freq`,
  `collapse_repeats`) and the data of each group are taken from it.

- groups:

  The two groups to compare, as a character vector of length 2. May be
  omitted when the model has exactly two groups.

- n_perm:

  Integer \>= 2. Label permutations. Default `1000`.

- alpha:

  Significance level for the `significant` flag on the adjusted
  p-values. Default `0.05`.

- parallel, n_cores, seed:

  As in
  [`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md).

## Value

An object of class `net_hon_compare`: a list with `edges` (one row per
pooled rule: `from`, `to`, `order`, `count`, the two groups' counts and
probabilities (columns named after the groups), `diff` (probability
difference, first minus second), `p_value`, `p_adj` (BH), `significant`,
`n_perm_used`), `global` (`statistic`, the pooled-count-weighted mean
absolute difference, and `p_value`), `names`, `n_perm`, `alpha`,
`max_order`, `min_freq`, `n_trajectories` (per group) and `seed`. Has
`print`, `summary` and `plot` methods;
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
returns the rule table (`significant = TRUE` restricts it).

## Details

Per-sequence counts are precomputed once and every permutation is a
weighted aggregation; permutations are drawn before any parallel work,
so `parallel = TRUE` reproduces the serial result under the same `seed`.
Rules whose context is unobserved in a group under some permutation
contribute only their valid permutations (`n_perm_used`).

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances* 2(5),
e1600028.
[doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028)

Good, P. (2005). *Permutation, Parametric and Bootstrap Tests of
Hypotheses* (3rd ed.). Springer.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate. *Journal of the Royal Statistical Society B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

## See also

[`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
[`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md),
[`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md)

## Examples

``` r
set.seed(1)
first_order <- replicate(6, sample(c("a", "b", "c"), 12, replace = TRUE),
                         simplify = FALSE)
second_order <- replicate(6, rep(c("a", "b", "c", "b"), 3),
                          simplify = FALSE)
by_group <- hon(c(first_order, second_order),
                group = rep(c("first", "second"), each = 6), max_order = 2)
comparison <- hg_compare(by_group, n_perm = 99, seed = 1)
comparison
#> HON comparison: first (6 sequences) vs second (6 sequences)
#>   12 pooled rule edges, 99 permutations
#>   Global weighted |diff|: 0.4814, p = 0.01
#>   Significant edges (BH, alpha = 0.05): 11
#>    from to order count count_first count_second probability_first
#>       a  a     1     7           7            0         0.3500000
#>       a  b     1    23           5           18         0.2500000
#>       a  c     1     8           8            0         0.4000000
#>       b  a     1    18           6           12         0.2307692
#>       b  b     1    14          14            0         0.5384615
#>       b  c     1    24           6           18         0.2307692
#>       c  a     1     6           6            0         0.3000000
#>       c  b     1    26           8           18         0.4000000
#>       c  c     1     6           6            0         0.3000000
#>  c -> b  a     2    14           2           12         0.3333333
#>  probability_second       diff p_value      p_adj significant n_perm_used
#>                 0.0  0.3500000    0.02 0.02400000        TRUE          99
#>                 1.0 -0.7500000    0.01 0.01333333        TRUE          99
#>                 0.0  0.4000000    0.01 0.01333333        TRUE          99
#>                 0.4 -0.1692308    0.07 0.07000000       FALSE          99
#>                 0.0  0.5384615    0.01 0.01333333        TRUE          99
#>                 0.6 -0.3692308    0.01 0.01333333        TRUE          99
#>                 0.0  0.3000000    0.01 0.01333333        TRUE          99
#>                 1.0 -0.6000000    0.01 0.01333333        TRUE          99
#>                 0.0  0.3000000    0.01 0.01333333        TRUE          99
#>                 1.0 -0.6666667    0.01 0.01333333        TRUE          99
#> ... 2 more rows
hg_get(comparison, sort_by = "p_adj", top = 6)
#>   from to order count count_first count_second probability_first
#> 1    a  b     1    23           5           18         0.2500000
#> 2    a  c     1     8           8            0         0.4000000
#> 3    b  b     1    14          14            0         0.5384615
#> 4    b  c     1    24           6           18         0.2307692
#> 5    c  a     1     6           6            0         0.3000000
#> 6    c  b     1    26           8           18         0.4000000
#>   probability_second       diff p_value      p_adj significant n_perm_used
#> 1                1.0 -0.7500000    0.01 0.01333333        TRUE          99
#> 2                0.0  0.4000000    0.01 0.01333333        TRUE          99
#> 3                0.0  0.5384615    0.01 0.01333333        TRUE          99
#> 4                0.6 -0.3692308    0.01 0.01333333        TRUE          99
#> 5                0.0  0.3000000    0.01 0.01333333        TRUE          99
#> 6                1.0 -0.6000000    0.01 0.01333333        TRUE          99
```
