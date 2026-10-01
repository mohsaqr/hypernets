# Markov stability of the states of a transition network

Per-state persistence (self-transition probability), stationary
probability, return time, sojourn time and mean first-passage times to
and from the other states of a first-order Markov chain.

## Usage

``` r
hg_markov_stability(
  data,
  normalize = TRUE,
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

  A row-stochastic transition matrix, a network object (`netobject`,
  `tna`, `cograph_network`, `netobject_group`) whose weights are the
  transition matrix, or sequences in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md),
  from which the first-order transition probabilities are estimated.

- normalize:

  Logical. Normalise the rows to sum to one. Default `TRUE`.

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide, list or model input.

## Value

A `net_markov_stability` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
one row per state (`what = "states"`, the default), the mean
first-passage times (`"passage_time"`) or the stationary distribution
(`"stationary"`).

## Conditions

The quantities are defined only for an irreducible chain, one in which
every state can reach every other. A chain that is not irreducible (a
transient or absorbing state, or several closed classes, as a pruned
higher-order network can have) raises `hypernets_not_ergodic`, a
`hypernets_bad_input`, naming those states.

## References

Kemeny, J. G., & Snell, J. L. (1976). *Finite Markov Chains*. Springer.

## See also

[`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md)

## Examples

``` r
P <- matrix(c(0.5, 0.3, 0.2,
              0.2, 0.5, 0.3,
              0.1, 0.4, 0.5), 3, 3, byrow = TRUE,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
ms <- hg_markov_stability(P)
hg_get(ms, sort_by = "stationary")
#>   state persistence stationary return_time sojourn_time avg_time_to_others
#> 1     B         0.5     0.4182        2.39            2               4.92
#> 2     C         0.5     0.3455        2.89            2               4.77
#> 3     A         0.5     0.2364        4.23            2               3.63
#>   avg_time_from_others
#> 1                 2.83
#> 2                 3.95
#> 3                 6.54
hg_get(ms, what = "passage_time", from = "A")
#>   from to    steps
#> 1    A  A 4.230769
#> 2    A  B 3.043478
#> 3    A  C 4.210526
```
