# Tables of a Markov stability analysis

Tables of a Markov stability analysis

## Usage

``` r
# S3 method for class 'net_markov_stability'
hg_get(
  x,
  what = c("states", "passage_time", "stationary"),
  ...,
  sort_by = NULL,
  decreasing = TRUE,
  from = NULL,
  to = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_markov_stability` object from
  [`hg_markov_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_markov_stability.md).

- what:

  Which table: `"states"` (default, one row per state), `"passage_time"`
  (one row per ordered pair of states) or `"stationary"` (one row per
  state, distribution only).

- ...:

  Ignored.

- sort_by:

  `NULL` (construction order, the default) or a column name of the
  selected table.

- decreasing:

  Logical. Sort `sort_by` from largest to smallest (`TRUE`, the default)
  or smallest to largest (`FALSE`). Ties are broken by the row key,
  ascending in both directions: `state`, or `from` then `to` for the
  passage table.

- from, to:

  `NULL` (default, every state) or a character vector of states.
  Restrict the `what = "passage_time"` table to passages leaving the
  `from` states and arriving at the `to` states. Only that table has
  `from`/`to` columns; using either with another `what`, or naming a
  state the chain does not have, raises `hypernets_bad_input`.

- top:

  Integer or `NULL`. Keep only the first `top` rows, applied last: after
  `what`, `from`/`to` and `sort_by`.

## Value

A base `data.frame`. For `what = "states"`, one row per state with
`state`, `persistence`, `stationary`, `return_time`, `sojourn_time`,
`avg_time_to_others` and `avg_time_from_others`, as the estimator stores
them (probabilities rounded to 4 decimals, times to 2). For
`what = "passage_time"`, one row per ordered pair with `from`, `to` and
`steps`. For `what = "stationary"`, one row per state with `state`,
`stationary` and `return_time`; the passage and stationary tables are
unrounded.

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
hg_get(ms, what = "stationary")
#>   state stationary return_time
#> 1     A  0.2363636    4.230769
#> 2     B  0.4181818    2.391304
#> 3     C  0.3454545    2.894737
```
