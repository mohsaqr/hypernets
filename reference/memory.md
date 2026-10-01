# Per-context path dependence

For every context of `order - 1` states, compares the observed
next-state distribution with the first-order prediction from the most
recent state alone: the Kullback-Leibler divergence and entropy drop say
how much the longer history adds, and `flips` marks the contexts where
it changes the most likely next state.

## Usage

``` r
memory(
  data,
  order = 2L,
  min_count = 5L,
  base = 2,
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
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md):
  a long event table (with `action`), a wide data.frame, a list of
  vectors, or a model object carrying its sequences.

- order:

  Integer \>= 2. Length of the conditioning context plus one: `2`
  compares two-step memory with one-step. Default `2`.

- min_count:

  Integer. Contexts observed fewer times are dropped. Default `5`.

- base:

  Logarithm base of the entropies and divergences. Default `2` (bits).

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide, list or model input.

## Value

A `net_path_dependence` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
one row per context (`what = "contexts"`), sorted by divergence, largest
first.

## References

Kullback, S., & Leibler, R. A. (1951). On information and sufficiency.
*The Annals of Mathematical Statistics*, 22(1), 79–86.
[doi:10.1214/aoms/1177729694](https://doi.org/10.1214/aoms/1177729694)

Cover, T. M., & Thomas, J. A. (2006). *Elements of Information Theory*
(2nd ed.). Wiley.

## See also

[`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
[`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md)

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
pd <- memory(seqs, order = 2, min_count = 1)
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
```
