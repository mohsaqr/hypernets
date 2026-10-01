# Multi-order generative model (MOGen) of sequences

Fits Markov chains of every order from 0 to `max_order` as layers of one
multi-order model (higher-order De Bruijn graphs) and selects the order
that best balances fit and parsimony.

## Usage

``` r
mogen(
  data,
  max_order = 5L,
  criterion = c("aic", "bic", "lrt"),
  lrt_alpha = 0.01,
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

- max_order:

  Integer. Highest Markov order tested. Default `5`; capped at the
  longest sequence minus one, with a message.

- criterion:

  `"aic"` (default), `"bic"` or `"lrt"` (likelihood-ratio test).

- lrt_alpha:

  Significance level of the likelihood-ratio test. Default `0.01`.

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide, list or model input.

## Value

A `net_mogen` object (also a `cograph_network`). Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the order-selection table (`what = "orders"`, the default), the fitted
transitions of a layer (`"transitions"`), the path counts (`"paths"`,
with `k =`) and the pathways (`"pathways"`).

## References

Scholtes, I. (2017). When is a network a network? Multi-order graphical
model selection in pathways and temporal networks. *Proceedings of KDD
2017*, 1037–1046.
[doi:10.1145/3097983.3098145](https://doi.org/10.1145/3097983.3098145)

Gote, C., Casiraghi, G., Schweitzer, F., & Scholtes, I. (2023).
Predicting variable-length paths in networked systems using multi-order
generative models. *Applied Network Science*, 8, 68.
[doi:10.1007/s41109-023-00596-0](https://doi.org/10.1007/s41109-023-00596-0)

## See also

[`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
[`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md)

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
mg <- mogen(seqs, max_order = 2)
hg_get(mg)
#>   order log_likelihood      aic      bic df layer_df optimal
#> 1     0     -37.457050 82.91410 87.62631  4        4   FALSE
#> 2     1     -12.712215 35.42443 41.31470  5        1   FALSE
#> 3     2      -7.167038 24.33408 30.22434  5        0    TRUE
hg_get(mg, what = "paths", k = 3)
#>   k        path count proportion
#> 1 3 a -> b -> c     4      0.250
#> 2 3 x -> b -> d     4      0.250
#> 3 3 b -> c -> a     2      0.125
#> 4 3 b -> d -> x     2      0.125
#> 5 3 c -> a -> b     2      0.125
#> 6 3 d -> x -> b     2      0.125
```
