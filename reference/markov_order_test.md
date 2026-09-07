# Test the Markov order of a sequential process

Principled test of whether a categorical sequence is best described as a
\\k\\-th order Markov chain. At each order \\k = 1, \ldots,\\
`max_order`, the function computes the classical likelihood-ratio
statistic (\\G^2\\) for the conditional independence \\s \perp x \mid
w\\, where \\w\\ is the \\(k-1)\\-gram context, \\x\\ is the extra
(k-th-back) state added at order \\k\\, and \\s\\ is the next state.
Under \\H_0\\ (order-\\(k-1)\\ is correct), \\s\\ is independent of
\\x\\ given \\w\\.

The null distribution is obtained by an **exact within-\\w\\ permutation
test**: for each context \\w\\ the successor labels are exchangeable
under \\H_0\\, so shuffling \\s\\ within each \\w\\-group yields an
exact reference distribution for \\G^2\\. No plug-in MLE bias and no
refitting per replicate. An asymptotic \\\chi^2\\ p-value is reported
alongside for reference.

The optimal order is the smallest \\k\\ that is **not** significantly
better than \\k - 1\\ at level `alpha`: we keep raising the order while
the test rejects, and stop at the first non-rejection.

## Usage

``` r
markov_order_test(
  data,
  max_order = 3L,
  n_perm = 500L,
  alpha = 0.05,
  parallel = FALSE,
  n_cores = 2L,
  seed = NULL
)
```

## Arguments

- data:

  A data.frame (wide format, one sequence per row) or list of character
  vectors (one per trajectory). NAs are treated as end of sequence.

- max_order:

  Integer. Highest Markov order to test. Default 3.

- n_perm:

  Integer. Number of within-\\w\\ permutations per order. Default 500.

- alpha:

  Numeric. Significance level for order selection. Default 0.05.

- parallel:

  Logical. Use
  [`parallel::mclapply`](https://rdrr.io/r/parallel/mclapply.html) for
  permutations. Default `FALSE` (set `TRUE` only on Unix-like systems).

- n_cores:

  Integer. Cores for parallel execution. Default 2.

- seed:

  Optional integer seed for reproducibility.

## Value

An object of class `net_markov_order` with elements:

- optimal_order:

  Integer. Selected order via sequential permutation test.

- bic_order:

  Integer. Order minimising BIC (reported by `print`/`plot`; not used in
  the permutation selection).

- aic_order:

  Integer. Order minimising AIC (reported by `print`/`plot`; not used in
  the permutation selection).

- test_table:

  Tidy data.frame, one row per order tested with columns `order`,
  `loglik`, `AIC`, `BIC`, `df`, `g2`, `p_permutation`, `p_asymptotic`,
  `significant`. `AIC`/`BIC` are the information-criterion values used
  for the `ic` plot panel and to derive `aic_order`/`bic_order`; the
  order-0 row has `NA` for the test columns (`df`, `g2`,
  `p_permutation`, `p_asymptotic`, `significant`).

- permutation_null:

  List of numeric vectors (length `max_order`), one empirical null
  \\G^2\\ distribution per order.

- logliks:

  Named numeric vector of log-likelihoods per order (for AIC / BIC panel
  only, not used in the test).

- layer_dofs:

  Named integer vector of model degrees of freedom per order (free
  parameters added at each layer), used to compute the AIC / BIC
  columns.

- transition_matrices:

  List of fitted transition matrices.

- states:

  Character vector of observed state labels.

- n_sequences, n_observations:

  Data summary.

- n_perm, alpha, max_order:

  Call settings.

## Examples

``` r
# \donttest{
# First-order Markov data: test should select order 1
set.seed(1)
states <- letters[1:4]
tm <- matrix(runif(16), 4, 4, dimnames = list(states, states))
tm <- tm / rowSums(tm)
seqs <- lapply(1:30, function(.) {
  s <- character(50); s[1] <- sample(states, 1)
  for (i in 2:50) s[i] <- sample(states, 1, prob = tm[s[i - 1], ])
  s
})
res <- markov_order_test(seqs, max_order = 3, n_perm = 300, seed = 1)
as.data.frame(res)
#>   order    loglik      AIC      BIC  df        g2 p_permutation p_asymptotic
#> 1     0 -2011.025 4028.050 4043.990  NA        NA            NA           NA
#> 2     1 -1805.112 3640.224 3719.922   9 411.78019   0.003322259 4.193315e-83
#> 3     2 -1789.863 3703.726 4033.146  36  30.35171   0.803986711 7.337874e-01
#> 4     3 -1701.468 3808.936 4887.520 136 173.56575   0.099667774 1.633824e-02
#>   significant
#> 1          NA
#> 2        TRUE
#> 3       FALSE
#> 4       FALSE
summary(res)
plot(res)

# }
```
