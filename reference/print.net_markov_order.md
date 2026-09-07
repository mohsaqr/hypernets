# Print Method for net_markov_order

Print Method for net_markov_order

## Usage

``` r
# S3 method for class 'net_markov_order'
print(x, ...)
```

## Arguments

- x:

  A `net_markov_order` object.

- ...:

  Ignored.

## Value

The input object, invisibly.

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
