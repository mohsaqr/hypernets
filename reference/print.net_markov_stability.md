# Print Method for net_markov_stability

Print Method for net_markov_stability

## Usage

``` r
# S3 method for class 'net_markov_stability'
print(x, digits = 4L, top = NULL, ...)
```

## Arguments

- x:

  A `net_markov_stability` object.

- digits:

  Integer. Decimal places for the probability columns; time columns are
  shown at two fewer. Default `4`.

- top:

  Integer or `NULL`. Show only the `top` states with the largest
  stationary probability. Default `NULL` (all states).

- ...:

  Ignored.

## Value

`x` invisibly.

## Examples

``` r
# A first-order chain, straight from a transition matrix
P <- matrix(c(0.7, 0.2, 0.1,
              0.3, 0.5, 0.2,
              0.2, 0.3, 0.5),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
markov_stability(P)
#> Markov Stability -- 3 states
#> 
#>  state persistence stationary_prob return_time sojourn_time avg_time_to_others
#>      A         0.7          0.4634        2.16         3.33               6.20
#>      B         0.5          0.3171        3.15         2.00               5.18
#>      C         0.5          0.2195        4.56         2.00               4.03
#>  avg_time_from_others
#>                  3.95
#>                  4.23
#>                  7.22
#> 
#> as.data.frame(x) for the full table, as.data.frame(x, what = "passage_time") for first passage times.

# The higher-order walk: states that carry their memory
hon <- build_hon(split(human_long$code, human_long$session_id),
                 max_order = 2L, min_freq = 50L)
ms <- markov_stability(hon)
as.data.frame(ms, sort_by = "stationary_prob", top = 5L)
#>                  state persistence stationary_prob return_time sojourn_time
#> 1   Specify -> Command   0.0000000      0.11167480    8.954571     1.000000
#> 2   Command -> Specify   0.0000000      0.10021296    9.978749     1.000000
#> 3   Specify -> Specify   0.4842520      0.09193677   10.877041     1.938931
#> 4   Command -> Command   0.5784314      0.09161941   10.914718     2.372093
#> 5 Specify -> Interrupt   0.0000000      0.08410419   11.890014     1.000000
#>   avg_time_to_others avg_time_from_others
#> 1           81.10439             8.358104
#> 2           80.59760            10.552172
#> 3           80.48824            19.269256
#> 4           80.43635            23.299006
#> 5           79.03732            11.067830
as.data.frame(ms, what = "passage_time", sort_by = "steps", top = 5L)
#>                 from                 to    steps
#> 1  Specify -> Refine Specify -> Correct 361.1430
#> 2  Refine -> Specify Specify -> Correct 360.1430
#> 3             Refine Specify -> Correct 358.9391
#> 4 Specify -> Specify Specify -> Correct 357.5931
#> 5 Specify -> Inquire Specify -> Correct 357.1221
```
