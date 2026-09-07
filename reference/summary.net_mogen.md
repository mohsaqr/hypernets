# Summary Method for net_mogen

Summary Method for net_mogen

## Usage

``` r
# S3 method for class 'net_mogen'
summary(object, ...)
```

## Arguments

- object:

  A `net_mogen` object.

- ...:

  Additional arguments (ignored).

## Value

A per-order model-selection data.frame with columns `order`,
`layer_dof`, `cum_dof`, `loglik`, `aic`, `bic`, `best`
(`"AIC"`/`"BIC"`/`"AIC+BIC"` marker) and `selected` (`"<--"` on the
chosen order), returned visibly; the summary text is printed as a side
effect.

Returned **invisibly**: `summary(x)` prints the summary and nothing
else; assign the result to keep the table.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
mg <- build_mogen(seqs, max_order = 2)
summary(mg)
#> Multi-Order Generative Model (MOGen) Summary
#> 
#>   States: A, B, C, D
#>   Paths: 3 | Observations: 12
#>   Best by AIC: order 1  |  Best by BIC: order 1
#>   Selected:    order 1 (by aic)
#> 

# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
mog <- build_mogen(seqs, max_order = 2L)
summary(mog)
#> Multi-Order Generative Model (MOGen) Summary
#> 
#>   States: A, B, C
#>   Paths: 5 | Observations: 15
#>   Best by AIC: order 1  |  Best by BIC: order 1
#>   Selected:    order 1 (by aic)
#> 
# }
```
