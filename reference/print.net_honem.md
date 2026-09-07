# Print Method for net_honem

Print Method for net_honem

## Usage

``` r
# S3 method for class 'net_honem'
print(x, ...)
```

## Arguments

- x:

  A `net_honem` object.

- ...:

  Additional arguments (ignored).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon_2 <- build_hon(seqs, max_order = 2)
hem <- build_honem(hon_2, dim = 2)
print(hem)
#> HONEM: Higher-Order Network Embedding
#>   Nodes:      4
#>   Dimensions: 2
#>   Max power:  10
#>   Variance explained: 78.2%

# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
hon   <- build_hon(seqs, max_order = 2L)
honem <- build_honem(hon, dim = 2L)
print(honem)
#> HONEM: Higher-Order Network Embedding
#>   Nodes:      3
#>   Dimensions: 2
#>   Max power:  10
#>   Variance explained: 82.6%
# }
```
