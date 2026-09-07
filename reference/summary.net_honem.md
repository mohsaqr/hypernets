# Summary Method for net_honem

Summary Method for net_honem

## Usage

``` r
# S3 method for class 'net_honem'
summary(object, ...)
```

## Arguments

- object:

  A `net_honem` object.

- ...:

  Additional arguments (ignored).

## Value

A data.frame with one row per node: column `node` (node label) followed
by `dim1`, `dim2`, ..., `dim`*d* embedding coordinates, returned
visibly; the summary text is printed as a side effect.

Returned **invisibly**: `summary(x)` prints the summary and nothing
else; assign the result to keep the table.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon_2 <- build_hon(seqs, max_order = 2)
hem <- build_honem(hon_2, dim = 2)
summary(hem)
#> HONEM Summary
#> 
#>   Nodes: 4 | Dimensions: 2
#>   Variance explained: 78.2%
#>   Top singular values: 1.023, 0.6
#>   Embedding range: [-0.573, 0.430]

# \donttest{
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon <- build_hon(seqs, max_order = 3)
he <- build_honem(hon, dim = 2)
summary(he)
#> HONEM Summary
#> 
#>   Nodes: 4 | Dimensions: 2
#>   Variance explained: 78.2%
#>   Top singular values: 1.023, 0.6
#>   Embedding range: [-0.573, 0.430]
# }
```
