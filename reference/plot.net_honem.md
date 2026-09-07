# Plot Method for net_honem

Plot Method for net_honem

## Usage

``` r
# S3 method for class 'net_honem'
plot(x, dims = c(1L, 2L), ...)
```

## Arguments

- x:

  A `net_honem` object.

- dims:

  Integer vector of length 2. Dimensions to plot (default: `c(1, 2)`).

- ...:

  Additional arguments passed to
  [`plot`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The input object, invisibly.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon_2 <- build_hon(seqs, max_order = 2)
hem <- build_honem(hon_2, dim = 2)
plot(hem)


# \donttest{
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon <- build_hon(seqs, max_order = 3)
he <- build_honem(hon, dim = 2)
plot(he)
# }
```
