# Plot Persistence Landscape

Plot Persistence Landscape

## Usage

``` r
# S3 method for class 'net_persistence_landscape'
plot(x, ...)
```

## Arguments

- x:

  A `net_persistence_landscape` object.

- ...:

  Ignored.

## Value

A ggplot.

## Examples

``` r
# \donttest{
mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
rownames(mat) <- colnames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 5)
pl <- persistence_landscape(ph, k_max = 3, dimension = 0)
plot(pl)

# }
```
