# Print Persistence Landscape

Print Persistence Landscape

## Usage

``` r
# S3 method for class 'net_persistence_landscape'
print(x, ...)
```

## Arguments

- x:

  A `net_persistence_landscape` object.

- ...:

  Ignored.

## Value

The input, invisibly.

## Examples

``` r
mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
rownames(mat) <- colnames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 5)
pl <- persistence_landscape(ph, k_max = 3, dimension = 0)
print(pl)
#> Persistence Landscape (dimension 0, k_max = 3)
#>   Evaluated at 200 grid points in [0.5000, 0.6000]
#>   Sup norms by k: k1=0.0497  k2=0.0000  k3=0.0000 
```
