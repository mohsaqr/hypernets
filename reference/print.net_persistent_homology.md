# Print persistent homology results

Print persistent homology results

## Usage

``` r
# S3 method for class 'net_persistent_homology'
print(x, ...)
```

## Arguments

- x:

  A `net_persistent_homology` object.

- ...:

  Additional arguments (unused).

## Value

The input object, invisibly.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 10)
print(ph)
#> Persistent Homology
#>   10 filtration steps [0.6000 → 0.0060]
#>   Features: b0: 2 (1 persistent) 
#>   Longest-lived:
#>     b0: 0.6000 → 0.0000 (life: 0.6000)
#>     b0: 0.6000 → 0.5000 (life: 0.1000)
```
