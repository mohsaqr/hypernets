# Bottleneck distance between persistence diagrams

The largest displacement of the best matching between two persistence
diagrams, where a point may also match the diagonal. Essential classes
are matched only to essential classes; a dimension whose essential
counts differ is `Inf` apart.

## Usage

``` r
hg_bottleneck(d1, d2, dimension = NULL, tol = .Machine$double.eps^0.5)
```

## Arguments

- d1, d2:

  [`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md)
  results, or data.frames with columns `dimension`, `birth`, `death`.

- dimension:

  Integer vector of dimensions to compare. `NULL` compares every
  dimension present in either diagram.

- tol:

  Tolerance of the binary search. Default `.Machine$double.eps^0.5`.

## Value

A named numeric vector, one distance per dimension, named `"dim_<k>"`.

## References

Cohen-Steiner, D., Edelsbrunner, H., & Harer, J. (2007). Stability of
persistence diagrams. *Discrete & Computational Geometry*, 37, 103–120.
[doi:10.1007/s00454-006-1276-5](https://doi.org/10.1007/s00454-006-1276-5)

Edelsbrunner, H., & Harer, J. (2010). *Computational Topology: An
Introduction*. American Mathematical Society.

## See also

[`hg_wasserstein()`](https://mohsaqr.github.io/hypernets/reference/hg_wasserstein.md)

## Examples

``` r
d1 <- data.frame(dimension = 0L, birth = 0, death = 2)
d2 <- data.frame(dimension = 0L, birth = 0, death = 3)
hg_bottleneck(d1, d2)
#> dim_0 
#>     1 
```
