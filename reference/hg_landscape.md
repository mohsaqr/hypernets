# Persistence landscape of a persistence diagram

The persistence landscape: at every filtration value t, the k-th largest
of the tent functions \\\max(0, \min(t - b, d - t))\\ of the diagram's
(birth, death) pairs in one homology dimension. Essential classes, which
never die, are left out.

## Usage

``` r
hg_landscape(ph, k_max = 5L, dimension = 1L, t_grid = NULL)
```

## Arguments

- ph:

  A
  [`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md)
  result, or a data.frame with columns `dimension`, `birth`, `death`.

- k_max:

  Highest landscape level. Default `5`.

- dimension:

  Homology dimension. Default `1`.

- t_grid:

  Evaluation points. `NULL` uses 200 points spanning the pairs.

## Value

A `persistence_landscape` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
one row per level and grid point (`what = "landscape"`).

## References

Bubenik, P. (2015). Statistical topological data analysis using
persistence landscapes. *Journal of Machine Learning Research*, 16,
77–102.

## See also

[`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md)

## Examples

``` r
set.seed(1)
w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
w <- (w + t(w)) / 2
diag(w) <- 0
pl <- hg_landscape(hg_homology(w, n_steps = 10), dimension = 0)
hg_get(pl, k = 1, top = 5)
#>   k         t        value
#> 1 1 0.5838612 0.0000000000
#> 2 1 0.5846595 0.0007983042
#> 3 1 0.5854578 0.0015966085
#> 4 1 0.5862561 0.0023949127
#> 5 1 0.5870544 0.0031932170
```
