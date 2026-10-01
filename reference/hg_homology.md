# Persistent homology of a weighted network or a distance matrix

Tracks the homology classes (components, loops, voids) of a filtration:
a weighted clique filtration of a similarity matrix (`type = "clique"`,
strongest edges first) or a Vietoris-Rips filtration of a distance
matrix (`type = "vr"`). Every class is paired with the simplex that
creates it (birth) and the one that destroys it (death), by exact
boundary-matrix reduction over Z/2.

## Usage

``` r
hg_homology(x, n_steps = 20L, max_dim = 3L, type = "clique", max_scale = NULL)
```

## Arguments

- x:

  A square matrix, a network object carrying one (`netobject`, `tna`),
  points for `type = "vr"` (see
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)),
  a `simplicial_complex` built with `type = "vr"` (its stored filtration
  is used as is), or a window complex from
  `simplicial(type = "window")`. For a window complex the filtration is
  the count filtration (Petri et al. 2013): the complex at count `t`
  holds the sets of actions observed in at least `t` windows and their
  faces, so a class is born at the count where it appears and dies at a
  lower count; `birth`, `death` and the Betti-curve `threshold` are
  counts, and `n_steps` is not used (every observed count is a
  threshold).

- n_steps:

  Grid points of the reported Betti curve. Default `20`; the persistence
  diagram itself is exact.

- max_dim:

  Highest simplex dimension tracked. Default `3`.

- type:

  `"clique"` (default) or `"vr"`.

- max_scale:

  For `type = "vr"`: the longest edge admitted. `NULL` uses the largest
  distance.

## Value

A `persistent_homology` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the persistence diagram (`what = "persistence"`, the default) or the
Betti curves (`"betti"`).

## References

Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2002). Topological
persistence and simplification. *Discrete & Computational Geometry*, 28,
511–533.
[doi:10.1007/s00454-002-2885-2](https://doi.org/10.1007/s00454-002-2885-2)

Petri, G., Scolamiero, M., Donato, I., & Vaccarino, F. (2013).
Topological strata of weighted complex networks. *PLoS ONE*, 8(6),
e66506.
[doi:10.1371/journal.pone.0066506](https://doi.org/10.1371/journal.pone.0066506)

## See also

[`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md),
[`hg_bottleneck()`](https://mohsaqr.github.io/hypernets/reference/hg_bottleneck.md),
[`hg_wasserstein()`](https://mohsaqr.github.io/hypernets/reference/hg_wasserstein.md)

## Examples

``` r
set.seed(1)
w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
w <- (w + t(w)) / 2
diag(w) <- 0
ph <- hg_homology(w, n_steps = 10)
hg_get(ph, sort_by = "persistence")
#>    dimension     birth     death persistence
#> 1          0 0.7427237 0.0000000  0.74272370
#> 2          3 0.2344513 0.0000000  0.23445130
#> 3          3 0.2344513 0.0000000  0.23445130
#> 4          0 0.7427237 0.5838612  0.15886254
#> 5          3 0.1558863 0.0000000  0.15588635
#> 6          3 0.1558863 0.0000000  0.15588635
#> 7          3 0.1558863 0.0000000  0.15588635
#> 8          1 0.5170309 0.3655044  0.15152644
#> 9          0 0.7427237 0.6583996  0.08432412
#> 10         0 0.7427237 0.6902349  0.05248880
#> 11         0 0.7427237 0.7162022  0.02652146
#> 12         1 0.6441215 0.6299381  0.01418338
```
