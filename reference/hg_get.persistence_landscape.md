# Table of a persistence landscape

Table of a persistence landscape

## Usage

``` r
# S3 method for class 'persistence_landscape'
hg_get(x, what = "landscape", ..., k = NULL, top = NULL)
```

## Arguments

- x:

  A `persistence_landscape` object from
  [`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md).

- what:

  `"landscape"`, the only table.

- ...:

  Additional arguments (ignored).

- k:

  Integer or `NULL`. Keep only this landscape level.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per (level, grid point): `k` (landscape level),
`t` (filtration value), `value` (landscape height).

## Examples

``` r
set.seed(1)
w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
w <- (w + t(w)) / 2
diag(w) <- 0
pl <- hg_landscape(hg_homology(w, n_steps = 10),
                            dimension = 0)
hg_get(pl, k = 1, top = 5)
#>   k         t        value
#> 1 1 0.5838612 0.0000000000
#> 2 1 0.5846595 0.0007983042
#> 3 1 0.5854578 0.0015966085
#> 4 1 0.5862561 0.0023949127
#> 5 1 0.5870544 0.0031932170
```
