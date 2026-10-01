# Tables of a persistent homology result

Tables of a persistent homology result

## Usage

``` r
# S3 method for class 'persistent_homology'
hg_get(
  x,
  what = c("persistence", "betti"),
  ...,
  dimension = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `persistent_homology` object from
  [`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md).

- what:

  `"persistence"` (default) for the persistence diagram, or `"betti"`
  for the Betti curves across the filtration.

- ...:

  Additional arguments (ignored).

- dimension:

  Integer or `NULL`. Keep only this homological dimension.

- sort_by:

  `NULL` (construction order, default) or `"persistence"` -
  longest-lived features first. Only for `what = "persistence"`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "persistence"`, one row per feature:
`dimension`, `birth`, `death`, `persistence`. For `what = "betti"`, one
row per (threshold, dimension): `threshold`, `dimension`, `betti`.

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
hg_get(ph, what = "betti", dimension = 0)
#>      threshold dimension betti
#> 1  0.742723701         0     5
#> 2  0.661024094         0     3
#> 3  0.579324487         0     1
#> 4  0.497624880         0     1
#> 5  0.415925273         0     1
#> 6  0.334225665         0     1
#> 7  0.252526058         0     1
#> 8  0.170826451         0     1
#> 9  0.089126844         0     1
#> 10 0.007427237         0     1
```
