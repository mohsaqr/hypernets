# Seed stability of a hypergraph clustering across resolutions

Fits
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
twice per requested `k` with different k-means seeds and reports whether
the partitions agree – the reproducibility criterion for choosing a
resolution: a partition that changes with the seed is not estimable from
the data, whatever its eigengap looks like.

## Usage

``` r
hg_stability(
  hg,
  k,
  type = c("zhou", "random_walk"),
  seeds = c(1L, 99L),
  nstart = 25L
)

hypergraph_stability(
  hg,
  k,
  type = c("zhou", "random_walk"),
  seeds = c(1L, 99L),
  nstart = 25L
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hypergraph`).

- k:

  Vector of cluster counts to test, each at least 2.

- type:

  `"zhou"` or `"random_walk"`, as in
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md).

- seeds:

  Two distinct k-means seeds (default `c(1L, 99L)`).

- nstart:

  Number of k-means starts per fit (default `25L`).

## Value

A base `data.frame`, one row per resolution, with columns `k`,
`identical_partition` (are the two partitions literally identical,
labels included) and `ari` (their adjusted Rand index; 1 means the same
partition up to label names).

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
hg_stability(hg, k = 2, type = "random_walk")
#>   k identical_partition ari
#> 1 2                TRUE   1
```
