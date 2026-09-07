# Seed labels from a clustering, for spreading or training

Turns an unsupervised partition into the labeled examples that
[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md),
[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md)
and
[`hg_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md)
take: from each cluster, the `n` nodes with the highest stationary
probability `pi` (the cluster's most representative members under the
random walk). Ties are broken alphabetically by node name so the
selection is deterministic.

## Usage

``` r
hg_seeds(embedding, n = 5L)

hypergraph_seeds(embedding, n = 5L)
```

## Arguments

- embedding:

  The data.frame from `hg_cluster(what = "embedding")` – columns `node`,
  `cluster`, `pi`.

- n:

  Seeds per cluster (default `5L`); clusters smaller than `n` contribute
  all their nodes.

## Value

A named character vector – names are node identifiers, values their
cluster labels – ready to pass as the `labels` argument of
[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md),
[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md)
or
[`hg_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md).

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
pool <- hg_cluster(hg, k = 2, seed = 1, what = "embedding")
hg_seeds(pool, n = 1)
#>   cooking_2     space_1 
#> "Cluster 1" "Cluster 2" 
```
