# Topic sizes

Counts the documents of each cluster of a clustered hypergraph, as a
number and as a share of the clustered documents, and, when a weight per
document is given (a repeat count, a citation count), the weighted size
and share as well.

## Usage

``` r
hg_topic_sizes(hg, clusters, weights = NULL)

# S3 method for class 'hypernets_topic_sizes'
plot(x, ...)

hypergraph_topic_sizes(hg, clusters, weights = NULL)
```

## Arguments

- hg:

  The document hypergraph the clustering was computed on.

- clusters:

  The tidy table returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  (columns `node`, `cluster`), or a named vector of cluster labels.

- weights:

  `NULL` (default), or a numeric vector named by document giving each
  document's weight.

- x:

  A table returned by the verb.

- ...:

  Unused; for S3 consistency.

## Value

A base `data.frame` of class `hypernets_topic_sizes`, one row per topic
in natural order: `topic`, `n`, `share`, and with `weights` also
`weighted_n` and `weighted_share`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
shares as horizontal bars, weighted beside unweighted when both exist.
Raises `hypernets_bad_input` for unknown node names or weights that do
not name every clustered document.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_topic_sizes(hg, topics)
#>       topic n share
#> 1 Cluster 1 2   0.5
#> 2 Cluster 2 2   0.5
hg_topic_sizes(hg, topics, weights = c(cooking_1 = 3, cooking_2 = 1,
                                       space_1 = 1, space_2 = 1))
#>       topic n share weighted_n weighted_share
#> 1 Cluster 1 2   0.5          4      0.6666667
#> 2 Cluster 2 2   0.5          2      0.3333333
```
