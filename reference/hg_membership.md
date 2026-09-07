# Soft topic membership from the spectral embedding

Turns a hard partition into a membership of every document in every
topic. The documents are placed in the spectral embedding of the
hypergraph Laplacian that
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
cuts (the same `type` and `edge_weights`), each topic's centre is the
mean position of its documents, and the membership of a document in a
topic is the fuzzy c-means weight with fuzziness 2: the inverse squared
distance to that centre, normalised over the topics. A document at a
centre has membership 1 there; a document halfway between two centres
has 0.5 in each. The values sum to one over the topics, so their level
depends on the number of topics `k`: the uniform value is `1 / k`, and
with many topics even a clearly assigned document has a modest
membership. Read them as ratios between topics (a document with 0.15 and
0.11 in two topics is 1.4 times closer to the first), not as
probabilities of belonging.

## Usage

``` r
hg_membership(
  hg,
  clusters,
  type = c("zhou", "random_walk"),
  edge_weights = NULL
)

# S3 method for class 'honets_membership'
plot(x, ...)

hypergraph_membership(
  hg,
  clusters,
  type = c("zhou", "random_walk"),
  edge_weights = NULL
)
```

## Arguments

- hg:

  The document hypergraph the clustering was computed on.

- clusters:

  The tidy table returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  (columns `node`, `cluster`), or a named vector of cluster labels.

- type, edge_weights:

  Passed to
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  to reproduce the embedding the partition was cut in (defaults as
  there).

- x:

  A table returned by the verb.

- ...:

  Unused; for S3 consistency.

## Value

A base `data.frame` of class `honets_membership`, one row per document
and topic: `node`, `cluster` (the hard label), `topic`, `membership`
(rows of one document sum to one).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws, per
topic, the distribution of its documents' membership in it: a topic
whose documents sit near 1 is compact, one whose documents spread
towards 0.5 overlaps its neighbours. Raises `honets_bad_input` for
unknown node names.

## References

Bezdek, J. C. (1981). *Pattern Recognition with Fuzzy Objective Function
Algorithms*. Plenum.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_membership(hg, topics)
#>        node   cluster     topic  membership
#> 1 cooking_1 Cluster 1 Cluster 1 0.981830203
#> 2 cooking_1 Cluster 1 Cluster 2 0.018169797
#> 3 cooking_2 Cluster 1 Cluster 1 0.971963064
#> 4 cooking_2 Cluster 1 Cluster 2 0.028036936
#> 5   space_1 Cluster 2 Cluster 1 0.005819465
#> 6   space_1 Cluster 2 Cluster 2 0.994180535
#> 7   space_2 Cluster 2 Cluster 1 0.007346161
#> 8   space_2 Cluster 2 Cluster 2 0.992653839
```
