# Relations between topics: a weighted network of shared vocabulary

For a clustered document hypergraph, builds the topic-by-topic
co-occurrence network of the bibliometric kind: the strength between two
topics is the sum, over all words, of the product of the number of
documents in each topic that contain the word (full counting, the
aggregation bibnets uses for keyword co-occurrence). The raw sum can be
normalised by the same similarity measures as bibnets' `normalize()`,
with the diagonal of the co-occurrence matrix as each topic's total. The
result is an edge list with `source`, `target` and `weight`, the three
columns Gephi reads, or a `cograph_network` for plotting.

## Usage

``` r
hg_relations(
  hg,
  clusters,
  similarity = c("none", "association", "cosine", "jaccard", "inclusion", "equivalence"),
  what = c("edges", "network")
)

hypergraph_relations(
  hg,
  clusters,
  similarity = c("none", "association", "cosine", "jaccard", "inclusion", "equivalence"),
  what = c("edges", "network")
)
```

## Arguments

- hg:

  The document hypergraph the clustering was computed on.

- clusters:

  The tidy table returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  (columns `node`, `cluster`), or a named vector of cluster labels.

- similarity:

  `"none"` (default, the raw sum), `"association"`, `"cosine"`,
  `"jaccard"`, `"inclusion"` or `"equivalence"`.

- what:

  `"edges"` (default) for the edge list, `"network"` for a
  `cograph_network` built from it with
  [`cograph::as_cograph()`](https://sonsoles.me/cograph/reference/as_cograph.html),
  whose node table carries each topic's `size`.

## Value

For `what = "edges"`: a base `data.frame`, one row per pair of topics
with a positive weight, columns `source`, `target`, `weight`, pairs in
the topics' natural order. For `what = "network"`: a `cograph_network`
with one node per topic (`label`, `name`, `size`) and one undirected
weighted edge per pair. Raises `honets_bad_input` for unknown node
names.

## Details

With `D_i` the diagonal entry of topic i and `A_ij` the raw sum:
`"association"` is `A_ij / (D_i D_j)`, `"cosine"` is
`A_ij / sqrt(D_i D_j)`, `"jaccard"` is `A_ij / (D_i + D_j - A_ij)`,
`"inclusion"` is `A_ij / min(D_i, D_j)` and `"equivalence"` is
`A_ij^2 / (D_i D_j)`.

## References

van Eck, N. J., & Waltman, L. (2009). How to normalize cooccurrence
data? An analysis of some well-known similarity measures. *Journal of
the American Society for Information Science and Technology*, 60(8),
1635–1651.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_relations(hg, topics)
#>      source    target weight
#> 1 Cluster 1 Cluster 2      1
hg_relations(hg, topics, similarity = "cosine")
#>      source    target     weight
#> 1 Cluster 1 Cluster 2 0.07715167
```
