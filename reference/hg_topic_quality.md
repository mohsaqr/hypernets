# Topic coherence and exclusivity

Scores each cluster's owned vocabulary on two axes. Coherence is the
mean normalised pointwise mutual information (NPMI) between pairs of the
cluster's top words, computed on document co-occurrence within the
cluster: 1 when two words always occur together, 0 at independence,
negative when they avoid each other. Exclusivity is the mean share of
the same words, the fraction of each word's corpus count that falls in
the cluster: 1 when the cluster owns its vocabulary outright. The top
words are the ones
[`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
ranks by share with the given support floor, so the two scores describe
the vocabulary a reader is shown.

## Usage

``` r
hg_topic_quality(hg, clusters, n = 10L, min_docs = 1L)

# S3 method for class 'honets_topic_quality'
plot(x, ...)

hypergraph_topic_quality(hg, clusters, n = 10L, min_docs = 1L)
```

## Arguments

- hg:

  The document hypergraph the clustering was computed on (bag-of-words,
  documents as nodes).

- clusters:

  The tidy table returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  (columns `node`, `cluster`), or a named vector of cluster labels.

- n:

  Top words per cluster to score (default `10`).

- min_docs:

  Support floor passed to
  [`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  (default `1`).

- x:

  A table returned by the verb.

- ...:

  Unused; for S3 consistency.

## Value

A base `data.frame` of class `honets_topic_quality`, one row per topic
in natural order: `topic`, `size`, `n_words` (top words scored),
`coherence`, `exclusivity`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
exclusivity against coherence with one labelled point per topic. Raises
`honets_bad_input` for unknown node names.

## References

Bouma, G. (2009). Normalized (pointwise) mutual information in
collocation extraction. *Proceedings of GSCL*, 31–40.

Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
for structural topic models. *Journal of Statistical Software*, 91(2).

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_topic_quality(hg, topics, n = 3)
#>       topic size n_words  coherence exclusivity
#> 1 Cluster 1    2       3 -0.3333333           1
#> 2 Cluster 2    2       3  0.3333333           1
```
