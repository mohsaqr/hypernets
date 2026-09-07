# Hypergraph PageRank with edge-dependent vertex weights

Ranks the vertices of a weighted hypergraph by the stationary importance
of the Chitra & Raphael (2019) random walk: from a vertex, pick an
incident hyperedge with probability proportional to its edge weight,
then land on a member vertex with probability proportional to that
vertex's *edge-dependent* incidence weight (for a text hypergraph, e.g.,
the tf-idf of the word in the document). Damping and an optional
personalization vector turn the stationary distribution into PageRank
(Page et al. 1999): with probability `damping` the walker follows the
hypergraph walk, otherwise it teleports.

## Usage

``` r
hg_pagerank(
  hg,
  damping = 0.85,
  personalized = NULL,
  edge_weights = NULL,
  sort_by = NULL,
  n = Inf,
  max_iter = 1000L,
  tol = 1e-12
)

hypergraph_pagerank(
  hg,
  damping = 0.85,
  personalized = NULL,
  edge_weights = NULL,
  sort_by = NULL,
  n = Inf,
  max_iter = 1000L,
  tol = 1e-12
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md),
  or any hypernets `net_hypergraph` (connected when `damping = 1`).

- damping:

  Probability of following the hypergraph walk (default `0.85`);
  `1 - damping` is the teleport probability. Must be in `(0, 1]`;
  `damping = 1` gives the pure stationary distribution and requires a
  connected hypergraph to converge.

- personalized:

  Optional restart preference: a character vector of vertex names
  (uniform teleport over exactly those vertices), or a named
  non-negative vector of teleport weights over (a subset of) the vertex
  names; unnamed vertices get teleport probability 0. `NULL` (default)
  teleports uniformly.

- edge_weights:

  Positive hyperedge weights (one per hyperedge), or `NULL` (default)
  for the Hayashi et al. heuristic used by the Laplacian engines: the
  population standard deviation of each edge's non-zero vertex weights
  plus one, which reduces to unit weights on a binary incidence.

- sort_by:

  `NULL` (default, vertex order) or `"pagerank"` to sort descending
  (ties broken by vertex name).

- n:

  Return only the first `n` rows after sorting (default all).

- max_iter, tol:

  Power-iteration cap and L1 convergence tolerance.

## Value

A base `data.frame`, one row per vertex, with columns `node` and
`pagerank` (non-negative, summing to 1).

## Details

Chitra & Raphael prove that with edge-*independent* vertex weights the
walk collapses to a random walk on a weighted clique-expansion graph —
so the hypergraph brings genuinely new information exactly when the
incidence weights differ across the hyperedges a vertex belongs to
(tested against the closed-form graph stationary distribution). With
`damping = 1` and default `edge_weights`, the result equals the
stationary distribution of the Hayashi et al. (2020) EDVW walk used by
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
(tested at `1e-12` against
[`hypergraph_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_laplacian.md)),
and with uniform teleportation it equals
`hypergraph_centrality(type = "pagerank")` (tested).

## Conditions

Raises `hypernets_bad_input` for broken contracts and warns with
`hypernets_no_converge` (returning the last iterate) when `max_iter` is
reached before `tol`.

## References

Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
edge-dependent vertex weights. *ICML 2019*.

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Page, L., Brin, S., Motwani, R., & Winograd, T. (1999). The PageRank
citation ranking: Bringing order to the web. Stanford InfoLab.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"),
   weight = "tfidf")
hg_pagerank(hg, sort_by = "pagerank")
#>        node  pagerank
#> 1 cooking_2 0.2819453
#> 2   space_1 0.2478410
#> 3   space_2 0.2459361
#> 4 cooking_1 0.2242776
hg_pagerank(hg, personalized = c(cooking_1 = 1), sort_by = "pagerank")
#>        node   pagerank
#> 1 cooking_1 0.65606771
#> 2 cooking_2 0.25577437
#> 3   space_2 0.05757252
#> 4   space_1 0.03058540
```
