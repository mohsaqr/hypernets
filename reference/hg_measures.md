# Structural measures of a hypergraph, as tidy tables

Delegates to
[`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md)
and returns the requested slice as a tidy data.frame.

## Usage

``` r
hg_measures(
  hg,
  what = c("nodes", "edges", "overlap", "summary", "distribution", "components"),
  measure = c("hyperdegree", "strength", "n_neighbors", "size")
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hypergraph`).

- what:

  Which table: `"nodes"` (default; one row per node with `hyperdegree`,
  `strength`, `max_edge_size` and `n_neighbors`, the distinct nodes it
  shares a hyperedge with), `"edges"` (one row per hyperedge with its
  `size`), `"overlap"` (one row per hyperedge pair with `overlap`,
  `overlap_coefficient`, `jaccard`), `"summary"` (one row per scalar
  measure), `"distribution"` (the empirical distribution of `measure`,
  as a `hypernets_distribution` table whose
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
  CCDF), or `"components"` (one row per connected component through
  shared hyperedges, with its `n_nodes`, `n_edges`, `share` of nodes and
  `diameter`).

- measure:

  For `what = "distribution"`: `"hyperdegree"` (default), `"strength"`,
  `"n_neighbors"` (node measures) or `"size"` (hyperedge cardinality).

## Value

A base `data.frame`, one row per node, edge, edge pair, measure,
distinct value, or component according to `what`.

## Examples

``` r
hg <- text_hypergraph(c(
  a = "salt and soup and onions",
  b = "soup and salt",
  c = "stars and sky"
))
hg_measures(hg)
#>   node hyperdegree strength max_edge_size n_neighbors
#> 1    a           4        8             3           2
#> 2    b           3        7             3           2
#> 3    c           3        5             3           2
hg_measures(hg, what = "summary")
#>                  measure     value
#> 1                n_nodes 3.0000000
#> 2           n_hyperedges 6.0000000
#> 3                density 0.5555556
#> 4          avg_edge_size 1.6666667
#> 5 pairwise_participation 1.0000000
hg_measures(hg, what = "components")
#>   component n_nodes n_edges share diameter
#> 1         1       3       6     1        1
```
