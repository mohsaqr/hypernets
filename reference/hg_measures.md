# Structural measures for a hypergraph

Computes a comprehensive structural-statistics suite for a
[net_hg](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md):
node-level, hyperedge-level, and global measures. All measures are
derived in a few BLAS calls on the incidence matrix.

## Usage

``` r
# S3 method for class 'net_hg_measures'
print(x, ...)

hg_measures(
  hg,
  what = c("nodes", "edges", "overlap", "summary", "distribution", "components"),
  measure = c("hyperdegree", "strength", "n_neighbors", "size")
)
```

## Arguments

- x:

  A `hg_measures` object.

- ...:

  Additional arguments (ignored).

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hg`).

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

The input `x` invisibly.

A base `data.frame`, one row per node, edge, edge pair, measure,
distinct value, or component according to `what`.

## Details

`what =` selects which table is returned, as a tidy data.frame.

All measures are computed via standard matrix operations on the binary
incidence \\B = (b\_{ij})\\ where \\b\_{ij} = 1\\ iff node \\i\\ is in
hyperedge \\j\\:

- `hyperdegree = rowSums(B)`, `edge_sizes = colSums(B)`

- `co_degree = tcrossprod(B)` (with zero diagonal)

- `edge_pairwise_overlap = crossprod(B)` (with zero diagonal)

- `overlap_coefficient[i, j] = overlap[i, j] / min(edge_sizes[i], edge_sizes[j])`

- `jaccard[i, j] = overlap[i, j] / (edge_sizes[i] + edge_sizes[j] - overlap[i, j])`

Empty hypergraph (`n_hyperedges == 0`) returns trivial zeros and empty
matrices.

## References

Lee, G., Bu, F., Eliassi-Rad, T., & Shin, K. (2025). A survey on
hypergraph mining: patterns, tools, and generators. *ACM Computing
Surveys*, 57(8), 203.
[doi:10.1145/3719002](https://doi.org/10.1145/3719002)

Do, M. T., Yoon, S., Hooi, B., & Shin, K. (2020). Structural patterns
and generative models of real-world hypergraphs. arXiv:2006.07060.

## See also

[`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
[`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md).

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
