# The s-line graph of a hypergraph

Turns the hyperedges into vertices: two hyperedges are adjacent when
they share at least `s` vertices, and the edge weight is how many they
share. `s = 1` is the ordinary line graph (any overlap connects);
raising `s` keeps only substantial overlaps, which is how walk-based
measures on hyperedges are parametrised (Aksoy et al. 2020).

## Usage

``` r
hg_line_graph(hg, s = 1, what = c("edges", "matrix"))

hypergraph_line_graph(hg, s = 1, what = c("edges", "matrix"))
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md),
  or any hypernets `net_hypergraph`.

- s:

  Minimum number of shared vertices for two hyperedges to be adjacent. A
  single integer of at least 1; `1` (default) is the ordinary line
  graph.

- what:

  `"edges"` (default) for the tidy edge list, or `"matrix"` for the
  symmetric overlap matrix to hand to a graph engine.

## Value

With `what = "edges"`, a base data.frame with one row per unordered
hyperedge pair sharing at least `s` vertices, sorted by `from` then
`to`, with columns `from`, `to` (hyperedge names) and `weight` (the
number of shared vertices). A zero-row data.frame with those columns
when no pair reaches `s`. With `what = "matrix"`, the symmetric
`n_hyperedges` x `n_hyperedges` overlap matrix, zero on the diagonal and
wherever the overlap is below `s`, sparse if `hg`'s incidence is sparse.

## Details

This is the projection of the dual: `hg_line_graph(hg, s = 1)` returns
the same graph as `hg_project(dual_hypergraph(hg), weighted = FALSE)`,
which the tests assert.

## References

Aksoy, S. G., Joslyn, C., Ortiz Marrero, C., Praggastis, B., & Purvine,
E. (2020). Hypernetwork science via high-order hypergraph walks. *EPJ
Data Science*, 9(1), 16.
[doi:10.1140/epjds/s13688-020-00231-0](https://doi.org/10.1140/epjds/s13688-020-00231-0)

## See also

[`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md)
for the projection onto vertices,
[`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md)
for the role swap itself.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
                        c = "stars and salt"))
hg_line_graph(hg)
#>   from    to weight
#> 1  and  salt      2
#> 2  and  soup      2
#> 3  and stars      2
#> 4 salt  soup      1
#> 5 salt stars      1
#> 6 soup stars      1
hg_line_graph(hg, s = 2)
#>   from    to weight
#> 1  and  salt      2
#> 2  and  soup      2
#> 3  and stars      2
```
