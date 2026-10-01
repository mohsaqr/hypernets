# Clique expansion of a hypergraph

Projects a
[net_hg](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md)
to a standard pairwise netobject (the *clique expansion* - also called
the "downgrade" of a hypergraph to a dyadic graph). Each hyperedge of
size k contributes 1 (or its weight) to every pair of its members. The
resulting edge weight `W[i, j]` equals the number of hyperedges
containing both `i` and `j` (binary incidence) or the sum of incidence
products (weighted incidence).

## Usage

``` r
hg_clique_expansion(hg, weighted = TRUE)
```

## Arguments

- hg:

  A `net_hg` object as returned by
  [`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md)
  or
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md).

- weighted:

  Logical. If `TRUE` (default), use the hypergraph's incidence values
  directly (so weighted hypergraphs from
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
  produce weighted projections). If `FALSE`, binarise the incidence
  first so `W[i, j]` is just the count of shared hyperedges.

## Value

A `netobject` (also `cograph_network`) with
`method = "hg_clique_expansion"`, undirected, with weighted symmetric
adjacency `W = incidence %*% t(incidence)` and zero diagonal.

## Details

The clique expansion is the standard "loss-y but lossless-on-pairwise"
projection: it preserves *which pairs co-occurred* and *how often* but
discards the higher-order grouping. Comparing `hg_clique_expansion(hg)`
to a directly-estimated pairwise co-occurrence network of the same data
quantifies how much information was carried by the hyperedge structure.

Computed in one BLAS call via `tcrossprod(incidence)`; runs in
`O(n_nodes^2 * n_hyperedges)` time, fast for typical sizes.

Closes the I/O cycle: event data -\>
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
-\> `hg_clique_expansion()` -\> any function that accepts a `netobject`
(centrality, bootstrap, clustering, plotting via cograph).

## Note

(experimental) Validated against `tcrossprod(incidence)` with zero
diagonal. No external R package exposes clique expansion as a primitive;
the implementation is a direct one-line restatement of the definition.

## References

Tian, Y., & Zafarani, R. (2024). Higher-order network analysis methods.
*SIGKDD Explorations* 26(1), Section 5.1.5.

## See also

[`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md).

## Examples

``` r
df <- data.frame(
  person  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
)
hg  <- group_hypergraph(df, actor = "person", group = "session")
net <- hg_clique_expansion(hg)
net
#> Network (method: hg_clique_expansion) [undirected]
#>   Weights: [1.000, 2.000]  |  mean: 1.167
#> 
#>   Weight matrix:
#>     A B C D E
#>   A 0 2 1 0 0
#>   B 2 0 1 0 0
#>   C 1 1 0 1 1
#>   D 0 0 1 0 1
#>   E 0 0 1 1 0 
```
