# Graph and hypergraph representations of the same data

Reports, side by side, the four ways Coupette et al. (2024, Table 2)
model one dataset: the hypergraph as it is (`mh`, repeated member sets
counted separately), the binary hypergraph of distinct member sets
(`bh`), and the graph the hypergraph reduces to, with (`mg`) or without
(`bg`) edge multiplicities. The graph is the clique expansion for
collaboration data (every pair of members of a tribunal is an edge) or
the citation graph for data with hyperedge sources (the citing decision
is joined to every decision in the block). The degree statistics show
how much the choice of representation changes what a node looks like.

## Usage

``` r
hg_representations(hg, graph = c("clique", "citation"), edge_source = NULL)

hypergraph_representations(
  hg,
  graph = c("clique", "citation"),
  edge_source = NULL
)
```

## Arguments

- hg:

  A static `net_hypergraph`, typically a snapshot of a
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- graph:

  `"clique"` (default) or `"citation"`, the graph projection compared;
  see
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md).
  `"citation"` needs hyperedge sources.

- edge_source:

  Hyperedge sources for `graph = "citation"`, as in
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md);
  omitted when the hypergraph carries them.

## Value

A base data.frame with one row per representation, in the order `bg`,
`mg`, `bh`, `mh`, and columns `representation`, `type` (`"graph"` or
`"hypergraph"`), `n_nodes`, `n_edges`, `mean_degree`, `median_degree`. A
directed citation graph adds `median_out_degree` and `median_in_degree`,
and its `mean_degree` is the mean out-degree, which equals the mean
in-degree. Hypergraph degree is the number of hyperedges containing a
node; every node of the hypergraph counts, including nodes in no
hyperedge.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
tribunals <- data.frame(
  member = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a1", "a2"),
  case = rep(c("A", "B", "C"), each = 3)
)
hg <- group_hypergraph(tribunals, "member", "case")
hg_representations(hg)
#>   representation       type n_nodes n_edges mean_degree median_degree
#> 1             bg      graph       5       6         2.4             2
#> 2             mg      graph       5       9         3.6             4
#> 3             bh hypergraph       5       2         1.2             1
#> 4             mh hypergraph       5       3         1.8             2
```
