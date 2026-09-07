# Sub-hypergraph by hyperedges, nodes or hyperedge attributes

Keeps part of a hypergraph and returns it as a `net_hypergraph` with the
same incidence weights, edge metadata and multiplicities. Three
selectors combine by intersection. `edges` names the hyperedges to keep.
`nodes` keeps the hyperedges whose members all lie in the set, the
induced sub-hypergraph of HypergraphX. `where` keeps the hyperedges
whose attributes take given values, so `where = c(citing = "153-001")`
on a citation-block hypergraph is the hypergraph of one citing decision
(Coupette et al. 2024, Figure 3).

## Usage

``` r
hg_subset(hg, edges = NULL, nodes = NULL, where = NULL, drop_isolated = TRUE)

hypergraph_subset(
  hg,
  edges = NULL,
  nodes = NULL,
  where = NULL,
  drop_isolated = TRUE
)
```

## Arguments

- hg:

  A `net_hypergraph` (from
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  [`hypergraph_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md),
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  ...).

- edges:

  Hyperedge names to keep: a character vector, or a data.frame with an
  `edge` column such as the table
  [`hg_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md)
  or
  [`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
  returns, so a ranking can be passed straight through.

- nodes:

  Character vector of node names. Only hyperedges contained in this set
  are kept, and the node set of the result is exactly `nodes`.

- where:

  Hyperedge attribute values to keep: a named vector or list, one
  element per attribute column of the edge metadata (the columns
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md)
  keeps because they are constant within a hyperedge), each holding the
  value or values to keep.

- drop_isolated:

  Drop nodes that belong to no retained hyperedge? Default `TRUE`.
  Ignored when `nodes` is given.

## Value

A `net_hypergraph` whose incidence matrix is the selected sub-matrix of
the input, sparse if the input is sparse. Edge metadata (`edge_data`)
and duplicate multiplicities (`edge_multiplicity`) are subset alongside.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "b", "c", "d", "d", "e"),
  event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3"),
  kind = c("x", "x", "x", "x", "x", "x", "y", "y")
)
hg <- group_hypergraph(dat, actor = "member", group = "event")
hg_subset(hg, edges = c("e1", "e2"))
#> Hypergraph: 4 nodes, 2 hyperedges
#> Size distribution:
#>   size_3   : 2
#> Source: group membership (member = member, group = event)
hg_subset(hg, nodes = c("b", "c", "d", "e"))
#> Hypergraph: 4 nodes, 2 hyperedges
#> Size distribution:
#>   size_2   : 1
#>   size_3   : 1
#> Source: group membership (member = member, group = event)
hg_subset(hg, where = c(kind = "y"))
#> Hypergraph: 2 nodes, 1 hyperedges
#> Size distribution:
#>   size_2   : 1
#> Source: group membership (member = member, group = event)
```
