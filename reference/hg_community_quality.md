# Quality of a projected-hypergraph partition

Reports graph-partition diagnostics on the association projection:
unweighted coverage, weighted coverage, performance, weighted
modularity, and conductance. Conductance is the maximum (worst)
community conductance \\cut(S, \bar S) / min(vol(S), vol(\bar S))\\;
lower is better.

## Usage

``` r
hg_community_quality(
  hg,
  partition,
  method = c("association", "citation"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL
)

hypergraph_community_quality(
  hg,
  partition,
  method = c("association", "citation"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL
)
```

## Arguments

- hg:

  A static `net_hypergraph`.

- partition:

  An
  [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  result, a tidy node/label table, or a named label vector.

- method:

  The projection the partition is scored on: the `"association"` graph
  (default) or the undirected `"citation"` graph.

- duplicate_edges, self_association, edge_source:

  Projection controls passed to
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md).
  Together these reproduce the paper's binary/multi and self-association
  representations.

## Value

A one-row data frame.
