# Tidy tables of a temporal hypergraph

Tidy tables of a temporal hypergraph

## Usage

``` r
# S3 method for class 'net_temporal_hypergraph'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("memberships", "edges", "nodes"),
  ...
)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- row.names, optional:

  Ignored; present for compatibility.

- what:

  `"memberships"` (default) for one row per node-in-hyperedge spell
  (`member`, `edge`, `start`, `end`, `weight`), `"edges"` for one row
  per hyperedge with its clock and attributes, `"nodes"` for the node
  universe with entry times. Times are on the hypergraph's clock.

- ...:

  Ignored.

## Value

A base `data.frame`.
