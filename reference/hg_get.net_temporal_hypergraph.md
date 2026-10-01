# Tidy tables of a temporal hypergraph

Tidy tables of a temporal hypergraph

## Usage

``` r
# S3 method for class 'net_temporal_hypergraph'
hg_get(x, what = c("memberships", "edges", "nodes"), ...)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- what:

  `"memberships"` (default) for one row per node-in-hyperedge spell
  (`member`, `edge`, `start`, `end`, `weight`), `"edges"` for one row
  per hyperedge with its clock and attributes, `"nodes"` for the node
  universe with entry times. Times are on the hypergraph's clock.

- ...:

  Ignored.

## Value

A base `data.frame`.
