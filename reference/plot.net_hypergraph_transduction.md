# Plot method for net_hypergraph_transduction

Heatmap of the full node-by-class score matrix: rows are nodes (grouped
by predicted class), columns are classes, tile shading and printed
values are the spreading scores. Seed nodes (given labels) carry a black
tile border, and each node's winning class is marked with a dot, so
agreement between seeds, scores, and decisions is visible in one panel.
Rows whose winning and runner-up scores are close (small `margin`) are
the assignments to distrust.

## Usage

``` r
# S3 method for class 'net_hypergraph_transduction'
plot(x, ...)
```

## Arguments

- x:

  A `net_hypergraph_transduction` object.

- ...:

  Additional arguments (ignored).

## Value

A ggplot object, invisibly printable.
