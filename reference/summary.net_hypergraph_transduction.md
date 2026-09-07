# Summary method for net_hypergraph_transduction

Summary method for net_hypergraph_transduction

## Usage

``` r
# S3 method for class 'net_hypergraph_transduction'
summary(object, ...)
```

## Arguments

- object:

  A `net_hypergraph_transduction` object.

- ...:

  Additional arguments (ignored).

## Value

A data.frame, one row per class: `class`, `n_labeled`, `n_predicted`,
`mean_margin` (mean winning margin among the nodes predicted into the
class). Returned **invisibly**: `summary(x)` prints the summary and
nothing else; assign the result to keep the table.
