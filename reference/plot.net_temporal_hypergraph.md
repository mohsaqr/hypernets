# Plot a temporal-hypergraph snapshot through cograph

Plot a temporal-hypergraph snapshot through cograph

## Usage

``` r
# S3 method for class 'net_temporal_hypergraph'
plot(
  x,
  at = NULL,
  mode = c("active", "cumulative"),
  method = c("association", "clique"),
  ...
)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- at:

  Snapshot time; defaults to the end of observation.

- mode:

  Snapshot mode passed to
  [`hypergraph_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md).

- method:

  Projection weighting passed to
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md).

- ...:

  Additional arguments passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## Value

The cograph plot object, invisibly when rendered interactively.
