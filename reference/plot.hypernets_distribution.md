# Empirical distribution of a hypergraph measure

The table returned by
[`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md)
and
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md)
with `what = "distribution"`: one row per distinct observed value with
its count, proportion and complementary cumulative share, the CCDF that
the descriptive figures of Coupette et al. (2024) draw on a logarithmic
axis. A grouping column (`time`, `s`) present in the table gives one
curve per group.

## Usage

``` r
# S3 method for class 'hypernets_distribution'
plot(x, log = TRUE, ...)
```

## Arguments

- x:

  A `hypernets_distribution` table.

- log:

  Draw the CCDF on a logarithmic y axis (default `TRUE`).

- ...:

  Unused; for S3 consistency.

## Value

A ggplot object.
