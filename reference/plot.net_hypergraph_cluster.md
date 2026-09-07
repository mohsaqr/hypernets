# Plot method for net_hypergraph_cluster

Two diagnostic panels. `"spectrum"`: scree plot of the Laplacian
spectrum with the k used for clustering marked - the eigengap after k
supports (or questions) the choice of k. `"embedding"`: the nodes in the
first two spectral-embedding dimensions, labelled, coloured and shaped
by cluster, sized by stationary probability - the geometry k-means
actually clustered. `"both"` (default) arranges the two side by side
(via gridExtra when available, base grid viewports otherwise).

## Usage

``` r
# S3 method for class 'net_hypergraph_cluster'
plot(x, what = c("both", "spectrum", "embedding"), n_values = NULL, ...)
```

## Arguments

- x:

  A `net_hypergraph_cluster` object.

- what:

  Character. `"both"` (default), `"spectrum"`, or `"embedding"`.

- n_values:

  Integer. How many smallest eigenvalues to show in the spectrum panel
  (default: `min(3 * k, n_nodes)`).

- ...:

  Additional arguments (ignored).

## Value

For `"spectrum"`/`"embedding"`, the ggplot object. For `"both"`, the
arranged gtable when gridExtra is installed (drawn on the current
device), otherwise the two panels are drawn via grid viewports and the
list of the two ggplots is returned invisibly.
