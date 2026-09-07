# Compare community structure across representations

Sets several
[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
fits of the same data side by side, as the paper does for its eight GFCC
representations (Coupette et al. 2024, Figure 8): the cluster-size
distribution of each AMI medoid, the pairwise AMI, ARI and NMI between
medoids, and per-fit summaries (number of communities, singletons, the
balance between the two largest clusters).

## Usage

``` r
hg_compare_communities(..., hg = NULL, edge_source = NULL)

# S3 method for class 'honets_community_comparison'
print(x, ...)

# S3 method for class 'honets_community_comparison'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("summary", "similarity", "sizes", "quality", "matrix"),
  ...
)

# S3 method for class 'honets_community_comparison'
summary(object, ...)

# S3 method for class 'honets_community_comparison'
plot(x, what = c("sizes", "similarity"), ...)
```

## Arguments

- ...:

  For `plot`, unused.

- hg:

  Optional: the static `net_hypergraph` the fits were computed on. When
  given, every medoid is scored with
  [`hg_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md)
  on the projection its own fit used, and the scores are available as
  `what = "quality"`.

- edge_source:

  Hyperedge sources for the citation and self-association projections
  when scoring, as in
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md).

- x:

  A `honets_community_comparison` object.

- row.names, optional:

  Unused; present for the base S3 contract.

- what:

  Which table: `"summary"` (default), `"similarity"`, `"sizes"`,
  `"quality"`, or `"matrix"` for the similarity as a square matrix with
  AMI below and ARI above the diagonal (the layout of
  [`plot()`](https://rdrr.io/r/graphics/plot.default.html)).

- object:

  A `honets_community_comparison` object.

## Value

A `honets_community_comparison` object.
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
its `"summary"` (default; one row per fit with `model`, `medoid_seed`,
`n_communities`, `n_singletons`, `n_nontrivial`, `largest`, `second` and
`balance` = second / largest), `"similarity"` (one row per pair of fits
with `model_a`, `model_b`, `ami`, `ari`, `nmi`, on the nodes the two
medoids share), `"sizes"` (one row per community of every medoid with
`model`, `rank`, `n_nodes`) or, when `hg` was given, `"quality"` (one
row per fit with the columns of
[`hg_community_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_community_quality.md)).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
cluster-size distributions (`what = "sizes"`, the number of communities
at least as large as each size, on logarithmic axes) or the similarity
matrix (`what = "similarity"`, AMI below and ARI above the diagonal,
values printed in the cells).

For `plot`, a ggplot object.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
  edge = rep(paste0("e", 1:4), each = 3)
)
h <- group_hypergraph(dat, "member", "edge")
if (requireNamespace("igraph", quietly = TRUE)) {
  multi <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
  binary <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2,
                           duplicate_edges = "collapse")
  comparison <- hg_compare_communities(mh = multi, bh = binary)
  as.data.frame(comparison)
  as.data.frame(comparison, what = "similarity")
}
#>   model_a model_b n_nodes ami ari nmi
#> 1      mh      bh       6   1   1   1
```
