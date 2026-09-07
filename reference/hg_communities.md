# Stable Infomap communities of a hypergraph projection

Builds the normalized association graph of Coupette et al. (2024), runs
Infomap repeatedly, compares every pair of partitions with AMI, ARI and
NMI, and returns the run with the largest summed AMI as the medoid. The
paper uses 50 seeds and 100 Infomap trials per seed, which are the
defaults.

## Usage

``` r
hg_communities(
  hg,
  n_runs = 50L,
  trials = 100L,
  seeds = NULL,
  method = c("association", "citation"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL,
  directed = FALSE
)

# S3 method for class 'hg_communities'
print(x, ...)

# S3 method for class 'hg_communities'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("medoid", "partitions", "runs", "sizes", "ami", "ari", "nmi"),
  ...
)

# S3 method for class 'hg_communities'
plot(x, ...)

hypergraph_communities(
  hg,
  n_runs = 50L,
  trials = 100L,
  seeds = NULL,
  method = c("association", "citation"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL,
  directed = FALSE
)
```

## Arguments

- hg:

  A static `net_hypergraph`.

- n_runs:

  Number of independent seeded Infomap runs (default 50).

- trials:

  Infomap trials within each run (default 100).

- seeds:

  Integer seeds. `NULL` uses `seq_len(n_runs)`.

- method:

  Which graph Infomap runs on: the `"association"` projection (default;
  the paper's hypergraph-derived representations `bh`, `bhs`, `mh`,
  `mhs`) or the `"citation"` projection (its classic graph
  representations `bg`, `mg`, and with `directed = FALSE` their
  undirected variants `bgu`, `mgu`). See
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md).

- duplicate_edges, self_association, edge_source:

  Projection controls passed to
  [`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md).
  Together these reproduce the paper's binary/multi and self-association
  representations.

- directed:

  For `method = "citation"`: run Infomap with directed flow on the
  source-to-member graph? Default `FALSE`.

- x:

  An `hg_communities` object.

- ...:

  Additional arguments passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

- row.names, optional:

  Unused; present for the base S3 contract.

- what:

  Component to return: `"medoid"`, `"partitions"`, `"runs"`, `"sizes"`,
  `"ami"`, `"ari"`, or `"nmi"`.

## Value

An `hg_communities` object containing `medoid` (a tidy node/community
table), all `partitions`, AMI/ARI/NMI similarity matrices, run metadata,
community sizes, and the graph `projection` Infomap ran on.

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
  fit <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
  as.data.frame(fit)
}
#>   node community
#> 1    a         1
#> 2    b         1
#> 3    c         1
#> 4    x         2
#> 5    y         2
#> 6    z         2
```
