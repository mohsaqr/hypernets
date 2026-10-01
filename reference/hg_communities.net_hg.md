# Stable Infomap communities of a hypergraph projection

The hypergraph method of
[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md).
Builds the normalized association graph of Coupette et al. (2024), runs
Infomap repeatedly, compares every pair of partitions with AMI, ARI and
NMI, and returns the run with the largest summed AMI as the medoid. The
paper uses 50 seeds and 100 Infomap trials per seed, which are the
defaults.

## Usage

``` r
# S3 method for class 'net_hg'
hg_communities(
  x,
  n_runs = 50L,
  trials = 100L,
  seeds = NULL,
  method = c("association", "citation"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL,
  directed = FALSE,
  type = c("infomap", "irmm"),
  delta = 0.01,
  max_iter = 50L,
  edge_weights = NULL,
  ...
)
```

## Arguments

- x:

  A static `net_hg`.

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

- type:

  The community algorithm: `"infomap"` (default; Infomap on the
  projection chosen by `method`, as in Coupette et al. 2024) or `"irmm"`
  (iteratively reweighted modularity maximisation, Kumar et al. 2020:
  Louvain on the random-walk clique reduction \\A = H W (D_e - I)^{-1}
  H^T\\, then every hyperedge is reweighted to \\w'(e) =
  \frac{1}{m}\sum\_{i=1}^{c} \frac{\delta(e) + c}{k_i(e) + 1}\\ and
  averaged with its previous weight, until the largest weight change is
  at most `delta`). IRMM reads the hypergraph directly, so `trials`,
  `method`, `duplicate_edges`, `self_association`, `edge_source` and
  `directed` do not apply to it and raise `hypernets_bad_input` when
  supplied. Each of the `n_runs` runs is one full IRMM fit whose Louvain
  steps draw from R's RNG under that run's seed (the caller's RNG state
  is restored); igraph is required.

- delta:

  For `type = "irmm"`: stop when no hyperedge weight changes by more
  than `delta` in a pass (default 0.01, the paper's threshold).

- max_iter:

  For `type = "irmm"`: maximum reweighting passes (default 50). A run
  that hits it is flagged `converged = FALSE` in the `"runs"` table and
  raises one `hypernets_no_converge` warning.

- edge_weights:

  For `type = "irmm"`: initial positive hyperedge weights (one per
  hyperedge, or one value recycled). `NULL` uses the window counts of a
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  else unit weights.

- ...:

  Must be empty: an argument that only the memory-network method takes
  raises `hypernets_bad_input`.

## Value

An `hg_communities` object containing `medoid` (a tidy node/community
table), all `partitions`, AMI/ARI/NMI similarity matrices, run metadata,
community sizes, and the graph `projection` Infomap ran on. For
`type = "irmm"` the `"runs"` table has `run`, `seed`, `n_communities`,
`iterations`, `converged`, `max_weight_change` and `modularity` (the
linear hypergraph modularity of
[`hg_modularity()`](https://mohsaqr.github.io/hypernets/reference/hg_modularity.md)),
`projection` is the reweighted clique reduction of the medoid run, and
`hg_get(fit, what = "weights")` gives one row per hyperedge with `edge`,
`size`, `initial_weight` and the medoid run's final `weight`.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

Kumar, T., Vaidyanathan, S., Ananthapadmanabhan, H., Parthasarathy, S.,
& Ravindran, B. (2020). Hypergraph clustering by iteratively reweighted
modularity maximization. *Applied Network Science*, 5, 52.
[doi:10.1007/s41109-020-00300-3](https://doi.org/10.1007/s41109-020-00300-3)

Blondel, V. D., Guillaume, J.-L., Lambiotte, R., & Lefebvre, E. (2008).
Fast unfolding of communities in large networks. *Journal of Statistical
Mechanics*, 2008(10), P10008.
[doi:10.1088/1742-5468/2008/10/P10008](https://doi.org/10.1088/1742-5468/2008/10/P10008)

## See also

[`hg_modularity()`](https://mohsaqr.github.io/hypernets/reference/hg_modularity.md)
to score any partition.

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
  edge = rep(paste0("e", 1:4), each = 3)
)
h <- group_hypergraph(dat, "member", "edge")
if (requireNamespace("igraph", quietly = TRUE)) {
  fit <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
  hg_get(fit)
  irmm <- hg_communities(h, type = "irmm", n_runs = 3, seeds = 1:3)
  hg_get(irmm)
  hg_get(irmm, what = "runs")
}
#>   run seed n_communities iterations converged max_weight_change modularity
#> 1   1    1             2          6      TRUE       0.008789062       0.25
#> 2   2    2             2          6      TRUE       0.008789062       0.25
#> 3   3    3             2          6      TRUE       0.008789062       0.25
```
