# Y, T and O motif census with a hypergraph configuration null

Counts the three possible induced four-node motifs in a simple 3-uniform
hypergraph: Y contains two of the four possible triples, T contains
three, and O contains all four. Repeated identical hyperedges are
collapsed for the motif census, matching HypergraphX 1.5. The null
repeatedly reshuffles pairs of equal-size hyperedges while preserving
every node degree and edge cardinality before duplicate-edge collapse.

## Usage

``` r
hg_motifs(
  hg,
  n = 1000L,
  seed = NULL,
  what = c("test", "counts", "draws"),
  alternative = c("two_sided", "greater", "less"),
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  at = NULL,
  snapshot_mode = c("active", "cumulative"),
  multiedges = TRUE
)

hypergraph_motifs(
  hg,
  n = 1000L,
  seed = NULL,
  what = c("test", "counts", "draws"),
  alternative = c("two_sided", "greater", "less"),
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  at = NULL,
  snapshot_mode = c("active", "cumulative"),
  multiedges = TRUE
)

# S3 method for class 'honets_motifs'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("test", "draws"),
  ...
)

# S3 method for class 'honets_motifs'
plot(x, motif = c("Y", "T", "O"), ...)
```

## Arguments

- hg:

  A 3-uniform `net_hypergraph` or a
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- n:

  Number of configuration-model draws. The paper uses 1000.

- seed:

  Optional reproducibility seed; the caller's RNG state is restored on
  exit.

- what:

  `"test"` for observed/null summaries, `"counts"` for observed counts
  only, or `"draws"` for every null count.

- alternative:

  Empirical permutation-test direction.

- start, end, step, window, at:

  Measurement grid passed to
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  when `hg` is temporal.

- snapshot_mode, multiedges:

  Snapshot `mode` (`"active"` or `"cumulative"`) and multi-edge handling
  passed to
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  when `hg` is temporal.

- x:

  A `honets_motifs` test table.

- row.names, optional:

  Unused; present for the base S3 contract.

- ...:

  Unused; for S3 consistency.

- motif:

  Which motif's null distribution to draw: `"Y"` (default), `"T"` or
  `"O"`.

## Value

A tidy data frame. Test output includes observed count, null mean and
standard deviation, z-score, empirical p-value, relative abundance
`delta`, and the normalized motif profile used by HypergraphX.

For `as.data.frame`, the test table (`what = "test"`) or every null
count (`what = "draws"`, columns `run`, `motif`, `count`) as a plain
data.frame.

For `plot`, a ggplot object: the null distribution of the motif count as
a histogram with the observed count as a vertical line and the z-score
annotated (the paper's Figure 7).

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
dat <- data.frame(
  member = c(1, 2, 3, 1, 2, 4), edge = rep(c("e1", "e2"), each = 3)
)
h <- group_hypergraph(dat, "member", "edge")
hg_motifs(h, what = "counts")
#>   motif count
#> 1     Y     1
#> 2     T     0
#> 3     O     0
```
