# s-betweenness and s-closeness centrality of hyperedges

Hyperedges become vertices in the unweighted s-line graph: two are
adjacent when they share at least `s` nodes. Betweenness and closeness
are then ordinary shortest-path vertex centralities on that graph. At
`s = 1` these are the paper's hyperedge 1-betweenness and 1-closeness.

## Usage

``` r
hg_edge_centrality(
  hg,
  s = 1L,
  measure = c("betweenness", "closeness"),
  normalized = TRUE,
  top = NULL,
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  at = NULL,
  snapshot_mode = c("active", "cumulative"),
  multiedges = TRUE
)

hypergraph_edge_centrality(
  hg,
  s = 1L,
  measure = c("betweenness", "closeness"),
  normalized = TRUE,
  top = NULL,
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  at = NULL,
  snapshot_mode = c("active", "cumulative"),
  multiedges = TRUE
)
```

## Arguments

- hg:

  A static `net_hypergraph` or a
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- s:

  One or more positive integer intersection thresholds.

- measure:

  Any of `"betweenness"` and `"closeness"`.

- normalized:

  Use NetworkX/HypergraphX normalization: undirected betweenness is
  scaled by `2 / ((n - 1) * (n - 2))`, and closeness uses the
  Wasserman–Faust correction for disconnected line graphs.

- top:

  Optional number of highest-scoring hyperedges to retain per
  `(time, s, measure)` group.

- start, end, step, window, at:

  Measurement grid passed to
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  when `hg` is temporal.

- snapshot_mode, multiedges:

  Snapshot `mode` (`"active"` or `"cumulative"`) and multi-edge handling
  passed to
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  when `hg` is temporal.

## Value

A tidy data frame with `edge`, `s`, `measure`, and `value`; temporal
input adds a leading `time` column.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
h <- group_hypergraph(
  data.frame(member = c("a", "b", "c", "b", "c", "d"),
             event = rep(c("e1", "e2"), each = 3)),
  "member", "event"
)
hg_edge_centrality(h, s = 1)
#>   edge s     measure value
#> 1   e1 1 betweenness     0
#> 2   e2 1 betweenness     0
#> 3   e1 1   closeness     1
#> 4   e2 1   closeness     1
```
