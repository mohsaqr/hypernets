# Hyperedge-level measures, as a tidy table

One row per hyperedge. `hg_measures(what = "edges")` reports cardinality
alone; this adds the incident weight, how many other hyperedges a
hyperedge meets, and how large its vertex neighbourhood is — the
hyperedge-side counterparts of degree and neighbourhood size, and the
descriptive layer used to characterise higher-order structure (Coupette
et al. 2024). A temporal hypergraph is evaluated snapshot by snapshot,
with a leading `time` column, so the distribution of hyperedge sizes at
several dates (the paper's Figure 4b) or the neighbourhood size of a
tribunal over time (Figure 5c) come from one call.

## Usage

``` r
hg_edges(
  hg,
  what = c("edges", "distribution", "summary"),
  measure = c("size", "weight", "n_incident_edges", "n_neighbors"),
  s = 1L,
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  at = NULL,
  snapshot_mode = c("active", "cumulative"),
  multiedges = TRUE
)

hypergraph_edges(
  hg,
  what = c("edges", "distribution", "summary"),
  measure = c("size", "weight", "n_incident_edges", "n_neighbors"),
  s = 1L,
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

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md),
  any hypernets `net_hypergraph`, or a
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- what:

  `"edges"` (default) for one row per hyperedge, `"distribution"` for
  the empirical distribution of `measure` across hyperedges, or
  `"summary"` for its mean, standard deviation and quartiles.

- measure:

  Which column `what = "distribution"` and `what = "summary"` describe:
  `"size"` (default), `"weight"`, `"n_incident_edges"` or
  `"n_neighbors"`.

- s:

  Minimum number of shared vertices for another hyperedge to count as
  incident. The default `1` is ordinary incidence; larger values match
  the thresholds used by
  [`hg_edge_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_edge_centrality.md).
  Several values give one block of rows each, with an `s` column.

- start, end, step, window, at:

  Measurement grid passed to
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  when `hg` is temporal: the bounds of the period, how often to look,
  how much time each look covers, or the instants themselves.

- snapshot_mode, multiedges:

  Snapshot `mode` (`"active"` or `"cumulative"`) and multi-edge handling
  passed to
  [`hypergraph_snapshots()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshots.md)
  when `hg` is temporal.

## Value

With `what = "edges"`, a base data.frame with one row per hyperedge and
columns `edge` (name), `size` (integer, vertices it contains), `weight`
(numeric, its incidence weights summed), `n_incident_edges` (integer,
other hyperedges sharing at least `s` vertices) and `n_neighbors`
(integer, vertices adjacent to a member without being one). With
`what = "distribution"`, a `hypernets_distribution` table with one row
per distinct observed value of `measure`, ascending, with columns
`value`, `n`, `proportion` and `ccdf` — the complementary cumulative
distribution \\P(X \ge value)\\, so the first row's `ccdf` is always 1;
its [`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
CCDF. With `what = "summary"`, one row with `n_edges`, `mean`, `sd`,
`min`, `q25`, `median`, `q75` and `max`. Several `s` values add an `s`
column; temporal input adds a leading `time` column, and the summary is
then a `hypernets_series` whose
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws each
statistic against time.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## See also

[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md)
for vertex-level measures,
[`hg_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md)
for the graph whose degrees `n_incident_edges` reports.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
                        c = "stars and salt"))
hg_edges(hg)
#>    edge size weight n_incident_edges n_neighbors
#> 1   and    3      3                3           0
#> 2  salt    2      2                3           1
#> 3  soup    2      2                3           1
#> 4 stars    2      2                3           1
hg_edges(hg, what = "distribution")
#>   value n proportion ccdf
#> 1     2 3       0.75 1.00
#> 2     3 1       0.25 0.25
hg_edges(hg, what = "summary", measure = "n_neighbors")
#>   n_edges mean  sd min  q25 median q75 max
#> 1       4 0.75 0.5   0 0.75      1   1   1
```
