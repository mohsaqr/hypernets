# Hypergeometric anomaly detection for co-occurring node pairs

Scores every co-occurring pair of nodes against a hypergeometric null in
which a pair's propensity to share hyperedges is the product of the two
nodes' hyperdegrees. A pair that co-occurs far more often than that
product predicts is over-represented; far less often, under-represented.
This is the hypergraph counterpart of
[`build_hypa()`](https://mohsaqr.github.io/hypernets/reference/build_hypa.md),
which applies the same null to the transitions of a higher-order
network, and it returns the same columns so the two read alike.

## Usage

``` r
hg_hypa(hg, min_count = 5L, alpha = 0.05, p_adjust = "BH", top = NULL)
```

## Arguments

- hg:

  A `net_hypergraph`, or a snapshot of a temporal one.

- min_count:

  Minimum observed co-occurrence count for a pair to be scored (default
  `5L`). Pairs below the floor are not tested and do not enter the
  multiplicity correction.

- alpha:

  Significance threshold for the `anomaly` column (default `0.05`),
  applied to the adjusted p-values.

- p_adjust:

  Multiplicity correction, any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method or
  `"none"` (default `"BH"`).

- top:

  Optional number of rows to keep, the most significantly
  over-represented first. `NULL` (default) keeps every scored pair.

## Value

A base `data.frame`, one row per scored pair, most significantly
over-represented first (by `p_adjusted_over`, ties broken by `ratio`),
with columns `from`, `to`, `observed` (co-occurrence count), `expected`,
`ratio` (observed / expected), `p_under`, `p_over`, `p_adjusted_under`,
`p_adjusted_over` and `anomaly` (`"over"`, `"under"` or `"normal"`).

## Details

The test is analytic, so it needs no resampling. That is what makes it
usable where
[`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)
is not:
[`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)
answers whether the hypergraph as a whole carries more repetition than
chance, by Monte Carlo, and returns one row; this returns one row per
pair.

For a pair \\(u, v)\\ the propensity is \\\Xi\_{uv} = d_u d_v\\, the
product of the hyperdegrees. With \\N = \sum \Xi\\ over scored pairs,
\\K = \Xi\_{uv}\\ and \\n\\ the total co-occurrence mass, the observed
count is referred to a \\\mathrm{Hypergeometric}(N, K, n)\\ law:
`p_under` is \\P(X \le f)\\ and `p_over` is \\P(X \ge f)\\.

`min_count` is not a convenience. A pair seen once or twice cannot reach
significance after correction, so scoring it spends multiplicity budget
that the pairs which can reach significance then have to pay for.
Raising the floor raises the power of the whole test. The floor selects
the hypotheses by count, never by p-value, so restricting the correction
to them is legitimate. \\N\\ and the draw count are still taken over
every co-occurring pair, so the null describes the whole hypergraph.

Time is handled by scoring a cumulative snapshot rather than by slicing
into periods. Disjoint periods split a pair's evidence across cells
while multiplying the number of tests, so both terms move the wrong way;
a sequence of `hypergraph_snapshot(hg, at = ..., mode = "cumulative")`
accumulates evidence instead and shows when a pair becomes anomalous.

## References

LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies
in time series data on networks. *Proceedings of the 2020 SIAM
International Conference on Data Mining*, 460-468.

## See also

[`build_hypa()`](https://mohsaqr.github.io/hypernets/reference/build_hypa.md)
for the same null on a higher-order network,
[`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)
for the Monte Carlo whole-hypergraph tests.

## Examples

``` r
memberships <- data.frame(
  actor = c("a","b","c", "a","b","d", "a","b","e", "c","d","e", "a","b","f"),
  group = rep(c("g1","g2","g3","g4","g5"), each = 3)
)
hg <- group_hypergraph(memberships, actor = "actor", group = "group")
hg_hypa(hg, min_count = 2L)
#>   from to observed expected ratio   p_under    p_over p_adjusted_under
#> 1    a  b        4 2.857143   1.4 0.8809635 0.3074243        0.8809635
#>   p_adjusted_over anomaly
#> 1       0.3074243  normal
```
