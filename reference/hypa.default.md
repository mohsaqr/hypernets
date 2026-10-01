# Path anomalies (HYPA) in sequences

The sequence method of
[`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md).
Builds the De Bruijn graph of order `order` from the sequences and
scores every path of `order + 1` states against a hypergeometric null
whose propensities are fitted to the observed path counts: paths
observed far more often than the null predicts are over-represented, far
less often under-represented.

## Usage

``` r
# Default S3 method
hypa(
  x,
  order = 2L,
  alpha = 0.05,
  min_count = 5L,
  p_adjust = "BH",
  type = c("all", "over", "under"),
  order_by = c("sig", "ratio", "freq", "path"),
  n = 10L,
  action = NULL,
  actor = NULL,
  time = NULL,
  session = NULL,
  time_threshold = 900,
  timezone = "UTC",
  ...
)

# S3 method for class 'hypernets_hypa'
print(x, n = NULL, ...)
```

## Arguments

- x:

  Sequences in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md):
  a long event table (with `action`), a wide data.frame, a list of
  vectors, or a model object carrying its sequences.

- order:

  Integer (or integer vector) of De Bruijn orders. Default `2`.

- alpha:

  Significance level of the anomaly labels, in (0, 0.5). Default `0.05`.

- min_count:

  Integer. Paths observed fewer times are always labelled `"normal"`.
  Default `5`.

- p_adjust:

  Multiplicity correction: any
  [stats::p.adjust](https://rdrr.io/r/stats/p.adjust.html) method or
  `"none"`. Default `"BH"`; the two tails are adjusted separately.

- type:

  Which anomalous paths the printed result lists: `"all"` (default),
  `"over"` (paths observed more often than expected) or `"under"` (less
  often).

- order_by:

  How the printed paths are ordered: `"sig"` (default, smallest p-value
  of the path's own direction first), `"ratio"` (largest deviation
  first: highest observed/expected ratio for over-represented paths,
  lowest for under-represented ones), `"freq"` (most observed first) or
  `"path"` (alphabetical).

- n:

  Number of paths printed per direction. Default `10`.

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide, list or model input.

- ...:

  Must be empty: an argument that only the hypergraph method takes
  (`top`) raises `hypernets_bad_input`.

## Value

A `net_hypa` object (also a `cograph_network`) that prints the selected
anomalous paths. `type`, `order_by` and `n` change only what is printed
and the default of
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md);
every path is scored. Read the tables with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the anomalous paths (`what = "anomalies"`, the default), every scored
path (`"scores"`), the over- or under-represented ones (`"over"`,
`"under"`), or the anomalous paths as pathways (`"pathways"`).

## References

LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies
in time series data on networks. *Proceedings of the 2020 SIAM
International Conference on Data Mining*, 460-468.
[doi:10.1137/1.9781611976236.52](https://doi.org/10.1137/1.9781611976236.52)

## See also

[`hypa.net_hg()`](https://mohsaqr.github.io/hypernets/reference/hypa.net_hg.md)
for the same null on hypergraph co-occurrence.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
fit <- hypa(seqs, order = 2, min_count = 1, type = "over")
fit
#> Path anomalies (HYPA, order 2): 0 of 6 paths anomalous at alpha = 0.05 (0 over, 0 under, BH)
#> (no rows)
hg_get(fit)
#> [1] path      count     expected  ratio     p_adj     direction
#> <0 rows> (or 0-length row.names)
```
