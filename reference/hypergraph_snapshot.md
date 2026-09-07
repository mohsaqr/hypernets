# Extract one snapshot from a temporal hypergraph

A snapshot is a static `net_hypergraph` holding the hyperedges a window
measures. `mode = "active"` follows the format: an interval hyperedge is
active on its closed interval, a contact hyperedge at its instant.
`mode = "cumulative"` keeps every hyperedge begun by the end of the
window, the growing view of a contact log; with `at` unset it is the
static aggregate of everything observed.

## Usage

``` r
hypergraph_snapshot(
  x,
  at = NULL,
  window = 0,
  mode = c("active", "cumulative"),
  multiedges = TRUE
)

hg_snapshot(
  x,
  at = NULL,
  window = 0,
  mode = c("active", "cumulative"),
  multiedges = TRUE
)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- at:

  One time value, as a number on the hypergraph's clock or a date for a
  calendar hypergraph. `NULL` (default) is the end of observation.

- window:

  How much time the snapshot covers, in the hypergraph's time unit: `0`
  (default) evaluates the instant `at`; a positive width covers
  `[at, at + window)`; `"all"` covers the whole observation period.

- mode:

  `"active"` (default) or `"cumulative"`, see Description. The former
  `"all"` is deprecated: it is `"cumulative"` with `at` unset.

- multiedges:

  Keep distinct edge identities with identical member sets? `TRUE`
  matches a multi-hypergraph; `FALSE` collapses them to one edge and
  records their counts in `edge_multiplicity`.

## Value

A static `net_hypergraph` usable by every hypernets hypergraph verb; its
`params` record `at`, `window`, `temporal_mode` and the clock.
