# Extract a sequence of temporal-hypergraph snapshots

Snapshots on a measurement grid with Dynet's meaning: `start` and `end`
bound the measured period (default: the observation window), `step` is
how often to look, `window` how much time each look covers, and `at`
names instants directly. Without `step` and `at`, the grid is the event
times of the hypergraph, each looked at as a point.

## Usage

``` r
hypergraph_snapshots(
  x,
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  mode = c("active", "cumulative"),
  at = NULL,
  multiedges = TRUE
)

hg_snapshots(
  x,
  start = NULL,
  end = NULL,
  step = NULL,
  window = NULL,
  mode = c("active", "cumulative"),
  at = NULL,
  multiedges = TRUE
)
```

## Arguments

- x:

  A
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md).

- start, end:

  Bounds of the measured period, as numbers on the hypergraph's clock or
  dates for a calendar hypergraph. Default: the observation window.

- step:

  How often to measure, in the hypergraph's time unit. `NULL` (default)
  measures at every event time.

- window:

  How much time each measurement covers. Defaults to `step` when a step
  is given (a partition of the period) and to `0` otherwise (a point
  sample); larger than `step` gives a rolling window; `"all"` measures
  the whole period as one window.

- mode:

  `"active"` (default) or `"cumulative"`, see Description. The former
  `"all"` is deprecated: it is `"cumulative"` with `at` unset.

- at:

  Instants to measure instead of a grid; may be a vector.

- multiedges:

  Keep distinct edge identities with identical member sets? `TRUE`
  matches a multi-hypergraph; `FALSE` collapses them to one edge and
  records their counts in `edge_multiplicity`.

## Value

A named list of `net_hypergraph` objects with class
`net_hypergraph_snapshots`, named by the window starts on the
hypergraph's clock.

## Examples

``` r
seats <- data.frame(
  case = rep(c("A", "B", "C"), each = 3),
  arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3", "p1", "a4", "a5"),
  constituted = rep(c(1, 2, 4), each = 3), concluded = rep(c(4, 3, 6), each = 3)
)
thg <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                           start = "constituted", end = "concluded")
yearly <- hypergraph_snapshots(thg, step = 2)
names(yearly)
#> [1] "1" "3" "5"
```
