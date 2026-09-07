# Temporal hypergraph from an edge list or co-occurrence data

Builds a temporal hypergraph the way a temporal network is defined from
data in this ecosystem (Dynet's `dynet()`): one row per relation, with
the columns that name its ends and its clock. Two input shapes are
accepted. An **edge list** names `from` and `to`, and every row is a
hyperedge of size two, so an ordinary temporal network is the same
object. **Co-presence data** name an `actor` and a `group`: every actor
sharing one value of `group` (a case, a citation block, a seminar)
belongs to one hyperedge. Two clocks are understood, and the one you
name selects the format:

## Usage

``` r
temporal_hypergraph(
  data,
  from = NULL,
  to = NULL,
  actor = NULL,
  group = NULL,
  time = NULL,
  start = NULL,
  end = NULL,
  weight = NULL,
  nodes = NULL,
  time_unit = "auto",
  observation_start = NULL,
  observation_end = NULL,
  sparse = FALSE,
  cooccur_by = NULL
)
```

## Arguments

- data:

  A data frame with one row per relation, or a sequence table (a wide
  data frame of states, a list of character vectors, or a `tna` /
  `netobject` model); see Details.

- from, to:

  Column names of a pairwise edge list. Detected from the alias table
  when neither is given and `actor`/`group` are not named.

- actor, group:

  Column names of co-presence data: the node, and the grouping whose
  shared values bind nodes into one hyperedge. Naming either selects the
  co-presence format; the other is then detected by alias if not given.

- time:

  Column with the instant of a contact hyperedge.

- start, end:

  Columns with the interval on which a hyperedge is active. `start`
  without an `end` column is read as a contact clock.

- weight:

  Optional membership-weight column.

- nodes:

  Optional node universe: a character vector of node names, or a data
  frame whose first column holds the names and whose `start`, `time` or
  `date` column, if present, gives the time each node enters (on the
  same clock as `data`). Nodes that never appear in a hyperedge are then
  kept as zero-degree nodes of every snapshot, and
  [`hg_growth()`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  counts a node from its own start. Every observed node must be in the
  universe.

- time_unit:

  Unit for calendar times: `"auto"` (default), `"seconds"`, `"minutes"`,
  `"hours"`, `"days"` or `"weeks"`. Numeric times are left alone and
  reported as `"step"`.

- observation_start, observation_end:

  Optional bounds of the observation window, as numbers on the stored
  clock or as dates for a calendar hypergraph. Either may be omitted;
  the corresponding limit of the data is then used.

- sparse:

  Store every snapshot's incidence as a sparse `Matrix`? Default
  `FALSE`.

- cooccur_by:

  Deprecated name of `group`; using it warns with a `honets_deprecated`
  condition.

## Value

A `net_temporal_hypergraph` holding the membership table (`node`,
`edge`, `start`, `end`, `weight`), the edge metadata (`edge`, `start`,
`end` and the hyperedge attributes), the node universe with entry times,
the sorted event times, `format` (`"interval"` or `"contact"`),
`time_unit`, `origin` and the `observation` bounds.
`as.data.frame(x, what = "memberships" | "edges" | "nodes")` returns the
three tables; [`summary()`](https://rdrr.io/r/base/summary.html) the
one-row description.

## Details

- interval:

  `start` and `end`: the hyperedge is active on the closed interval
  between them. A missing `end` means it stays active through the end of
  observation.

- contact:

  `time`: the hyperedge is an instantaneous event, as in a contact log –
  a citation, a message, a meeting on one day. It is present at that
  instant only. The cumulative view in which every hyperedge stays once
  it has appeared (the point-aggregation model of Coupette et al. 2024)
  is a snapshot `mode`, not a property of the data: ask for it with
  `mode = "cumulative"` in
  [`hypergraph_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_snapshot.md),
  [`hg_growth()`](https://mohsaqr.github.io/hypernets/reference/hg_growth.md)
  and
  [`hg_edges()`](https://mohsaqr.github.io/hypernets/reference/hg_edges.md).

A third shape is a **sequence table**: one row per session and one
column per position holding the state at that step (the wide format of
tna and TraMineR), a list of character vectors, or a `tna` / `netobject`
model built from one. It is recognised when no relational column is
named or detected and every column is categorical. Each session is one
hyperedge whose members are the states it contains, and a state's
membership is a contact at its position, `1` to the session's length, on
a `"step"` clock: the simple co-occurrence reading in which time is
order.

In any shape a membership may carry its own time, as when a log has one
row per attendance rather than one time per group. The hyperedge then
spans from its first to its last membership, each membership is present
on its own spell only, and a snapshot keeps the memberships present in
its window: `mode = "cumulative"` at step `t` is what each session had
shown by `t`, and `window = 3` at `t` is what it showed on steps `t` to
`t + 2`. When every membership carries its hyperedge's time, as a
tribunal or a citation block does, nothing changes.

Column names are resolved case-insensitively from the same alias table
Dynet uses, so `Sender`/`Receiver`, `source`/`target`,
`onset`/`terminus` and `timestamp` are understood without being spelled
out; a name you give explicitly must exist as written. Times may be
numeric, `Date`, `POSIXct` or character date-time strings. Numeric times
are kept as they are and reported in `"step"` units; calendar input is
converted to elapsed time since the earliest time in the data, in a unit
chosen for the span (`"days"` beyond three days) or given as
`time_unit`, and the unit and origin are stored in the object and shown
by [`print()`](https://rdrr.io/r/base/print.html). Every time you pass
later – `at`, `start`, `end`, the observation bounds – may be a date for
a calendar hypergraph or a number on that clock, and every `time` column
in a result is on that clock.

`observation_start` and `observation_end` declare the study window when
it is known independently of the log, with Dynet's meaning: they bound
the snapshot times and the measurement grid, and an open-ended hyperedge
is active through `observation_end`, but the stored memberships are
never rewritten and
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the original spells. Without them the window is the span of the data.

Every other column that is constant within a hyperedge is kept as a
hyperedge attribute in the edge metadata, where
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) can colour by
it and where
[`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md)
finds the source a citation block belongs to. Columns that vary within a
hyperedge, such as the seat an arbitrator held, are not attributes of
the hyperedge and are left out.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
seats <- data.frame(
  case = c("A", "A", "A", "B", "B", "B"),
  arbitrator = c("p1", "a1", "a2", "p2", "a1", "a3"),
  constituted = c(1, 1, 1, 2, 2, 2), concluded = c(4, 4, 4, 5, 5, 5),
  sector = c("oil", "oil", "oil", "gas", "gas", "gas")
)
thg <- temporal_hypergraph(seats, actor = "arbitrator", group = "case",
                           start = "constituted", end = "concluded")
thg
#> Temporal hypergraph: 5 nodes, 2 hyperedges, 4 event times
#> Format: interval (a hyperedge is active from its start to its end)
#> Time: time (numeric steps); observed from 1 to 5
hypergraph_snapshot(thg, at = 3)
#> Hypergraph: 5 nodes, 2 hyperedges
#> Size distribution:
#>   size_3   : 2
#> Source: group membership (member = member, group = edge)

# a contact log on a calendar: instants at each date, cumulative on request
contacts <- data.frame(from = c("a", "b", "c"), to = c("b", "c", "a"),
                       date = as.Date(c("2024-01-01", "2024-01-05", "2024-01-09")))
calls <- temporal_hypergraph(contacts, from = "from", to = "to", time = "date")
calls
#> Temporal hypergraph: 3 nodes, 3 hyperedges, 3 event times
#> Format: contact (a hyperedge is an instantaneous event; mode = "cumulative" keeps every hyperedge once it appears)
#> Time: days since 2024-01-01; observed from 0 to 8
hypergraph_snapshot(calls, at = as.Date("2024-01-05"))
#> Hypergraph: 2 nodes, 1 hyperedges
#> Size distribution:
#>   size_2   : 1
#> Source: group membership (member = member, group = edge)
hypergraph_snapshot(calls, at = as.Date("2024-01-05"), mode = "cumulative")
#> Hypergraph: 3 nodes, 2 hyperedges
#> Size distribution:
#>   size_2   : 2
#> Source: group membership (member = member, group = edge)
```
