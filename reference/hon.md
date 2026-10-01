# Higher-order network (HON) from sequences

Builds a higher-order network: a state is split into memory nodes
(`"a -> b"`: state `b` reached from `a`) wherever the next-state
distribution depends on where the sequence came from, as decided by the
BuildHON / BuildHON+ rule extraction.

## Usage

``` r
hon(
  data,
  max_order = 5L,
  min_freq = 1L,
  collapse_repeats = FALSE,
  method = "hon+",
  action = NULL,
  actor = NULL,
  time = NULL,
  session = NULL,
  time_threshold = 900,
  timezone = "UTC",
  group = NULL
)
```

## Arguments

- data:

  Sequences in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md):
  a long event table (with `action`), a wide data.frame, a list of
  vectors, or a model object carrying its sequences.

- max_order:

  Integer. Highest order considered. Default `5`.

- min_freq:

  Integer. Transitions observed fewer times are treated as zero. Default
  `1`.

- collapse_repeats:

  Logical. Collapse adjacent repeated states before extraction. Default
  `FALSE`.

- method:

  `"hon+"` (default, parameter-free BuildHON+) or `"hon"` (the original
  eager BuildHON).

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide, list or model input.

- group:

  `NULL` (default) for one network, or the grouping of the sequences for
  a group model: the name of a column of `data` (long or wide format),
  or a vector with one label per sequence (per row of a wide data frame,
  per element of a list). A group model holds one network per group and
  the data of each group; it is the input of
  [`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md)
  and can be passed to
  [`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md)
  and
  [`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md).

## Value

A `net_hon` object (also a `cograph_network`, so it plots directly).
Read its tables with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the rules (`what = "rules"`, the default), the memory nodes (`"nodes"`)
and the higher-order pathways (`"pathways"`).

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances*, 2(5),
e1600028.
[doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028)

Saebi, M., Xu, J., Kaplan, L. M., Ribeiro, B., & Chawla, N. V. (2020).
Efficient modeling of higher-order dependencies in networks: from
algorithm to application for anomaly detection. *EPJ Data Science*, 9,
15.
[doi:10.1140/epjds/s13688-020-00233-y](https://doi.org/10.1140/epjds/s13688-020-00233-y)

## See also

[`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md),
[`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md),
[`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md),
[`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md)

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
net <- hon(seqs, max_order = 2)
hg_get(net)
#>          path   from to count probability from_order to_order
#> 1      a -> b      a  b     4         1.0          1        2
#> 2 a -> b -> c a -> b  c     4         1.0          2        1
#> 3      b -> c      b  c     4         0.5          1        1
#> 4      b -> d      b  d     4         0.5          1        1
#> 5      c -> a      c  a     2         1.0          1        1
#> 6      d -> x      d  x     2         1.0          1        1
#> 7      x -> b      x  b     4         1.0          1        2
#> 8 x -> b -> d x -> b  d     4         1.0          2        1

# the same network from a long event table
events <- data.frame(
  who  = rep(c("s1", "s2", "s3", "s4"), each = 6),
  step = rep(1:6, 4),
  what = unlist(seqs)
)
hg_get(hon(events, action = "what", actor = "who", time = "step",
           max_order = 2), what = "nodes")
#> Error in Nestimate::prepare(data[columns], actor = actor, action = action,     time = time, session = session, time_threshold = time_threshold,     timezone = timezone): unused argument (timezone = timezone)
```
