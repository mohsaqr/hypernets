# Windowed Sequence Hyperedges

Builds a hypergraph from categorical sequence data: a fixed-size window
moves over each sequence and the distinct states observed inside one
window form one hyperedge (Ding et al. 2020). Windows whose distinct
state sets coincide collapse into a single hyperedge whose weight is the
number of such windows (`window_counts`); the incidence cells hold the
total within-window occurrences of each state, so the incidence matrix
carries edge-dependent vertex weights (Chitra & Raphael 2019) that
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
with `type = "random_walk"` uses directly. The window counts are the
default hyperedge weights of the whole Laplacian family
([`hg_laplacian()`](https://mohsaqr.github.io/hypernets/reference/hg_laplacian.md),
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md)).

## Usage

``` r
window_hypergraph(
  data,
  window = 3L,
  step = 1L,
  action = NULL,
  actor = NULL,
  time = NULL,
  session = NULL,
  time_threshold = 900,
  timezone = "UTC",
  min_size = 1L,
  min_weight = 1L
)
```

## Arguments

- data:

  Sequences in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md):
  a long event table, a wide data.frame or character matrix (one
  sequence per row, trailing `NA`s stripped), or a list of character
  vectors.

- window:

  Integer \>= 2. Window size in sequence positions.

- step:

  Integer \>= 1. Offset between consecutive window starts: `1` slides,
  `window` tumbles.

- action, actor, time, session, time_threshold, timezone:

  Long-format arguments, as in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md).

- min_size:

  Integer \>= 1. Drop hyperedges with fewer distinct states after
  collapsing. The default `1` keeps everything.

- min_weight:

  Integer \>= 1. Drop hyperedges observed in fewer than `min_weight`
  windows, keeping only recurrent state combinations (the role
  `min_freq` plays in rule extraction). The default `1` keeps
  everything. The total dropped by `min_size` and `min_weight` together
  is recorded in `params$n_dropped`; with both at their defaults,
  `sum(window_counts)` equals the number of non-empty full windows.

## Value

A `net_hg` object (as from
[`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md)
and
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)):
a list with `hyperedges` (list of sorted node index vectors),
`incidence` (numeric node x hyperedge matrix of within-window occurrence
totals), `nodes`, `n_nodes`, `n_hyperedges`, `window_counts` (integer,
one weight per hyperedge: the number of windows collapsed into it),
`size_distribution`, and `params` (`source = "window_hypergraph"`,
`window`, `step`, `min_size`, `n_sequences`, `n_short_sequences`,
`n_windows`, `n_empty_windows`, `n_dropped`). Use
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
for the tidy one-row-per-hyperedge table.

## Details

`step = 1` (default) gives sliding windows; `step = window` gives
tumbling (non-overlapping) windows. Only full windows are formed: a
sequence shorter than `window` contributes none and is counted in
`params$n_short_sequences`. `NA` states are excluded from a window's
set; a window containing only `NA`s is skipped and counted in
`params$n_empty_windows`.

Whole-sequence hyperedges (each complete sequence as one hyperedge) are
the special case already covered by
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
on long-format data with `member = action` and `group = actor`; use this
verb when the hyperedges should be local in time.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification.
*Proceedings of EMNLP 2020*, 4927-4936.
[doi:10.18653/v1/2020.emnlp-main.399](https://doi.org/10.18653/v1/2020.emnlp-main.399)

Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
edge-dependent vertex weights. *Proceedings of the 36th International
Conference on Machine Learning*, PMLR 97, 1172-1181.

## See also

[`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
[`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md)

## Examples

``` r
hg <- window_hypergraph(human_long, action = "code",
                        actor = "session_id", time = "timestamp",
                        window = 3L)
#> Error in Nestimate::prepare(data[columns], actor = actor, action = action,     time = time, session = session, time_threshold = time_threshold,     timezone = timezone): unused argument (timezone = timezone)
hg
#> Error: object 'hg' not found
edges <- hg_get(hg)
#> Error: object 'hg' not found
head(edges)
#> Error: object 'edges' not found

# Tumbling windows over wide-format sequences
wide <- data.frame(
  t1 = c("plan", "code"), t2 = c("code", "test"),
  t3 = c("test", "code"), t4 = c("plan", "debug")
)
window_hypergraph(wide, window = 2L, step = 2L)
#> Hypergraph: 4 nodes, 4 hyperedges (sizes 2: 4)
#> Source: windowed sequences, window = 2, step = 2 (4 windows from 2 sequences)
#>  hyperedge size     members weight
#>         h1    2 code, debug      1
#>         h2    2  code, plan      1
#>         h3    2  code, test      1
#>         h4    2  plan, test      1
```
