# Hypergraph modularity of a partition

Scores a partition of the nodes of a hypergraph with the hypergraph
modularity of Kaminski et al. (2019): the weighted share of hyperedges
that a community "owns" (the edge contribution) minus the share expected
under the generalised Chung-Lu null model that preserves weighted
degrees and hyperedge sizes (the degree tax). The hypergraph is read as
a family of node sets: cell values of the incidence are ignored, only
membership and one weight per hyperedge count.

## Usage

``` r
hg_modularity(
  hg,
  partition,
  type = c("linear", "majority", "strict"),
  edge_weights = NULL,
  what = c("score", "communities")
)
```

## Arguments

- hg:

  A `net_hg`.

- partition:

  The partition to score: a data.frame with `node` and a label column
  (`community`, `cluster`, `label` or `predicted`), a named label
  vector, an unnamed vector in node order, or an
  [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  fit (its AMI medoid is scored). Every node must be labelled.

- type:

  Hyperedge weighting \\\omega(d, c)\\: `"linear"` (default),
  `"majority"` or `"strict"`.

- edge_weights:

  Positive hyperedge weights (one per hyperedge, or one value recycled).
  `NULL` uses the window counts of a
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  else unit weights.

- what:

  `"score"` (default; one row) or `"communities"` (one row per
  community; its `modularity` column sums to the score).

## Value

A base data.frame. For `what = "score"`: one row with `type`,
`modularity`, `edge_contribution`, `degree_tax` (both already divided by
`total_weight`, so `modularity = edge_contribution - degree_tax`),
`n_communities`, `n_edges` (non-empty hyperedges) and `total_weight`.
For `what = "communities"`: one row per community with `community`,
`n_nodes`, `volume` (summed weighted degree), `edge_contribution`,
`degree_tax` and `modularity`, sorted by decreasing `modularity`. Raises
`hypernets_bad_input` for a non-`net_hg` input, a partition that misses
a node or has `NA` labels, or a hypergraph with no non-empty hyperedge;
invalid `edge_weights` fail the shared weight check.

## Details

A hyperedge of size \\d\\ whose largest part holds \\c\\ of its nodes
contributes \\w_e\\\omega(d, c)\\ to the community holding that part
when \\c \> d/2\\. The three `type`s are HyperNetX's `wdc` functions:
`"strict"` (\\\omega = 1\\ only when the whole edge lies in one
community, the strict modularity of Kaminski et al. 2019), `"majority"`
(\\\omega = 1\\ when a strict majority does), and `"linear"` (\\\omega =
c/d\\ for a strict majority, the HyperNetX default; Kaminski, Pralat &
Theberge 2020). The degree tax of community \\A\\ is \\\sum_d W_d
\sum\_{c \> d/2} \omega(d, c)\\ \mathrm{Binom}(c; d,
\mathrm{vol}(A)/\mathrm{vol}(V))\\, with \\W_d\\ the total weight of
\\d\\-edges and \\\mathrm{vol}\\ the weighted degree. Both terms are
divided by the total hyperedge weight.

## References

Kaminski, B., Poulin, V., Pralat, P., Szufel, P., & Theberge, F. (2019).
Clustering via hypergraph modularity. *PLoS ONE*, 14(11), e0224307.
[doi:10.1371/journal.pone.0224307](https://doi.org/10.1371/journal.pone.0224307)

Kaminski, B., Pralat, P., & Theberge, F. (2020). Community detection
algorithm using hypergraph modularity. In *Complex Networks & Their
Applications IX*, Studies in Computational Intelligence 943, 152-163.
[doi:10.1007/978-3-030-65347-7_13](https://doi.org/10.1007/978-3-030-65347-7_13)

## See also

[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
with `type = "irmm"` to find a partition.

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z",
             "c", "x"),
  edge = c(rep(paste0("e", 1:4), each = 3), "e5", "e5")
)
h <- group_hypergraph(dat, "member", "edge")
split <- data.frame(node = c("a", "b", "c", "x", "y", "z"),
                    community = c(1, 1, 1, 2, 2, 2))
hg_modularity(h, split)
#>     type modularity edge_contribution degree_tax n_communities n_edges
#> 1 linear        0.1               0.8        0.7             2       5
#>   total_weight
#> 1            5
hg_modularity(h, split, type = "strict")
#>     type modularity edge_contribution degree_tax n_communities n_edges
#> 1 strict        0.5               0.8        0.3             2       5
#>   total_weight
#> 1            5
hg_modularity(h, split, what = "communities")
#>   community n_nodes volume edge_contribution degree_tax modularity
#> 1         1       3      7               0.4       0.35       0.05
#> 2         2       3      7               0.4       0.35       0.05
```
