# Degree assortativity of a hypergraph

Measures whether nodes of similar hyperdegree share hyperedges, with the
three choice functions of Chodrow (2020, Sec. 4.2). From every hyperedge
with at least two members a pair of nodes is chosen, and the coefficient
is the correlation of their scores over hyperedges, each hyperedge
counting once:

## Usage

``` r
hg_assortativity(hg, type = "uniform", scale = c("rank", "degree"))
```

## Arguments

- hg:

  A hypernets `net_hg`, dense or sparse.

- type:

  One or more of `"uniform"` (default), `"top_2"`, `"top_bottom"`.

- scale:

  `"rank"` (default; Chodrow's generalized Spearman) or `"degree"`
  (Pearson on raw hyperdegrees).

## Value

A base `data.frame`, one row per `type`, with columns `type`, `scale`,
`assortativity` (in `[-1, 1]`, `NA` when undefined) and `n_edges` (the
hyperedges of size at least two that entered).

## Details

- `"uniform"` – a uniformly random pair of distinct members (computed
  exactly, as the expectation over the random choice; the pair is
  unordered, so both orders enter with equal weight);

- `"top_2"` – the two members of largest degree (unordered);

- `"top_bottom"` – the ordered pair (largest-degree member,
  smallest-degree member), measuring whether the least-connected member
  of a group is better connected when its best-connected member is.

The score of a node is its hyperdegree (the number of hyperedges
containing it). With `scale = "rank"` (default) the score is the node's
rank by hyperdegree among all nodes (ties averaged), which gives
Chodrow's "generalized Spearman" coefficient; `scale = "degree"` uses
the raw hyperdegree (a Pearson coefficient, as XGI does).

On a graph (every hyperedge of size two) the three choice functions give
the same pairs, and `type = "uniform", scale = "degree"` is Newman's
degree assortativity (tested against
[`igraph::assortativity_degree()`](https://r.igraph.org/reference/assortativity.html)).

**Convention differences with XGI 0.10.2**
(`degree_assortativity(exact = TRUE)`), documented and tested: XGI's
exact `"uniform"` weights every hyperedge by its number of member pairs
\\m(m-1)\\ instead of counting each hyperedge once, so the two agree
only on uniform hypergraphs (XGI's Monte Carlo sampler, `exact = FALSE`,
samples Chodrow's edge-weighted quantity and agrees within sampling
error); XGI symmetrizes `"top-bottom"` (both orders), whereas Chodrow's
definition keeps the order (largest, smallest); XGI returns 0 where the
coefficient is undefined. `"top_2"` agrees with XGI exactly.

The coefficient is undefined (`NA`) when either chosen score has zero
variance across hyperedges, e.g. on a regular hypergraph.

**Zero is not the null value for `"top_2"` and `"top_bottom"`.** Both
pair order statistics of a hyperedge (its two largest scores, or its
largest and smallest), which are correlated even when members are drawn
at random, so random hypergraphs give positive values (Chodrow 2020
compares against a null model for this reason). Judge them against
coefficients on degree-matched random hypergraphs, not against 0;
`"uniform"` is centred near 0 under random mixing.

## Conditions

Raises `hypernets_bad_input` when `hg` is not a `net_hg` or has no
hyperedge with at least two members.

## References

Chodrow, P. S. (2020). Configuration models of random hypergraphs.
*Journal of Complex Networks*, 8(3), cnaa018.
[doi:10.1093/comnet/cnaa018](https://doi.org/10.1093/comnet/cnaa018)

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
[doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004)

## See also

[`hg_degree_correlation()`](https://mohsaqr.github.io/hypernets/reference/hg_degree_correlation.md),
[`hg_transitivity()`](https://mohsaqr.github.io/hypernets/reference/hg_transitivity.md).

## Examples

``` r
groups <- data.frame(
  member = c("a", "b", "c", "a", "d", "b", "d", "e", "c", "d", "e", "f"),
  group  = c("g1", "g1", "g1", "g2", "g2", "g3", "g3", "g3",
             "g4", "g4", "g4", "g4")
)
hg <- group_hypergraph(groups, actor = "member", group = "group")
hg_assortativity(hg, type = c("uniform", "top_2", "top_bottom"))
#>         type scale assortativity n_edges
#> 1    uniform  rank    -0.2934132       4
#> 2      top_2  rank    -0.6000000       4
#> 3 top_bottom  rank    -0.3333333       4
```
