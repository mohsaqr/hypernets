# Higher-order centralities with first-order projection

Computes centralities on the higher-order topology of a
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
network and, by default, projects them back onto the first-order states,
following Scholtes, Wider & Garas (2016). Because a higher-order node
carries the memory of how a state was reached, these centralities can
rank states differently from the same measures on the first-order
transition network - and that difference is the evidence that memory
matters for flow, not just for prediction.

## Usage

``` r
hon_centrality(
  hon,
  type = c("pagerank", "betweenness", "closeness"),
  project = TRUE,
  projection = c("scaled", "last", "first", "all"),
  damping = 0.85,
  weighted = FALSE,
  max_iter = 1000L,
  tol = 1e-12,
  max_paths = 1e+06,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- hon:

  A `net_hon` object from
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md).

- type:

  Character vector, any subset of
  `c("pagerank", "betweenness", "closeness")`. Default computes all
  three.

- project:

  Logical. `TRUE` (default) returns one row per first-order state;
  `FALSE` returns one row per higher-order node with the ordinary
  directed-graph centralities of the higher-order topology (useful to
  see which *contexts*, not states, dominate).

- projection:

  How a higher-order node's PageRank is distributed over the states of
  its path when `project = TRUE`: `"scaled"` (default, pathpy's
  default - each of the k states receives a k-th of the value, so the
  projected PageRank still sums to 1), `"last"` (all of it to the state
  actually occupied), `"first"`, or `"all"` (every state on the path
  receives the full value; the result no longer sums to 1). Ignored by
  `betweenness` and `closeness`, which are first-order by construction.

- damping:

  PageRank damping factor in (0, 1). Default `0.85`.

- weighted:

  Logical. `FALSE` (default, pathpy's convention) runs PageRank on the
  binary higher-order topology; `TRUE` uses the edge weights of
  `hon$matrix`.

- max_iter, tol:

  Power-iteration controls for PageRank.

- max_paths:

  Integer cap on the shortest paths enumerated for betweenness.
  Exceeding it raises a `hypernets_too_many_paths` error rather than
  silently truncating. Default `1e6`.

- sort_by:

  `NULL` (alphabetical by state or node, the default) or one of the
  requested `type`s - sort the table by that centrality, largest first,
  with ties broken by state/node so the order stays deterministic.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A `data.frame`, one row per first-order state (or per higher-order node
when `project = FALSE`), with the requested centralities as columns:

- state:

  First-order state (when `project = TRUE`).

- node, order:

  Higher-order node and its order (when `project = FALSE`).

- pagerank:

  Projected PageRank of the higher-order topology.

- betweenness:

  Share of first-order shortest paths (implied by the higher-order
  topology) passing through the state.

- closeness:

  Harmonic closeness on the implied first-order distances: the sum of
  inverse distances *into* the state.

Rows are ordered by state (or by node) so the result is deterministic.

## Details

Semantics match pathpy 2.2.0 (the reference implementation by the
method's author), generalized from fixed-order to the variable-order
networks
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
produces: a higher-order node's order is the number of states in the
path it represents, and the first-order distance implied by a
higher-order hop count adds that node's order minus one. With a uniform
order the generalization reduces exactly to pathpy's formulas.

## References

Scholtes, I., Wider, N., & Garas, A. (2016). Higher-order aggregate
networks in the analysis of temporal networks: path structures and
centralities. *The European Physical Journal B* 89, 61.
[doi:10.1140/epjb/e2016-60663-0](https://doi.org/10.1140/epjb/e2016-60663-0)

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances* 2(5),
e1600028.
[doi:10.1126/sciadv.1600028](https://doi.org/10.1126/sciadv.1600028)

## See also

[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md),
[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md),
[`path_dependence()`](https://mohsaqr.github.io/hypernets/reference/path_dependence.md)

## Examples

``` r
seqs <- c(
  replicate(6, rep(c("a", "b", "c"), 4), simplify = FALSE),
  replicate(6, rep(c("x", "b", "d"), 4), simplify = FALSE)
)
hon <- build_hon(seqs, max_order = 2)
hon_centrality(hon)
#>   state  pagerank betweenness closeness
#> 1     a 0.2428728           1       1.5
#> 2     b 0.1813376           2       3.0
#> 3     c 0.1664584           1       1.5
#> 4     d 0.1664584           1       1.5
#> 5     x 0.2428728           1       1.5

# Which contexts, rather than which states, carry the flow?
hon_centrality(hon, type = "pagerank", project = FALSE)
#>     node order   pagerank
#> 1      a     1 0.16291823
#> 2 a -> b     2 0.15990907
#> 3      b     1 0.02142857
#> 4      c     1 0.16645842
#> 5      d     1 0.16645842
#> 6      x     1 0.16291823
#> 7 x -> b     2 0.15990907
```
