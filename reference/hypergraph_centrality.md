# Hypergraph eigenvector centralities

Computes one or more eigenvector-style centralities on a
[net_hypergraph](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md):
*clique-motif* (CEC), *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC).
Each variant captures influence differently - CEC flattens group
structure via clique expansion, while ZEC and HEC propagate through the
higher-order groups directly.

## Usage

``` r
hypergraph_centrality(
  hg,
  type = c("clique", "Z", "H"),
  max_iter = 1000L,
  tol = 1e-08,
  normalize = TRUE,
  damping = 0.85,
  edge_weights = NULL,
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- hg:

  A `net_hypergraph` (from
  [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md),
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  or
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)).

- type:

  Character vector, any subset of
  `c("clique", "Z", "H", "pagerank", "subhypergraph")`. The default
  computes the three eigenvector variants; request `"pagerank"`
  explicitly.

- max_iter:

  Maximum number of power-iteration steps. Default `1000`.

- tol:

  Convergence tolerance on the L1 change between successive iterates.
  Default `1e-8`.

- normalize:

  Logical. If `TRUE` (default), each returned centrality vector is
  L2-normalized to unit norm (compatible with
  [`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)'s
  scale for type `"clique"`). Does not apply to `"pagerank"`, which
  always sums to 1, or `"subhypergraph"`, which is returned on its
  natural log scale.

- damping:

  Single numeric in (0, 1). PageRank damping factor (probability of
  following the walk rather than teleporting). Default `0.85`. Only used
  by `type = "pagerank"`.

- edge_weights:

  NULL, a single positive number (recycled to every hyperedge, e.g. `1`
  for unit weights), or a positive numeric vector, one per hyperedge.
  Hyperedge weights of the EDVW random walk behind `type = "pagerank"`;
  `NULL` defaults to the window counts for hypergraphs built by
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  else the Hayashi et al. dispersion heuristic (unit weights on a binary
  incidence). Only used by `type = "pagerank"`.

- sort_by:

  `NULL` (node order, default) or one of the requested `type`s - return
  the table ranked by that centrality, largest first, ties broken by
  node name so the order is deterministic.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame, one row per node: `node`, followed by one column per
requested `type`, in the order requested. Compare variants by reading
across the columns; rank by one with `sort_by =`.

## Details

**Clique-motif eigenvector centrality (CEC)**: forms the clique-expanded
pairwise graph \\W\\ where \\W\_{ij} = \|\\e : i, j \in e\\\|\\ and
returns the leading eigenvector of \\W\\. Equivalent to running
[`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)
on
[`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md)
output.

**Z-eigenvector centrality (ZEC)**: solves the linear eigen-equation on
the hyperedge tensor, \$\$\lambda\\ x_i \\=\\ \sum\_{e \ni i}\\
\prod\_{j \in e,\\ j \neq i} x_j,\$\$ via power iteration. Works for
hypergraphs with mixed edge sizes.

**H-eigenvector centrality (HEC)**: solves the power-k-1 eigen-equation,
\$\$\lambda\\ x_i^{k-1} \\=\\ \sum\_{e \ni i}\\ \prod\_{j \in e,\\ j
\neq i} x_j.\$\$ For uniform hypergraphs (all hyperedges of size \\k\\),
this is equivalent to normalizing the ZEC update by the geometric-mean
exponent \\1/(k-1)\\. For mixed sizes, the effective exponent is taken
from the largest hyperedge; expect slightly different rankings from ZEC
in the mixed case.

**Hypergraph PageRank** (`"pagerank"`): the stationary distribution of
the damped EDVW random walk of Chitra & Raphael (2019): from node \\v\\,
pick a hyperedge \\e \ni v\\ with probability proportional to its weight
\\w(e)\\, then a node \\u \in e\\ with probability proportional to its
edge-dependent vertex weight \\\gamma_e(u)\\ (the incidence cell, i.e.
occurrence totals for
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md));
with probability \\1 - damping\\ teleport uniformly. Their collapse
theorem: with edge-*independent* vertex weights (a binary incidence) the
walk is equivalent to PageRank on the weighted clique expansion with
edge weights \\\sum\_{e \ni u,v} w(e)/\delta(e)\\ – the hypergraph adds
information exactly when \\\gamma\\ is edge-dependent. Nodes left in no
hyperedge (possible after `min_weight`/`min_size` filtering) teleport
from every step and receive only teleportation mass. The undamped
stationary distribution of the same walk is the `pi` column reported by
[`hypergraph_cluster()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_cluster.md).

**Subhypergraph centrality** (`"subhypergraph"`): the logarithm of the
diagonal of the matrix exponential of the clique adjacency derived from
binary incidence (so entries count shared hyperedges),
\\\log\[\exp(W)\]\_{ii}\\. It counts closed walks based at each node
with a factorial penalty for length and matches the implementation used
by HypergraphX 1.5 in the legal-hypergraphs analysis.

## Note

The `"clique"`, `"Z"`, and `"H"` variants have an independent XGI
equivalence runner in `local_testing_and_equivalence/`; clique is also
validated against
[`igraph::eigen_centrality`](https://r.igraph.org/reference/eigen_centrality.html)
(cosine ~ 1). `"pagerank"` is checked against
[`igraph::page_rank`](https://r.igraph.org/reference/page_rank.html) on
the collapsed graph of the edge-independent case plus a dense
linear-system solve. Clean-room list-based tensor equations remain in
the unit suite for deterministic validation without Python.

## References

Benson, A. R. (2019). Three hypergraph eigenvector centralities. *SIAM
Journal on Mathematics of Data Science* 1(2), 293-312. arXiv:1807.09644.

Chitra, U., & Raphael, B. J. (2019). Random walks on hypergraphs with
edge-dependent vertex weights. *Proceedings of the 36th International
Conference on Machine Learning*, PMLR 97, 1172-1181.

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *Proceedings of CIKM 2020*,
495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Estrada, E., & Rodriguez-Velazquez, J. A. (2005). Complex networks as
hypergraphs. *arXiv preprint physics/0505137*.

## See also

[`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md),
[`clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/clique_expansion.md),
[`hypergraph_measures()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_measures.md).

## Examples

``` r
df <- data.frame(
  member  = c("A", "B", "C", "A", "B", "D", "C", "D", "E"),
  session = c("S1", "S1", "S1", "S2", "S2", "S3", "S3", "S3", "S3")
)
hg <- group_hypergraph(df, "member", "session")

# One row per node, one column per variant - compare across the columns
hypergraph_centrality(hg)
#>   node    clique            Z         H
#> 1    A 0.5345225 6.666667e-01 0.6147176
#> 2    B 0.5345225 6.666667e-01 0.6147176
#> 3    C 0.5345225 3.333333e-01 0.4216259
#> 4    D 0.2672612 3.983626e-09 0.1823130
#> 5    E 0.2672612 3.983626e-09 0.1823130

# Ranked by one of them
hypergraph_centrality(hg, sort_by = "clique")
#>   node    clique            Z         H
#> 1    A 0.5345225 6.666667e-01 0.6147176
#> 2    B 0.5345225 6.666667e-01 0.6147176
#> 3    C 0.5345225 3.333333e-01 0.4216259
#> 4    D 0.2672612 3.983626e-09 0.1823130
#> 5    E 0.2672612 3.983626e-09 0.1823130

# PageRank of the EDVW random walk
hypergraph_centrality(hg, type = "pagerank")
#>   node  pagerank
#> 1    A 0.2241942
#> 2    B 0.2241942
#> 3    C 0.2498116
#> 4    D 0.1509000
#> 5    E 0.1509000
```
