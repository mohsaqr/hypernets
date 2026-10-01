# Hypergraph eigenvector centralities

The hypergraph method of
[`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md).
Computes one or more eigenvector-style centralities on a
[net_hg](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md):
*clique-motif* (CEC), *Z-eigenvector* (ZEC), and *H-eigenvector* (HEC).
Each variant captures influence differently - CEC flattens group
structure via clique expansion, while ZEC and HEC propagate through the
higher-order groups directly.

## Usage

``` r
# S3 method for class 'net_hg'
hg_centrality(
  x,
  type = c("clique", "Z", "H"),
  sort_by = NULL,
  n = Inf,
  max_iter = 1000L,
  tol = 1e-08,
  normalize = TRUE,
  damping = 0.85,
  edge_weights = NULL,
  alpha = NULL,
  ...
)
```

## Arguments

- x:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hg`).

- type:

  Centralities to compute; any of `"clique"`, `"Z"`, `"H"` (default: all
  three), `"pagerank"`, `"subhypergraph"`, `"katz"`.

- sort_by:

  Optional centrality name to sort by, descending (ties broken by node
  name); default keeps node order.

- n:

  Return only the first `n` rows after sorting (default all) – e.g.
  `sort_by = "clique", n = 10` for the ten most central nodes.

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

- alpha:

  Attenuation factor of `type = "katz"`: a single number in \\(0,
  1/\lambda\_{max})\\, where \\\lambda\_{max}\\ is the largest
  eigenvalue of the hypergraph adjacency. There is no correct default;
  `type = "katz"` without `alpha` raises `hypernets_bad_input`, and so
  does an `alpha` at or above the bound (the message reports it).

- ...:

  Must be empty: an argument that only the memory-network method takes
  raises `hypernets_bad_input`.

## Value

A base `data.frame`, one row per node (or the `n` requested rows), with
one column per requested centrality.

## Details

`type` also accepts `"pagerank"` (EDVW hypergraph PageRank),
`"subhypergraph"` and `"katz"` (Katz centrality, which needs `alpha`);
the result is one row per node.

**Clique-motif eigenvector centrality (CEC)**: forms the clique-expanded
pairwise graph \\W\\ where \\W\_{ij} = \|\\e : i, j \in e\\\|\\ and
returns the leading eigenvector of \\W\\. Equivalent to running
[`igraph::eigen_centrality()`](https://r.igraph.org/reference/eigen_centrality.html)
on
[`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md)
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
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md).

**Subhypergraph centrality** (`"subhypergraph"`): the logarithm of the
diagonal of the matrix exponential of the clique adjacency derived from
binary incidence (so entries count shared hyperedges),
\\\log\[\exp(W)\]\_{ii}\\. It counts closed walks based at each node
with a factorial penalty for length and matches the implementation used
by HypergraphX 1.5 in the legal-hypergraphs analysis.

**Katz centrality** (`"katz"`): Katz's (1953) status index on the
hypergraph adjacency of Estrada & Rodriguez-Velazquez (2006),
\\A\_{ij}\\ = the number of hyperedges containing both \\i\\ and \\j\\
(the weighted clique expansion), \$\$x = \sum\_{k \ge 1} \alpha^k A^k
\mathbf{1} = (I - \alpha A)^{-1} \alpha A \mathbf{1},\qquad 0 \< \alpha
\< 1/\lambda\_{max}(A),\$\$ the \\\alpha\\-discounted number of walks
leaving each node (Battiston et al. 2020, Sec. III.B.3). Returned on
this natural scale, not normalized; the alpha centrality of Bonacich &
Lloyd (2001) with unit exogenous status, which also counts the length-0
walk, is `katz + 1` (tested against
[`igraph::alpha_centrality()`](https://r.igraph.org/reference/alpha_centrality.html)
on the clique expansion). XGI's `katz_centrality()` is the same vector
for its fixed \\\alpha = 2^{-n}\\, rescaled to sum to one (tested). It
is the only type that also runs on a sparse incidence.

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

Estrada, E., & Rodriguez-Velazquez, J. A. (2006). Subgraph centrality
and clustering in complex hyper-networks. *Physica A*, 364, 581-594.
[doi:10.1016/j.physa.2005.12.002](https://doi.org/10.1016/j.physa.2005.12.002)

Katz, L. (1953). A new status index derived from sociometric analysis.
*Psychometrika*, 18(1), 39-43.
[doi:10.1007/BF02289026](https://doi.org/10.1007/BF02289026)

Bonacich, P., & Lloyd, P. (2001). Eigenvector-like measures of
centrality for asymmetric relations. *Social Networks*, 23(3), 191-201.
[doi:10.1016/S0378-8733(01)00038-7](https://doi.org/10.1016/S0378-8733%2801%2900038-7)

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
[doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004)

## See also

[`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
[`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md),
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md).

## Examples

``` r
hg <- text_hypergraph(c(
  a = "salt and soup and onions",
  b = "soup and salt",
  c = "stars and salt"
))
hg_centrality(hg, type = "clique")
#>   node    clique
#> 1    a 0.6059128
#> 2    b 0.6059128
#> 3    c 0.5154991
hg_centrality(hg, type = c("clique", "katz"), alpha = 0.1)
#>   node    clique      katz
#> 1    a 0.6059128 0.9354839
#> 2    b 0.6059128 0.9354839
#> 3    c 0.5154991 0.7741935
```
