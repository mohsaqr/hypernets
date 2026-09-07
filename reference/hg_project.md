# Project a hypergraph onto a weighted graph

Collapses every hyperedge into pairwise vertex relations, giving the
weighted graph that a graph engine (paths, betweenness, communities) can
consume. Two weightings are available. `"clique"` is plain
co-occurrence, the clique expansion: a hyperedge of size \\\|e\|\\
contributes the same amount to each of its \\\|e\|(\|e\|-1)/2\\ pairs,
so a single large hyperedge can dominate every downstream measure.
`"association"` divides each hyperedge's contribution by \\\|e\|-1\\,
following Coupette et al. (2024):

## Usage

``` r
hg_project(
  hg,
  method = c("clique", "association", "citation"),
  weighted = TRUE,
  what = c("edges", "matrix"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL,
  directed = FALSE
)

hypergraph_project(
  hg,
  method = c("clique", "association", "citation"),
  weighted = TRUE,
  what = c("edges", "matrix"),
  duplicate_edges = c("count", "collapse"),
  self_association = FALSE,
  edge_source = NULL,
  directed = FALSE
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md),
  or any hypernets `net_hypergraph`.

- method:

  Weighting. `"clique"` (default) sums incidence products;
  `"association"` applies the \\1/(\|e\|-1)\\ normalisation above;
  `"citation"` joins each hyperedge's source to its members.

- weighted:

  `method = "clique"` only. `TRUE` (default) uses the incidence weights,
  `FALSE` their membership pattern. Setting it together with
  `method = "association"` is an error, because the association
  weighting is defined on hyperedge cardinality and never on the
  incidence weights.

- what:

  `"edges"` (default) for the tidy edge list, or `"matrix"` for the
  symmetric weight matrix to hand to a graph engine.

- duplicate_edges:

  For `method = "association"`, `"count"` (default) lets repeated
  hyperedges contribute repeatedly (the paper's multi-hypergraph
  representation); `"collapse"` lets each distinct member set contribute
  once (its binary-hypergraph representation).

- self_association:

  Add the paper's source-to-member association term? Requires one source
  vertex per hyperedge through `edge_source` or `hg$edge_data$source`.
  For source `u`, every membership occurrence of target `v` contributes
  `1 / sum_e |e|` over hyperedges sourced by `u`, so the added incident
  weight from all of `u`'s citations sums to one.

- edge_source:

  Which node each hyperedge comes from: the name of a column of the
  hyperedge attributes (`hg$edge_data`, e.g. `"citing"`), a vector of
  length `n_hyperedges`, a named vector keyed by hyperedge, or a
  two-column data frame named `edge` and `source`. `NULL` uses an
  attribute column named `source`. Used with `self_association = TRUE`
  and `method = "citation"`.

- directed:

  For `method = "citation"`: return the directed source-to-member matrix
  (rows cite columns)? Default `FALSE`.

## Value

With `what = "edges"`, a base data.frame with one row per unordered
vertex pair of non-zero weight, sorted by `from` then `to`, with columns
`from`, `to` (vertex names) and `weight` (numeric). Vertices sharing no
hyperedge do not appear, so an edgeless hypergraph yields a zero-row
data.frame with those columns. With `what = "matrix"`, the symmetric
`n_nodes` x `n_nodes` weight matrix with zero diagonal and vertex names
as dimnames, sparse if `hg`'s incidence is sparse.

## Details

\$\$w(\\u,v\\) = \sum\_{e \supseteq \\u,v\\} \frac{1}{\|e\| - 1}\$\$

so the total weight a hyperedge adds around any one of its members is
exactly 1, and each vertex's weighted degree in the projection equals
the number of hyperedges of size at least two that contain it. That
normalisation is what makes the projection a legitimate random-walk
operator rather than an arbitrary co-occurrence count.

A hyperedge of size one has no pairs, so it contributes nothing and is
skipped rather than dividing by zero. Such hyperedges are common in
text: a word used in exactly one document is a singleton hyperedge in
the document orientation, which is why the degree identity above counts
only hyperedges of size at least two.

A third projection, `"citation"`, is the classic graph a hypergraph with
sources reduces to: one edge from the source of every hyperedge to each
of its members, so a citation-block hypergraph becomes the ordinary
citation graph of the paper's Table 2. `duplicate_edges = "count"`
weights an edge by how many blocks repeat it (the multi-graph, `mg`);
`"collapse"` keeps the binary graph (`bg`). `directed = TRUE` keeps the
source-to-member orientation; the default symmetrises.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

Zhou, D., Huang, J., & Schoelkopf, B. (2006). Learning with hypergraphs:
clustering, classification, and embedding. *NeurIPS 19*, 1601-1608.

## See also

[`hg_line_graph()`](https://mohsaqr.github.io/hypernets/reference/hg_line_graph.md)
for the dual projection, onto hyperedges.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars",
                        c = "stars and salt"))
hg_project(hg)
#>   from to weight
#> 1    a  b      2
#> 2    a  c      2
#> 3    b  c      2
hg_project(hg, method = "association")
#>   from to weight
#> 1    a  b    1.5
#> 2    a  c    1.5
#> 3    b  c    1.5
```
