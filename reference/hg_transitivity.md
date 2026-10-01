# Local clustering coefficients of a hypergraph

Computes one or more published local clustering coefficients for every
node of a
[net_hg](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md).
Each variant is named after its definition:

## Usage

``` r
hg_transitivity(
  hg,
  type = "projection",
  what = c("nodes", "summary"),
  sort_by = NULL,
  n = Inf
)
```

## Arguments

- hg:

  A hypernets `net_hg`
  ([`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  ...), dense or sparse.

- type:

  One or more of `"projection"` (default), `"extra_overlap"`,
  `"two_node_union"`, `"two_node_min"`, `"two_node_max"`.

- what:

  `"nodes"` (default; one row per node) or `"summary"` (one row per
  `type` with the mean local coefficient over the nodes where it is
  defined, the average local clustering used by Chodrow 2020, and the
  number of nodes where it is undefined).

- sort_by:

  For `what = "nodes"`, optional `type` to sort by, descending (ties
  broken by node name, `NA` last).

- n:

  Return only the first `n` node rows after sorting.

## Value

A base `data.frame`. For `what = "nodes"`: one row per node with column
`node` and one numeric column per requested `type` (in `[0, 1]`, `NA`
where undefined). For `what = "summary"`: one row per `type` with
columns `type`, `mean_local` and `n_undefined`.

## Details

- `"projection"` – the Watts & Strogatz (1998) local clustering of the
  unweighted clique expansion: the share of pairs of a node's neighbours
  that are themselves neighbours. A path through hyperedges is closed by
  *any* hyperedge, so a single hyperedge of size three already closes a
  triangle (the path-based generalization reviewed by Battiston et al.
  2020, Sec. III.C).

- `"extra_overlap"` – Zhou & Nakhleh (2011), in the statement of Klimm,
  Deane & Reinert (2021): for every unordered pair of hyperedges \\e,
  f\\ containing the node, the *extra overlap* \$\$EO(e, f) =
  \frac{\|N(D\_{e,f}) \cap D\_{f,e}\| + \|N(D\_{f,e}) \cap
  D\_{e,f}\|}{\|D\_{e,f}\| + \|D\_{f,e}\|},\qquad D\_{e,f} = e \setminus
  f,\$\$ is the share of the non-shared members that are linked to the
  other side by a *third* hyperedge (\\N(S)\\ is the union of the
  hypergraph neighbourhoods of \\S\\); \\EO = 0\\ when \\e = f\\. The
  node's coefficient is the mean EO over its \\k(k-1)/2\\ pairs of
  hyperedges. It reduces to the ordinary clustering coefficient on a
  graph and does not count the triangles *inside* a hyperedge.

- `"two_node_union"`, `"two_node_min"`, `"two_node_max"` – the two-mode
  clustering of Latapy, Magnien & Del Vecchio (2008), applied to
  hypernetworks by Gallagher & Goldberg (2013): for nodes \\u, v\\ with
  hyperedge sets \\M(u), M(v)\\, \\cc(u, v) = \|M(u) \cap M(v)\| /
  \|M(u) \cup M(v)\|\\ (`union`), or with \\\min\\ / \\\max(\|M(u)\|,
  \|M(v)\|)\\ in the denominator; the node's coefficient averages
  \\cc(u, v)\\ over its hypergraph neighbours \\v\\ (Latapy et al.'s
  \\N(N(u))\\).

**Undefined values are `NA`.** A node with fewer than two projection
neighbours has no `"projection"` coefficient, a node in no hyperedge has
no `"extra_overlap"` coefficient (hyperdegree 1 gives 0, as Klimm et al.
define it), and a node without neighbours has no `"two_node_*"`
coefficient. XGI 0.10.2 reports 0 in all of these cases; the values
agree everywhere else (the oracle test maps `NA` to 0 before comparing).

The global clustering coefficient \\C_2(H)\\ of Estrada &
Rodriguez-Velazquez (2006) is not provided (no independent oracle).

Duplicate hyperedges are distinct hyperedges (as in XGI): they raise the
hyperdegree and, for `"extra_overlap"`, contribute pairs with EO = 0.
`"extra_overlap"` forms a dense edge-by-edge matrix and visits every
pair of intersecting hyperedges, so it is meant for hypergraphs with at
most a few thousand hyperedges; the other variants use sparse products
when the incidence is sparse.

## Conditions

Raises `hypernets_bad_input` when `hg` is not a `net_hg` or `n` is not a
count.

## References

Watts, D. J., & Strogatz, S. H. (1998). Collective dynamics of
'small-world' networks. *Nature*, 393(6684), 440-442.
[doi:10.1038/30918](https://doi.org/10.1038/30918)

Zhou, W., & Nakhleh, L. (2011). Properties of metabolic graphs:
biological organization or representation artifacts? *BMC
Bioinformatics*, 12, 132.
[doi:10.1186/1471-2105-12-132](https://doi.org/10.1186/1471-2105-12-132)

Klimm, F., Deane, C. M., & Reinert, G. (2021). Hypergraphs for
predicting essential genes using multiprotein complex data. *Journal of
Complex Networks*, 9(2), cnaa028.
[doi:10.1093/comnet/cnaa028](https://doi.org/10.1093/comnet/cnaa028)

Latapy, M., Magnien, C., & Del Vecchio, N. (2008). Basic notions for the
analysis of large two-mode networks. *Social Networks*, 30(1), 31-48.
[doi:10.1016/j.socnet.2007.04.006](https://doi.org/10.1016/j.socnet.2007.04.006)

Gallagher, S. R., & Goldberg, D. S. (2013). Clustering coefficients in
protein interaction hypernetworks. *Proceedings of the ACM Conference on
Bioinformatics, Computational Biology and Biomedical Informatics (BCB
'13)*, 552-560.
[doi:10.1145/2506583.2506635](https://doi.org/10.1145/2506583.2506635)

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1-92.
[doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004)

Chodrow, P. S. (2020). Configuration models of random hypergraphs.
*Journal of Complex Networks*, 8(3), cnaa018.
[doi:10.1093/comnet/cnaa018](https://doi.org/10.1093/comnet/cnaa018)

## See also

[`hg_assortativity()`](https://mohsaqr.github.io/hypernets/reference/hg_assortativity.md),
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
[`hg_clique_expansion()`](https://mohsaqr.github.io/hypernets/reference/hg_clique_expansion.md).

## Examples

``` r
groups <- data.frame(
  member = c("a", "b", "c", "a", "d", "b", "d", "c", "d", "e"),
  group  = c("g1", "g1", "g1", "g2", "g2", "g3", "g3", "g4", "g4", "g4")
)
hg <- group_hypergraph(groups, actor = "member", group = "group")
hg_transitivity(hg, type = c("projection", "extra_overlap"))
#>   node projection extra_overlap
#> 1    a  1.0000000     1.0000000
#> 2    b  1.0000000     1.0000000
#> 3    c  0.6666667     0.7500000
#> 4    d  0.6666667     0.7777778
#> 5    e  1.0000000     0.0000000
hg_transitivity(hg, type = "two_node_union", what = "summary")
#>             type mean_local n_undefined
#> 1 two_node_union  0.3305556           0
```
