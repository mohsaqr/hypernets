# Correlation of hyperdegrees across interaction orders

For every pair of hyperedge sizes \\s \< t\\ present in the hypergraph,
the Pearson correlation across nodes between the size-\\s\\ hyperdegree
(the number of size-\\s\\ hyperedges containing the node) and the
size-\\t\\ hyperdegree, as reported by Lotito et al. (2023) for the
high-school contact hypergraph (\\\rho\_{2,3}\\, \\\rho\_{2,4}\\,
\\\rho\_{3,4}\\). A high value means that nodes active in pairs are also
active in larger groups.

## Usage

``` r
hg_degree_correlation(hg)
```

## Arguments

- hg:

  A hypernets `net_hg`, dense or sparse.

## Value

A base `data.frame`, one row per pair of present hyperedge sizes
`size_1 < size_2`, with columns `size_1`, `size_2`, `correlation`
(Pearson, `NA` when undefined) and `n_nodes`.

## Details

Every node of `hg` enters every correlation (with degree 0 at the sizes
it does not take part in). Sizes below two are ignored, as in
hypergraphx. A pair is `NA` when one of the two sequences is constant.
hypergraphx 1.8.0 `degree_correlation()` returns the full
`(max_size - 1)`-square matrix, including sizes with no hyperedge (all
`NaN`); hypernets returns one row per pair of sizes that occur. Tested
against hypergraphx at `1e-12`.

## Conditions

Raises `hypernets_bad_input` when `hg` is not a `net_hg` or has fewer
than two distinct hyperedge sizes of at least two.

## References

Lotito, Q. F., Contisciani, M., De Bacco, C., Di Gaetano, L., Gallo, L.,
Montresor, A., Musciotto, F., Ruggeri, N., & Battiston, F. (2023).
Hypergraphx: a library for higher-order network analysis. *Journal of
Complex Networks*, 11(3), cnad019.
[doi:10.1093/comnet/cnad019](https://doi.org/10.1093/comnet/cnad019)

## See also

[`hg_assortativity()`](https://mohsaqr.github.io/hypernets/reference/hg_assortativity.md).

## Examples

``` r
groups <- data.frame(
  member = c("a", "b", "a", "c", "b", "c", "a", "d",
             "a", "b", "c", "a", "b", "d", "c", "d", "e"),
  group  = c("p1", "p1", "p2", "p2", "p3", "p3", "p4", "p4",
             "t1", "t1", "t1", "t2", "t2", "t2", "t3", "t3", "t3")
)
hg <- group_hypergraph(groups, actor = "member", group = "group")
hg_degree_correlation(hg)
#>   size_1 size_2 correlation n_nodes
#> 1      2      3   0.7844645       5
```
