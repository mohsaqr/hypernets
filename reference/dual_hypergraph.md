# Dual of a hypergraph

Returns the dual hypergraph: every hyperedge becomes a vertex and every
vertex becomes a hyperedge, with the transposed weighted incidence. For
a bag-construction
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
the dual is identical to rebuilding with the opposite `nodes`
orientation (tested), so document-level and word-level analyses can
share one constructed object. Duals of windowed and kNN hypergraphs are
returned as plain `net_hypergraph` objects (their corpus bookkeeping
does not transpose meaningfully).

## Usage

``` r
dual_hypergraph(hg)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md),
  or any hypernets `net_hypergraph`.

## Value

A hypergraph whose incidence is the transpose of `hg`'s: a
`text_hypergraph` with flipped `nodes` for bag constructions, otherwise
a `net_hypergraph`. Accepted by all `hg_*` verbs.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
dual <- dual_hypergraph(hg)
dual
#> Text hypergraph: 2 documents, 4 words (words as nodes, weight = n)
#> Hyperedges: 2 (documents); sizes 3-3, median 3
hg_measures(dual, what = "edges")
#>   edge size
#> 1    a    3
#> 2    b    3
```
