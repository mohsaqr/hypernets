# Centralities of a hypergraph or a memory network

One verb, two estimators, chosen by the class of `x`:

- a hypergraph (`net_hg`, e.g. from
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)):

  eigenvector-style hypergraph centralities – clique-motif, Z- and
  H-eigenvector, EDVW PageRank, subhypergraph and Katz; see
  [`hg_centrality.net_hg()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.net_hg.md).

- a memory network (`net_hon`, from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md)):

  PageRank, betweenness and closeness of the higher-order topology,
  projected onto the first-order states; see
  [`hg_centrality.net_hon()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.net_hon.md).

Each method keeps its own arguments; passing an argument that only the
other method takes raises `hypernets_bad_input`.

## Usage

``` r
hg_centrality(x, ...)

# Default S3 method
hg_centrality(x, ...)
```

## Arguments

- x:

  A `net_hg` or a `net_hon`.

- ...:

  Arguments of the method for `class(x)`.

## Value

A base `data.frame`, one row per node (or per state), with one column
per requested centrality. Any other input raises `hypernets_bad_input`.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
hg_centrality(hg, type = "clique")
#>   node    clique
#> 1    a 0.7071068
#> 2    b 0.7071068

seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
hg_centrality(hon(seqs, max_order = 2), type = "pagerank")
#>   state  pagerank
#> 1     a 0.1719145
#> 2     b 0.3222546
#> 3     c 0.1669582
#> 4     d 0.1669582
#> 5     x 0.1719145
```
