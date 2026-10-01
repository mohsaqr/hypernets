# Communities of a hypergraph or a memory network

One verb, two estimators, chosen by the class of `x`:

- a hypergraph (`net_hg`):

  an ensemble of Infomap runs on a projection (or IRMM with
  `type = "irmm"`), with the adjusted-mutual- information medoid; see
  [`hg_communities.net_hg()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.net_hg.md).

- a memory network (`net_hon`, from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md)):

  the map equation for memory networks: overlapping modules of the
  physical states, compared with the first-order map; see
  [`hg_communities.net_hon()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.net_hon.md).

Each method keeps its own arguments; passing an argument that only the
other method takes raises `hypernets_bad_input`.

## Usage

``` r
hg_communities(x, ...)

# Default S3 method
hg_communities(x, ...)
```

## Arguments

- x:

  A `net_hg` or a `net_hon`.

- ...:

  Arguments of the method for `class(x)`.

## Value

An `hg_communities` object (hypergraph) or a `net_hon_communities`
object (memory network); read either with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md).
Any other input raises `hypernets_bad_input`.

## Examples

``` r
seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
             c("c", "h", "d", "c", "h", "d", "c"))
comm <- hg_communities(hon(seqs, max_order = 2L), trials = 2L)
hg_get(comm)
#>   node state community      flow
#> 1    a     a         1 0.1669582
#> 2    b     b         1 0.1611273
#> 3    c     c         1 0.1669582
#> 4    d     d         1 0.1611273
#> 5    h     h         1 0.3438290
```
