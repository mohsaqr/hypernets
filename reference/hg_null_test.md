# Degree-preserving null-model test of hypergraph structure

Tests whether an observed structural statistic exceeds what the degree
sequences alone imply. The null model fixes both margins of the binary
membership (every vertex keeps its hyperdegree, every hyperedge its
size) and randomizes the memberships by checkerboard swaps (Gotelli
2000) – a sequential MCMC with burn-in `10 * nnz` swap attempts and
thinning `nnz` between samples, `nnz` being the number of memberships.
Statistics are evaluated on the binarized hypergraph (weights carry no
meaning under this null), through the same delegated measures as
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md).

## Usage

``` r
hg_null_test(
  hg,
  statistic = c("pairwise_participation", "density", "avg_edge_size", "avg_jaccard",
    "repeated_edges", "repeated_pairs"),
  method = c("swap", "configuration", "assignment"),
  n = 199L,
  seed = NULL,
  alternative = c("two_sided", "greater", "less")
)

hypergraph_null_test(
  hg,
  statistic = c("pairwise_participation", "density", "avg_edge_size", "avg_jaccard",
    "repeated_edges", "repeated_pairs"),
  method = c("swap", "configuration", "assignment"),
  n = 199L,
  seed = NULL,
  alternative = c("two_sided", "greater", "less")
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md),
  or any honets `net_hypergraph`.

- statistic:

  Statistics to test; any of `"pairwise_participation"`, `"density"`,
  `"avg_edge_size"`, `"avg_jaccard"` (mean pairwise edge Jaccard),
  `"repeated_edges"` (hyperedges whose member set already occurred: the
  paper's repeated collaboration trios) and `"repeated_pairs"` (member
  pairs counted with multiplicity minus the distinct co-occurring pairs:
  its repeated collaboration pairs). Several allowed.

- method:

  Which null. `"swap"` (default) is the checkerboard swap chain above:
  both margins are preserved *exactly*, and no vertex may repeat within
  a hyperedge. `"configuration"` draws independent stub-matchings
  (Chodrow 2020), the standard higher-order configuration model: vertex
  stubs are matched uniformly at random to hyperedge slots, so margins
  hold only up to collapse (a vertex whose stubs land twice in one
  hyperedge loses a degree). `"assignment"` is the degree-ordered
  assignment Coupette et al. (2024) use for their repeated-collaboration
  test: vertices in decreasing degree each choose their hyperedges
  without replacement among those with a free slot, so both margins hold
  exactly and no vertex repeats within a hyperedge, at the price of a
  draw that is not uniform over configurations. `"swap"` is the stricter
  null; use `"configuration"` to compare against the higher-order
  network literature, which reports it.

- n:

  Number of null samples (default `199L`).

- seed:

  Seed for the null draws; set it for a reproducible test. The global
  RNG state is restored on exit.

- alternative:

  `"two_sided"` (default), `"greater"`, or `"less"`.

## Value

A base `data.frame`, one row per statistic: `statistic`, `observed`,
`null_mean`, `null_lo`, `null_hi`, `z`, `p_value`, `n`, `method`.

## Details

The permutation p-value is `(1 + extreme) / (n + 1)` (Phipson & Smyth
2010), two-sided by default around the null mean; `null_lo`/`null_hi`
are the 2.5% and 97.5% null quantiles – report them with the observed
value, not the p-value alone.

## Conditions

Raises `honets_bad_input` for broken contracts.
`method = "configuration"` signals a `honets_configuration_collapse`
warning when stub matching collapses more than 1% of memberships on
average, which happens whenever hyperedges are large relative to the
vertex set – the usual case in the document orientation.

## References

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. *Ecology*, 81(9).

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero. *Statistical Applications in Genetics and Molecular Biology*,
9(1).

Chodrow, P. S. (2020). Configuration models of random hypergraphs.
*Journal of Complex Networks*, 8(3), cnaa018.
[doi:10.1093/comnet/cnaa018](https://doi.org/10.1093/comnet/cnaa018)

## Examples

``` r
hg <- text_hypergraph(c(
  a = "salt and soup and onions",
  b = "soup and salt",
  c = "stars and sky and salt",
  d = "stars and sky"
))
hg_null_test(hg, statistic = "pairwise_participation", n = 49, seed = 1)
#>                statistic observed null_mean null_lo null_hi  z p_value  n
#> 1 pairwise_participation        1         1       1       1 NA       1 49
#>   method
#> 1   swap
```
