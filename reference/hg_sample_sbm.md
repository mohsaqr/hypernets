# Sample Stochastic-Block-Model Hypergraphs

Samples a graph stochastic block model and augments every sampled dyad
to a hyperedge. Additional members come from the endpoint blocks;
`impurity` members can then be replaced by nodes outside those blocks.
This is the hypergraph SBM construction used by HyperG.

## Usage

``` r
hg_sample_sbm(
  n = NULL,
  P,
  block_sizes,
  d,
  impurity = 0L,
  variable_size = FALSE,
  absolute_purity = TRUE,
  seed = NULL
)

hypergraph_sample_sbm(
  n = NULL,
  P,
  block_sizes,
  d,
  impurity = 0L,
  variable_size = FALSE,
  absolute_purity = TRUE,
  seed = NULL
)
```

## Arguments

- n:

  Total node count; `NULL` uses `sum(block_sizes)`.

- P:

  Symmetric block-to-block dyad probability matrix.

- block_sizes:

  Positive integer community sizes.

- d:

  Hyperedge size, recycled across sampled dyads. With
  `variable_size = TRUE`, it is instead the Poisson mean and two is
  added.

- impurity:

  Non-negative integer number of augmented members replaced by nodes
  outside the endpoint blocks when possible.

- variable_size:

  Draw sizes as `2 + Poisson(d)`.

- absolute_purity:

  If true, impurity replacements must come from outside the endpoint
  blocks; otherwise any nonmember may be used.

- seed:

  Optional reproducibility seed. The caller's RNG state is restored on
  exit.

## Value

A `net_hypergraph`; `$blocks` records each node's planted block.

## Examples

``` r
P <- matrix(c(.5, .05, .05, .5), 2, 2)
h <- hg_sample_sbm(P = P, block_sizes = c(10, 10), d = 3, seed = 1)
as.data.frame(h, what = "nodes")
#>    node degree block
#> 1    V1      9     1
#> 2    V2      9     1
#> 3    V3     12     1
#> 4    V4      5     1
#> 5    V5      5     1
#> 6    V6      5     1
#> 7    V7      8     1
#> 8    V8      6     1
#> 9    V9      8     1
#> 10  V10      5     1
#> 11  V11      5     2
#> 12  V12      5     2
#> 13  V13      6     2
#> 14  V14     10     2
#> 15  V15     11     2
#> 16  V16     11     2
#> 17  V17      9     2
#> 18  V18      7     2
#> 19  V19     10     2
#> 20  V20     10     2
```
