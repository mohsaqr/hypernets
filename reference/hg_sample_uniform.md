# Sample Uniform or Regular Random Hypergraphs

`hg_sample_uniform()` makes every hyperedge contain exactly `k` distinct
nodes. `hg_sample_regular()` makes every node incident to exactly `k`
distinct hyperedges. Sampling probabilities can be unequal.

## Usage

``` r
hg_sample_uniform(n, m, k, prob = NULL, seed = NULL)

hg_sample_regular(n, m, k, prob = NULL, seed = NULL)

hypergraph_sample_uniform(n, m, k, prob = NULL, seed = NULL)

hypergraph_sample_regular(n, m, k, prob = NULL, seed = NULL)
```

## Arguments

- n:

  Number of nodes.

- m:

  Number of hyperedges. For
  [`hg_sample_gnp()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_gnp.md),
  `NULL` draws it from a Poisson distribution. For the uniform and
  regular models it is required.

- k:

  Hyperedge size for the uniform model; node degree for the regular
  model.

- prob:

  Sampling weights: length `n` for the uniform model and length `m` for
  the regular model. `NULL` is uniform.

- seed:

  Optional reproducibility seed. The caller's RNG state is restored on
  exit.

## Value

A `net_hypergraph`.

## Examples

``` r
u <- hg_sample_uniform(20, 8, k = 3, seed = 1)
r <- hg_sample_regular(20, 8, k = 2, seed = 1)
hg_measures(u, what = "distribution", measure = "size")
#>   value n proportion ccdf
#> 1     3 8          1    1
hg_measures(r, what = "distribution", measure = "hyperdegree")
#>   value  n proportion ccdf
#> 1     2 20          1    1
```
