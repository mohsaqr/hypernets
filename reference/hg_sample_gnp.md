# Sample Bernoulli-Incidence Random Hypergraphs

Generates an Erdos–Renyi-style random hypergraph with `n` nodes and `m`
hyperedges by drawing every incidence independently as Bernoulli(`p`).
If `m` is omitted it is Poisson with mean `lambda`, or `n * p` when
`lambda` is also omitted. Empty and singleton hyperedges are retained
because they are valid outcomes of this incidence model.

## Usage

``` r
hg_sample_gnp(n, m = NULL, p, lambda = NULL, seed = NULL)

hypergraph_sample_gnp(n, m = NULL, p, lambda = NULL, seed = NULL)
```

## Arguments

- n:

  Number of nodes.

- m:

  Number of hyperedges. For `hg_sample_gnp()`, `NULL` draws it from a
  Poisson distribution. For the uniform and regular models it is
  required.

- p:

  Incidence probability in `[0, 1]`.

- lambda:

  Optional Poisson mean for `m`.

- seed:

  Optional reproducibility seed. The caller's RNG state is restored on
  exit.

## Value

A `net_hypergraph` with binary incidence and model parameters in
`$params`.

## References

Marchette, D. J. (2021). HyperG: Hypergraphs in R. R package version
1.0.0.

## Examples

``` r
h <- hg_sample_gnp(n = 20, m = 8, p = 0.2, seed = 1)
summary(h)
#> Hypergraph summary
#>   Nodes:         20
#>   Hyperedges:    8
#>   Mean size:     3.00
#>   Max size:      5
```
