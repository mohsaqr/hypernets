# Q-Analysis

Computes Q-connectivity structure (Atkin 1974). Two maximal simplices
are q-connected if they share a face of dimension \\\geq q\\. Reports:

- **Q-vector**: number of connected components at each q-level

- **Structure vector**: highest simplex dimension per node

## Usage

``` r
q_analysis(sc)
```

## Arguments

- sc:

  A `net_simplicial` object.

## Value

A `net_q_analysis` object with `$q_vector`, `$structure_vector`, and
`$max_q`.

## References

Atkin, R. H. (1974). *Mathematical Structure in Human Affairs*.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
q_analysis(sc)
#> Q-Analysis (max q = 2)
#>   Components: q2:1 q1:1 q0:1
#>   Fully connected at all q levels
#>   Structure: A:2 B:2 C:2
```
