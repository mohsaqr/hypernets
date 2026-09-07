# Simplicial Degree

Counts how many simplices of each dimension contain each node.

## Usage

``` r
simplicial_degree(sc, normalized = FALSE, top = NULL)
```

## Arguments

- sc:

  A `net_simplicial` object.

- normalized:

  Divide by maximum possible count. Default `FALSE`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

Data frame with `node`, columns `d0` through `d_k`, and `total` (sum of
d1+). Sorted by total descending.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
simplicial_degree(sc)
#>   node d0 d1 d2 total
#> 1    A  1  2  1     3
#> 2    B  1  2  1     3
#> 3    C  1  2  1     3
```
