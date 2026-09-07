# Euler Characteristic

Computes \\\chi = \sum\_{k=0}^{d} (-1)^k f_k\\ where \\f_k\\ is the
number of k-simplices. By the Euler-Poincare theorem, \\\chi = \sum\_{k}
(-1)^k \beta_k\\.

## Usage

``` r
euler_characteristic(sc)
```

## Arguments

- sc:

  A `net_simplicial` object.

## Value

Integer.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
euler_characteristic(sc)
#> [1] 1
```
