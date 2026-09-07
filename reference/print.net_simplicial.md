# Print a simplicial complex

Print a simplicial complex

## Usage

``` r
# S3 method for class 'net_simplicial'
print(x, ...)
```

## Arguments

- x:

  A `net_simplicial` object.

- ...:

  Additional arguments (unused).

## Value

The input object, invisibly.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
print(sc)
#> Clique Complex 
#>   3 nodes, 7 simplices, dimension 2
#>   Density: 100.0%  |  Mean dim: 0.71  |  Euler: 1
#>   f-vector: (f0=3 f1=3 f2=1)
#>   Betti: b0=1
#>   Nodes: A, B, C 
```
