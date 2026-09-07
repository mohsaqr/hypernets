# Print Q-analysis results

Print Q-analysis results

## Usage

``` r
# S3 method for class 'net_q_analysis'
print(x, ...)
```

## Arguments

- x:

  A `net_q_analysis` object.

- ...:

  Additional arguments (unused).

## Value

The input object, invisibly.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
qa <- q_analysis(sc)
print(qa)
#> Q-Analysis (max q = 2)
#>   Components: q2:1 q1:1 q0:1
#>   Fully connected at all q levels
#>   Structure: A:2 B:2 C:2
```
