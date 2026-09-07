# Summary Method for net_hon

Summary Method for net_hon

## Usage

``` r
# S3 method for class 'net_hon'
summary(object, ...)
```

## Arguments

- object:

  A `net_hon` object.

- ...:

  Additional arguments (ignored).

## Value

The edge data.frame `object$edges` (columns `path`, `from`, `to`,
`count`, `probability`, `from_order`, `to_order`), returned visibly; the
summary text is printed as a side effect.

Returned **invisibly**: `summary(x)` prints the summary and nothing
else; assign the result to keep the table.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon <- build_hon(seqs, max_order = 2)
summary(hon)
#> Higher-Order Network (HON) Summary
#>   Nodes: 4 | Edges: 5 | Trajectories: 3
#>   First-order states: A, B, C, D
#>   Max order observed: 1 (requested: 2)
#>   Min frequency: 1
#>   Node order distribution:
#>     Order 1: 4 nodes

# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B"),
  V2 = c("B","C","A","B","C"),
  V3 = c("C","A","B","C","A")
)
hon <- build_hon(seqs, max_order = 2L)
summary(hon)
#> Higher-Order Network (HON) Summary
#>   Nodes: 3 | Edges: 3 | Trajectories: 5
#>   First-order states: A, B, C
#>   Max order observed: 1 (requested: 2)
#>   Min frequency: 1
#>   Node order distribution:
#>     Order 1: 3 nodes
# }
```
