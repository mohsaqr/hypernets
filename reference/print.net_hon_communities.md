# Print a memory-network community result

Print a memory-network community result

## Usage

``` r
# S3 method for class 'net_hon_communities'
print(x, n = 10L, ...)
```

## Arguments

- x:

  A `net_hon_communities` object.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  Ignored.

## Value

`x`, invisibly.

## Examples

``` r
seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
             c("c", "h", "d", "c", "h", "d", "c"))
print(hg_communities(hon(seqs, max_order = 2L), trials = 2L))
#> Memory-network communities (map equation)
#>   Memory nodes: 5 | states: 5 | communities: 1
#>   States in more than one community: 0
#>   Codelength memory:      2.2406 bits (one module 2.2406)
#>   Codelength first-order: 2.2505 bits (one module 2.2505; communities: 1)
#>   Runs: 2 (teleportation 0.15, seed 1)
#>  node state community      flow
#>     a     a         1 0.1669582
#>     b     b         1 0.1611273
#>     c     c         1 0.1669582
#>     d     d         1 0.1611273
#>     h     h         1 0.3438290
```
