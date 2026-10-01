# Euler characteristic of a simplicial complex

\\\chi = \sum_k (-1)^k f_k\\, the alternating sum of the number of
k-simplices; by the Euler-Poincare theorem it equals the alternating sum
of the Betti numbers.

## Usage

``` r
hg_euler(sc)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md).

## Value

A single integer.

## References

Hatcher, A. (2002). *Algebraic Topology*. Cambridge University Press.

## See also

[`hg_betti()`](https://mohsaqr.github.io/hypernets/reference/hg_betti.md)

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
hg_euler(simplicial(adj))
#> [1] 1
```
