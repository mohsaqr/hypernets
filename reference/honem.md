# Embedding of a higher-order network (HONEM)

Low-dimensional node embedding of a higher-order network that preserves
its higher-order dependencies: exponentially decaying powers of the
network's transition matrix, followed by a truncated singular value
decomposition.

## Usage

``` r
honem(hon, dim = 32L, max_power = 10L)
```

## Arguments

- hon:

  A `net_hon` from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md), or a
  square weighted adjacency matrix.

- dim:

  Integer. Embedding dimension, capped at the number of nodes minus one.
  Default `32`.

- max_power:

  Integer. Longest walk length in the neighbourhood matrix. Default
  `10`.

## Value

A `net_honem` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the node coordinates (`what = "embeddings"`, the default) or the
variance per dimension (`"variance"`).

## References

Saebi, M., Ciampaglia, G. L., Kaplan, L. M., & Chawla, N. V. (2020).
HONEM: Learning embedding for higher order networks. *Big Data*, 8(4),
255–269.
[doi:10.1089/big.2019.0169](https://doi.org/10.1089/big.2019.0169)

## See also

[`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md)

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
emb <- honem(hon(seqs, max_order = 2), dim = 2)
hg_get(emb)
#>     node      dim1          dim2
#> 1      a -0.366535 -4.082483e-01
#> 2 a -> b -0.428539 -4.082483e-01
#> 3      b -0.428539  2.578122e-15
#> 4      c -0.366535 -4.082483e-01
#> 5      d -0.366535  4.082483e-01
#> 6      x -0.366535  4.082483e-01
#> 7 x -> b -0.428539  4.082483e-01
```
