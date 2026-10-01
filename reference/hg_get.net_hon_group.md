# Tables of a group model of memory networks

Returns the requested table of every group's memory network, stacked,
with a `group` column first.

## Usage

``` r
# S3 method for class 'net_hon_group'
print(x, n = 10L, ...)

# S3 method for class 'net_hon_group'
hg_get(x, ...)

# S3 method for class 'net_hon_boot_group'
print(x, n = 10L, ...)

# S3 method for class 'net_hon_boot_group'
hg_get(x, ...)
```

## Arguments

- x:

  A `net_hon_group` from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md) with
  `group`.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  Passed to
  [`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
  on each group's network (`what`, `order_min`, `sort_by`, `top`, ...).

## Value

[`print()`](https://rdrr.io/r/base/print.html) returns `x` invisibly.

A data.frame: `group`, then the columns of the requested table.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b"), c("a", "b", "c", "a", "b"),
             c("x", "b", "d", "x", "b"), c("x", "b", "d", "x", "b"))
by_group <- hon(seqs, group = c("g1", "g1", "g2", "g2"), max_order = 2)
by_group
#> Memory networks by group: 2 groups (g1: 2 sequences, g2: 2 sequences)
#>  group   path from to count probability from_order to_order
#>     g1 a -> b    a  b     4           1          1        1
#>     g1 b -> c    b  c     2           1          1        1
#>     g1 c -> a    c  a     2           1          1        1
#>     g2 b -> d    b  d     2           1          1        1
#>     g2 d -> x    d  x     2           1          1        1
#>     g2 x -> b    x  b     4           1          1        1
hg_get(by_group, top = 3)
#>   group   path from to count probability from_order to_order
#> 1    g1 a -> b    a  b     4           1          1        1
#> 2    g1 b -> c    b  c     2           1          1        1
#> 3    g1 c -> a    c  a     2           1          1        1
#> 4    g2 b -> d    b  d     2           1          1        1
#> 5    g2 d -> x    d  x     2           1          1        1
#> 6    g2 x -> b    x  b     4           1          1        1
```
