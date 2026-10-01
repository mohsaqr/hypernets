# Q-analysis of a simplicial complex

Q-connectivity (Atkin 1974): two maximal simplices are q-connected when
they share a face of dimension at least q. Reports the number of
q-connected components at every level and each node's highest level.

## Usage

``` r
hg_qanalysis(sc)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md).

## Value

A `q_analysis` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
the components per level (`what = "q_levels"`, the default) or each
node's highest level (`"nodes"`).

## References

Atkin, R. H. (1974). *Mathematical Structure in Human Affairs*.
Heinemann.

## See also

[`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
hg_get(hg_qanalysis(simplicial(adj)))
#>   q components
#> 1 2          2
#> 2 1          1
#> 3 0          1
```
