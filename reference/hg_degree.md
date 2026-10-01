# Simplicial degree of every node

How many simplices of each dimension contain each node.

## Usage

``` r
hg_degree(sc, normalized = FALSE)
```

## Arguments

- sc:

  A `simplicial_complex` from
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md).

- normalized:

  Logical. Divide each count by its largest possible value. Default
  `FALSE`.

## Value

A base `data.frame`, one row per node, most connected first: `node`, one
column per dimension (`d0`, `d1`, ...), and `total` (the simplices of
dimension one and above). The same table is
`hg_get(sc, what = "degree")`.

## References

Serrano, D. H., & Gomez, D. S. (2020). Centrality measures in simplicial
complexes: Applications of topological data analysis to network science.
*Applied Mathematics and Computation*, 382, 125331.
[doi:10.1016/j.amc.2020.125331](https://doi.org/10.1016/j.amc.2020.125331)

## See also

[`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
hg_degree(simplicial(adj))
#>   node d0 d1 d2 total
#> 1    b  1  3  2     5
#> 2    c  1  3  2     5
#> 3    a  1  2  1     3
#> 4    d  1  2  1     3
```
