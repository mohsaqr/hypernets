# Tables of a simplicial complex

Tables of a simplicial complex

## Usage

``` r
# S3 method for class 'simplicial_complex'
hg_get(
  x,
  what = c("simplices", "f_vector", "betti", "degree", "validation"),
  ...,
  dimension = NULL,
  normalized = FALSE,
  top = NULL
)
```

## Arguments

- x:

  A `simplicial_complex` object from
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md).

- what:

  `"simplices"` (default) for one row per simplex, `"f_vector"` for the
  face counts by dimension, `"betti"` for the Betti numbers by
  dimension, `"degree"` for the simplicial degree of every node (the
  table of
  [`hg_degree()`](https://mohsaqr.github.io/hypernets/reference/hg_degree.md)),
  or `"validation"` for the tests of a window complex built with
  `validate = TRUE`.

- ...:

  Additional arguments (ignored).

- dimension:

  Integer or `NULL`. For `what = "simplices"`: keep only simplices of
  this dimension (`0` vertices, `1` edges, `2` triangles). For
  `what = "validation"`: keep only the sets of `dimension + 1` actions.

- normalized:

  Logical. Only for `what = "degree"`: divide each count by its largest
  possible value. Default `FALSE`.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "simplices"`, one row per simplex: `simplex`
(its number), `dimension`, `size`, `members` (the node names, comma
separated). For `what = "f_vector"`, one row per dimension: `dimension`,
`count`. For `what = "betti"`, one row per dimension: `dimension`,
`betti` (the number of connected parts, loops, voids, ...). For
`what = "degree"`, one row per node, most connected first: `node`, `d0`,
`d1`, ... and `total`. For `what = "validation"`, one row per tested set
of actions, by size, then p-value, then `z`: `members` (the actions,
comma separated), `size`, `count` (windows with exactly these actions),
`n_windows` (windows with this many distinct actions), `expected` (the
mean count under the null model), `z` (the count minus `expected`, in
null standard deviations), `p_value`, `p_adj` (Benjamini-Hochberg over
`n_tests` possible sets), `n_tests`, `significant`.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
sc <- simplicial(adj, type = "clique")
hg_get(sc)
#>    simplex dimension size members
#> 1        1         0    1       c
#> 2        2         0    1       b
#> 3        3         1    2    b, c
#> 4        4         0    1       d
#> 5        5         1    2    b, d
#> 6        6         2    3 b, c, d
#> 7        7         1    2    c, d
#> 8        8         0    1       a
#> 9        9         1    2    a, b
#> 10      10         2    3 a, b, c
#> 11      11         1    2    a, c
hg_get(sc, what = "f_vector")
#>   dimension count
#> 1         0     4
#> 2         1     5
#> 3         2     2
```
