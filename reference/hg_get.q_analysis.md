# Tables of a q-analysis

Tables of a q-analysis

## Usage

``` r
# S3 method for class 'q_analysis'
hg_get(x, what = c("q_levels", "nodes"), ..., top = NULL)
```

## Arguments

- x:

  A `q_analysis` object from
  [`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md).

- what:

  `"q_levels"` (default) for the number of q-connected components at
  each level, or `"nodes"` for each node's highest q-level.

- ...:

  Additional arguments (ignored).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "q_levels"`, one row per level: `q`,
`components`. For `what = "nodes"`, one row per node: `node`, `max_q`.

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
qa <- hg_qanalysis(simplicial(adj, type = "clique"))
hg_get(qa)
#>   q components
#> 1 2          2
#> 2 1          1
#> 3 0          1
hg_get(qa, what = "nodes")
#>   node max_q
#> 1    a     2
#> 2    b     2
#> 3    c     2
#> 4    d     2
```
