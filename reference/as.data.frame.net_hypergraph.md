# Coerce a net_hypergraph to a tidy hyperedge table

One row per hyperedge: its identifier, size (number of distinct member
states), the member states, and its weight (for hypergraphs built by
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
the number of windows collapsed into the hyperedge; `NA` for the other
constructors, whose hyperedges are unweighted).

## Usage

``` r
# S3 method for class 'net_hypergraph'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  ...,
  what = c("edges", "nodes", "memberships"),
  sort_by = NULL,
  top = NULL
)
```

## Arguments

- x:

  A `net_hypergraph` object.

- row.names:

  Ignored (present for S3 consistency with the generic).

- optional:

  Ignored (present for S3 consistency with the generic).

- ...:

  Additional arguments (ignored).

- what:

  `"edges"` (default) for one row per hyperedge, `"nodes"` for one row
  per node, or `"memberships"` for one row per node-in- hyperedge cell
  of the incidence matrix.

- sort_by:

  `NULL` (construction order, default), `"weight"`, or `"size"` - sort
  the table by that column, largest first (ties broken by hyperedge id
  so the order is deterministic).

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data.frame. For `what = "edges"`, one row per hyperedge with columns
`hyperedge` (character id), `size` (integer), `states` (comma-separated
member states), and `weight` (numeric window count, or `NA`). For
`what = "nodes"`, one row per node with columns `node` and `degree` (the
number of hyperedges it belongs to), plus `block` for a hypergraph with
planted blocks
([`hg_sample_sbm()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_sbm.md));
`sort_by = "degree"` orders it. For `what = "memberships"`, one row per
non-zero incidence cell with columns `node`, `hyperedge` and `weight`
(the incidence value: 1 for a binary hypergraph, the summed weight for a
weighted one), in hyperedge order; `sort_by = "weight"` orders it.

## Examples

``` r
hg <- window_hypergraph(list(s1 = c("a", "b", "a", "c")), window = 2L)
as.data.frame(hg)
#>   hyperedge size states weight
#> 1        h1    2   a, b      2
#> 2        h2    2   a, c      1
as.data.frame(hg, sort_by = "weight")
#>   hyperedge size states weight
#> 1        h1    2   a, b      2
#> 2        h2    2   a, c      1
as.data.frame(hg, sort_by = "weight", top = 3)
#>   hyperedge size states weight
#> 1        h1    2   a, b      2
#> 2        h2    2   a, c      1
as.data.frame(hg, what = "nodes", sort_by = "degree")
#>   node degree
#> 1    a      2
#> 2    b      1
#> 3    c      1
as.data.frame(hg, what = "memberships")
#>   node hyperedge weight
#> 1    a        h1      2
#> 2    b        h1      2
#> 3    a        h2      1
#> 4    c        h2      1
```
