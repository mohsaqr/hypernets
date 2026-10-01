# Summaries of results

[`summary()`](https://rdrr.io/r/base/summary.html) on a result returns
every table of the result as a list of data frames, named after the
`what` values of
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md),
so `summary(x)$validation` is the same data frame as
`hg_get(x, what = "validation")`. Results with overall figures add them
as further one-row or per-group tables. Tables that a particular result
does not hold are left out. The result itself is not changed.

## Usage

``` r
# S3 method for class 'hypernets_summary'
hg_get(x, what = names(x)[1L], ...)

# S3 method for class 'hypernets_summary'
print(x, n = 5L, ...)

# S3 method for class 'hypernets_result'
summary(object, ...)

# S3 method for class 'hypernets_hypa'
summary(object, ...)

# S3 method for class 'hypernets_community_comparison'
summary(object, ...)

# S3 method for class 'hg_communities'
summary(object, ...)

# S3 method for class 'net_hg_mmsbm'
summary(object, ...)

# S3 method for class 'hypernets_motifs'
summary(object, ...)

# S3 method for class 'net_hon_communities'
summary(object, ...)

# S3 method for class 'net_hon_boot_group'
summary(object, ...)

# S3 method for class 'net_hon_group'
summary(object, ...)

# S3 method for class 'net_hon_boot'
summary(object, ...)

# S3 method for class 'net_hon_compare'
summary(object, ...)

# S3 method for class 'hypernets_simplicial'
summary(object, ...)

# S3 method for class 'net_hg_topics'
summary(object, ...)
```

## Arguments

- x:

  A `hypernets_summary`.

- what:

  For
  [`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
  on a summary: the name of one of its tables. Default: the first.

- ...:

  Unused.

- n:

  For [`print()`](https://rdrr.io/r/base/print.html): number of rows of
  each table to show. Default `5`.

- object:

  A result: from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
  [`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md),
  [`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md),
  [`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
  [`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md),
  [`hg_markov_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_markov_stability.md),
  [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md) on
  sequences,
  [`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md),
  [`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md),
  [`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md),
  [`hg_compare_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_compare_communities.md),
  [`hg_mmsbm()`](https://mohsaqr.github.io/hypernets/reference/hg_mmsbm.md),
  [`hg_motifs()`](https://mohsaqr.github.io/hypernets/reference/hg_motifs.md),
  [`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md),
  [`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md),
  [`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md)
  or
  [`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md).

## Value

A `hypernets_summary`: a named list of data frames. Its print method
shows the first rows of each.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
fit_summary <- summary(memory(seqs, order = 2, min_count = 1))
fit_summary
#> Memory of order 2: 6 contexts, 1 change the most likely next state
#> Tables: contexts, overall
#> 
#> contexts (6 rows)
#>  context count entropy_first_order entropy_order_k entropy_drop kl
#>   a -> b     4                   1               0            1  1
#>   x -> b     4                   1               0            1  1
#>   b -> c     2                   0               0            0  0
#>   b -> d     2                   0               0            0  0
#>   c -> a     2                   0               0            0  0
#>  top_first_order top_order_k flips
#>                c           c FALSE
#>                c           d  TRUE
#>                a           a FALSE
#>                x           x FALSE
#>                b           b FALSE
#> 
#> overall (1 row)
#>  n_contexts n_flips kl_weighted entropy_drop_weighted
#>           6       1         0.5                   0.5
fit_summary$contexts
#>   context count entropy_first_order entropy_order_k entropy_drop kl
#> 1  a -> b     4                   1               0            1  1
#> 2  x -> b     4                   1               0            1  1
#> 3  b -> c     2                   0               0            0  0
#> 4  b -> d     2                   0               0            0  0
#> 5  c -> a     2                   0               0            0  0
#> 6  d -> x     2                   0               0            0  0
#>   top_first_order top_order_k flips
#> 1               c           c FALSE
#> 2               c           d  TRUE
#> 3               a           a FALSE
#> 4               x           x FALSE
#> 5               b           b FALSE
#> 6               b           b FALSE
```
