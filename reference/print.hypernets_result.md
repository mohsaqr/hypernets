# Print a result of an imported estimator

Prints one line naming the result and its main settings, then the first
rows of its default table, the table
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md)
returns.

## Usage

``` r
# S3 method for class 'hypernets_result'
hg_get(x, ...)

# S3 method for class 'hypernets_result'
print(x, n = 10L, ...)

# S3 method for class 'hypernets_result'
plot(x, y, what = c("states", "passage_time"), ...)
```

## Arguments

- x:

  A result of
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
  [`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md),
  [`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md),
  [`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md),
  [`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md),
  [`hg_markov_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_markov_stability.md),
  [`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md),
  [`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md)
  or
  [`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md).

- ...:

  Unused.

- n:

  Number of rows of the table to show. Default `10`.

- y:

  Unused.

- what:

  For a Markov stability result: `"states"` (default), one panel of bars
  per measure, or `"passage_time"`, the mean first-passage times as a
  heatmap (rows are the starting state).

## Value

`x`, invisibly.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
ggplot for a persistence landscape (the landscapes that are not zero
everywhere, one colour and line type each) and for a Markov stability
result; every other result is plotted by its estimator's method.

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
print(hon(seqs, max_order = 2), n = 3)
#> Higher-order network: 5 states, 5 nodes, 6 rules (highest order 1, min_freq 1, 2 sequences)
#>    path from to count probability from_order to_order
#>  a -> b    a  b     2         1.0          1        1
#>  b -> c    b  c     2         0.5          1        1
#>  b -> d    b  d     2         0.5          1        1
#> ... 3 more rows
```
