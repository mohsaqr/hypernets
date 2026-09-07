# Extract Pathways from Higher-Order Network Objects

Extracts higher-order pathway strings suitable for
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
Each pathway represents a multi-step dependency: source states lead to a
target state.

For `net_hon`: extracts edges where the source node is higher-order
(order \> 1), i.e., the transitions that differ from first-order Markov.

For `net_hypa`: extracts anomalous paths (over- or under-represented
relative to the hypergeometric null model).

For `net_mogen`: extracts all transitions at the optimal order (or a
specified order).

## Usage

``` r
pathways(x, ...)

# S3 method for class 'net_hon'
pathways(x, min_count = 1L, min_prob = 0, top = NULL, order = NULL, ...)

# S3 method for class 'net_hypa'
pathways(x, type = "all", ...)

# S3 method for class 'net_mogen'
pathways(x, order = NULL, min_count = 1L, min_prob = 0, top = NULL, ...)
```

## Arguments

- x:

  A higher-order network object (`net_hon`, `net_hypa`, or `net_mogen`).

- ...:

  Additional arguments.

- min_count:

  Integer. Minimum transition count to include (default: 1).

- min_prob:

  Numeric. Minimum transition probability to include (default: 0).

- top:

  Integer or NULL. Return only the top N pathways ranked by count
  (default: NULL = all).

- order:

  Integer or NULL. Markov order to extract. Default: optimal order from
  model selection.

- type:

  Character. Which anomalies to include: `"all"` (default), `"over"`, or
  `"under"`.

## Value

A character vector of pathway strings in arrow notation (e.g.
`"A B -> C"`), suitable for
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).

A character vector of pathway strings.

A character vector of pathway strings.

A character vector of pathway strings.

## Methods (by class)

- `pathways(net_hon)`: Extract higher-order pathways from HON

- `pathways(net_hypa)`: Extract anomalous pathways from HYPA

- `pathways(net_mogen)`: Extract transition pathways from MOGen

## Examples

``` r
# \donttest{
seqs <- list(c("A","B","C","D"), c("A","B","C","A"))
hon <- build_hon(seqs, max_order = 3)
pw <- pathways(hon)
# }
```
