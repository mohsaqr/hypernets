# Summary Method for net_hypa

Summary Method for net_hypa

## Usage

``` r
# S3 method for class 'net_hypa'
summary(
  object,
  n = 10L,
  type = c("all", "over", "under"),
  order_by = c("sig", "freq", "frequency", "ratio", "path"),
  ...
)
```

## Arguments

- object:

  A `net_hypa` object.

- n:

  Integer. Maximum number of paths to display per category (default:
  10).

- type:

  Character. Which anomalies to show: `"all"` (default), `"over"`, or
  `"under"`.

- order_by:

  Character. Ranking used within each anomaly direction: `"sig"` ranks
  by the active tail probability, `"freq"` by observed count, `"ratio"`
  by observed/expected ratio, and `"path"` alphabetically.

- ...:

  Additional arguments (ignored).

## Value

A data frame with path, observed, expected, ratio, p_tail, and direction
columns.

Returned **invisibly**: `summary(x)` prints the summary and nothing
else; assign the result to keep the table.

## Examples

``` r
seqs <- list(c("A","B","C"), c("B","C","A"), c("A","C","B"), c("A","B","C"))
hyp <- build_hypa(seqs, k = 2)
#> Warning: 'k' is deprecated; use 'order' instead.
summary(hyp)
#> HYPA Summary
#> 
#>   Order(s): 2 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0) | order_by: sig
#> 
#>   No anomalous paths detected.

# \donttest{
seqs <- data.frame(
  V1 = c("A","B","C","A","B","C","A","B","C","A"),
  V2 = c("B","C","A","B","C","A","B","C","A","B"),
  V3 = c("C","A","B","C","A","B","C","A","B","C"),
  V4 = c("A","B","C","A","B","C","A","B","C","A")
)
hypa <- build_hypa(seqs, k = 2L)
#> Warning: 'k' is deprecated; use 'order' instead.
summary(hypa)
#> HYPA Summary
#> 
#>   Order(s): 2 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0) | order_by: sig
#> 
#>   No anomalous paths detected.
summary(hypa, type = "over", n = 5)
#> HYPA Summary
#> 
#>   Order(s): 2 | Edges: 3
#>   Alpha: 0.05 | p_adjust: BH
#>   Anomalous: 0 (over: 0, under: 0) | order_by: sig
#> 
#>   No anomalous paths detected.
# }
```
