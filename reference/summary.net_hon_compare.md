# Summary method for net_hon_compare

Summary method for net_hon_compare

## Usage

``` r
# S3 method for class 'net_hon_compare'
summary(object, ...)
```

## Arguments

- object:

  A `net_hon_compare` object.

- ...:

  Additional arguments (ignored).

## Value

A data.frame, one row per rule order: `order`, `n_edges`,
`n_significant`, `max_abs_diff`. Returned **invisibly**: `summary(x)`
prints the summary and nothing else; assign the result to keep the
table.
