# Summary method for net_hon_boot

Summary method for net_hon_boot

## Usage

``` r
# S3 method for class 'net_hon_boot'
summary(object, ...)
```

## Arguments

- object:

  A `net_hon_boot` object.

- ...:

  Additional arguments (ignored).

## Value

A data.frame, one row per rule order: `order`, `n_edges`,
`mean_support`, `min_support`, `mean_ci_width`. Returned **invisibly**:
`summary(x)` prints the summary and nothing else; assign the result to
keep the table.
