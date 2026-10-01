# Tables of a hypergraph community ensemble

Reads, prints and plots the result of
[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
on a hypergraph.

## Usage

``` r
# S3 method for class 'hg_communities'
print(x, n = 10L, ...)

# S3 method for class 'hg_communities'
hg_get(
  x,
  what = c("medoid", "partitions", "runs", "sizes", "ami", "ari", "nmi", "weights"),
  ...
)

# S3 method for class 'hg_communities'
plot(x, ...)
```

## Arguments

- x:

  An `hg_communities` object.

- n:

  Number of rows of the default table to print. Default `10`.

- ...:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html),
  additional arguments passed to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html);
  otherwise unused.

- what:

  Table to return: `"medoid"` (default: one row per node with its
  community in the AMI-medoid run), `"partitions"`, `"runs"`, `"sizes"`,
  `"ami"`, `"ari"`, `"nmi"`, or (IRMM fits only) `"weights"`. The three
  similarity tables have one row per distinct pair of runs (`run_a`,
  `run_b`, and the similarity), without the diagonal.

## Value

[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md): a
base data.frame. [`print()`](https://rdrr.io/r/base/print.html): `x`,
invisibly. [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
cograph plot of the projection, coloured by the medoid communities.

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "a", "b", "c", "x", "y", "z", "x", "y", "z"),
  edge = rep(paste0("e", 1:4), each = 3)
)
h <- group_hypergraph(dat, "member", "edge")
if (requireNamespace("igraph", quietly = TRUE)) {
  fit <- hg_communities(h, n_runs = 2, trials = 2, seeds = 1:2)
  hg_get(fit, what = "runs")
}
#>   run seed n_communities codelength
#> 1   1    1             2   2.052397
#> 2   2    2             2   2.052397
```
