# Hypergeometric anomalies of paths or of co-occurring node pairs

One verb, two estimators, chosen by the class of `x`:

- a hypergraph (`net_hg`):

  co-occurring node pairs scored against a hypergeometric null built
  from the hyperdegrees; returns a table. See
  [`hypa.net_hg()`](https://mohsaqr.github.io/hypernets/reference/hypa.net_hg.md).

- sequences (a long, wide or list input, or a model object carrying
  sequences; see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md)):

  paths of a De Bruijn graph scored against the same null (HYPA);
  returns a `net_hypa` object. See
  [`hypa.default()`](https://mohsaqr.github.io/hypernets/reference/hypa.default.md).

Each method keeps its own arguments; passing an argument that only the
other method takes raises `hypernets_bad_input`.

## Usage

``` r
hypa(x, ...)
```

## Arguments

- x:

  A `net_hg`, or sequences.

- ...:

  Arguments of the method for `class(x)`.

## Value

See the methods.

## References

LaRock, T., Nanumyan, V., Scholtes, I., Casiraghi, G., Eliassi-Rad, T.,
and Schweitzer, F. (2020). HYPA: Efficient detection of path anomalies
in time series data on networks. *Proceedings of the 2020 SIAM
International Conference on Data Mining*, 460-468.
[doi:10.1137/1.9781611976236.52](https://doi.org/10.1137/1.9781611976236.52)

## Examples

``` r
seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"),
             c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
hg_get(hypa(seqs, order = 2, min_count = 1))
#> [1] path      count     expected  ratio     p_adj     direction
#> <0 rows> (or 0-length row.names)
```
