# Persistence Landscape

Computes the persistence landscape (Bubenik 2015) from a persistence
diagram. Each (birth, death) pair contributes a tent function
\$\$\Lambda\_{(b,d)}(t) = \max(0, \min(t - b, d - t)).\$\$ The \\k\\-th
landscape function \\\lambda^{(k)}(t)\\ is the \\k\\-th largest of
\\\\\Lambda\_{(b_i,d_i)}(t)\\\_i\\ at each \\t\\. Landscapes are stable
under bottleneck distance and form a Banach-space embedding of
persistence diagrams.

## Usage

``` r
persistence_landscape(ph, k_max = 5L, dimension = 1L, t_grid = NULL)
```

## Arguments

- ph:

  A `net_persistent_homology` object or a data.frame with columns
  `dimension`, `birth`, `death`.

- k_max:

  Maximum landscape index to compute (default 5). Must be a single
  positive integer.

- dimension:

  Integer scalar – which homology dimension to compute the landscape
  for. Default 1.

- t_grid:

  Numeric vector of evaluation points. `NULL` (default) uses an even
  grid of 200 points covering the union of pair intervals.

## Value

A `net_persistence_landscape` object with:

- landscape:

  Data frame: `k`, `t`, `value`.

- dimension:

  Integer scalar.

- k_max:

  Integer scalar.

- t_grid:

  Numeric vector.

## References

Bubenik, P. (2015). Statistical topological data analysis using
persistence landscapes. *Journal of Machine Learning Research* **16**,
77-102.

## Examples

``` r
mat <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
rownames(mat) <- colnames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 5)
pl <- persistence_landscape(ph, k_max = 3, dimension = 0)
```
