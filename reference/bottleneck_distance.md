# Bottleneck Distance Between Persistence Diagrams

Computes the bottleneck distance between two persistence diagrams. For
finite pairs, the bottleneck distance is \$\$W\_\infty(D_1, D_2) =
\inf\_{\gamma} \sup\_{p \in D_1} \\p - \gamma(p)\\\_\infty,\$\$ where
\\\gamma\\ ranges over bijections \\D_1 \cup \Delta \to D_2 \cup
\Delta\\ and \\\Delta = \\(x,x)\\\\ is the diagonal. Each point may
match a point in the other diagram or its projection onto the diagonal
at cost \\\|d - b\|/2\\. Computed via binary search on \\\varepsilon\\
plus a Kuhn bipartite-matching feasibility check.

Essential classes (death = Inf in VR mode, or death = 0 in clique mode)
are matched one-to-one within each dimension. If the diagrams have
different numbers of essential classes in some dimension, the bottleneck
distance for that dimension is `Inf`.

## Usage

``` r
bottleneck_distance(d1, d2, dimension = NULL, tol = .Machine$double.eps^0.5)
```

## Arguments

- d1, d2:

  `net_persistent_homology` objects, or data.frames with columns
  `dimension`, `birth`, `death`.

- dimension:

  Integer vector of dimensions to compare. `NULL` (default) compares all
  dimensions appearing in either diagram and returns a named numeric
  vector.

- tol:

  Numerical tolerance for binary search (default
  `.Machine$double.eps ^ 0.5`).

## Value

Named numeric vector. Names are `"dim_<k>"`. `Inf` indicates a
structural mismatch (different essential counts in that dimension); a
self-distance is always 0.

## References

Edelsbrunner, H. & Harer, J. (2010). *Computational Topology: An
Introduction*. AMS. Section VIII.

## Examples

``` r
mat1 <- matrix(c(0, .6, .5, .6, 0, .4, .5, .4, 0), 3, 3)
rownames(mat1) <- colnames(mat1) <- c("A","B","C")
ph1 <- persistent_homology(mat1, n_steps = 5)
bottleneck_distance(ph1, ph1)  # self-distance is 0
#> dim_0 
#>     0 
```
