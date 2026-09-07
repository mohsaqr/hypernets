# Persistent Homology

Computes persistent homology via full boundary-matrix reduction over
\\\mathbb{Z}/2\\ (Edelsbrunner, Letscher & Zomorodian 2000). The
returned persistence diagram pairs each k-dimensional homology class to
the simplex whose addition creates it (birth) and the simplex whose
addition destroys it (death). Essential classes - those never killed -
are reported with `death = 0` in clique mode (similarity scale,
descending) and `death = Inf` in VR mode (distance scale, ascending).

Two filtration modes are supported:

- `type = "clique"`:

  Weighted clique filtration. Input is treated as a similarity matrix;
  high-weight simplices appear early. For each k-simplex \\\sigma\\, the
  filtration value is \\\min\_{(i,j) \in \sigma}\\\|w(i,j)\|\\.
  Thresholds run high to low.

- `type = "vr"`:

  Vietoris-Rips filtration on a non-negative distance matrix. For each
  k-simplex \\\sigma\\, the filtration value is \\\max\_{(i,j) \in
  \sigma}\\d(i,j)\\. Thresholds run low to high. Use `max_scale` to cap
  the filtration diameter.

## Usage

``` r
persistent_homology(
  x,
  n_steps = 20L,
  max_dim = 3L,
  type = "clique",
  max_scale = NULL
)
```

## Arguments

- x:

  A square matrix, `tna`, or `netobject`. For `type = "vr"`, must be a
  non-negative distance matrix.

- n_steps:

  Number of grid points for the reported Betti curve (default 20). The
  persistence diagram itself is exact - it does not depend on `n_steps`.

- max_dim:

  Maximum simplex dimension to track (default 3).

- type:

  Filtration: `"clique"` (default, similarity-weighted) or `"vr"`
  (Vietoris-Rips on distances).

- max_scale:

  For `type = "vr"` only: cap on edge length. Edges with
  `d(i,j) > max_scale` are excluded. `NULL` (default) uses `max(d)`.

## Value

A `net_persistent_homology` object with:

- betti_curve:

  Data frame: `threshold`, `dimension`, `betti`.

- persistence:

  Data frame of birth-death pairs: `dimension`, `birth`, `death`,
  `persistence`. Sorted by descending persistence.

- thresholds:

  Numeric vector of grid thresholds.

- mode:

  Either `"clique"` or `"vr"`.

## References

Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2000). Topological
persistence and simplification. *Discrete & Computational Geometry*
**28**, 511-533.

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
ph <- persistent_homology(mat, n_steps = 10)
print(ph)
#> Persistent Homology
#>   10 filtration steps [0.6000 → 0.0060]
#>   Features: b0: 2 (1 persistent) 
#>   Longest-lived:
#>     b0: 0.6000 → 0.0000 (life: 0.6000)
#>     b0: 0.6000 → 0.5000 (life: 0.1000)
```
