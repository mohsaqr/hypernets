# Build a Simplicial Complex

Constructs a simplicial complex from a network or higher-order pathway
object. Two construction methods are available:

- **Clique complex** (`"clique"`): every clique in the thresholded
  non-zero graph becomes a simplex. Edges with absolute weight \\\geq\\
  `threshold` are retained. The standard bridge from graph theory to
  algebraic topology.

- **Pathway complex** (`"pathway"`): each higher-order pathway from a
  `net_hon` or `net_hypa` becomes a simplex.

For `type = "vr"` (or alias `"rips"`), the input is treated as a
non-negative distance / dissimilarity matrix and a Vietoris-Rips
filtration is constructed: each k-simplex \\\sigma\\ enters at
\\\max\_{(i,j) \in \sigma} d(i,j)\\. Use `max_scale` to cap the
filtration diameter; edges with `d(i,j) > max_scale` are excluded.
Filtration values are attached as `$filtration` on the returned object
so
[`persistent_homology()`](https://mohsaqr.github.io/hypernets/reference/persistent_homology.md)
can read them directly.

## Usage

``` r
build_simplicial(
  x,
  type = "clique",
  threshold = 0,
  max_dim = 10L,
  max_pathways = NULL,
  anomaly = c("all", "over", "under"),
  max_scale = NULL,
  ...
)
```

## Arguments

- x:

  A square matrix, `tna`, `netobject`, `net_hon`, `net_hypa`, or
  `net_mogen`.

- type:

  Construction type: `"clique"` (default), `"pathway"`, or `"vr"` (alias
  `"rips"`).

- threshold:

  For `type = "clique"`: minimum non-zero absolute edge weight to
  include an edge (default 0). Edges below this are ignored; zero-weight
  non-edges are never included. Ignored for `type = "vr"` - use
  `max_scale` instead.

- max_dim:

  Maximum simplex dimension (default 10). Must be a single non-negative
  integer. A k-simplex has k+1 nodes.

- max_pathways:

  For `type = "pathway"`: maximum number of pathways to include, ranked
  by count (HON) or ratio (HYPA). `NULL` includes all. Default `NULL`.

- anomaly:

  For HYPA pathway complexes, which anomaly direction to include:
  `"all"` (default), `"over"`, or `"under"`. Under-represented HYPA
  paths are ranked by smallest observed/expected ratio; over-represented
  paths are ranked by largest ratio.

- max_scale:

  For `type = "vr"`: maximum edge length to include in the filtration.
  `NULL` (default) uses `max(d)`.

- ...:

  Additional arguments passed to
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
  when `x` is a `tna`/`netobject` with `type = "pathway"`.

## Value

A `net_simplicial` object. For `type = "vr"` an additional `$filtration`
numeric vector is attached (parallel to `$simplices`).

## See also

[`betti_numbers`](https://mohsaqr.github.io/hypernets/reference/betti_numbers.md),
[`persistent_homology`](https://mohsaqr.github.io/hypernets/reference/persistent_homology.md),
[`simplicial_degree`](https://mohsaqr.github.io/hypernets/reference/simplicial_degree.md),
[`q_analysis`](https://mohsaqr.github.io/hypernets/reference/q_analysis.md)

## Examples

``` r
mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
colnames(mat) <- rownames(mat) <- c("A","B","C")
sc <- build_simplicial(mat, threshold = 0.3)
print(sc)
#> Clique Complex 
#>   3 nodes, 7 simplices, dimension 2
#>   Density: 100.0%  |  Mean dim: 0.71  |  Euler: 1
#>   f-vector: (f0=3 f1=3 f2=1)
#>   Betti: b0=1
#>   Nodes: A, B, C 
betti_numbers(sc)
#> b0 b1 b2 
#>  1  0  0 

# Vietoris-Rips on a distance matrix:
d <- 1 - mat
diag(d) <- 0
sc_vr <- build_simplicial(d, type = "vr", max_scale = 0.6)
```
