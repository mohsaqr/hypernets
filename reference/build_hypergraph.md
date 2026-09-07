# Higher-order hypergraph from a network's clique structure

Takes a network and produces a hypergraph by promoting k-cliques (k \>=
3) to k-hyperedges. Each k-clique is independently included as a
k-hyperedge with probability `p`. Optionally retains the underlying
pairwise edges as 2-hyperedges. Foundation for higher-order analyses.

## Usage

``` r
build_hypergraph(
  net,
  p = 1,
  type = c("clique", "vr", "rips"),
  include_pairwise = TRUE,
  max_size = 3L,
  threshold = 0,
  seed = NULL
)

# S3 method for class 'net_hypergraph'
print(x, ...)

# S3 method for class 'net_hypergraph'
summary(object, ...)
```

## Arguments

- net:

  A `netobject`, `cograph_network`, `net_simplicial`, or numeric
  adjacency / weight matrix. Directed inputs are symmetrised by the
  underlying clique enumerator.

- p:

  Probability in `[0, 1]` that each k-clique with k \>= 3 becomes a
  k-hyperedge. Default `1` (deterministic - every found clique is
  included).

- type:

  Hyperedge enumeration. `"clique"` (default) promotes k-cliques in the
  binarised adjacency to k-hyperedges. A metric Vietoris-Rips
  construction is **not implemented**; `"vr"` / `"rips"` are accepted by
  `match.arg` but raise an error rather than silently aliasing
  `"clique"`.

- include_pairwise:

  Logical. Include 2-edges from the input network as 2-hyperedges.
  Default `TRUE`. Set `FALSE` for a "fully higher-order" hypergraph
  containing only k-hyperedges with k \>= 3.

- max_size:

  Integer \>= 2. Maximum hyperedge size to extract. Default `3L`
  (triangles only). `4L` also includes 4-cliques as 4-hyperedges, etc.

- threshold:

  Numeric. Edge weight cutoff used to binarise the adjacency for clique
  enumeration. Default `0` (any non-zero weight is an edge).

- seed:

  Optional integer for reproducible Bernoulli sampling when `0 < p < 1`.

- x:

  A `net_hypergraph` object (for `print`).

- ...:

  Additional arguments (ignored).

- object:

  A `net_hypergraph` object (for `summary`).

## Value

A `net_hypergraph` object: a list with components

- `hyperedges`:

  List of integer vectors. Each entry is a hyperedge given as the sorted
  node indices it spans.

- `incidence`:

  Numeric matrix of size `n_nodes` x `n_hyperedges`.
  `incidence[i, j] = 1` iff node i belongs to hyperedge j. Row names are
  node names; column names are `h1`, `h2`, ...

- `nodes`:

  Character vector of node names.

- `n_nodes`, `n_hyperedges`:

  Scalar counts.

- `size_distribution`:

  Named integer vector: number of hyperedges of each size, named
  `size_2`, `size_3`, ...

- `params`:

  Recorded call parameters: `source`, `type`, `p`, `include_pairwise`,
  `max_size`, `threshold`, `seed`.

The input `x` invisibly.

A data.frame, one row per node: `node`, `degree` (the number of
hyperedges the node belongs to). Returned **invisibly**: `summary(x)`
prints the summary and nothing else; assign the result to keep the
table.

## Details

The construction follows Burgio, Matamalas, Gomez & Arenas (2020) on
simplicial / hypergraph contagion. For each k-clique with k \>= 3 found
in the underlying graph (via
[`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md)),
an independent Bernoulli(`p`) trial decides whether that clique becomes
a k-hyperedge. Underlying pairwise edges are always retained when
`include_pairwise = TRUE`, so the resulting hypergraph contains both the
original 2-edge structure and the sampled higher-order interactions.

At the limits:

- `p = 0` with `include_pairwise = TRUE` reproduces the input pairwise
  network as a hypergraph of size-2 edges.

- `p = 1` with `include_pairwise = FALSE` returns a fully higher-order
  hypergraph containing only the k-hyperedges (k \>= 3) found in the
  network's clique complex.

## References

Burgio, G., Matamalas, J. T., Gomez, S., & Arenas, A. (2020). Evolution
of cooperation in the presence of higher-order interactions: from
networks to hypergraphs. *Entropy* 22(7), 744.
[doi:10.3390/e22070744](https://doi.org/10.3390/e22070744)

## See also

[`build_simplicial()`](https://mohsaqr.github.io/hypernets/reference/build_simplicial.md)
(underlying clique enumeration), `Nestimate::build_network()`.

## Examples

``` r
set.seed(1)
n <- 8
adj <- matrix(stats::rbinom(n * n, 1, 0.5), n, n)
diag(adj) <- 0
adj <- (adj + t(adj)) > 0
rownames(adj) <- colnames(adj) <- LETTERS[seq_len(n)]
hg <- build_hypergraph(adj, p = 1, max_size = 3L)
print(hg)
#> Hypergraph: 8 nodes, 51 hyperedges
#> Size distribution:
#>   size_2   : 23
#>   size_3   : 28
#> Source: network cliques (p = 1.00, include_pairwise = TRUE, max_size = 3)
summary(hg)
#> Hypergraph summary
#>   Nodes:         8
#>   Hyperedges:    51
#>   Mean size:     2.55
#>   Max size:      3
```
