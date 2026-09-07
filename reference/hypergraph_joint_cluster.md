# Joint NMF clustering with an auxiliary vertex relation

Implements the two multi-view objectives in Section 6.3.2 of Hayashi et
al. (2020). `method = "joint"` is Eq. 18 (J-NMF), jointly factorizing
the edge-dependent incidence `X` and a vertex relation `S`.
`method = "joint_symmetric"` is Eq. 19 (JS-NMF), jointly factorizing the
representative-digraph similarity `T = I - L` and `S`. The paper uses a
symmetric patent-citation indicator for `S`; any finite non-negative
node-by-node relation matrix can be supplied here.

## Usage

``` r
hypergraph_joint_cluster(
  hg,
  relations,
  k,
  method = c("joint", "joint_symmetric"),
  alpha = 1,
  beta = 1,
  gamma = 1,
  type = c("random_walk", "zhou"),
  edge_weights = NULL,
  nstart = 10L,
  seed = NULL,
  max_iter = 500L,
  tol = 1e-06
)
```

## Arguments

- hg:

  A connected `net_hypergraph`.

- relations:

  A non-negative `n_nodes` by `n_nodes` numeric matrix. If it has
  dimnames, rows and columns are reordered to `hg$nodes`.

- k:

  Number of clusters, between 2 and `n_nodes - 1`.

- method:

  `"joint"` for J-NMF (Eq. 18) or `"joint_symmetric"` for JS-NMF (Eq.
  19).

- alpha:

  Non-negative symmetry penalty in JS-NMF.

- beta:

  Non-negative agreement penalty between the content and relation vertex
  factors.

- gamma:

  Non-negative weight of the auxiliary-relation loss.

- type, edge_weights:

  Laplacian inputs for JS-NMF.

- nstart:

  Number of random initializations.

- seed:

  Optional integer initialization seed.

- max_iter:

  Maximum multiplicative-update iterations.

- tol:

  Relative objective tolerance.

## Value

A `net_hypergraph_cluster` object. Its `$embedding` is the fitted vertex
factor `M`; `$params` contains all fitted factors, the objective trace
and convergence diagnostics.

## References

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*, 495-504.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)
