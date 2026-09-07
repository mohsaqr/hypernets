# AllSet hypergraph node classifier

Implements the AllSet framework of Chien et al. (2022) as one complete
node-to-edge-to-node propagation layer. `model = "deepsets"` uses Eq. 7,
an MLP before and after each multiset sum. `model = "transformer"` uses
Eq. 8: pooling by multihead attention from a learned query, followed by
residual MLPs and layer normalization. Both directions have independent
learned multiset functions and are permutation invariant.

## Usage

``` r
hypergraph_allset(
  hg,
  labels,
  features = "incidence",
  model = c("deepsets", "transformer"),
  hidden = 128L,
  heads = 4L,
  epochs = 200L,
  lr = 0.001,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  seed = 1L,
  verbose = FALSE
)

hg_allset(
  hg,
  labels,
  features = "incidence",
  model = c("deepsets", "transformer"),
  hidden = 128L,
  heads = 4L,
  epochs = 200L,
  lr = 0.001,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  seed = 1L,
  verbose = FALSE
)

hypergraph_alldeepsets(...)

hypergraph_allset_transformer(...)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hypergraph`), dense or sparse.

- labels:

  The known labels: a named character vector (names are node
  identifiers, values class labels) or a tidy data.frame with a `node`
  column and a `label`, `cluster` or `predicted` column. At least two
  classes.

- features:

  Vertex feature matrix, one row per node in `hg` node order (rownames
  must match the node names), e.g. sbert embeddings. The default
  `"incidence"` uses the hypergraph's own weighted incidence rows
  (tf-idf bag-of-words features for a `nodes = "doc"` text hypergraph).

- model:

  `"deepsets"` for AllDeepSets or `"transformer"` for AllSetTransformer.

- hidden:

  Hidden width (default 128).

- heads:

  Number of attention heads for AllSetTransformer. `hidden` must be
  divisible by `heads`.

- epochs:

  Training epochs (default 200).

- lr, weight_decay:

  Adam learning rate and L2 penalty.

- dropout:

  Hidden-layer dropout.

- validation:

  Fraction of labeled nodes held out to select an epoch.

- seed:

  Reproducible R and torch seed.

- verbose:

  Message loss every 20 epochs.

- ...:

  Arguments passed to `hypergraph_allset()`.

## Value

A prediction data.frame in the same format as
[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md),
with training history and the AllSet model name attached.

## References

Chien, E., Pan, C., Peng, J., & Milenkovic, O. (2022). You are AllSet: A
multiset function framework for hypergraph neural networks. *ICLR 2022*.
