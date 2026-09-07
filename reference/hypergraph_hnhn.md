# HNHN node classifier

Trains Hypergraph Networks with Hyperedge Neurons (Dong et al. 2020). A
full layer sends node features to learned hyperedge representations and
then sends those representations back to nodes. The paper's
normalization is implemented directly: `alpha` changes the relative
contribution of large hyperedges in edge-to-node messages, while `beta`
changes the contribution of high-degree nodes in node-to-edge messages.
Zero recovers ordinary neighborhood means in each direction.

## Usage

``` r
hypergraph_hnhn(
  hg,
  labels,
  features = "incidence",
  hidden = 128L,
  alpha = 0,
  beta = 0,
  epochs = 200L,
  lr = 0.01,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  seed = 1L,
  verbose = FALSE
)

hg_hnhn(
  hg,
  labels,
  features = "incidence",
  hidden = 128L,
  alpha = 0,
  beta = 0,
  epochs = 200L,
  lr = 0.01,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any honets `net_hypergraph`), dense or sparse.

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

- hidden:

  Hidden width (default 128).

- alpha:

  Hyperedge-cardinality normalization exponent.

- beta:

  Vertex-degree normalization exponent.

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

## Value

A prediction data.frame in the same format as
[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md),
with training history and the normalization exponents attached.

## References

Dong, Y., Sawin, W., & Bengio, Y. (2020). HNHN: Hypergraph networks with
hyperedge neurons. *ICML Graph Representation Learning and Beyond
Workshop*.
