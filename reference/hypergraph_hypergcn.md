# HyperGCN node classifier

Implements Yadati et al.'s (2019) mediator expansion and two-layer GCN.
Every hyperedge selects the pair of vertices farthest apart in the
current transformed representation. With `method = "hypergcn"`, the
expansion is recomputed at both layers and every epoch (Algorithm 1).
`"fast"` fixes it once from the input features (Algorithm 2), and
`"one"` uses only the farthest pair without mediators (Algorithm 3).
Exact distance ties are resolved lexicographically for reproducibility.

## Usage

``` r
hypergraph_hypergcn(
  hg,
  labels,
  features = "incidence",
  method = c("hypergcn", "fast", "one"),
  hidden = 128L,
  epochs = 200L,
  lr = 0.01,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  edge_weights = NULL,
  seed = 1L,
  verbose = FALSE
)

hg_hypergcn(
  hg,
  labels,
  features = "incidence",
  method = c("hypergcn", "fast", "one"),
  hidden = 128L,
  epochs = 200L,
  lr = 0.01,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  edge_weights = NULL,
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

- method:

  `"hypergcn"`, `"fast"`, or `"one"` for Algorithms 1–3.

- hidden:

  Hidden width (default 128).

- epochs:

  Training epochs (default 200).

- lr, weight_decay:

  Adam learning rate and L2 penalty.

- dropout:

  Hidden-layer dropout.

- validation:

  Fraction of labeled nodes held out to select an epoch.

- edge_weights:

  Optional positive hyperedge weights.

- seed:

  Reproducible R and torch seed.

- verbose:

  Message loss every 20 epochs.

## Value

A prediction data.frame in the same format as
[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md),
with training history attached.

## References

Yadati, N., Nimishakavi, M., Yadav, P., Nitin, V., Louis, A., &
Talukdar, P. (2019). HyperGCN: A new method of training graph
convolutional networks on hypergraphs. *NeurIPS 32*.
