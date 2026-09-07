# Hypergraph neural network classifier (HGNN)

Trains the two-layer hypergraph convolutional network of Feng et al.
(2019) on a hypergraph: each layer propagates linearly transformed
vertex features through \\G = D_v^{-1/2} H W D_e^{-1} H^T D_v^{-1/2}\\
(the Zhou similarity operator), with ReLU and dropout between layers and
a cross-entropy loss on the labeled vertices. Needs the suggested torch
package.

## Usage

``` r
hg_neural(
  hg,
  labels,
  features = "incidence",
  hidden = 128L,
  epochs = 600L,
  lr = 0.01,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  edge_weights = NULL,
  seed = 1L,
  verbose = FALSE
)

hypergraph_neural(
  hg,
  labels,
  features = "incidence",
  hidden = 128L,
  epochs = 600L,
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

- hidden:

  Width of the hidden layer.

- epochs:

  Training epochs (Adam).

- lr, weight_decay:

  Adam learning rate and L2 penalty. The defaults (`lr = 0.01`,
  `epochs = 600`) were selected by validation accuracy on held-out seeds
  for high-dimensional sparse text features; the original paper's
  `lr = 0.001`, 200 epochs underfits that regime.

- dropout:

  Dropout rate after the hidden layer.

- validation:

  Fraction of the labeled seeds held out (stratified) to pick the best
  epoch; `0` trains on all seeds for `epochs` and keeps the final
  weights.

- edge_weights:

  Optional positive hyperedge weights `W` (one per hyperedge); default
  all 1, the paper's initialization.

- seed:

  Integer seed for weight initialization, dropout and the validation
  split (torch and R RNGs); results are deterministic given a seed.

- verbose:

  If `TRUE`, message the loss every 20 epochs.

## Value

A base `data.frame`, one row per node, with columns `node`, `label` (the
given label or `NA`), `predicted`, `score` (softmax probability of the
winning class) and `margin` (winner minus runner-up probability). The
training history is attached as attribute `"history"` (data.frame:
`epoch`, `loss`, `val_accuracy`).

## References

Feng, Y., You, H., Zhang, Z., Ji, R., & Gao, Y. (2019). Hypergraph
neural networks. *AAAI 33*.

## Examples

``` r
# \donttest{
if (requireNamespace("torch", quietly = TRUE)) {
  hg <- text_hypergraph(c(
    cooking_1 = "simmer the soup with onions and carrots",
    cooking_2 = "this soup recipe needs salt on a cold night",
    space_1 = "the telescope revealed a distant galaxy and stars",
    space_2 = "astronomers aimed the telescope at the stars all night"
  ), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
  hg_neural(hg, labels = c(cooking_1 = "cooking", space_1 = "space"),
            hidden = 8, epochs = 50, validation = 0)
}
# }
```
