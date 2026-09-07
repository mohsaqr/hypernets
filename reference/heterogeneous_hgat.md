# Heterogeneous graph attention classifier (HGAT)

Implements Linmei et al.'s (2019) heterogeneous graph convolution and
dual-level attention. Node types retain separate input feature spaces
and transformation matrices. For every target node, type-level attention
first weighs the available neighboring information types; node-level
attention then weighs individual neighbors before their type-specific
transformed representations are combined. This is HGAT for heterogeneous
information networks, distinct from Ding et al.'s document-hypergraph
[`text_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md).

## Usage

``` r
heterogeneous_hgat(
  adjacency,
  node_types,
  features,
  labels,
  hidden = 128L,
  layers = 2L,
  epochs = 200L,
  lr = 0.005,
  weight_decay = 5e-04,
  dropout = 0.5,
  validation = 0.1,
  seed = 1L,
  verbose = FALSE
)
```

## Arguments

- adjacency:

  Symmetric non-negative adjacency matrix with identical unique row and
  column names.

- node_types:

  Named character vector assigning every node a type.

- features:

  Named list with one numeric feature matrix per type. Rows are named
  nodes of that type; feature widths may differ across types.

- labels:

  Named character labels for a subset of nodes.

- hidden:

  Common hidden dimension.

- layers:

  Number of dual-attention layers (paper default 2).

- epochs, lr, weight_decay, dropout, validation, seed, verbose:

  Training controls as in
  [`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md).

## Value

A prediction data.frame for every heterogeneous-network node, with
training history and node types attached as attributes.

## References

Linmei, H., Yang, T., Shi, C., Ji, H., & Li, X. (2019). Heterogeneous
graph attention networks for semi-supervised short text classification.
*EMNLP-IJCNLP 2019*, 4821–4830.
