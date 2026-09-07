# HyperGAT document classifier

Trains the dual-attention hypergraph network of Ding et al. (2020) –
every document becomes its own hypergraph (its unique words as vertices,
its sentences as hyperedges), two attention layers aggregate words into
sentence and optional LDA-topic hyperedges and those hyperedges back
into words, and a masked mean pool feeds a linear classifier. Inductive:
only labeled documents are trained on, and every document (labeled or
not) is scored. Needs the suggested torch package. Semantics follow the
official implementation. `semantic = "none"` reproduces its
sentence-only ("w/o semantic") ablation. `semantic = "lda"` adds the
full paper path: online variational-Bayes LDA is fitted to labeled
documents only, with the topic count defaulting to the number of
classes, and each document receives one edge per topic containing the
topic's top words present in that document.

## Usage

``` r
hg_hypergat(
  x,
  labels,
  column = NULL,
  id = NULL,
  stop_words = stop_words_en(),
  min_count = 1L,
  lowercase = TRUE,
  semantic = c("none", "lda"),
  lda_keywords = NULL,
  lda_topics = NULL,
  lda_top_n = 10L,
  lda_max_iter = 10L,
  lda_batch_size = 128L,
  lda_offset = 50,
  lda_decay = 0.7,
  lda_seed = 0L,
  embed_dim = 300L,
  hidden = 100L,
  epochs = 10L,
  lr = 0.001,
  dropout = 0.3,
  batch_size = 8L,
  weight_decay = 1e-06,
  lr_decay = 0.1,
  lr_step = 3L,
  validation = 0.1,
  class_weights = c("balanced", "none"),
  embeddings = NULL,
  seed = 1L,
  verbose = FALSE,
  what = c("predictions", "attention")
)

text_hypergat(
  x,
  labels,
  column = NULL,
  id = NULL,
  stop_words = stop_words_en(),
  min_count = 1L,
  lowercase = TRUE,
  semantic = c("none", "lda"),
  lda_keywords = NULL,
  lda_topics = NULL,
  lda_top_n = 10L,
  lda_max_iter = 10L,
  lda_batch_size = 128L,
  lda_offset = 50,
  lda_decay = 0.7,
  lda_seed = 0L,
  embed_dim = 300L,
  hidden = 100L,
  epochs = 10L,
  lr = 0.001,
  dropout = 0.3,
  batch_size = 8L,
  weight_decay = 1e-06,
  lr_decay = 0.1,
  lr_step = 3L,
  validation = 0.1,
  class_weights = c("balanced", "none"),
  embeddings = NULL,
  seed = 1L,
  verbose = FALSE,
  what = c("predictions", "attention")
)
```

## Arguments

- x:

  A character vector of documents (names become ids) or a data.frame
  with a text column.

- labels:

  The known labels: a named character vector (names are document ids,
  values class labels) or a tidy data.frame with a `node` column and a
  `label`, `cluster` or `predicted` column. At least two classes.

- column, id:

  When `x` is a data.frame: the text column and the optional id column,
  as in
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md).

- stop_words:

  Words removed before building sentences (default
  [`stop_words_en()`](https://mohsaqr.github.io/hypernets/reference/stop_words_en.md)).

- min_count:

  Minimum corpus frequency for a word to become a vertex.

- lowercase:

  Lowercase the text first.

- semantic:

  `"none"` (default) for sentence hyperedges only, or `"lda"` to append
  the paper's semantic topic hyperedges.

- lda_keywords:

  Optional precomputed topic keywords: a list of character vectors, one
  per topic. With `semantic = "lda"`, `NULL` fits LDA natively;
  supplying the list bypasses fitting and reproduces a saved official
  preprocessing run.

- lda_topics:

  Number of LDA topics. `NULL` (default) uses the number of classes, as
  in the paper. Ignored when `lda_keywords` is supplied.

- lda_top_n:

  Number of highest-probability words per fitted topic (paper default
  10).

- lda_max_iter, lda_batch_size:

  Online variational-Bayes passes and minibatch size (scikit-learn
  defaults used by the official code: 10 and 128).

- lda_offset, lda_decay:

  Online learning schedule (official values 50 and 0.7).

- lda_seed:

  LDA initialization seed (official value 0), separate from the neural
  training `seed`.

- embed_dim, hidden:

  Embedding and hidden width (official defaults 300 and 100).

- epochs, lr, dropout, batch_size, weight_decay:

  Training hyperparameters (official defaults: 10, 0.001, 0.3, 8, 1e-6).

- lr_decay, lr_step:

  StepLR schedule: multiply the learning rate by `lr_decay` every
  `lr_step` epochs (official 0.1 every 3).

- validation:

  Fraction of the labeled documents held out (stratified) to pick the
  best epoch; `0` keeps the final weights.

- class_weights:

  `"balanced"` (default, inverse-frequency loss weights as in the
  official run script) or `"none"`.

- embeddings:

  Optional pretrained word-vector matrix (rownames are words,
  `embed_dim` columns); words not covered keep their random
  initialization.

- seed:

  Integer seed (R and torch); results are deterministic given a seed.

- verbose:

  Message the loss each epoch.

- what:

  `"predictions"` (default) returns the per-document classification
  table; `"attention"` returns the trained network's node-level
  attention per document and word, the input
  `hg_keywords(type = "attention")` takes.

## Value

With `what = "predictions"`, a base `data.frame`, one row per (kept)
document: `node`, `label` (the given label or `NA`), `predicted`,
`score` (softmax probability of the winning class), `margin` (winner
minus runner-up). The training history is attached as attribute
`"history"` (`epoch`, `loss`, `val_accuracy`). Attribute `"semantic"`
records the topic keywords and LDA settings.

With `what = "attention"`, a base `data.frame`, one row per (document,
word) pair: `node`, `word`, `attention` (the word's node-level attention
weights from the **first** attention layer, the one that attends over
the word embeddings, in evaluation mode, summed over the hyperedges it
belongs to – each hyperedge's weights sum to one, so a document's
`attention` column sums to its hyperedge count), `attention_2` (the same
from the second layer) and `n_edges` (how many of the document's
hyperedges contain the word). Ordered by `node`, then `word`. The second
layer attends over first-layer outputs, which are identical for every
word of a document that has a single hyperedge (one sentence), so
`attention_2` is uniform within such documents by construction;
`hg_keywords(type = "attention")` uses `attention`.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*.

## Examples

``` r
# \donttest{
if (requireNamespace("torch", quietly = TRUE)) {
  docs <- c(
    cooking_1 = "Simmer the soup. Add onions and carrots.",
    cooking_2 = "This soup recipe needs salt. Serve on a cold night.",
    space_1 = "The telescope revealed a galaxy. Stars everywhere.",
    space_2 = "Astronomers aimed the telescope. The stars were sharp."
  )
  hg_hypergat(docs, labels = c(cooking_1 = "cooking", space_1 = "space"),
              embed_dim = 16, hidden = 8, epochs = 5, validation = 0)
}
# }
```
