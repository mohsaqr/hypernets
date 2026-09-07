# Tidy tables of a text hypergraph

Tidy tables of a text hypergraph

## Usage

``` r
# S3 method for class 'text_hypergraph'
as.data.frame(
  x,
  row.names = NULL,
  optional = FALSE,
  what = c("weights", "documents", "vocabulary", "sentences"),
  ...
)
```

## Arguments

- x:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  object.

- row.names, optional:

  Ignored; present for S3 consistency.

- what:

  Which table: `"weights"` (default) – for the bag construction one row
  per document-word pair (`doc`, `word`, `n`, `weight`); for the window
  construction one row per window-content/word membership (`edge`,
  `word`, `weight` = window count); for the sentence construction one
  row per sentence/word membership (`edge`, `word`, `weight` =
  in-sentence count); for the knn construction one row per hyperedge
  membership (`doc`, `edge`, `weight` = cosine similarity).
  `"sentences"` (sentence construction only) gives one row per sentence
  hyperedge (`edge`, `doc`, `sentence`, `n_tokens`). `"documents"` gives
  one row per document (with `n_tokens`/`n_types` for token-based
  constructions, plus any metadata columns carried from the input).
  `"vocabulary"` gives one row per word (`word`, `n`, `doc_freq`, and
  `idf` under tf-idf weighting; empty for the knn construction, which
  has no token layer).

- ...:

  Unused.

## Value

A base `data.frame` as described under `what`.

## Examples

``` r
hg <- text_hypergraph(c(a = "salt and soup", b = "soup and stars"))
as.data.frame(hg)
#>   doc  word n weight
#> 1   a   and 1      1
#> 2   a  salt 1      1
#> 3   a  soup 1      1
#> 4   b   and 1      1
#> 5   b  soup 1      1
#> 6   b stars 1      1
as.data.frame(hg, what = "documents")
#>   doc n_tokens n_types
#> 1   a        3       3
#> 2   b        3       3
```
