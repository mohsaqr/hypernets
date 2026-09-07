# Build a weighted hypergraph from a text corpus

Tokenizes a corpus (base R, deterministic) and builds a weighted
hypergraph (a
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
object with a text layer), in one of three constructions:

## Usage

``` r
text_hypergraph(
  x,
  column = NULL,
  id = NULL,
  construction = c("bag", "window", "sentence", "knn"),
  nodes = c("doc", "word"),
  weight = c("n", "tfidf"),
  stop_words = NULL,
  min_count = 1L,
  lowercase = TRUE,
  window = 3L,
  window_mode = c("sliding", "tumbling"),
  k = 10L,
  embeddings = NULL,
  model = NULL,
  sparse = NULL
)
```

## Arguments

- x:

  A character vector of documents, or a `data.frame` containing a text
  column.

- column:

  Name of the text column when `x` is a `data.frame`.

- id:

  Optional name of an ID column when `x` is a `data.frame`; its values
  (unique, non-missing) become the document identifiers. Defaults to the
  names of `x` when it is a named character vector, otherwise `"doc_1"`,
  `"doc_2"`, ...

- construction:

  `"bag"` (default), `"window"`, `"sentence"` or `"knn"` – see Details.

- nodes:

  Which entity is the vertex set for the bag construction: `"doc"`
  (default) or `"word"`. Ignored by `"window"` (vertices are words) and
  `"knn"` (vertices are documents).

- weight:

  Term weighting for the bag construction: `"n"` (raw count, default) or
  `"tfidf"`. The window construction always weights by window counts,
  the sentence construction by in-sentence counts, the knn construction
  by cosine similarity.

- stop_words:

  Optional character vector of words to drop after tokenization
  (compared after lowercasing when `lowercase = TRUE`); see
  [`stop_words_en()`](https://mohsaqr.github.io/hypernets/reference/stop_words_en.md).
  Not applicable to `"knn"`.

- min_count:

  Minimum total corpus count for a word to be kept (default `1L`, keep
  everything). Not applicable to `"knn"`.

- lowercase:

  Lowercase the text before tokenization (default `TRUE`).

- window:

  Window size in tokens for `construction = "window"` (default `3L`).

- window_mode:

  `"sliding"` (default) or `"tumbling"`, for `construction = "window"`.

- k:

  Number of nearest neighbors per hyperedge for `construction = "knn"`
  (default `10L`).

- embeddings:

  Numeric matrix of document embeddings for `construction = "knn"`: one
  row per document, either rownames matching the document IDs or rows in
  document order. `NULL` (default) encodes the text with `sbert`.

- model:

  Passed to `sbert::encode()` when embeddings are computed (`NULL` =
  sbert's default model).

- sparse:

  Store the incidence as a `Matrix::dgCMatrix` (bag and sentence
  constructions). Default `NULL` chooses: sparse when the
  document-by-word incidence would have a million cells or more, dense
  otherwise. Sparse hypergraphs scale to tens of thousands of documents;
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
  [`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md),
  [`hg_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md),
  and
  [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md)
  use sparse operator paths that agree with the dense engines (tested),
  while tensor centralities and the null test currently require the
  dense representation.

## Value

An object of class `c("text_hypergraph", "net_hypergraph")` – a
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md)
hypergraph accepted by every hypernets hypergraph verb and by
[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md),
[`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md),
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md),
and
[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md)
– with a `text` field recording the corpus tables. Use
[`as.data.frame.text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.text_hypergraph.md)
for the tidy weight table, and its `what` argument for the document and
vocabulary tables.

## Details

- `construction = "bag"` (default): the document-word incidence. With
  `nodes = "doc"`, documents are vertices and each word is a hyperedge
  over the documents containing it (the orientation for document
  clustering and transductive classification, Hayashi et al. 2020); with
  `nodes = "word"` the orientation reverses (the HyperGAT
  sentence-as-hyperedge setup generalized to whole documents).

- `construction = "window"`: words are vertices and every token window
  is a hyperedge – the sliding-window sequential construction of Ding et
  al. (2020, Sec. 3.2), extended with weights: each hyperedge is a
  distinct window content (the sorted set of words co-occurring in a
  window) and its incidence weight is the number of windows with that
  content. Windows never cross document boundaries.
  `window_mode = "sliding"` moves one token at a time (a document
  shorter than `window` forms one whole-document window); `"tumbling"`
  uses consecutive chunks, trailing partial chunk included. This is the
  text front end of
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md):
  on token sequences with no repeated word inside a window, the two
  constructions produce the same incidence matrix (tested), with the
  hyperedges here named by their word content. They differ only at the
  edges of the definition – a document shorter than `window`, a trailing
  partial chunk, or a word repeated inside one window, which a
  set-valued hyperedge collapses.

- `construction = "sentence"`: words are vertices and every sentence is
  a hyperedge – HyperGAT's sentence-as-hyperedge construction (Ding et
  al. 2020, Sec. 3.2) over the whole corpus. Documents are split at `.`,
  `!`, `?` and `;`; each sentence with at least one surviving token is
  one hyperedge named `<doc>#<sentence>`, and its incidence weight is
  the word's occurrence count in that sentence. Two words are then
  related when they share a sentence, not merely a document – the scope
  [`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  uses for `type = "sentence_centrality"`.

- `construction = "knn"`: documents are vertices and each document plus
  its `k` nearest neighbors in an embedding space is one hyperedge,
  weighted by cosine similarity (see
  [`knn_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/knn_hypergraph.md)).
  Pass a precomputed `embeddings` matrix, or leave it `NULL` to encode
  the text with the `sbert` package (if installed; models download only
  on explicit user confirmation, per sbert's policy).

Tokens are maximal runs of letters (with internal apostrophes); digits
and punctuation are separators. `weight = "tfidf"` (bag construction
only) uses the smoothed formula `tf * (log((1 + N) / (1 + df)) + 1)`
with `tf` the raw count, `N` the number of documents, and `df` the
word's document frequency – the `smooth_idf` variant of Manning,
Raghavan and Schutze (2008) as popularized by scikit-learn, which is
never zero and so never silently disconnects a vertex. `stop_words` and
`min_count` filtering happen before windowing, so removing a token
closes the gap it leaves in the sequence.

## Conditions

Raises `hypernets_bad_input` (broken argument contract, including
bag-only arguments passed to other constructions),
`hypernets_empty_corpus` (no document survives tokenization and
filtering), `hypernets_missing_embeddings` (`construction = "knn"` with
neither `embeddings` nor the sbert package), and warns with
`hypernets_dropped_documents` when some documents end up empty.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*.
[doi:10.18653/v1/2020.emnlp-main.399](https://doi.org/10.18653/v1/2020.emnlp-main.399)

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*.
[doi:10.1145/3340531.3412034](https://doi.org/10.1145/3340531.3412034)

Manning, C. D., Raghavan, P., & Schutze, H. (2008). *Introduction to
Information Retrieval*. Cambridge University Press.

## Examples

``` r
corpus <- c(
  cooking_1 = "Simmer the soup with onions and carrots",
  cooking_2 = "This soup recipe needs a pinch of salt",
  space_1   = "The telescope revealed a distant galaxy",
  space_2   = "Astronomers aimed the telescope at the night sky"
)
hg <- text_hypergraph(corpus, weight = "tfidf")
hg
#> Text hypergraph: 4 documents, 23 words (documents as nodes, weight = tfidf)
#> Hyperedges: 23 (words); sizes 1-3, median 1
as.data.frame(hg)
#>          doc        word n   weight
#> 1  cooking_1         and 1 1.916291
#> 2  cooking_1     carrots 1 1.916291
#> 3  cooking_1      onions 1 1.916291
#> 4  cooking_1      simmer 1 1.916291
#> 5  cooking_1        soup 1 1.510826
#> 6  cooking_1         the 1 1.223144
#> 7  cooking_1        with 1 1.916291
#> 8  cooking_2           a 1 1.510826
#> 9  cooking_2       needs 1 1.916291
#> 10 cooking_2          of 1 1.916291
#> 11 cooking_2       pinch 1 1.916291
#> 12 cooking_2      recipe 1 1.916291
#> 13 cooking_2        salt 1 1.916291
#> 14 cooking_2        soup 1 1.510826
#> 15 cooking_2        this 1 1.916291
#> 16   space_1           a 1 1.510826
#> 17   space_1     distant 1 1.916291
#> 18   space_1      galaxy 1 1.916291
#> 19   space_1    revealed 1 1.916291
#> 20   space_1   telescope 1 1.510826
#> 21   space_1         the 1 1.223144
#> 22   space_2       aimed 1 1.916291
#> 23   space_2 astronomers 1 1.916291
#> 24   space_2          at 1 1.916291
#> 25   space_2       night 1 1.916291
#> 26   space_2         sky 1 1.916291
#> 27   space_2   telescope 1 1.510826
#> 28   space_2         the 2 2.446287

win <- text_hypergraph(corpus, construction = "window", window = 3)
win
#> Text hypergraph: 4 documents, 23 words (windowed hyperedges: w = 3, sliding, 21 windows)
#> Hyperedges: 20 (distinct windows); sizes 3-3, median 3

# knn from precomputed embeddings (sbert-free, offline)
emb <- matrix(c(1, 0,  0.9, 0.1,  0, 1,  0.1, 0.9),
              nrow = 4, byrow = TRUE,
              dimnames = list(names(corpus), NULL))
text_hypergraph(corpus, construction = "knn", k = 1, embeddings = emb)
#> Text hypergraph: 4 documents (kNN embedding hyperedges: k = 1, cosine)
#> Hyperedges: 4 (kNN neighborhoods); sizes 2-2, median 2
```
