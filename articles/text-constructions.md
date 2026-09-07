# Windowed and embedding hypergraphs

The bag-of-words hypergraph of
[`vignette("text-hypergraphs")`](https://mohsaqr.github.io/hypernets/articles/text-hypergraphs.md)
ignores two things a corpus knows: the *order* of words within a
document and the *position* of documents in an embedding space. Two
further constructions use them. The question for each is the same: does
the added structure change a conclusion the bag construction would have
reached?

## Windowed hyperedges: word order as structure

`construction = "window"` slides a token window over each document
(never crossing document boundaries); each distinct window content — the
set of words co-occurring within `window` tokens of each other — is one
hyperedge, weighted by its window count. This is the sequential
construction of Ding et al. (2020, Sec. 3.2), extended with weights —
the text front end of
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
with which it shares its incidence matrix (a shipped package test).

The 2020 abstracts make a compact testbed:

``` r

stops <- c(
  stop_words_en(),
  "using", "used", "use", "based", "results", "study", "research",
  "paper", "article", "findings", "data"
)
early <- subset(covid_abstracts, year == 2020)

windowed <- text_hypergraph(
  early,
  column = "abstract",
  id = "doc",
  construction = "window",
  window = 3,
  stop_words = stops,
  min_count = 3L
)
windowed
#> Text hypergraph: 40 documents, 349 words (windowed hyperedges: w = 3, sliding, 2000 windows)
#> Hyperedges: 1880 (distinct windows); sizes 2-3, median 3
```

``` r

hg_centrality(windowed, type = "clique", sort_by = "clique", n = 10)
#>          node    clique
#> 1   education 0.4174717
#> 2       covid 0.2953234
#> 3    distance 0.2901158
#> 4    pandemic 0.2835453
#> 5    learning 0.2741282
#> 6    students 0.2405080
#> 7      online 0.2318112
#> 8    teachers 0.1884582
#> 9    teaching 0.1434927
#> 10 university 0.1113114
```

Compare the same corpus as a bag-of-words hypergraph over words:

``` r

bag_words <- text_hypergraph(
  early,
  column = "abstract",
  id = "doc",
  nodes = "word",
  weight = "tfidf",
  stop_words = stops,
  min_count = 3L
)
hg_centrality(bag_words, type = "clique", sort_by = "clique", n = 10)
#>           node    clique
#> 1        covid 0.2792762
#> 2     pandemic 0.2660941
#> 3    education 0.2569350
#> 4     students 0.2173556
#> 5     learning 0.1904566
#> 6     teachers 0.1463876
#> 7     teaching 0.1459636
#> 8       online 0.1355434
#> 9   university 0.1237287
#> 10 educational 0.1141029
```

The rankings agree on the theme words — but the windowed construction
puts **distance** third, and the bag construction does not rank it in
the top ten at all. Adjacency is doing real work: *distance* matters in
this corpus almost exclusively inside the phrases *distance learning*
and *distance education*, which a window sees and a bag cannot. The
window size is an arbitrary choice, so it gets a sensitivity check:

``` r

windowed_4 <- text_hypergraph(
  early,
  column = "abstract",
  id = "doc",
  construction = "window",
  window = 4,
  stop_words = stops,
  min_count = 3L
)
hg_centrality(windowed_4, type = "clique", sort_by = "clique", n = 10)
#>          node    clique
#> 1   education 0.4314610
#> 2    distance 0.3052559
#> 3       covid 0.2971267
#> 4    learning 0.2805257
#> 5    pandemic 0.2689772
#> 6    students 0.2455921
#> 7      online 0.2170722
#> 8    teachers 0.2044145
#> 9    teaching 0.1417211
#> 10 university 0.1085232
```

At `window = 4`, *distance* rises to second. The conclusion — adjacency
promotes phrase vocabulary the bag misses — survives the window choice.

## kNN embedding hyperedges: semantic space as structure

`construction = "knn"` places each document in an embedding space and
forms one hyperedge per document: itself plus its `k` nearest neighbors,
weighted by cosine similarity. The bundled `covid_embeddings` matrix
holds sbert embeddings of the `covid_abstracts` corpus (pinned
`all-MiniLM-L6-v2`, L2-normalized, built by
`data-raw/covid_embeddings.R`), so the pipeline runs offline; with the
sbert package installed, omitting `embeddings` encodes the text
directly.

``` r

knn_hg <- text_hypergraph(
  covid_abstracts,
  column = "abstract",
  id = "doc",
  construction = "knn",
  k = 10,
  embeddings = covid_embeddings
)
knn_hg
#> Text hypergraph: 165 documents (kNN embedding hyperedges: k = 10, cosine)
#> Hyperedges: 165 (kNN neighborhoods); sizes 11-11, median 11
```

``` r

knn_clusters <- hg_cluster(knn_hg, k = 4, type = "random_walk", seed = 1)
aggregate(node ~ cluster, data = knn_clusters, FUN = length)
#>     cluster node
#> 1 Cluster 1   41
#> 2 Cluster 2   41
#> 3 Cluster 3   37
#> 4 Cluster 4   46
```

The embedding hypergraph partitions the corpus into balanced groups
(41/41/37/46), where the bag-of-words clustering of
[`vignette("text-hypergraphs")`](https://mohsaqr.github.io/hypernets/articles/text-hypergraphs.md)
produced 16/44/29/76 — one dominant residual cluster. Semantic
neighborhoods spread structure more evenly than shared vocabulary does.
The two views agree only partially:

``` r

bag <- text_hypergraph(
  covid_abstracts,
  column = "abstract",
  id = "doc",
  weight = "tfidf",
  stop_words = stops,
  min_count = 3L
)
bag_clusters <- hg_cluster(bag, k = 4, type = "random_walk", seed = 1)
hg_agreement(knn_clusters, bag_clusters, what = "table")
#>      label_x   label_y  n
#> 1  Cluster 1 Cluster 1  9
#> 2  Cluster 1 Cluster 2  4
#> 3  Cluster 1 Cluster 3  8
#> 4  Cluster 1 Cluster 4 20
#> 5  Cluster 2 Cluster 1  1
#> 6  Cluster 2 Cluster 2 32
#> 7  Cluster 2 Cluster 3  2
#> 8  Cluster 2 Cluster 4  6
#> 9  Cluster 3 Cluster 1  2
#> 10 Cluster 3 Cluster 2  5
#> 11 Cluster 3 Cluster 3 17
#> 12 Cluster 3 Cluster 4 13
#> 13 Cluster 4 Cluster 1  4
#> 14 Cluster 4 Cluster 2  3
#> 15 Cluster 4 Cluster 3  2
#> 16 Cluster 4 Cluster 4 37
```

Two themes are found by both constructions (32 and 37 documents stay
together); the rest reshuffles. Neither view is “the” truth — they
answer different questions (shared vocabulary vs. semantic proximity) —
but a grouping stable across both is worth trusting.

`k` is the consequential arbitrary choice here, so it gets the
sensitivity check:

``` r

knn_5_hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
                            construction = "knn", k = 5,
                            embeddings = covid_embeddings)
knn_5 <- hg_cluster(knn_5_hg, k = 4, type = "random_walk", seed = 1)
knn_15_hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
                             construction = "knn", k = 15,
                             embeddings = covid_embeddings)
knn_15 <- hg_cluster(knn_15_hg, k = 4, type = "random_walk", seed = 1)
hg_agreement(knn_clusters, knn_15)
#>     n agreement aligned     ari
#> 1 165 0.8666667     143 0.67443
hg_agreement(knn_clusters, knn_5)
#>     n agreement aligned       ari
#> 1 165 0.7575758     125 0.5305385
```

The `aligned` column counts the documents that stay with the majority of
their cluster. Against `k = 15`, 143 of 165 documents keep their
cluster’s majority alignment; against `k = 5`, 125 do. The partition is
robust upward and softer downward — small neighborhoods fragment one of
the four groups. Report conclusions at the k where they are stable, and
say so.

One honesty note: everything in this section is conditional on the
encoder. `covid_embeddings` comes from one pinned model; a different
encoder yields a different geometry, and any claim that matters should
be re-checked with a second model
(`text_hypergraph(construction = "knn", model = )` makes that a
one-argument change when sbert is installed).

## When to use which construction

- **`"bag"`** — the default; shared vocabulary, interpretable hyperedges
  (words or documents), tf-idf weighting, the full
  clustering/transduction story of
  [`vignette("text-hypergraphs")`](https://mohsaqr.github.io/hypernets/articles/text-hypergraphs.md).
- **`"window"`** — when word order carries the signal: phrase
  vocabulary, collocations, sequence-like corpora. Vertices are words;
  use it with
  [`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md)
  and
  [`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md).
  Prefer modest corpora or `min_count` filtering — distinct windows
  multiply quickly.
- **`"knn"`** — when meaning matters more than vocabulary overlap:
  paraphrase-heavy or multilingual-adjacent corpora, or as the explicit,
  spectral alternative to UMAP+HDBSCAN pipelines. Vertices are
  documents; every hyperedge has `k + 1` members; the geometry is the
  encoder’s.

## References

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*. <doi:10.18653/v1/2020.emnlp-main.399>.

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*.
<doi:10.1145/3340531.3412034>.
