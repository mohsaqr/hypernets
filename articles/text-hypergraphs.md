# Hypergraph analysis of a text corpus

A corpus of abstracts shares vocabulary unevenly: a word like
*telescope* binds a few documents tightly, while *challenges* touches
half the corpus weakly. A weighted hypergraph represents that structure
directly — documents are vertices and every word is a hyperedge over the
documents it occurs in, weighted by tf-idf. Two questions drive this
vignette, asked of the bundled `covid_abstracts` corpus (165 COVID-19
education-research abstracts, 2020–2024): do the abstracts organize into
coherent groups purely through shared vocabulary, and can a handful of
labeled abstracts classify the rest?

The methods are Zhou, Huang & Schölkopf’s (2006) hypergraph Laplacian
and Hayashi et al.’s (2020) random-walk Laplacian with edge-dependent
vertex weights, the spectral engines of honets’ hypergraph family
(HyperNetX-parity verified), applied to a text hypergraph built by
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md).

## From corpus to hypergraph

[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
tokenizes deterministically in base R, counts document–word occurrences,
and hands the weighted incidence structure to
[`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md).
Function words and a few corpus-boilerplate terms are excluded, and
words in fewer than three abstracts are dropped — a word seen once
cannot connect anything.

``` r

stops <- c(
  stop_words_en(),
  "using", "used", "use", "based", "results", "study", "research",
  "paper", "article", "findings", "data"
)

hg <- text_hypergraph(
  covid_abstracts,
  column = "abstract",
  id = "doc",
  weight = "tfidf",
  stop_words = stops,
  min_count = 3L
)
hg
#> Text hypergraph: 165 documents, 1453 words (documents as nodes, weight = tfidf)
#> Hyperedges: 1453 (words); sizes 1-140, median 4
```

The corpus becomes 165 vertices connected by 1,453 word-hyperedges. The
median word touches 4 abstracts; the largest hyperedge spans 140 of the
165 — a near-universal term that tf-idf weighting correspondingly
down-weights.

``` r

hg_measures(hg, what = "summary")
#>                  measure        value
#> 1                n_nodes 1.650000e+02
#> 2           n_hyperedges 1.453000e+03
#> 3                density 4.330434e-02
#> 4          avg_edge_size 7.145217e+00
#> 5 pairwise_participation 9.957132e-01
```

A pairwise participation of 0.996 says nearly every pair of abstracts
shares at least one word: the hypergraph is one dense component, so
spectral methods apply to the whole corpus at once.

## Which words hold the corpus together?

Reversing the orientation makes words the vertices — two words are close
when they co-occur in many documents — and centrality then ranks the
corpus’s organizing vocabulary.

``` r

word_hg <- text_hypergraph(
  covid_abstracts,
  column = "abstract",
  id = "doc",
  nodes = "word",
  weight = "tfidf",
  stop_words = stops,
  min_count = 3L
)
hg_centrality(word_hg, type = "clique", sort_by = "clique", n = 10)
#>           node    clique
#> 1        covid 0.2798710
#> 2     pandemic 0.2707863
#> 3    education 0.2491869
#> 4     students 0.1920792
#> 5     learning 0.1906239
#> 6       online 0.1705769
#> 7     teaching 0.1651923
#> 8  educational 0.1158119
#> 9     teachers 0.1084348
#> 10  challenges 0.1079189
```

The ranking is the corpus in miniature: *covid*, *pandemic*,
*education*, *students*, *learning*, *online*, *teaching*. This is
remote-education research, and the centrality ranking recovers that
without reading a single abstract.

## Do the abstracts cluster?

[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
with `type = "random_walk"` is the Hayashi et al. (2020) pipeline: the
tf-idf weights become edge-dependent vertex weights of a random walk,
whose normalized Laplacian is embedded and clustered. The number of
clusters is explicit by design; four is the working choice here and is
checked below.

``` r

clusters <- hg_cluster(hg, k = 4, type = "random_walk", seed = 1)
aggregate(node ~ cluster, data = clusters, FUN = length)
#>     cluster node
#> 1 Cluster 1   16
#> 2 Cluster 2   44
#> 3 Cluster 3   29
#> 4 Cluster 4   76
```

``` r

merged <- merge(clusters, covid_abstracts, by.x = "node", by.y = "doc")
aggregate(year ~ cluster, data = merged, FUN = mean)
#>     cluster     year
#> 1 Cluster 1 2021.062
#> 2 Cluster 2 2021.682
#> 3 Cluster 3 2021.345
#> 4 Cluster 4 2021.711
```

The partition is 16/44/29/76, and every cluster’s mean publication year
sits between 2021.1 and 2021.7: the clusters are thematic, not temporal.
The vocabulary of this literature separates by topic while staying
stable across the pandemic years.

Is `k = 4` doing real work? Refitting at `k = 3` answers the sensitivity
question:

``` r

clusters_3 <- hg_cluster(hg, k = 3, type = "random_walk", seed = 1)
hg_agreement(clusters, clusters_3, what = "table")
#>     label_x   label_y  n
#> 1 Cluster 1 Cluster 1 15
#> 2 Cluster 1 Cluster 3  1
#> 3 Cluster 2 Cluster 1  1
#> 4 Cluster 2 Cluster 2 37
#> 5 Cluster 2 Cluster 3  6
#> 6 Cluster 3 Cluster 3 29
#> 7 Cluster 4 Cluster 1 22
#> 8 Cluster 4 Cluster 2 36
#> 9 Cluster 4 Cluster 3 18
```

The two compact themes survive the change almost intact (15 of 16, and
all 29, stay together), while the largest cluster redistributes. The
small clusters are stable structure worth interpreting; the large one is
the corpus’s residual mass, and conclusions should not lean on its
boundary.

## Classifying from six labels

Transductive label spreading (Zhou et al. 2006) takes labels for a few
vertices and propagates them over the hyperedges. Labeling the three
alphabetically first abstracts of 2020 as `"early"` and of 2024 as
`"late"` asks whether pandemic-era vocabulary dates an abstract:

``` r

labels <- c(
  "2-s2.0-85085897904" = "early",
  "2-s2.0-85087124472" = "early",
  "2-s2.0-85087646998" = "early",
  "2-s2.0-85107489479" = "late",
  "2-s2.0-85112557152" = "late",
  "2-s2.0-85113355756" = "late"
)
predictions <- hg_classify(hg, labels = labels)
predicted_years <- merge(predictions, covid_abstracts,
                         by.x = "node", by.y = "doc")
aggregate(year ~ predicted, data = predicted_years, FUN = mean)
#>   predicted     year
#> 1     early 2020.000
#> 2      late 2021.605
```

``` r

hg_agreement(predictions, covid_abstracts,
             node = c("node", "doc"), label = c("predicted", "year"),
             what = "table")
#>   label_x label_y  n
#> 1   early    2020  3
#> 2    late    2020 37
#> 3    late    2021 40
#> 4    late    2022 40
#> 5    late    2023 40
#> 6    late    2024  5
```

The classifier is conservative and correct where it commits: every
abstract it calls `"early"` is genuinely from 2020, but it finds only
three of them — the rest of the corpus, 2020 included, reads as
`"late"`. The decisive reading: this literature’s vocabulary converged
almost immediately, so publication year is weakly encoded in word
choice. Label spreading here is trustworthy for what it asserts and
honest about how little the structure supports.

## Ranking, dualizing, and testing the structure

Three verbs complete the analysis layer.
[`hg_pagerank()`](https://mohsaqr.github.io/hypernets/reference/hg_pagerank.md)
ranks vertices by the stationary importance of the Chitra & Raphael
(2019) random walk with edge-dependent vertex weights – the walk that
actually uses the tf-idf weights (verified against the HyperNetX
reference implementation at machine precision):

``` r

hg_pagerank(hg, sort_by = "pagerank", n = 3)
#>                 node   pagerank
#> 1 2-s2.0-85124537994 0.01648602
#> 2 2-s2.0-85116687295 0.01258984
#> 3 2-s2.0-85115059168 0.01231774
```

The top-ranked abstracts are the corpus’s connective tissue – documents
whose vocabulary reaches everywhere. A `personalized =` vector turns the
same verb into “importance from the perspective of these documents”.
[`dual_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/dual_hypergraph.md)
swaps the vertex and hyperedge roles of a fitted object, so the
word-level view never requires re-tokenizing:

``` r

dual_hypergraph(hg)
#> Text hypergraph: 165 documents, 1453 words (words as nodes, weight = tfidf)
#> Hyperedges: 165 (documents); sizes 22-121, median 63
```

Finally,
[`hg_null_test()`](https://mohsaqr.github.io/hypernets/reference/hg_null_test.md)
asks whether observed structure exceeds what the degree sequences alone
imply: it rewires the binary membership by degree-preserving
checkerboard swaps (Gotelli 2000) and compares. On the 2020 abstracts:

``` r

abstracts_2020 <- text_hypergraph(
  subset(covid_abstracts, year == 2020),
  column = "abstract",
  id = "doc",
  weight = "tfidf",
  stop_words = stops,
  min_count = 3L
)
hg_null_test(
  abstracts_2020,
  statistic = c("pairwise_participation", "avg_jaccard"),
  n = 99,
  seed = 1
)
#>                statistic   observed  null_mean    null_lo    null_hi         z
#> 1 pairwise_participation 0.98846154 0.99109039 0.98647436 0.99487179 -1.045800
#> 2            avg_jaccard 0.05659397 0.05583251 0.05519231 0.05634021  2.483364
#>   p_value  n method
#> 1    0.36 99   swap
#> 2    0.01 99   swap
```

Read the two rows together: near-universal pairwise participation
(0.988) is *not* evidence of structure – the null reproduces it (p =
0.36), because vocabulary-rich abstracts connect to everything under any
degree-preserving rewiring. The mean edge Jaccard *is* significant (z =
2.48, p = 0.01): words repeat each other’s document combinations more
than their frequencies require, which is exactly the topical
co-occurrence the clustering exploits. The decisive practice: never
report a structural statistic without its null interval.

## Using an existing quanteda or tidytext pipeline

[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
needs no NLP dependency, but any long document–word table is already the
incidence structure, so existing pipelines plug into the same engines
directly:

``` r

# quanteda: tidytext::tidy() turns a dfm into a long table with columns
# document, term, count (verified against quanteda 4.x / tidytext 0.4.x)
group_hypergraph(
  tidytext::tidy(my_dfm),
  actor = "term", group = "document", weight = "count"
)

# tidytext: word counts are already long (document, word, n)
group_hypergraph(
  my_word_counts,
  actor = "word", group = "document", weight = "n"
)
```

## When to use which

- **[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)**
  — no labels, discover thematic groups; prefer `type = "random_walk"`
  whenever the weights are tf-idf (edge-dependent vertex weights are
  exactly its assumption).
- **[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md)**
  — a few trusted labels, classify the rest; read `margin` as the
  confidence and expect conservative behavior when the vocabulary signal
  is weak.
- **[`hg_centrality()`](https://mohsaqr.github.io/hypernets/reference/hg_centrality.md)**
  — rank organizing vocabulary (words as nodes) or hub documents
  (documents as nodes).
- **[`hg_measures()`](https://mohsaqr.github.io/hypernets/reference/hg_measures.md)**
  — the structural facts (sizes, overlaps, density) that say whether
  spectral analysis is even applicable; check connectivity-adjacent
  quantities before clustering.

## References

Hayashi, K., Aksoy, S. G., Park, C. H., & Park, H. (2020). Hypergraph
random walks, Laplacians, and clustering. *CIKM 2020*.
<doi:10.1145/3340531.3412034>.

Zhou, D., Huang, J., & Schölkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 2006*.
