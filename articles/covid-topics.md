# Topic structure of the COVID-19 education literature as a document-word hypergraph

The text hypergraph pipeline in honets represents a corpus as a
document-word hypergraph and partitions it with a spectral method. This
vignette applies the pipeline to the `covid` corpus of the sbert
package, a Scopus export of 4,170 records on COVID-19 and education
published between 2020 and 2024. Each record has an abstract and a
publication year. The analysis proceeds in eight steps: text cleaning,
hypergraph construction, selection of the number of topics,
partitioning, topic description, topic sizes, topic quality, soft
membership, and the network of relations between topics.

## Text cleaning

[`clean_text()`](https://mohsaqr.github.io/hypernets/reference/clean_text.md)
removes the non-content that a bibliographic export carries in its
abstracts. It decodes HTML entities, strips tags, normalises typographic
characters, and removes bracketed citation numbers, list markers, URLs,
DOIs, trailing copyright notices and bare numbers. It takes a character
vector or a data frame with a text column and returns an object of the
same length or the same rows. A text whose share of alphabetic
characters falls below `min_content` is replaced by the empty string.
The `remove` argument takes additional patterns; here it names the
placeholder that Scopus inserts for a missing abstract.

``` r

abstracts <- clean_text(sbert::covid, column = "Abstract",
                        remove = "no abstract available", min_content = 0.5)
abstracts <- unique(abstracts)
```

The result is a data frame with the same columns as the input. The 323
placeholder rows and four low-content rows now hold empty abstracts, and
[`unique()`](https://rdrr.io/r/base/unique.html) has collapsed the
records that appear twice in the export. Empty abstracts are retained at
this stage so that the year column stays aligned with the text.

## Hypergraph construction

[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
builds a weighted hypergraph from a corpus. With the default
`construction = "bag"` and `nodes = "doc"`, each document is a node and
each word is a hyperedge that contains every document using the word.
With the default `weight = "n"`, the incidence cell of a document and a
word is the count of the word in the document. Tokens are lower-cased
alphabetic strings, so numbers do not enter the vocabulary. Documents
without any remaining token are dropped with a warning.

``` r

stops <- c(stop_words_en(), "covid", "pandemic")
hg <- text_hypergraph(abstracts, column = "Abstract", stop_words = stops,
                      min_count = 5L)
hg
```

The stop list contains the English function words of
[`stop_words_en()`](https://mohsaqr.github.io/hypernets/reference/stop_words_en.md)
and the two words that name the subject of every abstract; a word
present in every document does not separate any two documents. Words
used in fewer than five abstracts are excluded by `min_count`. The
incidence matrix has more than a million cells, so the default
`sparse = NULL` stores it as a sparse matrix.

## Number of topics

[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
with `what = "eigenvalues"` returns the leading eigenvalues of the
normalised hypergraph Laplacian (Zhou, Huang and Schölkopf, 2007) and
the gap after each. A large gap between consecutive eigenvalues
indicates a partition that the structure supports.

``` r

hg_cluster(hg, k = 18, seed = 1, what = "eigenvalues")
```

The gap after the second eigenvalue is the largest, followed by the gaps
after the third and fourth. The gaps after the fifth and later
eigenvalues are an order of magnitude smaller. The spectrum therefore
supports two or four clusters. The partition below uses sixteen
clusters, the number of topics in the structural topic model of this
literature reported by Saqr et al. (2023), so that the two sets of
topics can be compared. Sixteen is a finer partition than the spectrum
supports on its own; its justification is the vocabulary each cluster
owns, reported in the last section.

## Partition

[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
embeds the documents in the leading eigenvectors of the Laplacian and
partitions the embedding with k-means. With the default `type = "zhou"`
the Laplacian is built from hyperedge membership alone: each hyperedge
has unit weight and the incidence cells are not read. The function takes
the hypergraph, the number of clusters and a seed, and returns a data
frame with one row per document and its cluster label.

``` r

topics <- hg_cluster(hg, k = 16, seed = 1)
head(topics)
```

## Topic description

[`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
ranks the words of each cluster. It takes the hypergraph and the
partition and returns a data frame with one row per cluster and word,
containing the score, the word’s share and the number of cluster
documents that contain the word. The print method shows one row per
cluster with its words joined;
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) returns
the long form. The default score, `type = "mass"`, is the summed
incidence weight that a cluster’s documents place on the word, which
under `weight = "n"` is the word’s count within the cluster.

``` r

hg_keywords(hg, topics)
```

Under the default ranking every cluster lists the same words: learning,
students, education, online and teaching. These are the most frequent
words of the corpus, and they are the most frequent words of every
cluster as well. The `share` column, the fraction of a word’s corpus
count that falls in the cluster, separates the clusters.
`sort_by = "share"` ranks by that column, and `min_docs` sets a support
floor so that a word owned by a cluster only because it occurs in a few
of its abstracts does not reach the top. The `type` argument accepts
several scores at once and returns one block of rows per score.

``` r

words <- hg_keywords(hg, topics, type = c("frequency", "ctfidf", "centrality"),
                     sort_by = "share", min_docs = 15)
words
```

The three scores are the raw count of the word in the cluster
(`frequency`), the class-based tf-idf of Grootendorst (2022) (`ctfidf`),
and the PageRank of the word in the cluster’s own word hypergraph, in
which the cluster’s abstracts are the hyperedges (`centrality`; the
default measure). Ten words per cluster are returned by default. Ranked
by share, the three scores return the same owned vocabulary for each
cluster and differ in the order of the words. The sixteen clusters
correspond to distinct literatures: parents and childhood, infection
prevention, systematic reviews, lockdown policy, narrative essays,
messaging tools in Indonesian higher education, industry and recovery,
laboratory teaching, teacher agency, questionnaire validation, Spanish
methodological studies, technology acceptance models, dental and
surgical assessment, disadvantaged households, the spring campus
transition, and residency training.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) on the keyword
table draws one panel per cluster and score, with one bar per word.

``` r

plot(words, value = "share")
```

## Topic sizes

[`hg_topic_sizes()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_sizes.md)
counts the documents of each topic and their share of the clustered
documents. With `weights`, a numeric vector named by document, it also
reports the weighted size and share; a citation count or a repeat count
are typical weights. One row per topic.

``` r

sizes <- hg_topic_sizes(hg, topics)
sizes
plot(sizes)
```

## Topic coherence and exclusivity

[`hg_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
scores each topic’s owned vocabulary on two axes. Coherence is the mean
normalised pointwise mutual information between pairs of the topic’s top
words, computed on their co-occurrence within the topic’s documents: 1
when the words always occur together, 0 at independence, negative when
they avoid each other. Exclusivity is the mean share of the same words,
the fraction of each word’s corpus count that the topic holds. The top
words are the ten that
[`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
ranks by share with the support floor given.

``` r

quality <- hg_topic_quality(hg, topics, min_docs = 15)
quality
plot(quality)
```

Exclusivity is high for every topic because the words are selected by
share. Coherence separates them: the topics whose top words co-occur in
most of their documents describe a single kind of study; the topics with
low coherence collect abstracts that share a vocabulary but not a
subject.

## Soft membership

[`hg_membership()`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
turns the hard partition into a membership of every document in every
topic. Documents are placed in the spectral embedding that
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
cut, each topic’s centre is the mean position of its documents, and a
document’s membership in a topic is the fuzzy c-means weight with
fuzziness 2, the inverse squared distance to the centre normalised over
topics. A document at a centre has membership 1 there; one halfway
between two centres has 0.5 in each. One row per document and topic.

``` r

membership <- hg_membership(hg, topics)
head(membership, 16)
plot(membership)
```

The memberships sum to one over the sixteen topics, so the uniform value
is 0.0625 and a document’s own-topic membership of 0.15 is 2.3 times
that. They are read as ratios between topics: the first document above
is 1.3 times closer to topic 1 than to topic 2 and four times closer
than to topic 3. For the median document the nearest centre is 17%
closer than the second nearest, so most abstracts lie near a boundary
between topics; the hard label records the side they fell on. The plot
shows, per topic, the distribution of its own documents’ membership in
it. A compact topic has its documents well above the uniform value; a
diffuse topic sits close to it.

## Relations between topics

[`hg_relations()`](https://mohsaqr.github.io/hypernets/reference/hg_relations.md)
builds the topic-by-topic co-occurrence network through shared
vocabulary. For each pair of topics it sums, over all words, the product
of the number of documents in each topic that contain the word, the
full-counting aggregation of bibliometric keyword networks. The
`similarity` argument normalises that sum by the measures of van Eck and
Waltman (2009); `"cosine"` divides it by the geometric mean of the two
topics’ totals and lies in \[0, 1\]. The function returns an edge list
with `source`, `target` and `weight`, the columns Gephi reads, and with
`what = "network"` a `cograph_network`.

``` r

relations <- hg_relations(hg, topics, similarity = "cosine")
head(relations)
```

Every pair of topics shares the corpus vocabulary, so all 120 cosine
weights lie between 0.58 and 0.90. The plot keeps the upper quartile of
the edges (`minimum = 0.82`), places the topics with a spring layout and
scales line width by weight.

``` r

network <- hg_relations(hg, topics, similarity = "cosine", what = "network")
cograph::splot(network, layout = "spring", minimum = 0.82,
               edge_width_range = c(0.5, 8))
```

## Summary

The pipeline produced sixteen document clusters from 3,666 abstracts
without embeddings or truncation. The spectrum of the hypergraph
Laplacian supports two or four clusters; the finer partition was chosen
to match a published topic model. Ranking words by the share of their
corpus count that falls in a cluster, with a support floor, yields a
distinct vocabulary for each cluster, whereas ranking by frequency
yields the corpus vocabulary for all of them.

## References

Grootendorst, M. (2022). BERTopic: Neural topic modeling with a
class-based TF-IDF procedure. arXiv:2203.05794.

Saqr, M., Raspopovic Milic, M., Pancheva, K., Jovic, J., Peltekova, E.
V., and Conde, M. Á. (2023). A multimethod synthesis of Covid-19
education research: the tightrope between covidization and
meaningfulness. Universal Access in the Information Society.

van Eck, N. J., and Waltman, L. (2009). How to normalize cooccurrence
data? An analysis of some well-known similarity measures. Journal of the
American Society for Information Science and Technology, 60(8),
1635-1651.

Bezdek, J. C. (1981). Pattern Recognition with Fuzzy Objective Function
Algorithms. Plenum.

Bouma, G. (2009). Normalized (pointwise) mutual information in
collocation extraction. Proceedings of GSCL, 31-40.

Zhou, D., Huang, J., and Schölkopf, B. (2007). Learning with
hypergraphs: clustering, classification, and embedding. Advances in
Neural Information Processing Systems 19.
