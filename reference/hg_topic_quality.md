# Topic coherence and exclusivity

Scores the top words of each topic on two axes, with the published
definitions that topic-model software reports, so hypergraph topics and
topics from another model (LDA, STM) can be compared on the same corpus.
The corpus is the document hypergraph: its documents are the reference
collection and a word's documents are the documents it occurs in.

## Usage

``` r
hg_topic_quality(
  hg,
  clusters = NULL,
  words = NULL,
  n = 10L,
  min_docs = 1L,
  coherence = c("umass", "npmi", "npmi_cluster"),
  exclusivity = c("frex", "share", "none"),
  sort_by = c("probability", "share"),
  frexw = 0.7,
  topics = NULL
)

# S3 method for class 'hypernets_topic_quality'
plot(x, ...)
```

## Arguments

- hg:

  The document hypergraph the topics describe: a bag-of-words
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  with documents as nodes (the default construction).

- clusters:

  The tidy table returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
  (columns `node`, `cluster`), or a named vector of cluster labels.
  `NULL` when only `words` are scored.

- words:

  `NULL` (default); a data.frame of topic words with columns `word` and
  `topic` (or `cluster`, so
  [`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
  output is passed as it is) and optionally `rank` (the order within a
  topic; row order otherwise); or a character matrix with one column per
  topic and rows in rank order, as `topicmodels::terms()` returns it
  (column names are the topic labels). Every word must be in the
  vocabulary of `hg`. With `clusters` too, every topic must be a cluster
  label.

- n:

  Top words per topic to score (default `10`).

- min_docs:

  Support floor: with `clusters` and no `words`, keep only words found
  in at least this many of the cluster's documents (default `1`).

- coherence:

  `"umass"` (default), `"npmi"` or `"npmi_cluster"`; see Details.

- exclusivity:

  `"frex"` (default), `"share"` or `"none"`; see Details.

- sort_by:

  How a cluster's top words are chosen when `words` is not given:
  `"probability"` (default) or `"share"`; see Details.

- frexw:

  Weight on exclusivity in FREX, in (0, 1) (default `0.7`, as in `stm`).

- topics:

  A mixed-membership topic model of `hg` fitted by
  [`hg_topics()`](https://mohsaqr.github.io/hypernets/reference/hg_topics.md).
  Its topics are scored on their most probable words, and FREX is
  computed from its word distributions \\P(w \mid z)\\, as `stm`
  computes it from an STM's. `size` is then the topic's expected number
  of documents (its summed shares). Give `topics` instead of `clusters`
  and `words`.

- x:

  A table returned by the verb.

- ...:

  Unused; for S3 consistency.

## Value

A base `data.frame` of class `hypernets_topic_quality`, one row per
topic (cluster labels in natural order, or topics of `words` in natural
order): `topic`, `size` (the cluster's documents; `NA` without
`clusters`), `n_words` (top words scored), `coherence`, `exclusivity`
(`NA` for fewer than two words, or `exclusivity = "none"`),
`coherence_type` and `exclusivity_type`.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws
exclusivity against coherence with one labelled point per topic. Raises
`hypernets_bad_input` for unknown node names, a malformed `words` table,
words outside the vocabulary, topics of `words` that are not cluster
labels, a measure that needs `clusters` without them, or a hypergraph
whose nodes are not documents.

## Details

**Topics.** A topic is either a cluster of documents (`clusters`, from
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
or any labelling) or a list of words (`words`, e.g. the top terms of a
fitted LDA or STM model), or both. For a cluster, the word distribution
is the plug-in estimate \\\beta\_{kw} = c\_{kw} / \sum_v c\_{kv}\\,
where \\c\_{kw}\\ is the token count of word \\w\\ in the cluster's
documents, over the words that occur in some clustered document; its top
words are the `n` most probable (`sort_by = "probability"`, ties in
vocabulary order, as `stm` ranks \\\beta\\) or the `n` with the largest
share of their corpus count in the cluster (`sort_by = "share"`, the
ranking of
[`hg_keywords()`](https://mohsaqr.github.io/hypernets/reference/hg_keywords.md)
with `sort_by = "share"`), keeping words found in at least `min_docs` of
the cluster's documents. Given `words`, those words are scored in their
given order.

**Coherence** (`coherence =`), with \\D(v)\\ the number of documents
containing \\v\\, \\D(v, u)\\ the number containing both, \\N\\ the
number of documents and \\v_1, \dots, v_M\\ the top words in rank order:

- `"umass"` (default): the UMass coherence of Mimno et al. (2011),
  \\\sum\_{m=2}^{M} \sum\_{l=1}^{m-1} \log\\(D(v_m, v_l) + \epsilon) /
  (D(v_l) + \epsilon)\\\\, a sum over the \\M(M-1)/2\\ pairs that
  conditions on the higher-ranked word, counted on the whole corpus,
  with \\\epsilon = 0.01\\ in numerator and denominator as in
  `stm::semanticCoherence()` (Roberts et al. 2019); Mimno et al. add 1
  to the numerator only. At most 0; closer to 0 is more coherent.

- `"npmi"`: the mean normalised pointwise mutual information over the
  word pairs (Bouma 2009; the topic-coherence measure of Lau et al.
  2014), with document co-occurrence probabilities \\p(v) = D(v)/N\\ and
  \\p(v, u) = D(v, u)/N\\ on the whole corpus: \\\log\\p(v,u) /
  (p(v)p(u))\\ / -\log p(v,u)\\; -1 for a pair that never co-occurs (the
  limit), 1 for a pair present in every document.

- `"npmi_cluster"`: the same NPMI counted only on the cluster's own
  documents (the definition this verb used before 0.6.0; needs
  `clusters`). It rewards words that co-occur inside the topic whatever
  they do elsewhere, so it is not comparable with published values.

**Exclusivity** (`exclusivity =`), which needs the cluster word
distributions and so `clusters`:

- `"frex"` (default): the sum over the top words of FREX (Bischof &
  Airoldi 2012), the weighted harmonic mean of a word's rank in the
  topic's distribution and its rank in exclusivity \\\beta\_{kw} /
  \sum_j \beta\_{jw}\\, both as empirical CDFs over the vocabulary: \\1
  / (w / \mathrm{ex} + (1 - w) / \mathrm{fr})\\ with weight `frexw` on
  exclusivity, exactly as `stm::exclusivity()` computes it. Between 0
  and `n`; larger is more exclusive.

- `"share"`: the mean over the top words of the share of the word's
  corpus count that falls in the cluster, \\c\_{kw} / \sum_j c\_{jw}\\
  (the definition used before 0.6.0); 1 when the cluster owns its words.

- `"none"`: no exclusivity (`NA`); the choice when only `words` are
  given.

The measures before 0.6.0 are
`coherence = "npmi_cluster", exclusivity = "share", sort_by = "share"`.

## References

Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A.
(2011). Optimizing semantic coherence in topic models. *Proceedings of
the 2011 Conference on Empirical Methods in Natural Language Processing
(EMNLP)*, 262–272. <https://aclanthology.org/D11-1024/>

Bouma, G. (2009). Normalized (pointwise) mutual information in
collocation extraction. *Proceedings of GSCL*, 31–40.

Lau, J. H., Newman, D., & Baldwin, T. (2014). Machine reading tea
leaves: Automatically evaluating topic coherence and topic model
quality. *Proceedings of the 14th Conference of the European Chapter of
the Association for Computational Linguistics (EACL)*, 530–539.
[doi:10.3115/v1/E14-1056](https://doi.org/10.3115/v1/E14-1056)

Bischof, J. M., & Airoldi, E. M. (2012). Summarizing topical content
with word frequency and exclusivity. *Proceedings of the 29th
International Conference on Machine Learning (ICML)*. arXiv:1206.4631.

Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
for structural topic models. *Journal of Statistical Software*, 91(2),
1–40. [doi:10.18637/jss.v091.i02](https://doi.org/10.18637/jss.v091.i02)

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_topic_quality(hg, topics, n = 3)
#>       topic size n_words coherence exclusivity coherence_type exclusivity_type
#> 1 Cluster 1    2       3 -5.991489    2.358896          umass             frex
#> 2 Cluster 2    2       3 -1.376369    2.489531          umass             frex
hg_topic_quality(hg, topics, n = 3, coherence = "npmi")
#>       topic size n_words coherence exclusivity coherence_type exclusivity_type
#> 1 Cluster 1    2       3 0.0000000    2.358896           npmi             frex
#> 2 Cluster 2    2       3 0.6666667    2.489531           npmi             frex
# topic words from another model, scored on the same corpus
lda_words <- data.frame(topic = c("t1", "t1", "t2", "t2"),
                        word = c("soup", "salt", "telescope", "stars"))
hg_topic_quality(hg, words = lda_words, exclusivity = "none")
#>   topic size n_words  coherence exclusivity coherence_type exclusivity_type
#> 1    t1   NA       2 -0.6881844          NA          umass             none
#> 2    t2   NA       2  0.0000000          NA          umass             none
```
