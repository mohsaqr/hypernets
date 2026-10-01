# Mixed-membership topic model of a hypergraph

Factorizes the document-word counts of a hypergraph into topics, so that
every document is a mixture of topics and every topic a distribution
over words. The counts \\X\\ are approximated by \\W H\\, with
non-negative \\W\\ (documents by topics) and \\H\\ (topics by words), by
the multiplicative updates that minimise the generalised
Kullback-Leibler divergence (Lee & Seung 1999, 2001). With this
divergence the fit is the maximum-likelihood solution of probabilistic
latent semantic analysis (Hofmann 1999; Gaussier & Goutte 2005): the
normalised rows of \\W\\ are the documents' topic shares \\P(z \mid d)\\
and the normalised rows of \\H\\ the topics' word probabilities \\P(w
\mid z)\\. LDA is the same model with Dirichlet priors on these
distributions.

## Usage

``` r
hg_topics(
  hg,
  k,
  nstart = 10L,
  max_iter = 1000L,
  tol = 1e-05,
  depth = 10L,
  seed = 1L,
  parallel = FALSE,
  n_cores = 2L
)

# S3 method for class 'net_hg_topics'
hg_get(
  x,
  what = c("topics", "shares", "words", "documents", "restarts"),
  ...,
  n = 10L,
  topic = NULL,
  top = NULL
)

# S3 method for class 'net_hg_topics'
print(x, n = 10L, ...)

# S3 method for class 'net_hg_topics'
plot(x, y, n = 8L, ...)
```

## Arguments

- hg:

  A `net_hg` or `text_hypergraph`.

- k:

  Number of topics, a whole number of at least 2 and smaller than the
  numbers of documents and words.

- nstart:

  Number of random starts (default `10`).

- max_iter:

  Maximum multiplicative updates per start (default `1000`).

- tol:

  Relative change of the divergence between two checks, ten updates
  apart, below which a start has converged (default `1e-5`); `0` runs
  every start for `max_iter` updates.

- depth:

  Number of top words compared between starts (default `10`).

- seed:

  Seed of the first start; start `r` uses `seed + r - 1`. The caller's
  random-number stream is restored on exit. Default `1`.

- parallel:

  Run the starts with
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (not on Windows). Every start seeds itself, so the result equals the
  serial one. Default `FALSE`.

- n_cores:

  Cores when `parallel = TRUE` (default `2`).

- x:

  A `net_hg_topics` object.

- what:

  Which table: `"topics"` (default), `"shares"`, `"words"`,
  `"documents"` or `"restarts"`.

- ...:

  Unused.

- n:

  For `what = "words"`: the number of words per topic (default `10`).
  For [`print()`](https://rdrr.io/r/base/print.html): rows of the topic
  table shown.

- topic:

  For `what = "shares"` and `"words"`: keep only these topics (labels
  such as `"Topic 3"`).

- top:

  Keep only the first `top` rows of the returned table.

- y:

  Unused.

## Value

A `net_hg_topics` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
one row per topic (`what = "topics"`, the default: `topic`, `prevalence`
(its mean share over the documents), `documents` (the summed shares, the
expected number of documents), `agreement` (its mean matched average
Jaccard over the other starts) and `top_words`); one row per document
and topic (`"shares"`: `node`, `topic`, `share`, summing to one over a
document's topics); one row per topic and word (`"words"`: `topic`,
`rank`, `word`, `probability`, `n` words per topic); one row per
document with its largest share (`"documents"`: `node`, `topic`,
`share`); or one row per start (`"restarts"`: `run`, `divergence`,
`iterations`, `converged`, `best`, `agreement`).
[`print()`](https://rdrr.io/r/base/print.html) shows the topics,
[`summary()`](https://rdrr.io/r/base/summary.html) every table, and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the top words
of every topic. A start that reaches `max_iter` without converging
raises the warning `hypernets_no_converge`.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns a
ggplot of the `n` most probable words of every topic, one panel per
topic.

## Details

The fit depends on its random start, so it is run from `nstart` starts
and the start with the lowest divergence is kept. Every other start is
compared with it by the agreement of the topics' top words (Greene,
O'Callaghan & Cunningham 2014): the topics of the two fits are matched
one to one by the Hungarian method (Kuhn 1955) on the average Jaccard
coefficient of their ranked top `depth` words, and the matched
coefficients are averaged. A topic that the other starts find again has
a high `agreement`.

A text hypergraph from
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
is factorized on its word counts, whatever its node orientation and
incidence weighting; any other hypergraph on its incidence, with its
nodes as documents and its hyperedges as words.

## Conditions

`hypernets_bad_input` for an input that is not a hypergraph, a `k`
outside its range, or invalid control arguments.

## References

Lee, D. D., & Seung, H. S. (1999). Learning the parts of objects by
non-negative matrix factorization. *Nature*, 401, 788-791.
[doi:10.1038/44565](https://doi.org/10.1038/44565)

Lee, D. D., & Seung, H. S. (2001). Algorithms for non-negative matrix
factorization. *Advances in Neural Information Processing Systems*, 13,
556-562.

Hofmann, T. (1999). Probabilistic latent semantic indexing. *Proceedings
of SIGIR 1999*, 50-57.
[doi:10.1145/312624.312649](https://doi.org/10.1145/312624.312649)

Gaussier, E., & Goutte, C. (2005). Relation between PLSA and NMF and
implications. *Proceedings of SIGIR 2005*, 601-602.
[doi:10.1145/1076034.1076148](https://doi.org/10.1145/1076034.1076148)

Greene, D., O'Callaghan, D., & Cunningham, P. (2014). How many topics?
Stability analysis for topic models. *Machine Learning and Knowledge
Discovery in Databases (ECML PKDD 2014)*, 498-513.
[doi:10.1007/978-3-662-44848-9_32](https://doi.org/10.1007/978-3-662-44848-9_32)

Kuhn, H. W. (1955). The Hungarian method for the assignment problem.
*Naval Research Logistics Quarterly*, 2, 83-97.
[doi:10.1002/nav.3800020109](https://doi.org/10.1002/nav.3800020109)

## See also

[`hg_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
scores the topics' words,
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
for a partition of the documents.

## Examples

``` r
corpus <- c(
  a = "soup salt onion soup broth", b = "salt soup broth onion",
  c = "stars sky moon night", d = "sky stars night moon moon",
  e = "soup stars salt sky night broth")
fit <- hg_topics(text_hypergraph(corpus), k = 2, nstart = 3)
fit
#> Topic model (KL factorization): 2 topics, 5 documents, 8 words; best of 3 starts (divergence 1.92745, converged)
#>    topic prevalence documents agreement                      top_words
#>  Topic 1        0.5       2.5 0.8809524  soup, salt, broth, onion, sky
#>  Topic 2        0.5       2.5 0.8020833 moon, night, stars, sky, broth
hg_get(fit, what = "shares")
#>    node   topic         share
#> 1     a Topic 1  1.000000e+00
#> 2     a Topic 2 3.864335e-202
#> 3     b Topic 1  1.000000e+00
#> 4     b Topic 2 1.851001e-200
#> 5     c Topic 1 4.254015e-196
#> 6     c Topic 2  1.000000e+00
#> 7     d Topic 1 7.006930e-199
#> 8     d Topic 2  1.000000e+00
#> 9     e Topic 1  5.000000e-01
#> 10    e Topic 2  5.000000e-01
hg_get(fit, what = "words", n = 4)
#>     topic rank  word probability
#> 1 Topic 1    1  soup   0.3333333
#> 2 Topic 1    2  salt   0.2500000
#> 3 Topic 1    3 broth   0.2500000
#> 4 Topic 1    4 onion   0.1666667
#> 5 Topic 2    1  moon   0.2500000
#> 6 Topic 2    2 night   0.2500000
#> 7 Topic 2    3 stars   0.2500000
#> 8 Topic 2    4   sky   0.2500000
```
