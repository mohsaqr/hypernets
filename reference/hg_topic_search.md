# Choose the number of topics by coherence and exclusivity

Fits
[`hg_topics()`](https://mohsaqr.github.io/hypernets/reference/hg_topics.md)
for every number of topics in `k` and scores each model by the mean
semantic coherence (UMass; Mimno et al. 2011) and the mean exclusivity
(FREX; Bischof & Airoldi 2012) of its topics, computed by
[`hg_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md)
on the `n` most probable words. This is the diagnostic of
`stm::searchK()` (Roberts, Stewart & Tingley 2019). Adding topics
usually makes them more exclusive and less coherent, so no single number
is best; Roberts et al. (2014) choose among the models on the frontier
of the two, those that no other model beats on both. A model is on the
`frontier` when no other number of topics has both a higher mean
coherence and a higher mean exclusivity.

## Usage

``` r
hg_topic_search(
  hg,
  k,
  nstart = 1L,
  n = 10L,
  frexw = 0.7,
  max_iter = 1000L,
  tol = 1e-05,
  seed = 1L,
  parallel = FALSE,
  n_cores = 2L
)

# S3 method for class 'hypernets_topic_search'
plot(x, ...)
```

## Arguments

- hg:

  A `net_hg` or `text_hypergraph`.

- k:

  Numbers of topics to compare, at least two whole numbers.

- nstart:

  Number of random starts (default `10`).

- n:

  Top words per topic scored (default `10`).

- frexw:

  Weight on exclusivity in FREX, in (0, 1) (default `0.7`).

- max_iter:

  Maximum multiplicative updates per start (default `1000`).

- tol:

  Relative change of the divergence between two checks, ten updates
  apart, below which a start has converged (default `1e-5`); `0` runs
  every start for `max_iter` updates.

- seed:

  Seed of the first start; start `r` uses `seed + r - 1`. The caller's
  random-number stream is restored on exit. Default `1`.

- parallel:

  Fit the numbers of topics with
  [`parallel::mclapply()`](https://rdrr.io/r/parallel/mclapply.html)
  (not on Windows); every fit seeds itself, so the result equals the
  serial one. Default `FALSE`.

- n_cores:

  Cores when `parallel = TRUE` (default `2`).

- x:

  A table returned by the verb.

- ...:

  Unused.

## Value

A base `data.frame` of class `hypernets_topic_search`, one row per
number of topics: `k`, `coherence` and `exclusivity` (means over the
topics), `divergence` (of the kept start), `agreement` (mean agreement
of the topics across the starts; `NA` with `nstart = 1`), `converged`
and `frontier`. [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
draws exclusivity against coherence, one labelled point per number of
topics, the frontier joined by a line and drawn as filled points.
Conditions are those of
[`hg_topics()`](https://mohsaqr.github.io/hypernets/reference/hg_topics.md);
a start that does not converge raises `hypernets_no_converge`.

## References

Roberts, M. E., Stewart, B. M., & Tingley, D. (2019). stm: An R package
for structural topic models. *Journal of Statistical Software*, 91(2),
1-40. [doi:10.18637/jss.v091.i02](https://doi.org/10.18637/jss.v091.i02)

Roberts, M. E., Stewart, B. M., Tingley, D., Lucas, C., Leder-Luis, J.,
Gadarian, S. K., Albertson, B., & Rand, D. G. (2014). Structural topic
models for open-ended survey responses. *American Journal of Political
Science*, 58(4), 1064-1082.
[doi:10.1111/ajps.12103](https://doi.org/10.1111/ajps.12103)

Mimno, D., Wallach, H. M., Talley, E., Leenders, M., & McCallum, A.
(2011). Optimizing semantic coherence in topic models. *Proceedings of
EMNLP 2011*, 262-272.

Bischof, J. M., & Airoldi, E. M. (2012). Summarizing topical content
with word frequency and exclusivity. *Proceedings of ICML 2012*,
201-208.

## See also

[`hg_topics()`](https://mohsaqr.github.io/hypernets/reference/hg_topics.md),
[`hg_topic_quality()`](https://mohsaqr.github.io/hypernets/reference/hg_topic_quality.md).

## Examples

``` r
corpus <- c(
  a = "soup salt onion soup broth", b = "salt soup broth onion",
  c = "stars sky moon night", d = "sky stars night moon moon",
  e = "soup stars salt sky night broth", f = "onion salt stars broth")
search <- hg_topic_search(text_hypergraph(corpus), k = 2:3, nstart = 2,
                          n = 4)
search
#>   k  coherence exclusivity divergence agreement converged frontier
#> 1 2 -0.4887538    3.240670   4.241503 0.5833333      TRUE     TRUE
#> 2 3 -1.3415043    3.242266   1.927450 0.6888889      TRUE     TRUE
```
