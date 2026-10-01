# English function-word stop list

A small, fixed list of English function words (articles, prepositions,
conjunctions, pronouns, auxiliaries) for the `stop_words` argument of
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md).
Deliberately minimal and versioned with the package: corpus-specific
boilerplate (e.g. "study", "results" in an abstract corpus) should be
added by the caller, as in `c(stop_words_en(), "study", "results")`.

## Usage

``` r
stop_words_en()
```

## Value

A sorted character vector of lowercase English function words.

## References

Manning, C. D., Raghavan, P., & Schütze, H. (2008). *Introduction to
Information Retrieval*. Cambridge University Press.
[doi:10.1017/CBO9780511809071](https://doi.org/10.1017/CBO9780511809071)

## Examples

``` r
hg <- text_hypergraph(
  c(a = "the salt and the soup", b = "the soup and the stars"),
  stop_words = stop_words_en()
)
hg_get(hg, what = "vocabulary")
#>    word count doc_freq
#> 1  salt     1        1
#> 2  soup     2        2
#> 3 stars     1        1
```
