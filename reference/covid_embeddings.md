# Sentence embeddings of the COVID-19 abstracts

Sentence embeddings of
[covid_abstracts](https://mohsaqr.github.io/hypernets/reference/covid_abstracts.md)'
abstract texts, computed with the `sbert` package's pinned
`all-MiniLM-L6-v2` model (L2-normalized rows). Bundled so that
`text_hypergraph(construction = "knn")` runs offline; rebuilt by
`data-raw/covid_embeddings.R`.

## Usage

``` r
covid_embeddings
```

## Format

A numeric matrix with 165 rows (rownames = `covid_abstracts$doc`) and
384 columns.

## Source

Computed from
[covid_abstracts](https://mohsaqr.github.io/hypernets/reference/covid_abstracts.md)
with `sbert::encode()`.

## Examples

``` r
hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
                      construction = "knn", k = 10,
                      embeddings = covid_embeddings)
hg
#> Text hypergraph: 165 documents (kNN embedding hyperedges: k = 10, cosine)
#> Hyperedges: 165 (kNN neighborhoods); sizes 11-11, median 11
```
