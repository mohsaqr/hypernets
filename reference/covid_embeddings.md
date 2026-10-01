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
#>                 doc               edge    weight
#>  2-s2.0-85085897904 2-s2.0-85085897904 1.0000000
#>  2-s2.0-85092100063 2-s2.0-85085897904 0.5784754
#>  2-s2.0-85092577169 2-s2.0-85085897904 0.5940401
#>  2-s2.0-85096991345 2-s2.0-85085897904 0.6156138
#>  2-s2.0-85096994563 2-s2.0-85085897904 0.5771233
#>  2-s2.0-85108879636 2-s2.0-85085897904 0.6228458
#>  2-s2.0-85109406490 2-s2.0-85085897904 0.5927937
#>  2-s2.0-85119298190 2-s2.0-85085897904 0.6315365
#>  2-s2.0-85124102536 2-s2.0-85085897904 0.5749891
#>  2-s2.0-85124535349 2-s2.0-85085897904 0.7546323
#> ... 1805 more rows
```
