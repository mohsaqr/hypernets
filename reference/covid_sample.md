# A random sample of 1,000 COVID-19 education research abstracts

A simple random sample of 1,000 abstracts from the Scopus export of
COVID-19 education research that the `sbert` package also uses (4,187
records, 2020–2024). Records without an abstract of at least 400
characters are excluded before sampling, and an abstract that appears
twice is kept once. The years are represented in the proportions of the
source. Rebuilt by `data-raw/covid_sample.R`.

## Usage

``` r
covid_sample
```

## Format

A data frame with 1,000 rows and 4 columns:

- doc:

  Scopus EID, the unique document identifier.

- title:

  Article title.

- abstract:

  Abstract text (each at least 400 characters).

- year:

  Publication year (integer, 2020–2024).

## Source

Scopus export of COVID-19 education research, 2020–2024.

## Examples

``` r
hg <- text_hypergraph(covid_sample, column = "abstract", id = "doc",
                      stop_words = stop_words_en(), min_count = 5)
#> vocabulary pruned to 3325 of 10009 words, retaining 90.79% of tokens
hg
#> Text hypergraph: 1000 documents, 3325 words (documents as nodes, weight = n)
#> Hyperedges: 3325 (words); sizes 1-870, median 10
#>                 doc       word count weight
#>  2-s2.0-85082857029   addition     1      1
#>  2-s2.0-85082857029   although     1      1
#>  2-s2.0-85082857029    attempt     1      1
#>  2-s2.0-85082857029 background     1      1
#>  2-s2.0-85082857029      basic     1      1
#>  2-s2.0-85082857029 beneficial     1      1
#>  2-s2.0-85082857029     bridge     1      1
#>  2-s2.0-85082857029       care     1      1
#>  2-s2.0-85082857029    centers     1      1
#>  2-s2.0-85082857029  challenge     1      1
#> ... 80985 more rows
```
