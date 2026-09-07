# COVID-19 education research abstracts

A deterministic stratified sample of 165 abstracts (up to 40 per year,
2020–2024) from a Scopus export of COVID-19 education research, the same
source corpus used by the `sbert` package's topic-modeling articles.
Rebuilt by `data-raw/covid_abstracts.R`.

## Usage

``` r
covid_abstracts
```

## Format

A data frame with 165 rows and 4 columns:

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
hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc")
hg
#> Text hypergraph: 165 documents, 4114 words (documents as nodes, weight = n)
#> Hyperedges: 4114 (words); sizes 1-165, median 1
```
