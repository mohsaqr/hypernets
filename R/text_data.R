#' COVID-19 education research abstracts
#'
#' A deterministic stratified sample of 165 abstracts (up to 40 per year,
#' 2020--2024) from a Scopus export of COVID-19 education research, the same
#' source corpus used by the `sbert` package's topic-modeling articles.
#' Rebuilt by `data-raw/covid_abstracts.R`.
#'
#' @format A data frame with 165 rows and 4 columns:
#' \describe{
#'   \item{doc}{Scopus EID, the unique document identifier.}
#'   \item{title}{Article title.}
#'   \item{abstract}{Abstract text (each at least 400 characters).}
#'   \item{year}{Publication year (integer, 2020--2024).}
#' }
#' @source Scopus export of COVID-19 education research, 2020--2024.
#' @examples
#' hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc")
#' hg
"covid_abstracts"

#' A random sample of 1,000 COVID-19 education research abstracts
#'
#' A simple random sample of 1,000 abstracts from the Scopus export of
#' COVID-19 education research that the `sbert` package also uses
#' (4,187 records, 2020--2024). Records without an abstract of at least 400
#' characters are excluded before sampling, and an abstract that appears
#' twice is kept once. The years are represented in the proportions of the
#' source. Rebuilt by `data-raw/covid_sample.R`.
#'
#' @format A data frame with 1,000 rows and 4 columns:
#' \describe{
#'   \item{doc}{Scopus EID, the unique document identifier.}
#'   \item{title}{Article title.}
#'   \item{abstract}{Abstract text (each at least 400 characters).}
#'   \item{year}{Publication year (integer, 2020--2024).}
#' }
#' @source Scopus export of COVID-19 education research, 2020--2024.
#' @examples
#' hg <- text_hypergraph(covid_sample, column = "abstract", id = "doc",
#'                       stop_words = stop_words_en(), min_count = 5)
#' hg
"covid_sample"

#' Sentence embeddings of the COVID-19 abstracts
#'
#' Sentence embeddings of [covid_abstracts]' abstract texts, computed with
#' the `sbert` package's pinned `all-MiniLM-L6-v2` model (L2-normalized
#' rows). Bundled so that `text_hypergraph(construction = "knn")` runs
#' offline; rebuilt by `data-raw/covid_embeddings.R`.
#'
#' @format A numeric matrix with 165 rows (rownames = `covid_abstracts$doc`)
#'   and 384 columns.
#' @source Computed from [covid_abstracts] with `sbert::encode()`.
#' @examples
#' hg <- text_hypergraph(covid_abstracts, column = "abstract", id = "doc",
#'                       construction = "knn", k = 10,
#'                       embeddings = covid_embeddings)
#' hg
"covid_embeddings"
