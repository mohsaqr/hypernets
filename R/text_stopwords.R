#' English function-word stop list
#'
#' A small, fixed list of English function words (articles, prepositions,
#' conjunctions, pronouns, auxiliaries) for the `stop_words` argument of
#' [text_hypergraph()]. Deliberately minimal and versioned with the package:
#' corpus-specific boilerplate (e.g. "study", "results" in an abstract
#' corpus) should be added by the caller, as in
#' `c(stop_words_en(), "study", "results")`.
#'
#' @return A sorted character vector of lowercase English function words.
#' @references
#' Manning, C. D., Raghavan, P., & Schütze, H. (2008). *Introduction to
#' Information Retrieval*. Cambridge University Press.
#' \doi{10.1017/CBO9780511809071}
#' @examples
#' hg <- text_hypergraph(
#'   c(a = "the salt and the soup", b = "the soup and the stars"),
#'   stop_words = stop_words_en()
#' )
#' hg_get(hg, what = "vocabulary")
#' @export
stop_words_en <- function() {
  sort(c(
    "the", "a", "an", "and", "or", "of", "to", "in", "on", "for", "with",
    "by", "from", "as", "at", "this", "that", "these", "those", "is",
    "are", "was", "were", "be", "been", "being", "it", "its", "we", "our",
    "their", "they", "has", "have", "had", "not", "no", "but", "which",
    "who", "during", "into", "through", "between", "among", "also", "can",
    "could", "may", "will", "would", "than", "then", "there", "here",
    "such", "more", "most", "other", "both", "each", "all", "some", "any",
    "how", "what", "when", "where", "while", "because", "about"
  ))
}
