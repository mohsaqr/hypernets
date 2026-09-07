# Corpus repair before tokenization: HTML, encoding damage, references,
# URLs, copyright tails, numbers and stop words. Deterministic base regex;
# every step is a named switch so a document can say exactly what was
# removed. Returns the input shape (same length, same rows) so nothing
# silently drops out of alignment: rows that fail `min_content` become ""
# and text_hypergraph() reports them as dropped documents. Every non-ASCII
# character in this file is written as a \u escape so the source is
# portable.

#' Clean a corpus before building a text hypergraph
#'
#' Repairs and strips the non-content a bibliographic or web export leaves in
#' a text: HTML tags and entities, mojibake and typographic characters,
#' bracketed citation numbers and list or section numbering, URLs and DOIs,
#' trailing copyright notices, bare numbers, and optionally stop words. Each
#' step is a switch. The result has the same length (or the same rows) as the
#' input; a text that falls below `min_content` becomes the empty string,
#' which [text_hypergraph()] then drops with a `hypernets_dropped_documents`
#' warning, so nothing leaves the corpus silently.
#'
#' @param x A character vector, or a data.frame with a text column.
#' @param column When `x` is a data.frame, the name of its text column.
#' @param html Decode HTML entities (`&amp;`, `&nbsp;`, `&#8217;`, ...) and
#'   strip tags (default `TRUE`).
#' @param encoding Repair common UTF-8-read-as-Latin-1 mojibake (the
#'   three-character garble of a curly apostrophe) and normalise typographic
#'   quotes, dashes and non-breaking spaces to their ASCII forms (default
#'   `TRUE`).
#' @param citations Remove bracketed reference numbers (`[12]`, `[3, 4]`,
#'   `[1-5]`), leading list markers (`1.`, `(3)`, `a)`) and dotted section
#'   numbers (`3.2.1`) (default `TRUE`).
#' @param urls Remove URLs, `www.` addresses, `doi:` strings and bare DOIs
#'   (default `TRUE`).
#' @param copyright Remove a trailing copyright notice: everything from the
#'   last copyright sign, `(c)` or `Copyright` to the end when that tail is
#'   under `copyright_max` characters, plus "All rights reserved" (default
#'   `TRUE`).
#' @param copyright_max Longest tail treated as a notice (default `300`).
#' @param numbers Remove bare numbers, percentages and years (default
#'   `TRUE`). Numbers never enter a text hypergraph's vocabulary anyway (the
#'   tokeniser keeps alphabetic tokens), so this matters for display and for
#'   `min_content`.
#' @param remove Extra patterns to remove (Perl regular expressions,
#'   case-insensitive), applied after the repairs and before the content
#'   floor: a database's placeholder for a missing abstract, a publisher's
#'   own boilerplate. Default `NULL`.
#' @param stop_words Words to remove, case-insensitively at word boundaries
#'   (default `NULL`, none). Usually left to [text_hypergraph()]'s own
#'   `stop_words`; use it here when the cleaned text itself is shown.
#' @param min_content Minimum share of alphabetic characters among
#'   non-space characters, judged on the repaired text before numbers are
#'   removed (default `0`, keep everything). A reference list or a table
#'   fragment scores low; `0.5` is a sensible floor for abstracts.
#' @return `x` with the text cleaned: a character vector of the same length
#'   (names kept), or the same data.frame with `column` replaced. `NA` and
#'   low-content texts become `""`. Raises `hypernets_bad_input` for a
#'   malformed `x`, `column` or switch.
#' @examples
#' messy <- c(
#'   "1. The programme raised attainment in 2020, see section 3.2.1.",
#'   "Peer&nbsp;support <b>improved</b> outcomes [12] across schools.",
#'   "Full report at https://example.org/study.pdf and doi:10.1/x.",
#'   "Students' views shifted. (c) 2021 Informa UK Limited.",
#'   "(3) OJ No L 297, 24.11.1979, p. 1."
#' )
#' clean_text(messy)
#' clean_text(messy, min_content = 0.5)
#' clean_text(c("[No abstract available]", "A real abstract."),
#'            remove = "no abstract available")
#' @export
clean_text <- function(x, column = NULL, html = TRUE, encoding = TRUE,
                       citations = TRUE, urls = TRUE, copyright = TRUE,
                       copyright_max = 300L, numbers = TRUE,
                       remove = NULL, stop_words = NULL, min_content = 0) {
  flags <- list(html = html, encoding = encoding, citations = citations,
                urls = urls, copyright = copyright, numbers = numbers)
  bad_flag <- names(flags)[!vapply(flags, \(f) isTRUE(f) || isFALSE(f),
                                   logical(1))]
  if (length(bad_flag) > 0L) {
    stop(errorCondition(
      sprintf("`%s` must be TRUE or FALSE", bad_flag[[1L]]),
      class = "hypernets_bad_input", call = NULL
    ))
  }
  stopifnot(
    "`copyright_max` must be a single positive number" =
      length(copyright_max) == 1L && is.numeric(copyright_max) &&
      copyright_max >= 1,
    "`min_content` must be a single number in [0, 1]" =
      length(min_content) == 1L && is.numeric(min_content) &&
      min_content >= 0 && min_content <= 1,
    "`stop_words` must be NULL or a character vector" =
      is.null(stop_words) || is.character(stop_words),
    "`remove` must be NULL or a character vector of patterns" =
      is.null(remove) || (is.character(remove) && !anyNA(remove))
  )
  if (is.data.frame(x)) {
    if (!is.character(column) || length(column) != 1L ||
        !column %in% names(x)) {
      stop(errorCondition(
        "when `x` is a data.frame, `column` must name one of its columns",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    text <- as.character(x[[column]])
  } else {
    if (!is.character(x)) {
      stop(errorCondition(
        "`x` must be a character vector or a data.frame",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    if (!is.null(column)) {
      stop(errorCondition(
        "`column` applies only when `x` is a data.frame",
        class = "hypernets_bad_input", call = NULL
      ))
    }
    text <- x
  }
  text[is.na(text)] <- ""

  if (isTRUE(html)) text <- .thg_clean_html(text)
  if (isTRUE(encoding)) text <- .thg_clean_encoding(text)
  if (isTRUE(urls)) text <- .thg_clean_urls(text)
  if (isTRUE(copyright)) text <- .thg_clean_copyright(text, copyright_max)
  if (isTRUE(citations)) text <- .thg_clean_citations(text)
  for (pattern in remove) {
    text <- gsub(paste0("(?i)", pattern), " ", text, perl = TRUE)
  }
  # the content floor is judged before numbers go: stripping a reference's
  # dates and pages would otherwise raise its letter share
  content_ok <- if (min_content > 0) {
    .thg_content_ratio(text) >= min_content
  } else {
    rep(TRUE, length(text))
  }
  if (isTRUE(numbers)) text <- .thg_clean_numbers(text)
  if (!is.null(stop_words) && length(stop_words) > 0L) {
    text <- .thg_clean_stop_words(text, stop_words)
  }
  text <- .thg_clean_tidy(text)
  text[!content_ok] <- ""

  if (is.data.frame(x)) {
    x[[column]] <- text
    x
  } else {
    text
  }
}

.thg_clean_html <- function(text) {
  text <- gsub("<[^>]+>", " ", text, perl = TRUE)
  named <- c(nbsp = " ", amp = "&", lt = "<", gt = ">", quot = "\"",
             apos = "'", ndash = "-", mdash = "-", hellip = "...",
             lsquo = "'", rsquo = "'", ldquo = "\"", rdquo = "\"")
  for (entity in names(named)) {
    text <- gsub(paste0("&", entity, ";"), named[[entity]], text,
                 fixed = TRUE)
  }
  decode <- function(m, base) {
    code <- strtoi(sub("^&#[xX]?", "", sub(";$", "", m)), base = base)
    ok <- !is.na(code) & code > 0 & code < 1114112
    out <- m
    out[ok] <- vapply(code[ok], intToUtf8, character(1))
    out
  }
  text <- .thg_gsub_fn(text, "&#[0-9]{1,7};", \(m) decode(m, 10L))
  .thg_gsub_fn(text, "&#[xX][0-9A-Fa-f]{1,6};", \(m) decode(m, 16L))
}

# apply `fn` to every regex match (base R has no gsub with a callback)
.thg_gsub_fn <- function(text, pattern, fn) {
  hits <- gregexpr(pattern, text, perl = TRUE)
  regmatches(text, hits) <- lapply(regmatches(text, hits), \(m) {
    if (length(m) == 0L) m else fn(m)
  })
  text
}

.thg_clean_encoding <- function(text) {
  # UTF-8 bytes of curly quotes / dashes read as Latin-1 begin with
  # a-circumflex + euro sign; the third character tells which mark it was
  garble <- "\u00e2\u20ac"
  mojibake <- c("'", "'", "\"", "\"", "-", "-", "...")
  names(mojibake) <- paste0(garble, c("\u2122", "\u02dc", "\u0153",
                                      "\u009d", "\u201c", "\u201d",
                                      "\u00a6"))
  mojibake <- c(mojibake,
                "\u00c2\u00a0" = " ",        # non-breaking space garbled
                "\u00c3\u00a9" = "\u00e9",   # e-acute garbled
                "\u00ef\u00bf\u00bd" = "")   # replacement character garbled
  for (bad in names(mojibake)) {
    text <- gsub(bad, mojibake[[bad]], text, fixed = TRUE)
  }
  # typographic quotes, dashes, ellipsis, non-breaking / zero-width spaces
  typographic <- c("\u2018" = "'", "\u2019" = "'", "\u201a" = "'",
                   "\u201c" = "\"", "\u201d" = "\"", "\u201e" = "\"",
                   "\u2013" = "-", "\u2014" = "-", "\u2212" = "-",
                   "\u2026" = "...", "\u00a0" = " ", "\u202f" = " ",
                   "\u200b" = "", "\ufeff" = "")
  for (ch in names(typographic)) {
    text <- gsub(ch, typographic[[ch]], text, fixed = TRUE)
  }
  text
}

.thg_clean_urls <- function(text) {
  # a URL runs to the next space, but sentence punctuation right before
  # that space belongs to the sentence, not the address
  ending <- "\\S+?(?=[.,;:)]*(?:\\s|$))"
  text <- gsub(paste0("(?i)\\b(?:https?://|www\\.)", ending), " ", text,
               perl = TRUE)
  text <- gsub(paste0("(?i)\\bdoi:?\\s*10\\.\\d{1,9}/", ending), " ", text,
               perl = TRUE)
  gsub(paste0("\\b10\\.\\d{4,9}/", ending), " ", text, perl = TRUE)
}

.thg_clean_copyright <- function(text, max_chars) {
  sign <- "\u00a9"
  tail <- sprintf("(?i)(?:%s|\\(c\\)|\\bcopyright\\b)[^%s]{0,%d}$",
                  sign, sign, as.integer(max_chars))
  text <- gsub(tail, " ", text, perl = TRUE)
  gsub("(?i)\\ball rights reserved\\.?", " ", text, perl = TRUE)
}

.thg_clean_citations <- function(text) {
  # [12], [3, 4], [1-5], [12-14, 20]
  text <- gsub("\\[\\s*\\d+(?:\\s*[-,]\\s*\\d+)*\\s*\\]", " ", text,
               perl = TRUE)
  # leading list markers at the start or after sentence punctuation:
  # "1. ", "(3) ", "3) ", "a) ", "iv. "
  text <- gsub(
    "(^|[.!?;:]\\s+)(?:\\(?\\d{1,3}[.)]|\\(?[a-z][.)]|\\(?[ivx]{1,5}[.)])\\s+",
    "\\1", text, perl = TRUE
  )
  # dotted section numbers: 3.2, 3.2.1
  gsub("\\b\\d+(?:\\.\\d+)+\\b", " ", text, perl = TRUE)
}

.thg_content_ratio <- function(text) {
  compact <- gsub("[[:space:]]+", "", text, perl = TRUE)
  letters_n <- nchar(gsub("[^[:alpha:]]", "", compact, perl = TRUE))
  total <- nchar(compact)
  ifelse(total > 0L, letters_n / total, 0)
}

.thg_clean_numbers <- function(text) {
  # standalone numbers, decimals, thousands, percentages, ordinals
  gsub("(?<![[:alpha:]])[-+]?\\d+(?:[.,]\\d+)*(?:%|st|nd|rd|th)?(?![[:alpha:]])",
       " ", text, perl = TRUE)
}

.thg_clean_stop_words <- function(text, stop_words) {
  stop_words <- unique(tolower(stop_words))
  escaped <- gsub("([.\\\\+*?\\[^\\]$(){}=!<>|:-])", "\\\\\\1", stop_words,
                  perl = TRUE)
  pattern <- sprintf("(?i)(?<![[:alpha:]'])(?:%s)(?![[:alpha:]'])",
                     paste(escaped, collapse = "|"))
  gsub(pattern, " ", text, perl = TRUE)
}

.thg_clean_tidy <- function(text) {
  text <- gsub("\\s+([,.;:!?])", "\\1", text, perl = TRUE)
  text <- gsub("\\(\\s*\\)|\\[\\s*\\]", " ", text, perl = TRUE)
  text <- gsub("(?:\\s*[,;:]\\s*){2,}", ", ", text, perl = TRUE)
  text <- gsub(",\\s*\\.", ".", text, perl = TRUE)
  text <- gsub("[[:space:]]+", " ", text, perl = TRUE)
  text <- gsub("^[\\s,.;:]+", "", text, perl = TRUE)
  trimws(text)
}
