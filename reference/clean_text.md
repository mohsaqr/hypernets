# Clean a corpus before building a text hypergraph

Repairs and strips the non-content a bibliographic or web export leaves
in a text: HTML tags and entities, mojibake and typographic characters,
bracketed citation numbers and list or section numbering, URLs and DOIs,
trailing copyright notices, bare numbers, and optionally stop words.
Each step is a switch. The result has the same length (or the same rows)
as the input; a text that falls below `min_content` becomes the empty
string, which
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
then drops with a `hypernets_dropped_documents` warning, so nothing
leaves the corpus silently.

## Usage

``` r
clean_text(
  x,
  column = NULL,
  html = TRUE,
  encoding = TRUE,
  citations = TRUE,
  urls = TRUE,
  copyright = TRUE,
  copyright_max = 300L,
  numbers = TRUE,
  remove = NULL,
  stop_words = NULL,
  min_chars = 0L,
  min_content = 0
)
```

## Arguments

- x:

  A character vector, or a data.frame with a text column.

- column:

  When `x` is a data.frame, the name of its text column.

- html:

  Decode HTML entities (`&amp;`, `&nbsp;`, `&#8217;`, ...) and strip
  tags (default `TRUE`).

- encoding:

  Repair common UTF-8-read-as-Latin-1 mojibake (the three-character
  garble of a curly apostrophe) and normalise typographic quotes, dashes
  and non-breaking spaces to their ASCII forms (default `TRUE`).

- citations:

  Remove bracketed reference numbers (`[12]`, `[3, 4]`, `[1-5]`),
  leading list markers (`1.`, `(3)`, `a)`) and dotted section numbers
  (`3.2.1`) (default `TRUE`).

- urls:

  Remove URLs, `www.` addresses, `doi:` strings and bare DOIs (default
  `TRUE`).

- copyright:

  Remove a trailing copyright notice: everything from the last copyright
  sign, `(c)` or `Copyright` to the end when that tail is under
  `copyright_max` characters, plus "All rights reserved" (default
  `TRUE`).

- copyright_max:

  Longest tail treated as a notice (default `300`).

- numbers:

  Remove bare numbers, percentages and years (default `TRUE`). Numbers
  never enter a text hypergraph's vocabulary anyway (the tokeniser keeps
  alphabetic tokens), so this matters for display and for `min_content`.

- remove:

  Extra patterns to remove (Perl regular expressions, case-insensitive),
  applied after the repairs and before the content floor: a database's
  placeholder for a missing abstract, a publisher's own boilerplate.
  Default `NULL`.

- stop_words:

  Words to remove, case-insensitively at word boundaries (default
  `NULL`, none). Usually left to
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)'s
  own `stop_words`; use it here when the cleaned text itself is shown.

- min_chars:

  Minimum number of characters for a word to be kept (default `0L`, keep
  everything). Raise it to drop the single letters and short fragments
  left by initials, enumerations and hyphenated line breaks;
  `min_chars = 3` keeps words of three characters or more.

- min_content:

  Minimum share of alphabetic characters among non-space characters,
  judged on the repaired text before numbers are removed (default `0`,
  keep everything). A reference list or a table fragment scores low;
  `0.5` is a sensible floor for abstracts.

## Value

`x` with the text cleaned: a character vector of the same length (names
kept), or the same data.frame with `column` replaced. `NA` and
low-content texts become `""`. Raises `hypernets_bad_input` for a
malformed `x`, `column` or switch.

## Examples

``` r
messy <- c(
  "1. The programme raised attainment in 2020, see section 3.2.1.",
  "Peer&nbsp;support <b>improved</b> outcomes [12] across schools.",
  "Full report at https://example.org/study.pdf and doi:10.1/x.",
  "Students' views shifted. (c) 2021 Informa UK Limited.",
  "(3) OJ No L 297, 24.11.1979, p. 1."
)
clean_text(messy)
#> [1] "The programme raised attainment in, see section."
#> [2] "Peer support improved outcomes across schools."  
#> [3] "Full report at and."                             
#> [4] "Students' views shifted."                        
#> [5] "OJ No L, p.."                                    
clean_text(messy, min_content = 0.5)
#> [1] "The programme raised attainment in, see section."
#> [2] "Peer support improved outcomes across schools."  
#> [3] "Full report at and."                             
#> [4] "Students' views shifted."                        
#> [5] ""                                                
clean_text(c("[No abstract available]", "A real abstract."),
           remove = "no abstract available")
#> [1] ""                 "A real abstract."
```
