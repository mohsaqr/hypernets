skip_on_cran()

# non-ASCII input is written as \u escapes so this file is portable:
# \u2019 curly apostrophe, \u00a9 copyright sign, \u00e2\u20ac\u2122 the
# three-character mojibake of a curly apostrophe
messy <- c(
  list = "1. The programme raised attainment in 2020, see section 3.2.1.",
  journal = "(3) OJ No L 297, 24.11.1979, p. 1.",
  html = "Peer&nbsp;support <b>improved</b> outcomes [12] across schools.",
  url = "Full report at https://example.org/study.pdf and doi:10.1/x.",
  notice = paste0("Students\u2019 views shifted. \u00a9 2021 Informa UK ",
                  "Limited, trading as Taylor & Francis Group."),
  mojibake = "Teachers\u00e2\u20ac\u2122 workload rose by 25% [3, 4].",
  entity = "Learning &amp; teaching &#8217;online&#8217; in 2021.",
  rights = "Simulation training works. All rights reserved."
)

test_that("clean_text repairs and strips each class of non-content", {
  out <- clean_text(messy)
  expect_identical(length(out), length(messy))
  expect_identical(names(out), names(messy))
  expect_identical(out[["list"]],
                   "The programme raised attainment in, see section.")
  expect_identical(out[["html"]],
                   "Peer support improved outcomes across schools.")
  expect_identical(out[["url"]], "Full report at and.")
  expect_identical(out[["notice"]], "Students' views shifted.")
  expect_identical(out[["mojibake"]], "Teachers' workload rose by.")
  expect_identical(out[["entity"]], "Learning & teaching 'online' in.")
  expect_identical(out[["rights"]], "Simulation training works.")
  # every result is plain ASCII
  expect_true(all(vapply(out, \(s) all(utf8ToInt(s) < 128), logical(1))))
})

test_that("clean_text switches, stop words, min_content and NA behave", {
  keep_numbers <- clean_text(messy[["list"]], numbers = FALSE)
  expect_true(grepl("2020", keep_numbers, fixed = TRUE))
  expect_false(grepl("3.2.1", keep_numbers, fixed = TRUE))
  raw_html <- clean_text(messy[["html"]], html = FALSE)
  expect_true(grepl("<b>", raw_html, fixed = TRUE))
  kept_copyright <- clean_text("Copyright 2020 Elsevier. Great work.",
                               copyright = FALSE)
  expect_identical(kept_copyright, "Copyright Elsevier. Great work.")
  leading_notice <- clean_text("Copyright 2020 Elsevier. Great work.")
  expect_identical(leading_notice, "")
  trailing_notice <- clean_text("Great work. Copyright 2020 Elsevier.")
  expect_identical(trailing_notice, "Great work.")
  sw <- clean_text("The cat sat on the mat", stop_words = c("the", "on"))
  expect_identical(sw, "cat sat mat")
  # the content floor: the journal citation is nearly all digits/punctuation
  floored <- clean_text(messy, min_content = 0.5)
  expect_identical(unname(floored[["journal"]]), "")
  expect_true(nzchar(floored[["html"]]))
  with_na <- clean_text(c("a", NA_character_))
  expect_identical(with_na, c("a", ""))
  empty <- clean_text(character(0))
  expect_identical(empty, character(0))
  # idempotent: cleaning clean text changes nothing
  once <- clean_text(messy)
  twice <- clean_text(once)
  expect_identical(twice, once)
})

test_that("clean_text keeps data.frame rows aligned and feeds text_hypergraph", {
  df <- data.frame(id = c("a", "b", "c"),
                   text = c(messy[["html"]], messy[["journal"]],
                            messy[["notice"]]),
                   year = c(2020L, 1979L, 2021L))
  out <- clean_text(df, column = "text", min_content = 0.5)
  expect_identical(dim(out), dim(df))
  expect_identical(out$id, df$id)
  expect_identical(out$year, df$year)
  expect_identical(out$text[[2L]], "")
  # the empty row is dropped by the constructor, with its named warning
  expect_warning(hg <- text_hypergraph(out, column = "text", id = "id"),
                 class = "honets_dropped_documents")
  expect_setequal(hg$nodes, c("a", "c"))
  expect_error(clean_text(df), class = "honets_bad_input")
  expect_error(clean_text(df, column = "nope"), class = "honets_bad_input")
  expect_error(clean_text("x", column = "text"), class = "honets_bad_input")
  expect_error(clean_text(1), class = "honets_bad_input")
  expect_error(clean_text("x", html = NA), class = "honets_bad_input")
  expect_error(clean_text("x", min_content = 2), "min_content")
})

test_that("clean_text `remove` patterns empty placeholders before the floor", {
  x <- c("[No abstract available]", "A real abstract about schools.")
  out <- clean_text(x, remove = "no abstract available")
  expect_identical(out, c("", "A real abstract about schools."))
  # without the pattern the placeholder survives as ordinary text
  kept <- clean_text(x)
  expect_identical(kept[[1L]], "[No abstract available]")
  # several patterns, case-insensitive, regex
  out2 <- clean_text("Published by ELSEVIER Ltd; Data: Table 2.",
                     remove = c("published by \\w+ ltd", "table \\d"))
  expect_identical(out2, "Data:.")
  expect_error(clean_text("x", remove = NA_character_), "remove")
})
