# Packages in Suggests must never be touched while the namespace loads:
# R CMD INSTALL lazy-loads every top-level expression, so a top-level
# `torch::nn_module(...)` makes hypernets uninstallable wherever torch is
# absent (this broke CI on 2026-09-07). Every top-level expression that
# references a Suggests package must be a function definition.

.suggests_packages <- function() {
  desc <- read.dcf(system.file("DESCRIPTION", package = "hypernets"),
                   fields = "Suggests")
  pkgs <- unlist(strsplit(desc[[1L]], ",\\s*"))
  trimws(sub("\\s*\\(.*\\)$", "", pkgs))
}

.top_level_suggests_refs <- function(path, pkgs) {
  exprs <- parse(path, keep.source = FALSE)
  is_function_definition <- vapply(exprs, function(e) {
    is.call(e) && is.symbol(e[[1L]]) &&
      as.character(e[[1L]]) %in% c("<-", "=") &&
      is.call(e[[3L]]) && identical(e[[3L]][[1L]], as.symbol("function"))
  }, logical(1L))
  offenders <- vapply(exprs[!is_function_definition], function(e) {
    any(all.names(e) %in% pkgs)
  }, logical(1L))
  sum(offenders)
}

test_that("no top-level expression references a Suggests package", {
  src_dir <- testthat::test_path("..", "..", "R")
  skip_if_not(dir.exists(src_dir), "package source R/ not available")
  pkgs <- .suggests_packages()
  files <- list.files(src_dir, pattern = "\\.R$", full.names = TRUE)
  counts <- vapply(files, .top_level_suggests_refs, integer(1L), pkgs = pkgs)
  bad <- basename(files)[counts > 0L]
  expect_length(bad, 0L)
})

test_that("the detector itself catches a top-level torch call", {
  tmp <- tempfile(fileext = ".R")
  on.exit(unlink(tmp), add = TRUE, after = FALSE)
  writeLines(c(
    ".ok <- function() torch::nn_module()",
    ".bad <- torch::nn_module()"
  ), tmp)
  expect_identical(.top_level_suggests_refs(tmp, "torch"), 1L)
})
