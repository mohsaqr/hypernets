testthat::skip_on_cran()

test_that("the Legal Hypergraphs datasets have the documented shape", {
  expect_identical(dim(icsid_tribunals), c(2226L, 10L))
  expect_identical(names(icsid_tribunals),
                   c("case", "seat", "arbitrator", "registered", "constituted",
                     "concluded", "is_concluded", "economic_sector", "subject",
                     "respondent"))
  expect_s3_class(icsid_tribunals$constituted, "Date")
  expect_identical(length(unique(icsid_tribunals$case)), 742L)
  expect_identical(length(unique(icsid_tribunals$arbitrator)), 441L)
  expect_identical(sum(is.na(icsid_tribunals$concluded)), 207L * 3L)

  expect_identical(dim(gfcc_decisions), c(3618L, 5L))
  expect_identical(names(gfcc_decisions), c("decision", "date", "citation", "type", "title"))
  expect_false(anyDuplicated(gfcc_decisions$decision) > 0)

  expect_identical(dim(gfcc_citations), c(77284L, 5L))
  expect_identical(names(gfcc_citations),
                   c("citing", "cited", "block", "date_citing", "date_cited"))
  expect_identical(length(unique(gfcc_citations$block)), 46257L)
  expect_true(all(gfcc_citations$citing %in% gfcc_decisions$decision))
  expect_true(all(gfcc_citations$cited %in% gfcc_decisions$decision))
})

test_that("the datasets build the paper hypergraphs directly", {
  tribunals <- temporal_hypergraph(icsid_tribunals, actor = "arbitrator",
                                   cooccur_by = "case", start = "constituted",
                                   end = "concluded")
  expect_identical(length(tribunals$nodes), 441L)
  expect_identical(length(tribunals$edges), 742L)
  expect_identical(tribunals$evolution, "interval")
  expect_true("economic_sector" %in% names(tribunals$edge_data))
  expect_false("seat" %in% names(tribunals$edge_data))
  blocks <- temporal_hypergraph(gfcc_citations, actor = "cited", cooccur_by = "block",
                                time = "date_citing", nodes = gfcc_decisions,
                                sparse = TRUE)
  expect_identical(length(blocks$nodes), 3618L)
  expect_identical(length(blocks$edges), 46257L)
  expect_true("citing" %in% names(blocks$edge_data))
  aggregate <- hypergraph_snapshot(blocks, mode = "all")
  expect_identical(nrow(hg_subset(aggregate, where = c(citing = "153-001"))$edge_data), 68L)
})
