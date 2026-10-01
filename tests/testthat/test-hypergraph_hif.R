testthat::skip_on_cran()
testthat::skip_if_not_installed("jsonlite")

# g1 = {a, b, c} with weights 1, 2, 3; g2 = {a, d} with 0.5, 1; z isolated.
.hif_fixture <- function(sparse = FALSE, weighted = TRUE) {
  memberships <- data.frame(
    actor = c("a", "b", "c", "a", "d"),
    group = c("g1", "g1", "g1", "g2", "g2"),
    w = c(1, 2, 3, 0.5, 1),
    year = c(2001, 2001, 2001, 2002, 2002),
    stringsAsFactors = FALSE
  )
  group_hypergraph(memberships, actor = "actor", group = "group",
                   weight = if (weighted) "w" else NULL,
                   nodes = c("a", "b", "c", "d", "z"), sparse = sparse)
}

.hif_core <- c("hyperedges", "incidence", "nodes", "n_nodes",
               "n_hyperedges", "size_distribution")

test_that("a weighted group hypergraph round-trips exactly", {
  hg <- .hif_fixture()
  back <- hg_read_hif(hg_write_hif(hg))
  expect_s3_class(back, "net_hg")
  expect_identical(back[.hif_core], hg[.hif_core])
  expect_identical(back$edge_data, hg$edge_data)
  expect_identical(back$params$source, "hg_read_hif")
  expect_identical(back$params$metadata$source, "group_hypergraph")
})

test_that("a binary incidence stays integer and is written without weights", {
  hg <- .hif_fixture(weighted = FALSE)
  json <- hg_write_hif(hg)
  expect_false(grepl("\"weight\"", json, fixed = TRUE))
  back <- hg_read_hif(json)
  expect_identical(typeof(back$incidence), "integer")
  expect_identical(back[.hif_core], hg[.hif_core])
})

test_that("window, network, text and sparse hypergraphs round-trip", {
  hw <- window_hypergraph(list(s1 = c("a", "b", "a", "c", "b")), window = 2L)
  bw <- hg_read_hif(hg_write_hif(hw))
  expect_identical(bw[.hif_core], hw[.hif_core])
  expect_identical(bw$window_counts, as.numeric(hw$window_counts))

  adj <- matrix(c(0, 1, 1, 0, 1, 0, 1, 0, 1, 1, 0, 1, 0, 0, 1, 0), 4, 4,
                dimnames = list(letters[1:4], letters[1:4]))
  hn <- network_hypergraph(adj, max_size = 3L)
  expect_identical(hg_read_hif(hg_write_hif(hn))[.hif_core], hn[.hif_core])

  th <- text_hypergraph(c(d1 = "the cat sat at night",
                          d2 = "the dog ran at night"))
  bt <- hg_read_hif(hg_write_hif(th))
  expect_identical(bt[.hif_core], th[.hif_core])
  expect_equal(bt$edge_data, th$edge_data)   # integer counts come back double

  hs <- .hif_fixture(sparse = TRUE)
  bs <- hg_read_hif(hg_write_hif(hs), sparse = TRUE)
  expect_s4_class(bs$incidence, "dgCMatrix")
  expect_identical(bs[.hif_core], hs[.hif_core])
})

test_that("doubles, node weights and incidence attributes survive", {
  hg <- .hif_fixture()
  hg$incidence[hg$incidence != 0] <- c(1 / 3, 0.1 + 0.2, pi, exp(1), 1e-10)
  hg$node_data <- data.frame(node = hg$nodes, weight = c(1, 2, 3, 4, 0.25),
                             role = c("x", "y", NA, "x", "z"),
                             stringsAsFactors = FALSE)
  hg$incidence_data <- data.frame(node = c("a", "d"), edge = c("g2", "g2"),
                                  seat = c(1, 2), stringsAsFactors = FALSE)
  back <- hg_read_hif(hg_write_hif(hg))
  expect_identical(back$incidence, hg$incidence)
  expect_identical(back$node_data, hg$node_data)
  expect_identical(back$incidence_data, hg$incidence_data)
})

test_that("files work and read -> write is a fixed point", {
  hg <- .hif_fixture()
  path <- tempfile(fileext = ".json")
  on.exit(unlink(path), add = TRUE)
  expect_identical(hg_write_hif(hg, file = path), path)
  first <- hg_read_hif(path)
  json_1 <- hg_write_hif(first)
  json_2 <- hg_write_hif(hg_read_hif(json_1))
  expect_identical(json_1, json_2)
  expect_identical(first[.hif_core], hg[.hif_core])
})

test_that("reading is invariant to the order of records", {
  hg <- .hif_fixture()
  doc <- jsonlite::fromJSON(hg_write_hif(hg), simplifyVector = FALSE)
  doc$incidences <- rev(doc$incidences)
  doc$nodes <- rev(doc$nodes)
  doc$edges <- rev(doc$edges)
  shuffled <- hg_read_hif(as.character(jsonlite::toJSON(
    doc, auto_unbox = TRUE, digits = I(17))))
  expect_identical(shuffled$nodes, rev(hg$nodes))
  expect_identical(shuffled$incidence[hg$nodes, colnames(hg$incidence)],
                   hg$incidence)
  measures <- hg_measures(shuffled)
  measures <- measures[match(hg$nodes, measures$node), ]
  rownames(measures) <- NULL
  expect_equal(measures, hg_measures(hg))
})

test_that("a hand-written HIF file gives the hand-computed hypergraph", {
  # XGI-style: integer identifiers, no weights, an attribute on one node,
  # an isolated node and an empty edge (both only in the record arrays).
  json <- paste0(
    '{"network-type": "asc", "metadata": {"name": "toy"},',
    ' "nodes": [{"node": 0, "attrs": {"colour": "red"}}, {"node": 9}],',
    ' "edges": [{"edge": "e0", "weight": 2}, {"edge": "empty"}],',
    ' "incidences": [{"edge": "e0", "node": 0}, {"edge": "e0", "node": 1},',
    ' {"edge": "e1", "node": 1}, {"edge": "e1", "node": 2}]}')
  hg <- hg_read_hif(json)
  expected <- matrix(c(1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L),
                     4, 3, dimnames = list(c("0", "9", "1", "2"),
                                           c("e0", "empty", "e1")))
  expect_identical(hg$incidence, expected)
  expect_identical(hg$window_counts, c(2, 1, 1))
  expect_identical(hg$node_data$colour, c("red", NA, NA, NA))
  expect_identical(lengths(hg$hyperedges), c(2L, 0L, 2L))
  expect_identical(hg$params$network_type, "asc")
  expect_identical(hg$params$metadata, list(name = "toy"))
})

test_that("malformed or unsupported HIF raises hypernets_bad_input", {
  bad <- c(
    not_json = "{not json",
    no_incidences = '{"nodes": []}',
    unknown_top = '{"incidences": [], "extra": 1}',
    bad_type = '{"network-type": "mixed", "incidences": []}',
    directed = '{"network-type": "directed", "incidences": [{"edge": 1, "node": 1}]}',
    direction = '{"incidences": [{"edge": 1, "node": 1, "direction": "head"}]}',
    no_node = '{"incidences": [{"edge": 1}]}',
    fractional_id = '{"incidences": [{"edge": 1.5, "node": 1}]}',
    unknown_field = '{"incidences": [{"edge": 1, "node": 1, "colour": "red"}]}',
    text_weight = '{"incidences": [{"edge": 1, "node": 1, "weight": "heavy"}]}',
    zero_weight = '{"incidences": [{"edge": 1, "node": 1, "weight": 0}]}',
    duplicate_pair = '{"incidences": [{"edge": 1, "node": 1}, {"edge": 1, "node": 1}]}',
    duplicate_node = '{"incidences": [], "nodes": [{"node": 1}, {"node": 1}]}',
    attr_clash = '{"incidences": [], "nodes": [{"node": 1, "attrs": {"weight": 2}}]}'
  )
  lapply(names(bad), \(case) {
    expect_error(hg_read_hif(bad[[case]]), class = "hypernets_bad_input",
                 label = case)
  })
  expect_error(hg_read_hif(tempfile(fileext = ".json")),
               class = "hypernets_bad_input")
  expect_error(hg_write_hif(list(incidence = diag(2))),
               class = "hypernets_bad_input")
})

test_that("hg_get() reaches the HIF node and incidence tables", {
  hg <- .hif_fixture()
  hg$node_data <- data.frame(node = hg$nodes, weight = c(1, 2, 3, 4, 0.25),
                             stringsAsFactors = FALSE)
  hg$incidence_data <- data.frame(node = c("a", "d"), edge = c("g2", "g2"),
                                  seat = c(1, 2), stringsAsFactors = FALSE)
  back <- hg_read_hif(hg_write_hif(hg))
  expect_identical(hg_get(back, what = "node_data"), hg$node_data)
  expect_identical(hg_get(back, what = "incidence_data"), hg$incidence_data)
  expect_identical(nrow(hg_get(back, what = "node_data", top = 2L)), 2L)
  plain <- group_hypergraph(data.frame(p = c("a", "b"), g = c("x", "x")),
                            actor = "p", group = "g")
  expect_error(hg_get(plain, what = "node_data"), class = "hypernets_bad_input")
})

test_that("missing attribute values are omitted, not written as null", {
  hg <- .hif_fixture()
  hg$node_data <- data.frame(node = hg$nodes,
                             role = c("x", NA, NA, "y", NA),
                             stringsAsFactors = FALSE)
  json <- hg_write_hif(hg)
  expect_false(grepl("null", json, fixed = TRUE))
  back <- hg_read_hif(json)
  expect_identical(hg_get(back, what = "node_data"), hg$node_data)
  expect_identical(hg_write_hif(back), hg_write_hif(hg_read_hif(hg_write_hif(back))))
})

test_that("edge attributes are reachable and a read hypergraph says its source", {
  hg <- .hif_fixture()
  back <- hg_read_hif(hg_write_hif(hg))
  expect_identical(hg_get(back, what = "edge_data"), back$edge_data)
  expect_output(print(back), "Hypergraph Interchange Format")
})
