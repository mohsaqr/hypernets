# ---- Hypergraph Interchange Format (HIF) -----------------------------------
# Read and write the JSON interchange standard shared by XGI, HyperNetX,
# HypergraphX, HAT and SimpleHypergraphs.jl (Coll et al. 2025). A net_hg maps
# onto HIF as: one incidence record per non-zero cell of the incidence matrix,
# one node record per row, one edge record per column. jsonlite is a Suggests
# dependency, so both verbs check for it before touching JSON.

#' Read and write the Hypergraph Interchange Format (HIF)
#'
#' `hg_write_hif()` writes a `net_hg` as HIF JSON, and `hg_read_hif()` reads
#' a HIF file (or a HIF JSON string) back into a `net_hg`. HIF is the
#' interchange standard of the higher-order network libraries (XGI,
#' HyperNetX, HypergraphX, HAT, SimpleHypergraphs.jl), so these two verbs move
#' a hypergraph between hypernets and any of them.
#'
#' @section Mapping:
#' \describe{
#'   \item{incidences}{One record per non-zero cell of the incidence matrix,
#'     `{"node", "edge", "weight"}`. `weight` is the incidence value; it is
#'     left out only for an integer 0/1 incidence (a binary hypergraph). On
#'     reading, an absent weight is 1, and a file whose incidences carry no
#'     weight gives an integer 0/1 incidence matrix, otherwise a double one,
#'     so the storage mode round-trips.}
#'   \item{nodes}{One record per node, so isolated nodes survive. Node
#'     attributes are the columns of the hypergraph's node table (the
#'     `node_data` a previous `hg_read_hif()` stored) plus the planted
#'     `block` of [hg_sample_sbm()]. A column called `weight` is written as
#'     the node's HIF `weight`.}
#'   \item{edges}{One record per hyperedge, named by the incidence column
#'     name. The hyperedge weight (the window counts of
#'     [window_hypergraph()], which the Laplacian and modularity verbs read
#'     as hyperedge weights) is written as HIF `weight`; the hyperedge
#'     attribute table (`edge_data`: constant-within-group columns of
#'     [group_hypergraph()], word counts of [text_hypergraph()], a previous
#'     read) and the multiplicities of [hg_snapshot()] are written as
#'     `attrs`.}
#'   \item{metadata}{The scalar construction parameters and the writing
#'     package version; a hypergraph that came from `hg_read_hif()` writes its
#'     original metadata back unchanged.}
#' }
#'
#' @section What does not round-trip:
#' \itemize{
#'   \item Node and edge identifiers are strings: an integer identifier from
#'     another library (XGI's default `0, 1, 2, ...`) is read as `"0"`,
#'     `"1"`, `"2"`.
#'   \item JSON has one number type, so numeric attributes come back as
#'     double; factors, dates and times are written as their character form
#'     and come back as character.
#'   \item The construction record (`params`) is replaced by the source of
#'     the read; the constructor-specific fields that are not hypergraph data
#'     (the `$text` layer of [text_hypergraph()], the clustered-sequence
#'     tables of [group_hypergraph()]) are not written. `block` and
#'     `multiplicity` come back as node and edge attribute columns rather than
#'     as the `blocks` and `edge_multiplicity` fields.
#'   \item Directed hypergraphs (`"network-type": "directed"` or an incidence
#'     `direction`) have no `net_hg` counterpart and are refused.
#'   \item A zero incidence weight cannot be stored (membership in a `net_hg`
#'     is a non-zero cell) and is refused, as is a repeated (node, edge) pair,
#'     whose meaning HIF leaves open.
#' }
#' Weights are written with 17 significant digits, so doubles round-trip
#' exactly.
#'
#' @param hg A `net_hg` (from [group_hypergraph()], [window_hypergraph()],
#'   [text_hypergraph()], [network_hypergraph()], `hg_read_hif()`, ...).
#' @param file For `hg_write_hif()`: the path to write, or `NULL` (default) to
#'   return the JSON text. For `hg_read_hif()`: the path of a HIF file, or a
#'   HIF JSON string.
#' @param pretty Logical. Indent the JSON? Default `FALSE`.
#' @param sparse Logical. Store the incidence of the read hypergraph as a
#'   sparse `dgCMatrix`? Default `FALSE`.
#' @return `hg_write_hif()`: with `file = NULL`, the HIF JSON as a single
#'   character string; otherwise `file`, invisibly. `hg_read_hif()`: a
#'   `net_hg` with the usual fields (`incidence` node x hyperedge with node
#'   and hyperedge names as dimnames, `hyperedges`, `nodes`, `n_nodes`,
#'   `n_hyperedges`, `size_distribution`, `params`) plus, when the file has
#'   them, `window_counts` (hyperedge weights; an edge record without a
#'   weight counts 1), `edge_data` (one row per hyperedge: `edge` and one
#'   column per edge attribute), `node_data` (one row per node: `node`,
#'   `weight` when the file has node weights, and one column per node
#'   attribute) and `incidence_data` (one row per incidence with attributes:
#'   `node`, `edge` and the attribute columns). Nodes and hyperedges are in
#'   file order: the `nodes` / `edges` records first, then those first seen
#'   in `incidences`. `params` holds `source = "hg_read_hif"`,
#'   `network_type` and the file's `metadata`.
#'
#'   Both verbs raise `hypernets_missing_dependency` without jsonlite.
#'   `hg_read_hif()` raises `hypernets_bad_input` for text that is not JSON,
#'   a missing file, and a document that breaks the HIF schema (no
#'   `incidences`, a record without its identifier, an unknown field, a
#'   non-numeric weight, an unknown `network-type`), a duplicated node,
#'   edge or (node, edge) pair, a zero incidence weight, a directed
#'   hypergraph, or an attribute named like its record's identifier.
#'   `hg_write_hif()` raises `hypernets_bad_input` for a non-`net_hg` input.
#'
#' @examples
#' if (requireNamespace("jsonlite", quietly = TRUE)) {
#'   memberships <- data.frame(
#'     actor = c("a", "b", "c", "a", "d"),
#'     group = c("g1", "g1", "g1", "g2", "g2"),
#'     weight = c(1, 2, 3, 0.5, 1),
#'     year = c(2001, 2001, 2001, 2002, 2002)
#'   )
#'   hg <- group_hypergraph(memberships, actor = "actor", group = "group",
#'                          weight = "weight")
#'   path <- tempfile(fileext = ".json")
#'   hg_write_hif(hg, file = path)
#'   back <- hg_read_hif(path)
#'   hg_get(back, what = "memberships")
#'   hg_write_hif(hg, pretty = TRUE)
#' }
#'
#' @references
#' Coll, M., Joslyn, C. A., Landry, N. W., Lotito, Q. F., Myers, A.,
#' Pickard, J., Praggastis, B., & Szufel, P. (2025). HIF: The hypergraph
#' interchange format for higher-order networks. \emph{Network Science} 13,
#' e21. \doi{10.1017/nws.2025.10018}
#'
#' @export
hg_write_hif <- function(hg, file = NULL, pretty = FALSE) {
  .hif_need_jsonlite("hg_write_hif")
  .thg_check_hg(hg)
  stopifnot(
    "`file` must be NULL or a single path" =
      is.null(file) ||
        (is.character(file) && length(file) == 1L && !is.na(file) && nzchar(file)),
    "`pretty` must be TRUE or FALSE" =
      is.logical(pretty) && length(pretty) == 1L && !is.na(pretty)
  )
  node_names <- as.character(hg$nodes %||% rownames(hg$incidence) %||%
                               paste0("n", seq_len(hg$n_nodes)))
  edge_names <- as.character(colnames(hg$incidence) %||%
                               paste0("h", seq_len(hg$n_hyperedges)))

  # incidences: one record per non-zero cell, in hyperedge order
  cells <- Matrix::which(hg$incidence != 0, arr.ind = TRUE)
  cells <- cells[order(cells[, 2L], cells[, 1L]), , drop = FALSE]
  incidences <- data.frame(edge = edge_names[cells[, 2L]],
                           node = node_names[cells[, 1L]],
                           stringsAsFactors = FALSE)
  weights <- as.numeric(hg$incidence[cells])
  # an integer 0/1 incidence is a binary hypergraph and needs no weights;
  # a double one keeps them, so the storage mode round-trips
  if (!is.integer(hg$incidence) || any(weights != 1)) {
    incidences$weight <- weights
  }
  incidence_attrs <- .hif_attr_table(hg$incidence_data,
                                     keys = c("node", "edge"),
                                     ids = list(incidences$node, incidences$edge))
  if (!is.null(incidence_attrs)) incidences$attrs <- incidence_attrs

  # nodes: every node, so isolated nodes survive
  node_table <- hg$node_data
  if (!is.null(hg$blocks)) {
    node_table <- node_table %||% data.frame(node = node_names,
                                             stringsAsFactors = FALSE)
    node_table$block <- as.integer(hg$blocks)
  }
  nodes <- data.frame(node = node_names, stringsAsFactors = FALSE)
  if (!is.null(node_table) && "weight" %in% names(node_table)) {
    node_row <- match(node_names, as.character(node_table$node))
    nodes$weight <- as.numeric(node_table$weight[node_row])
    node_table$weight <- NULL
  }
  node_attrs <- .hif_attr_table(node_table, keys = "node", ids = list(node_names))
  if (!is.null(node_attrs)) nodes$attrs <- node_attrs

  # edges: every hyperedge, with its weight and attributes
  edges <- data.frame(edge = edge_names, stringsAsFactors = FALSE)
  if (!is.null(hg$window_counts)) edges$weight <- as.numeric(hg$window_counts)
  edge_table <- hg$edge_data
  if (!is.null(hg$edge_multiplicity)) {
    edge_table <- edge_table %||% data.frame(edge = edge_names,
                                             stringsAsFactors = FALSE)
    edge_table$multiplicity <- as.integer(hg$edge_multiplicity)
  }
  edge_attrs <- .hif_attr_table(edge_table, keys = "edge", ids = list(edge_names))
  if (!is.null(edge_attrs)) edges$attrs <- edge_attrs

  document <- list(
    `network-type` = "undirected",
    metadata = .hif_metadata(hg),
    incidences = incidences,
    nodes = nodes,
    edges = edges
  )
  json <- jsonlite::toJSON(document, dataframe = "rows", auto_unbox = TRUE,
                           digits = I(17), na = "null", null = "null",
                           pretty = pretty)
  json <- as.character(json)
  if (is.null(file)) return(json)
  writeLines(json, file, useBytes = TRUE)
  invisible(file)
}

#' @rdname hg_write_hif
#' @export
hg_read_hif <- function(file, sparse = FALSE) {
  .hif_need_jsonlite("hg_read_hif")
  stopifnot(
    "`file` must be a single path or a HIF JSON string" =
      is.character(file) && length(file) == 1L && !is.na(file),
    "`sparse` must be TRUE or FALSE" =
      is.logical(sparse) && length(sparse) == 1L && !is.na(sparse)
  )
  is_text <- grepl("^\\s*\\{", file)
  if (!is_text && !file.exists(file)) {
    .thg_bad_input(sprintf("`file` is neither an existing file nor HIF JSON text: %s",
                           file))
  }
  text <- if (is_text) file else paste(readLines(file, warn = FALSE,
                                                 encoding = "UTF-8"),
                                       collapse = "\n")
  document <- tryCatch(
    jsonlite::fromJSON(text, simplifyVector = TRUE, simplifyMatrix = FALSE),
    error = function(e) {
      .thg_bad_input(sprintf("`file` is not valid JSON: %s", conditionMessage(e)))
    }
  )
  .hif_validate_document(document)
  network_type <- document[["network-type"]] %||% "undirected"

  incidences <- .hif_records(document$incidences, "incidences",
                             ids = c("edge", "node"),
                             fields = c("edge", "node", "weight", "direction",
                                        "attrs"))
  if (identical(network_type, "directed") || "direction" %in% names(incidences)) {
    .thg_bad_input("directed HIF (network-type \"directed\" or an incidence `direction`) has no net_hg counterpart; only undirected hypergraphs can be read")
  }
  node_records <- .hif_records(document$nodes, "nodes", ids = "node",
                               fields = c("node", "weight", "attrs"))
  edge_records <- .hif_records(document$edges, "edges", ids = "edge",
                               fields = c("edge", "weight", "attrs"))
  .hif_unique(node_records$node, "node identifier in `nodes`")
  .hif_unique(edge_records$edge, "edge identifier in `edges`")
  .hif_unique(paste(incidences$node, incidences$edge, sep = "\x01"),
              "(node, edge) pair in `incidences`")

  node_names <- unique(c(node_records$node, incidences$node))
  edge_names <- unique(c(edge_records$edge, incidences$edge))
  n <- length(node_names)
  m <- length(edge_names)

  has_weight <- "weight" %in% names(incidences)
  weights <- rep(1L, nrow(incidences))
  if (has_weight) {
    weights <- as.numeric(incidences$weight)
    weights[is.na(weights)] <- 1
  }
  if (any(weights == 0)) {
    .thg_bad_input("an incidence has weight 0; a net_hg records membership as a non-zero cell, so a zero-weight incidence cannot be stored")
  }
  i <- match(incidences$node, node_names)
  j <- match(incidences$edge, edge_names)
  if (sparse) {
    incidence <- Matrix::sparseMatrix(i = i, j = j, x = as.numeric(weights),
                                      dims = c(n, m),
                                      dimnames = list(node_names, edge_names))
  } else {
    incidence <- matrix(if (has_weight) 0 else 0L, n, m,
                        dimnames = list(node_names, edge_names))
    incidence[cbind(i, j)] <- weights
  }
  if (m == 0L) colnames(incidence) <- NULL

  hg <- .thg_rebuild(list(), incidence, keep_edges = seq_len(m))
  hg <- hg[c("hyperedges", "incidence", "nodes", "n_nodes", "n_hyperedges",
             "size_distribution")]
  hg$params <- list(source = "hg_read_hif", network_type = network_type,
                    sparse = sparse, metadata = document$metadata %||% list())

  if ("weight" %in% names(edge_records)) {
    edge_weight <- edge_records$weight[match(edge_names, edge_records$edge)]
    edge_weight[is.na(edge_weight)] <- 1
    hg$window_counts <- as.numeric(edge_weight)
  }
  edge_data <- .hif_attr_frame(edge_records, key = "edge", ids = edge_names)
  if (!is.null(edge_data)) hg$edge_data <- edge_data
  node_data <- .hif_attr_frame(node_records, key = "node", ids = node_names)
  if ("weight" %in% names(node_records)) {
    node_weight <- node_records$weight[match(node_names, node_records$node)]
    node_weight[is.na(node_weight)] <- 1
    node_data <- node_data %||% data.frame(node = node_names,
                                           stringsAsFactors = FALSE)
    node_data <- data.frame(node_data["node"], weight = as.numeric(node_weight),
                            node_data[setdiff(names(node_data), "node")],
                            stringsAsFactors = FALSE, check.names = FALSE)
  }
  if (!is.null(node_data)) hg$node_data <- node_data
  if ("attrs" %in% names(incidences)) {
    incidence_data <- .hif_attr_frame(incidences, key = c("node", "edge"))
    if (!is.null(incidence_data)) {
      # one row per incidence that carries an attribute
      attr_columns <- setdiff(names(incidence_data), c("node", "edge"))
      carried <- Reduce(`|`, lapply(incidence_data[attr_columns], \(v) !is.na(v)))
      incidence_data <- incidence_data[carried, , drop = FALSE]
      rownames(incidence_data) <- NULL
      hg$incidence_data <- incidence_data
    }
  }
  structure(hg, class = "net_hg")
}

# ---- internal helpers --------------------------------------------------------

.hif_need_jsonlite <- function(fn) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(errorCondition(
      sprintf("`%s()` needs the suggested package `jsonlite`", fn),
      class = "hypernets_missing_dependency", call = NULL
    ))
  }
}

# Top level of the HIF schema: an object with `incidences`, an optional
# `network-type` from the enum, and no other fields.
.hif_validate_document <- function(document) {
  if (!is.list(document) || is.data.frame(document) ||
      is.null(names(document))) {
    .thg_bad_input("HIF must be a JSON object")
  }
  unknown <- setdiff(names(document),
                     c("network-type", "metadata", "incidences", "nodes", "edges"))
  if (length(unknown)) {
    .thg_bad_input(sprintf("HIF has unknown top-level field(s): %s",
                           paste(unknown, collapse = ", ")))
  }
  if (!"incidences" %in% names(document)) {
    .thg_bad_input("HIF requires an `incidences` array")
  }
  type <- document[["network-type"]]
  if (!is.null(type) && !(is.character(type) && length(type) == 1L &&
                          type %in% c("undirected", "directed", "asc"))) {
    .thg_bad_input("`network-type` must be \"undirected\", \"directed\" or \"asc\"")
  }
  if (!is.null(document$metadata) && !is.list(document$metadata)) {
    .thg_bad_input("HIF `metadata` must be a JSON object")
  }
  invisible(document)
}

# One HIF record array as a data.frame with character identifiers, checked
# against the schema's field list, identifier types and numeric weights.
.hif_records <- function(records, section, ids, fields) {
  if (is.null(records) || (is.list(records) && !is.data.frame(records) &&
                           length(records) == 0L)) {
    out <- as.data.frame(stats::setNames(
      replicate(length(ids), character(0L), simplify = FALSE), ids),
      stringsAsFactors = FALSE)
    return(out)
  }
  if (!is.data.frame(records)) {
    .thg_bad_input(sprintf("HIF `%s` must be an array of objects", section))
  }
  unknown <- setdiff(names(records), fields)
  if (length(unknown)) {
    .thg_bad_input(sprintf("HIF `%s` records have unknown field(s): %s",
                           section, paste(unknown, collapse = ", ")))
  }
  lapply(ids, \(id) {
    value <- records[[id]]
    if (is.null(value) || anyNA(value) || is.list(value) ||
        !(is.character(value) || is.numeric(value)) ||
        (is.numeric(value) && any(value != round(value)))) {
      .thg_bad_input(sprintf(
        "every HIF `%s` record needs `%s` as a string or an integer", section, id))
    }
  })
  records[ids] <- lapply(records[ids], \(value) {
    if (is.numeric(value)) format(value, scientific = FALSE, trim = TRUE) else value
  })
  if ("weight" %in% names(records) && !is.numeric(records$weight)) {
    .thg_bad_input(sprintf("HIF `%s` weights must be numbers", section))
  }
  if ("direction" %in% names(records) &&
      !all(records$direction[!is.na(records$direction)] %in% c("head", "tail"))) {
    .thg_bad_input("HIF incidence `direction` must be \"head\" or \"tail\"")
  }
  if ("attrs" %in% names(records) && !is.data.frame(records$attrs) &&
      !all(vapply(records$attrs, \(a) is.null(a) || is.list(a), logical(1L)))) {
    .thg_bad_input(sprintf("HIF `%s` attrs must be JSON objects", section))
  }
  records
}

.hif_unique <- function(values, what) {
  if (anyDuplicated(values)) {
    .thg_bad_input(sprintf("HIF has a duplicated %s", what))
  }
}

# The `attrs` of a record table as a data.frame keyed by `key` (NULL when no
# record carries an attribute). With `ids`, rows follow `ids` and records
# without a row get NA.
.hif_attr_frame <- function(records, key, ids = NULL) {
  attrs <- records$attrs
  if (is.null(attrs)) return(NULL)
  if (!is.data.frame(attrs)) {
    # a list of objects (records with and without attrs mixed as lists)
    attrs <- jsonlite::fromJSON(jsonlite::toJSON(
      lapply(attrs, \(a) a %||% stats::setNames(list(), character(0L))),
      auto_unbox = TRUE, digits = I(17)), simplifyMatrix = FALSE)
    if (!is.data.frame(attrs)) return(NULL)
  }
  if (!ncol(attrs)) return(NULL)
  clash <- intersect(names(attrs), c(key, if (identical(key, "node")) "weight"))
  if (length(clash)) {
    .thg_bad_input(sprintf(
      "a HIF attribute is named like its record's field (%s); rename it before reading",
      paste(clash, collapse = ", ")))
  }
  attrs[] <- lapply(attrs, \(column) {
    if (is.integer(column)) as.numeric(column) else column
  })
  out <- data.frame(records[key], attrs, stringsAsFactors = FALSE,
                    check.names = FALSE)
  if (!is.null(ids)) {
    out <- out[match(ids, out[[key]]), , drop = FALSE]
    out[[key]] <- ids
  }
  rownames(out) <- NULL
  out
}

# A key + attribute table as the nested `attrs` column of a record table,
# aligned with the record identifiers `ids` (NULL when the table has no
# attribute columns). NA cells are dropped from the JSON object.
.hif_attr_table <- function(table, keys, ids) {
  if (is.null(table)) return(NULL)
  columns <- setdiff(names(table), keys)
  if (!length(columns)) return(NULL)
  row_key <- do.call(paste, c(lapply(table[keys], as.character), sep = "\x01"))
  id_key <- do.call(paste, c(lapply(ids, as.character), sep = "\x01"))
  attrs <- table[match(id_key, row_key), columns, drop = FALSE]
  attrs[] <- lapply(attrs, \(column) {
    if (is.factor(column) || inherits(column, c("Date", "POSIXt"))) {
      as.character(column)
    } else {
      column
    }
  })
  rownames(attrs) <- NULL
  # one named list per record with its NA cells dropped, so a record without
  # a value omits the key instead of writing `null` (a read -> write of a
  # file where only some records carry an attribute is then a fixed point)
  lapply(seq_len(nrow(attrs)), \(i) {
    row <- lapply(attrs, \(column) column[[i]])
    kept <- row[!vapply(row, \(v) length(v) == 1L && is.na(v), logical(1L))]
    if (length(kept)) kept else stats::setNames(list(), character(0L))
  })
}

# Network-level metadata: the file's own metadata for a hypergraph that came
# from hg_read_hif(), otherwise the scalar construction parameters.
.hif_metadata <- function(hg) {
  params <- hg$params %||% list()
  if (identical(params$source, "hg_read_hif")) {
    metadata <- params$metadata %||% list()
    if (!length(metadata)) metadata <- stats::setNames(list(), character(0L))
    return(metadata)
  }
  scalar <- Filter(\(v) is.atomic(v) && length(v) == 1L && !is.na(v), params)
  c(scalar, list(generator = sprintf("hypernets %s",
                                     utils::packageVersion("hypernets"))))
}
