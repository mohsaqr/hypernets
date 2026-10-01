# Read and write the Hypergraph Interchange Format (HIF)

`hg_write_hif()` writes a `net_hg` as HIF JSON, and `hg_read_hif()`
reads a HIF file (or a HIF JSON string) back into a `net_hg`. HIF is the
interchange standard of the higher-order network libraries (XGI,
HyperNetX, HypergraphX, HAT, SimpleHypergraphs.jl), so these two verbs
move a hypergraph between hypernets and any of them.

## Usage

``` r
hg_write_hif(hg, file = NULL, pretty = FALSE)

hg_read_hif(file, sparse = FALSE)
```

## Arguments

- hg:

  A `net_hg` (from
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`network_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/network_hypergraph.md),
  `hg_read_hif()`, ...).

- file:

  For `hg_write_hif()`: the path to write, or `NULL` (default) to return
  the JSON text. For `hg_read_hif()`: the path of a HIF file, or a HIF
  JSON string.

- pretty:

  Logical. Indent the JSON? Default `FALSE`.

- sparse:

  Logical. Store the incidence of the read hypergraph as a sparse
  `dgCMatrix`? Default `FALSE`.

## Value

`hg_write_hif()`: with `file = NULL`, the HIF JSON as a single character
string; otherwise `file`, invisibly. `hg_read_hif()`: a `net_hg` with
the usual fields (`incidence` node x hyperedge with node and hyperedge
names as dimnames, `hyperedges`, `nodes`, `n_nodes`, `n_hyperedges`,
`size_distribution`, `params`) plus, when the file has them,
`window_counts` (hyperedge weights; an edge record without a weight
counts 1), `edge_data` (one row per hyperedge: `edge` and one column per
edge attribute), `node_data` (one row per node: `node`, `weight` when
the file has node weights, and one column per node attribute) and
`incidence_data` (one row per incidence with attributes: `node`, `edge`
and the attribute columns). Nodes and hyperedges are in file order: the
`nodes` / `edges` records first, then those first seen in `incidences`.
`params` holds `source = "hg_read_hif"`, `network_type` and the file's
`metadata`.

Both verbs raise `hypernets_missing_dependency` without jsonlite.
`hg_read_hif()` raises `hypernets_bad_input` for text that is not JSON,
a missing file, and a document that breaks the HIF schema (no
`incidences`, a record without its identifier, an unknown field, a
non-numeric weight, an unknown `network-type`), a duplicated node, edge
or (node, edge) pair, a zero incidence weight, a directed hypergraph, or
an attribute named like its record's identifier. `hg_write_hif()` raises
`hypernets_bad_input` for a non-`net_hg` input.

## Mapping

- incidences:

  One record per non-zero cell of the incidence matrix,
  `{"node", "edge", "weight"}`. `weight` is the incidence value; it is
  left out only for an integer 0/1 incidence (a binary hypergraph). On
  reading, an absent weight is 1, and a file whose incidences carry no
  weight gives an integer 0/1 incidence matrix, otherwise a double one,
  so the storage mode round-trips.

- nodes:

  One record per node, so isolated nodes survive. Node attributes are
  the columns of the hypergraph's node table (the `node_data` a previous
  `hg_read_hif()` stored) plus the planted `block` of
  [`hg_sample_sbm()`](https://mohsaqr.github.io/hypernets/reference/hg_sample_sbm.md).
  A column called `weight` is written as the node's HIF `weight`.

- edges:

  One record per hyperedge, named by the incidence column name. The
  hyperedge weight (the window counts of
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  which the Laplacian and modularity verbs read as hyperedge weights) is
  written as HIF `weight`; the hyperedge attribute table (`edge_data`:
  constant-within-group columns of
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  word counts of
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  a previous read) and the multiplicities of
  [`hg_snapshot()`](https://mohsaqr.github.io/hypernets/reference/hg_snapshot.md)
  are written as `attrs`.

- metadata:

  The scalar construction parameters and the writing package version; a
  hypergraph that came from `hg_read_hif()` writes its original metadata
  back unchanged.

## What does not round-trip

- Node and edge identifiers are strings: an integer identifier from
  another library (XGI's default `0, 1, 2, ...`) is read as `"0"`,
  `"1"`, `"2"`.

- JSON has one number type, so numeric attributes come back as double;
  factors, dates and times are written as their character form and come
  back as character.

- The construction record (`params`) is replaced by the source of the
  read; the constructor-specific fields that are not hypergraph data
  (the `$text` layer of
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  the clustered-sequence tables of
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md))
  are not written. `block` and `multiplicity` come back as node and edge
  attribute columns rather than as the `blocks` and `edge_multiplicity`
  fields.

- Directed hypergraphs (`"network-type": "directed"` or an incidence
  `direction`) have no `net_hg` counterpart and are refused.

- A zero incidence weight cannot be stored (membership in a `net_hg` is
  a non-zero cell) and is refused, as is a repeated (node, edge) pair,
  whose meaning HIF leaves open.

Weights are written with 17 significant digits, so doubles round-trip
exactly.

## References

Coll, M., Joslyn, C. A., Landry, N. W., Lotito, Q. F., Myers, A.,
Pickard, J., Praggastis, B., & Szufel, P. (2025). HIF: The hypergraph
interchange format for higher-order networks. *Network Science* 13, e21.
[doi:10.1017/nws.2025.10018](https://doi.org/10.1017/nws.2025.10018)

## Examples

``` r
if (requireNamespace("jsonlite", quietly = TRUE)) {
  memberships <- data.frame(
    actor = c("a", "b", "c", "a", "d"),
    group = c("g1", "g1", "g1", "g2", "g2"),
    weight = c(1, 2, 3, 0.5, 1),
    year = c(2001, 2001, 2001, 2002, 2002)
  )
  hg <- group_hypergraph(memberships, actor = "actor", group = "group",
                         weight = "weight")
  path <- tempfile(fileext = ".json")
  hg_write_hif(hg, file = path)
  back <- hg_read_hif(path)
  hg_get(back, what = "memberships")
  hg_write_hif(hg, pretty = TRUE)
}
#> [1] "{\n  \"network-type\": \"undirected\",\n  \"metadata\": {\n    \"source\": \"group_hypergraph\",\n    \"member\": \"actor\",\n    \"group\": \"group\",\n    \"weight\": \"weight\",\n    \"sparse\": false,\n    \"n_observations\": 5,\n    \"generator\": \"hypernets 0.6.0\"\n  },\n  \"incidences\": [\n    {\n      \"edge\": \"g1\",\n      \"node\": \"a\",\n      \"weight\": 1\n    },\n    {\n      \"edge\": \"g1\",\n      \"node\": \"b\",\n      \"weight\": 2\n    },\n    {\n      \"edge\": \"g1\",\n      \"node\": \"c\",\n      \"weight\": 3\n    },\n    {\n      \"edge\": \"g2\",\n      \"node\": \"a\",\n      \"weight\": 0.5\n    },\n    {\n      \"edge\": \"g2\",\n      \"node\": \"d\",\n      \"weight\": 1\n    }\n  ],\n  \"nodes\": [\n    {\n      \"node\": \"a\"\n    },\n    {\n      \"node\": \"b\"\n    },\n    {\n      \"node\": \"c\"\n    },\n    {\n      \"node\": \"d\"\n    }\n  ],\n  \"edges\": [\n    {\n      \"edge\": \"g1\",\n      \"attrs\": {\n        \"year\": 2001\n      }\n    },\n    {\n      \"edge\": \"g2\",\n      \"attrs\": {\n        \"year\": 2002\n      }\n    }\n  ]\n}"
```
