# Hypergraph from co-occurrence data or an edge list

Constructs a
[net_hypergraph](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
the way a network is defined from data. **Co-presence data** name an
`actor` and a `group`: every actor sharing one value of `group` (a
session, a team, a citation block) belongs to one hyperedge. An **edge
list** names `from` and `to`, and every row is a hyperedge of size two.
An optional `weight` column produces a weighted incidence matrix.

## Usage

``` r
group_hypergraph(
  data,
  actor = NULL,
  group = NULL,
  weight = NULL,
  nodes = NULL,
  sparse = FALSE,
  from = NULL,
  to = NULL,
  member = NULL,
  cooccur_by = NULL
)
```

## Arguments

- data:

  Data frame in long format, one row per actor-in-group or per edge.

- actor:

  Character. Name of the column whose values become the hypergraph's
  nodes (members, participants, actors).

- group:

  Character. Name of the column whose shared values bind actors into one
  hyperedge (groups, sessions, teams) – Dynet's co-presence vocabulary.
  When neither pair of columns is named, `from` and `to` are detected
  case-insensitively from the alias table
  [`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md)
  uses (`source`/`target`, `sender`/`receiver`, ...), and failing that
  `actor` and `group`.

- weight:

  Character or `NULL`. If supplied, the column is summed per
  `(actor, hyperedge)` pair to produce a weighted incidence matrix.
  Default `NULL` produces a 0/1 binary incidence matrix.

- nodes:

  Optional vector giving the complete node universe. This keeps nodes
  with no observed group memberships as zero-incidence rows, which is
  needed for representations such as citation hypergraphs where every
  decision is a node but some decisions are never cited.

- sparse:

  Logical. Store incidence as a sparse `Matrix`? Use this for large,
  sparse event data such as the full GFCC citation-block corpus.

- from, to:

  Column names of a pairwise edge list, as an alternative to `actor` and
  `cooccur_by`.

- member, cooccur_by:

  Deprecated names of `actor` and `group`; using them warns with a
  `honets_deprecated` condition.

## Value

A `net_hypergraph` object with the same structure produced by
[`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
(`hyperedges`, `incidence`, `nodes`, `n_nodes`, `n_hyperedges`,
`size_distribution`, `params`), plus `edge_data` when `data` carries
hyperedge attributes (see Details). The `params` list records
`source = "group_hypergraph"` and the original column names.

## Details

The bipartite representation preserves the full group structure without
projecting to a pairwise network. A group of three members A, B, C
produces a single 3-hyperedge containing all three, not three pairwise
edges AB, AC, BC. This avoids information loss when group interactions
are the primary unit of analysis (Perc et al. 2013).

Unlike
[`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
(which derives hyperedges from a network's clique structure),
`group_hypergraph()` takes group memberships directly. The two functions
are complementary:

- `group_hypergraph()` - when group membership is observed (sessions,
  transactions, co-authorships).

- [`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md) -
  when only pairwise interactions are observed and triadic structure
  must be inferred from triangles.

Rows with `NA` in the actor, hyperedge or weight column are dropped
silently.

Every other column of `data` that is constant within a hyperedge (a
session's date, a team's department) is kept as a hyperedge attribute in
the result's `edge_data` table, one row per hyperedge, where
[`hg_subset()`](https://mohsaqr.github.io/hypernets/reference/hg_subset.md)'s
`where` argument and
[`hg_project()`](https://mohsaqr.github.io/hypernets/reference/hg_project.md)'s
`edge_source` can use it. Columns that vary within a hyperedge describe
memberships, not hyperedges, and are left out.

## Note

Dense and sparse paths are tested for exact equality. The sparse path is
additionally exercised by the full GFCC reproduction from the Legal
Hypergraphs Zenodo archive (3,618 nodes, 46,165 hyperedges and 77,187
nonzero incidences).

## References

Perc, M., Gomez-Gardenes, J., Szolnoki, A., Floria, L. M., & Moreno, Y.
(2013). Evolutionary dynamics of group interactions on structured
populations: a review. *Journal of the Royal Society Interface* 10(80),
20120997.
[doi:10.1098/rsif.2012.0997](https://doi.org/10.1098/rsif.2012.0997)

## See also

[`build_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/build_hypergraph.md)
for the clique-based constructor,
[`temporal_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/temporal_hypergraph.md)
for the same inputs with a clock.

## Examples

``` r
df <- data.frame(
  person = c("Alice", "Bob", "Carol", "Alice", "Bob",
             "Dave", "Carol", "Dave", "Eve"),
  session = c("S1", "S1", "S1", "S2", "S2",
              "S3", "S3", "S3", "S3")
)
hg <- group_hypergraph(df, actor = "person", group = "session")
print(hg)
#> Hypergraph: 5 nodes, 3 hyperedges
#> Size distribution:
#>   size_2   : 1
#>   size_3   : 2
#> Source: group membership (member = person, group = session)
summary(hg)
#> Hypergraph summary
#>   Nodes:         5
#>   Hyperedges:    3
#>   Mean size:     2.67
#>   Max size:      3

contacts <- data.frame(from = c("a", "b"), to = c("b", "c"))
group_hypergraph(contacts, from = "from", to = "to")
#> Hypergraph: 3 nodes, 2 hyperedges
#> Size distribution:
#>   size_2   : 2
#> Source: group membership (member = actor, group = edge)
```
