# Per-actor state sequences from a clustered corpus

Turns document-level topic assignments into the canonical long sequence
table the memory family consumes. Each document of a
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
carries whatever metadata came with the corpus (an author, a turn, a
timestamp);
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
gives each document a topic. `hg_sequences()` joins the two, orders the
documents within each actor and returns one row per document as `actor`
/ `time` / `action` – the long shape
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md),
[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md),
[`compare_hon()`](https://mohsaqr.github.io/hypernets/reference/compare_hon.md)
and
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
read through their own `action` / `actor` / `time` arguments, so no
coercion is written at the call site:

## Usage

``` r
hg_sequences(hg, clusters = NULL, actor, order_by, state = NULL)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  – more generally, any `net_hypergraph` carrying a documents table, the
  one reached with `as.data.frame(hg, what = "documents")`.

- clusters:

  A partition of the documents, as returned by
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md):
  a data.frame with one row per node, a `node` column holding document
  identifiers and a state column (`cluster` by default, or the column
  named by `state`). `NULL` (the default) takes the state from the
  documents table instead, in which case `state` is required.

- actor:

  Name of the documents-table column identifying who produced each
  document. One sequence is returned per distinct value.

- order_by:

  Name of the documents-table column ordering the documents within an
  actor – a turn number, an index, a date or a timestamp. Anything
  [`order()`](https://rdrr.io/r/base/order.html) can sort is accepted,
  and the values are carried through to the `time` column unchanged.

- state:

  Name of the column holding each document's state. With `clusters`
  supplied it names a column of `clusters` and defaults to `"cluster"`;
  with `clusters = NULL` it names a column of the documents table and
  has no default.

## Value

A base `data.frame` with one row per document of `hg`, ordered by
`actor` and then `time`, with columns

- actor:

  character. The value of the `actor` column; one sequence per distinct
  value.

- time:

  The value of the `order_by` column, type unchanged, so a numeric turn
  stays numeric and a `Date` stays a `Date`.

- action:

  character. The document's state – its cluster label under `clusters`,
  or the `state` column of the documents table.

## Details

    seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
    build_hon(seqs, action = "action", actor = "actor", time = "time")

Those three column names are fixed, which is what makes the second call
the same every time. The memory-family verbs that do not (yet) take the
long form –
[`build_mogen()`](https://mohsaqr.github.io/hypernets/reference/build_mogen.md),
[`build_hypa()`](https://mohsaqr.github.io/hypernets/reference/build_hypa.md),
[`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md)
– read a list of trajectories instead; pass the result of
`hg_sequences()` to them only through one of the verbs above, since a
long table handed to them bare is read as a *wide* one (one trajectory
per row) and silently gives the wrong model.

This closes the one missing edge in the cross-family design:
[`pathways()`](https://mohsaqr.github.io/hypernets/reference/pathways.md)
bridges the memory family to the others and
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
bridges sequences to hypergraphs, but until now the text family
dead-ended at clustering and the caller had to assemble the sequence
table themselves.

## Conditions

A `hypernets_bad_input` error is raised when `hg` is not a hypergraph,
when it carries no documents table, when `actor`, `order_by` or `state`
does not name a column of the table it is looked up in, when `clusters`
lacks a `node` column or assigns a document more than one state, when
some document of `hg` has no row in `clusters`, or when the actor, time
or state values contain `NA` (which would silently shorten a sequence).

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances*, 2(5),
e1600028.

Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
clustering, classification, and embedding. *Advances in Neural
Information Processing Systems*, 19, 1601-1608.

## See also

[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
for the partition,
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
and
[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md)
for what to do with the sequences,
[`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
to read them back into a hypergraph.

## Examples

``` r
posts <- data.frame(
  student = rep(c("ana", "bo"), each = 4),
  turn = rep(1:4, times = 2),
  phase = rep(c("early", "early", "late", "late"), times = 2),
  text = c(
    "simmer the soup with onions and carrots",
    "this soup recipe needs salt on a cold night",
    "the telescope revealed a distant galaxy and stars",
    "astronomers aimed the telescope at the stars all night",
    "a galaxy of stars seen through the telescope",
    "salt and onions make the soup",
    "the soup simmered all night with carrots",
    "astronomers watched the distant galaxy"
  ),
  stringsAsFactors = FALSE
)
hg <- text_hypergraph(posts, column = "text",
                      stop_words = c("the", "with", "and", "a", "this",
                                     "at", "on", "all", "of", "through"))
topics <- hg_cluster(hg, k = 2, seed = 1)
hg_sequences(hg, topics, actor = "student", order_by = "turn")
#>   actor time    action
#> 1   ana    1 Cluster 1
#> 2   ana    2 Cluster 1
#> 3   ana    3 Cluster 2
#> 4   ana    4 Cluster 2
#> 5    bo    1 Cluster 2
#> 6    bo    2 Cluster 1
#> 7    bo    3 Cluster 1
#> 8    bo    4 Cluster 2

# A state that is already a column of the corpus needs no clustering.
hg_sequences(hg, actor = "student", order_by = "turn", state = "phase")
#>   actor time action
#> 1   ana    1  early
#> 2   ana    2  early
#> 3   ana    3   late
#> 4   ana    4   late
#> 5    bo    1  early
#> 6    bo    2  early
#> 7    bo    3   late
#> 8    bo    4   late

# Straight into the memory family, no further coercion.
seqs <- hg_sequences(hg, topics, actor = "student", order_by = "turn")
build_hon(seqs, action = "action", actor = "actor", time = "time",
          max_order = 2L)
#> Higher-Order Network (HON)
#>   Nodes:        2 (2 first-order states)
#>   Edges:        4
#>   Max order:    1 (requested 2)
#>   Min freq:     1
#>   Trajectories: 2
```
