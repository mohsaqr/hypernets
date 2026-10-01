# Read a table from any hypernets object

`hg_get()` is the single reader of the package. Every object a
constructor or a verb returns – a hypergraph, a memory network, a
simplicial complex, an inference or a community result – hands over its
tables through `hg_get()`, never through `$`. `what =` selects the
table; `what = NULL` (the default) returns the object's primary table.
Filters, orderings and `top =` are named arguments of the method for
that class.

## Usage

``` r
hg_get(x, what = NULL, ...)

# Default S3 method
hg_get(x, what = NULL, ...)
```

## Arguments

- x:

  A hypernets object.

- what:

  `NULL` (default) for the object's primary table, or the name of one of
  its tables (see *Tables by class*).

- ...:

  Arguments of the method for `class(x)`: filters (`order_min`,
  `min_count`, `dim`, `k`, `dimension`, `significant`, ...), `sort_by`,
  and `top` (the first `top` rows, applied last).

## Value

A base `data.frame`, one row per observation of the selected table. An
object of a class without a method raises `hypernets_bad_input`, naming
the class.

## Tables by class

- `net_hg` (hypergraphs):

  `"edges"` (default), `"nodes"`, `"memberships"`, `"sets"`,
  `"state_counts"`, `"node_data"`, `"edge_data"`, `"incidence_data"`;
  see
  [`hg_get.net_hg()`](https://mohsaqr.github.io/hypernets/reference/hg_get.net_hg.md).

- `text_hypergraph`:

  `"weights"` (default), `"documents"`, `"vocabulary"`, `"sentences"`;
  see
  [`hg_get.text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/hg_get.text_hypergraph.md).

- `net_temporal_hypergraph`:

  `"memberships"` (default), `"edges"`, `"nodes"`.

- `net_hon`
  ([`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md)):

  `"rules"` (default), `"nodes"`, `"pathways"`.

- `net_honem`
  ([`honem()`](https://mohsaqr.github.io/hypernets/reference/honem.md)):

  `"embeddings"` (default), `"variance"`.

- `net_hypa`
  ([`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md) on
  sequences):

  `"scores"` (default), `"over"`, `"under"`, `"pathways"`.

- `net_mogen`
  ([`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md)):

  `"orders"` (default), `"transitions"`, `"paths"`, `"pathways"`.

- `net_markov_order`
  ([`markov_order()`](https://mohsaqr.github.io/hypernets/reference/markov_order.md)):

  `"orders"` (default), `"null"`.

- `net_path_dependence`
  ([`memory()`](https://mohsaqr.github.io/hypernets/reference/memory.md)):

  `"contexts"`.

- `net_markov_stability`
  ([`hg_markov_stability()`](https://mohsaqr.github.io/hypernets/reference/hg_markov_stability.md)):

  `"states"` (default), `"passage_time"`, `"stationary"`.

- `net_hon_boot`
  ([`hg_bootstrap()`](https://mohsaqr.github.io/hypernets/reference/hg_bootstrap.md)):

  `"edges"`.

- `net_hon_compare`
  ([`hg_compare()`](https://mohsaqr.github.io/hypernets/reference/hg_compare.md)):

  `"edges"`.

- `net_hon_communities`
  ([`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  on a memory network):

  `"states"` (default), `"physical"`, `"modules"`, `"codelength"`,
  `"trials"`.

- `simplicial_complex`
  ([`simplicial()`](https://mohsaqr.github.io/hypernets/reference/simplicial.md)):

  `"simplices"` (default), `"f_vector"`, `"degree"`.

- `q_analysis`
  ([`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md)):

  `"q_levels"` (default), `"nodes"`.

- `persistent_homology`
  ([`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md)):

  `"persistence"` (default), `"betti"`.

- `persistence_landscape`
  ([`hg_landscape()`](https://mohsaqr.github.io/hypernets/reference/hg_landscape.md)):

  `"landscape"`.

- `hg_communities`
  ([`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md)
  on a hypergraph):

  see
  [`hg_get.hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_get.hg_communities.md).

- `hypernets_community_comparison`, `hypernets_motifs`, `net_hg_mmsbm`:

  see their methods.

## Examples

``` r
hg <- window_hypergraph(list(s1 = c("a", "b", "a", "c")), window = 2L)
hg_get(hg)
#>   hyperedge size members weight
#> 1        h1    2    a, b      2
#> 2        h2    2    a, c      1
hg_get(hg, what = "nodes")
#>   node degree
#> 1    a      2
#> 2    b      1
#> 3    c      1

seqs <- list(c("a", "b", "c", "a", "b", "c"), c("x", "b", "d", "x", "b", "d"))
hg_get(hon(seqs, max_order = 2), what = "nodes")
#>   id node
#> 1  1    a
#> 2  2    b
#> 3  3    c
#> 4  4    d
#> 5  5    x
```
