# Plot a hypergraph with hyperedges as blobs

Places the nodes with a graph layout of the clique projection and draws
every hyperedge as a smooth translucent blob around its members, through
[`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html).
Blobs can be coloured by hyperedge size, by a column of the edge
metadata, or by any value supplied per hyperedge, and outlined with a
line type by the same kinds of selector; honets adds the legend for that
mapping. This is the drawing convention used for the decision and
tribunal figures of Coupette et al. (2024). A hyperedge with a single
member is drawn as its node alone.

## Usage

``` r
# S3 method for class 'net_hypergraph'
plot(
  x,
  layout = c("spring", "circle"),
  seed = 1L,
  color_by = NULL,
  linetype_by = NULL,
  labels = TRUE,
  label_size = 3,
  node_size = 2.5,
  alpha = 0.45,
  padding = 0.06,
  legend_title = NULL,
  ...
)
```

## Arguments

- x:

  A `net_hypergraph` with at least one hyperedge of two or more members.
  Draw a part of a large hypergraph by passing
  [`hg_subset()`](https://mohsaqr.github.io/hypernets/reference/hg_subset.md)
  first.

- layout:

  `"spring"` (default;
  [`cograph::layout_spring()`](https://sonsoles.me/cograph/reference/layout_spring.html)
  on the clique projection), `"circle"`, or a data.frame with `node`,
  `x` and `y` columns to reuse coordinates across panels.

- seed:

  Seed for the spring layout (default `1`).

- color_by:

  Fill of the blobs: `NULL` (one colour), `"size"` (hyperedge
  cardinality, sequential scale), the name of a column in the edge
  metadata (`x$edge_data`), a vector named by hyperedge, or a vector
  with one value per hyperedge. Character or factor values get the
  Okabe-Ito palette; numeric values a sequential scale built from it.

- linetype_by:

  Outline line type of the blobs, resolved like `color_by`; discrete
  only.

- labels:

  `TRUE` (default) writes the node names, `FALSE` writes none, and a
  character vector named by node replaces the names shown.

- label_size, node_size:

  Text and point sizes.

- alpha:

  Fill transparency of the blobs (default `0.45`).

- padding:

  Blob radius around the member positions as a fraction of the layout
  extent (default `0.06`, the radius cograph draws when the layout fills
  its canvas; smaller values enlarge the canvas).

- legend_title:

  Legend title for `color_by`; defaults to the selector's name.

- ...:

  Unused; for S3 consistency.

## Value

A ggplot object.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
dat <- data.frame(
  member = c("a", "b", "c", "b", "c", "d", "d", "e", "f"),
  event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e4")
)
hg <- group_hypergraph(dat, "member", "event")
plot(hg, color_by = "size")
```
