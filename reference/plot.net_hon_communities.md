# Plot memory-network communities

- `type = "physical"` (default):

  every community is drawn as a pebble around the physical nodes it
  holds, through
  [`plot.net_hg()`](https://mohsaqr.github.io/hypernets/reference/plot.net_hg.md)
  on the community hypergraph (member = physical node, group =
  `"Community k"`, from `hg_get(x, what = "physical")`). A physical node
  shared by several communities lies inside each of their pebbles.
  Pebbles are coloured by community (Okabe-Ito, legend "Community").
  Each node is a black circle whose area follows its physical flow (the
  memory network's visit rate summed over its state nodes), and a
  triangle inside the circle points at the physical node its link flow
  most often goes to next. Each community has a title box with its name
  and flow (the visit rate of its state nodes) beside its pebble, drawn
  by
  [`plot.net_hg()`](https://mohsaqr.github.io/hypernets/reference/plot.net_hg.md)
  (haloed labels, legend below, wide margins); the boxes sit a little
  further out than there (`title_gap = 0.14`) to clear the labels of the
  nodes at a pebble's rim. A community holding a single physical node
  has no pebble of its own; its title box sits by the node and the
  caption names it.

- `type = "states"`:

  the state nodes of the higher-order network, drawn by
  [`cograph::overlay_communities()`](https://sonsoles.me/cograph/reference/overlay_communities.html),
  one ring per community
  ([`cograph::layout_groups()`](https://sonsoles.me/cograph/reference/layout_groups.html)),
  with community blobs, cograph's group legend, short state labels, and
  faded edges (transition probability below 0.05 hidden).

State nodes with zero flow (reached only by teleportation, see
[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md))
carry no physical flow and are never drawn in the physical view; the
state view leaves them and the singleton communities they form out
unless `show_zero_flow = TRUE`. A caption says how many were left out.
The tables returned by
[`hg_get.net_hon_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_get.net_hon_communities.md)
are unaffected.

## Usage

``` r
# S3 method for class 'net_hon_communities'
plot(x, type = c("physical", "states"), show_zero_flow = FALSE, ...)
```

## Arguments

- x:

  A `net_hon_communities` object.

- type:

  `"physical"` (default) or `"states"`.

- show_zero_flow:

  Draw zero-flow state nodes in the state view? Default `FALSE`.

- ...:

  For `type = "physical"`, passed to
  [`plot.net_hg()`](https://mohsaqr.github.io/hypernets/reference/plot.net_hg.md)
  (e.g. `seed`, `label_size`, `arrow_style = "outside"`, `layout`),
  overriding the defaults set here; for `type = "states"`, passed to
  [`cograph::overlay_communities()`](https://sonsoles.me/cograph/reference/overlay_communities.html)
  and on to
  [`cograph::splot()`](https://sonsoles.me/cograph/reference/splot.html).

## Value

For `type = "physical"`, a ggplot object (print it to draw). For
`type = "states"`, `x`, invisibly (cograph draws with base graphics).

## Conditions

`hypernets_bad_input` from
[`plot.net_hg()`](https://mohsaqr.github.io/hypernets/reference/plot.net_hg.md)
for arguments passed through `...` that it rejects.

## Examples

``` r
seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
             c("c", "h", "d", "c", "h", "d", "c"))
comm <- hg_communities(hon(seqs, max_order = 2L), trials = 2L)
plot(comm)

plot(comm, type = "states")
```
