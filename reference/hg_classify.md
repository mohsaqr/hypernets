# Transductive label spreading on a hypergraph, as a tidy table

Calls the in-package
[`hypergraph_transduction()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_transduction.md)
engine (Zhou et al. 2006): labels known for a few nodes spread over the
hypergraph structure to classify every node.

## Usage

``` r
hg_classify(
  hg,
  labels,
  xi = 0.99,
  type = c("zhou", "random_walk"),
  normalization = c("none", "class_mass")
)

hypergraph_classify(
  hg,
  labels,
  xi = 0.99,
  type = c("zhou", "random_walk"),
  normalization = c("none", "class_mass")
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any honets `net_hypergraph`).

- labels:

  The known labels: a named character vector (names are node identifiers
  – documents under `nodes = "doc"` – values their class labels), or a
  tidy data.frame with a `node` column and a `label`, `cluster` or
  `predicted` column.

- xi, type:

  Passed to
  [`hypergraph_transduction()`](https://mohsaqr.github.io/hypernets/reference/hypergraph_transduction.md).

- normalization:

  Decision rule for turning spread scores into predictions: `"none"`
  (default, the raw Zhou 2006 argmax) or `"class_mass"` (class-mass
  normalization, Zhu et al. 2003). Use `"class_mass"` when the labeled
  seeds are class-imbalanced – the raw rule can collapse every
  prediction onto the majority class.

## Value

A base `data.frame`, one row per node, with columns `node`, `label` (the
given label or `NA`), `predicted`, `score`, and `margin`.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
hg_classify(hg, labels = c(cooking_1 = "cooking", space_1 = "space"))
#>        node   label predicted     score      margin
#> 1 cooking_1 cooking   cooking 0.2665654 0.082276891
#> 2 cooking_2    <NA>   cooking 0.2538871 0.009941802
#> 3   space_1   space     space 0.3016056 0.117317071
#> 4   space_2    <NA>     space 0.2663331 0.072737125
```
