# Tidy tables of a memory-network community result

Tidy tables of a memory-network community result

## Usage

``` r
# S3 method for class 'net_hon_communities'
hg_get(
  x,
  what = c("states", "physical", "modules", "trials", "first_order", "codelength"),
  ...,
  community = NULL,
  overlapping = FALSE
)
```

## Arguments

- x:

  A `net_hon_communities` object.

- what:

  Which table:

  `"states"` (default)

  :   one row per memory node: `node` (the memory node, `"a -> b"`),
      `state` (the physical state it belongs to, `b`), `community`,
      `flow`.

  `"physical"`

  :   one row per state x community it appears in: `state`, `community`,
      `flow` (\\\pi\_{i \cap m}\\), `share` (of that state's flow),
      `n_communities` (\> 1 marks an overlapping state).

  `"modules"`

  :   one row per community: `community`, `n_nodes` (memory nodes),
      `n_states`, `flow`, `exit_flow`, `enter_flow`.

  `"trials"`

  :   one row per search run: `run`, `codelength`, `n_communities`,
      `ari_to_best` (adjusted Rand index with the kept partition),
      `best`.

  `"first_order"`

  :   one row per state: `state`, its `community` in the best
      first-order partition, and its first-order `flow`.

  `"codelength"`

  :   one row per model (`memory`, `first_order`): `codelength`,
      `index_codelength`, `module_codelength`, `one_level_codelength`,
      `savings_bits`, `savings_pct`, `n_communities`.

- ...:

  Ignored.

- community:

  Integer vector or `NULL`: keep only these communities (tables
  `"states"`, `"physical"`, `"modules"`).

- overlapping:

  Logical: for `what = "physical"`, keep only states that appear in more
  than one community. Default `FALSE`.

## Value

A base data.frame as described under `what`.

## Examples

``` r
seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
             c("c", "h", "d", "c", "h", "d", "c"))
comm <- hg_communities(hon(seqs, max_order = 2L), trials = 2L)
hg_get(comm, what = "modules")
#>   community n_nodes n_states flow exit_flow enter_flow
#> 1         1       5        5    1         0          0
hg_get(comm, what = "physical", overlapping = TRUE)
#> [1] state         community     flow          share         n_communities
#> <0 rows> (or 0-length row.names)
```
