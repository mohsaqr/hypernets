# Memory-network communities with the map equation

Finds modules in a higher-order network by minimising the map equation
for memory networks (Rosvall et al. 2014; Edler, Bohlin & Rosvall 2017).
The state nodes of the HON (`"a -> b"` is physical node `b` reached from
`a`) are clustered, but each module's codebook is written over
**physical** nodes: state nodes of one physical node that land in the
same module share a codeword. This is not first-order Infomap run on the
state graph, and a physical node whose state nodes fall in different
modules belongs to all of them, so modules of physical nodes can
overlap.

The same search is run on the first-order network obtained by projecting
the memory network's link flow onto physical nodes, so the result says
whether memory compresses the flow better than a first-order
description.

## Usage

``` r
# S3 method for class 'net_hon'
hg_communities(
  x,
  partition = NULL,
  trials = 10L,
  teleportation = 0.15,
  seed = 1L,
  ...
)
```

## Arguments

- x:

  A `net_hon` from
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md).

- partition:

  `NULL` (default) to search, or a given state-node partition to
  evaluate without searching: a named vector (names = state labels as in
  `hg_get(hon, what = "nodes")`) or a data.frame with columns `node` and
  `community`.

- trials:

  Number of independent search trials (default 10). Ignored when
  `partition` is given, except for the first-order search.

- teleportation:

  Teleportation probability \\\tau \in \[0, 1)\\ (default 0.15,
  Infomap's default).

- seed:

  Integer seed for the random node orders (default 1). The caller's
  random-number stream is restored on exit.

- ...:

  Must be empty: an argument that only the hypergraph method takes
  raises `hypernets_bad_input`.

## Value

A `net_hon_communities` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
`what = "states"` (one row per state node), `"physical"` (one row per
physical node x module), `"modules"`, `"trials"`, `"first_order"` and
`"codelength"`; see
[`hg_get.net_hon_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_get.net_hon_communities.md).

## Details

**Flow.** State node visit rates are the stationary distribution of the
walk on the HON transition matrix \\P\_{\alpha\beta} =
w\_{\alpha\beta}/w\_\alpha\\ (Edler et al. 2017, Eq. 8-9). With
probability `teleportation` (\\\tau\\), and always from a state with no
out-links, the walker restarts at the tail of a link drawn in proportion
to its weight (teleportation to links; Lambiotte & Rosvall 2012).
Teleportation is unrecorded: link flow is \\f\_{\alpha\beta} =
\pi\_\alpha P\_{\alpha\beta}\\ renormalised to sum to one, and a state's
flow is its link in-flow. This is Infomap's directed flow model; node
flows and codelengths match the Python `infomap` package (2.15.1).
`teleportation = 0` uses the plain stationary distribution and needs an
irreducible chain.

**Codelength** (Edler et al. 2017, Eq. 24-26). With module enter flow
\\q^{in}\_m\\, exit flow \\q^{out}\_m\\, and physical flow inside a
module \\\pi\_{i \cap m} = \sum\_{\beta_i \in m} \pi\_{\beta_i}\\,
\$\$L(M) = q^{in} H(Q) + \sum_m p_m H(P_m),\$\$ where \\q^{in} = \sum_m
q^{in}\_m\\, \\H(Q)\\ is the entropy of the enter rates \\q^{in}\_m /
q^{in}\\, \\p_m = q^{out}\_m + \sum_i \pi\_{i \cap m}\\, and \\H(P_m)\\
is the entropy of \\\\q^{out}\_m, \pi\_{i \cap m}\\\\ normalised by
\\p_m\\. The one-module codelength is the entropy of the physical-node
flow.

**Search** (Edler et al. 2017, Algorithms 2-5): every state node starts
in its own module; in random sequential order each node moves to the
neighbouring or new module that most decreases \\L\\; modules are then
aggregated into nodes and the moves repeat, until \\L\\ stops decreasing
(repeated node aggregation, the Louvain machinery). The result is
fine-tuned by restarting the aggregation from single state nodes placed
in their current modules. Coarse-tuning (Algorithm 6) is not
implemented. `trials` independent runs are made from one seeded random
stream and the shortest codelength is kept; `hg_get(x, what = "trials")`
reports every run and its agreement (adjusted Rand index) with the kept
one.

**First-order comparison.** The link flow of the memory network is
summed over state nodes into a physical-node network, \\w\_{ij} =
\sum\_{\alpha_i, \beta_j} f\_{\alpha\beta}\\, and searched with the same
flow model and the same number of trials. Rosvall et al. (2014) compare
memory and first-order maps of the same pathway data this way.

## Conditions

`hypernets_bad_input` for a non-`net_hon` input, invalid arguments, or a
partition that does not cover every state; `hypernets_not_ergodic` when
`teleportation = 0` and the walk has no unique stationary flow.

## References

Rosvall, M., Esquivel, A. V., Lancichinetti, A., West, J. D., &
Lambiotte, R. (2014). Memory in network flows and its effects on
spreading dynamics and community detection. *Nature Communications*, 5,
4630. [doi:10.1038/ncomms5630](https://doi.org/10.1038/ncomms5630)

Edler, D., Bohlin, L., & Rosvall, M. (2017). Mapping higher-order
network flows in memory and multilayer networks with Infomap.
*Algorithms*, 10(4), 112.
[doi:10.3390/a10040112](https://doi.org/10.3390/a10040112)

Rosvall, M., & Bergstrom, C. T. (2008). Maps of random walks on complex
networks reveal community structure. *Proceedings of the National
Academy of Sciences*, 105(4), 1118-1123.
[doi:10.1073/pnas.0706851105](https://doi.org/10.1073/pnas.0706851105)

Lambiotte, R., & Rosvall, M. (2012). Ranking and clustering of nodes in
networks with smart teleportation. *Physical Review E*, 85(5), 056107.
[doi:10.1103/PhysRevE.85.056107](https://doi.org/10.1103/PhysRevE.85.056107)

Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
Classification*, 2(1), 193-218.
[doi:10.1007/BF01908075](https://doi.org/10.1007/BF01908075)

## Examples

``` r
seqs <- list(c("a", "h", "b", "a", "h", "b", "a"),
             c("c", "h", "d", "c", "h", "d", "c"),
             c("a", "h", "b", "a", "h", "b"),
             c("c", "h", "d", "c", "h", "d"))
hon <- hon(seqs, max_order = 2L)
comm <- hg_communities(hon, trials = 3L)
comm
#> Memory-network communities (map equation)
#>   Memory nodes: 7 | states: 5 | communities: 3
#>   (1 community holds only nodes with zero flow, reached by teleportation alone)
#>   States in more than one community: 1
#>   Codelength memory:      1.6061 bits (one module 2.2589)
#>   Codelength first-order: 2.2525 bits (one module 2.2525; communities: 1)
#>   Runs: 3 (teleportation 0.15, seed 1)
#>    node state community      flow
#>       a     a         2 0.1664584
#>  a -> h     h         2 0.1629182
#>       b     b         2 0.1706234
#>       c     c         1 0.1664584
#>  c -> h     h         1 0.1629182
#>       d     d         1 0.1706234
#>       h     h         3 0.0000000
summary(comm)
#> Memory-network communities (map equation)
#> Tables: states, physical, modules, trials, first_order, codelength
#> 
#> states (7 rows)
#>    node state community      flow
#>       a     a         2 0.1664584
#>  a -> h     h         2 0.1629182
#>       b     b         2 0.1706234
#>       c     c         1 0.1664584
#>  c -> h     h         1 0.1629182
#> 
#> physical (6 rows)
#>  state community      flow share n_communities
#>      a         2 0.1664584   1.0             1
#>      b         2 0.1706234   1.0             1
#>      c         1 0.1664584   1.0             1
#>      d         1 0.1706234   1.0             1
#>      h         1 0.1629182   0.5             2
#> 
#> modules (3 rows)
#>  community n_nodes n_states flow  exit_flow enter_flow
#>          1       3        3  0.5 0.00000000 0.01071429
#>          2       3        3  0.5 0.00000000 0.01071429
#>          3       1        1  0.0 0.02142857 0.00000000
#> 
#> trials (3 rows)
#>  run codelength n_communities ari_to_best  best
#>    1   1.606134             3           1  TRUE
#>    2   1.606134             3           1 FALSE
#>    3   1.606134             3           1 FALSE
#> 
#> first_order (5 rows)
#>  state community      flow
#>      a         1 0.1668403
#>      b         1 0.1669077
#>      c         1 0.1668403
#>      d         1 0.1669077
#>      h         1 0.3325040
#> 
#> codelength (2 rows)
#>        model codelength index_codelength module_codelength one_level_codelength
#>       memory   1.606134       0.02142857          1.584705             2.258869
#>  first_order   2.252456       0.00000000          2.252456             2.252456
#>  savings_bits savings_pct n_communities
#>      0.652735    28.89654             3
#>      0.000000     0.00000             1
hg_get(comm)
#>     node state community      flow
#> 1      a     a         2 0.1664584
#> 2 a -> h     h         2 0.1629182
#> 3      b     b         2 0.1706234
#> 4      c     c         1 0.1664584
#> 5 c -> h     h         1 0.1629182
#> 6      d     d         1 0.1706234
#> 7      h     h         3 0.0000000
hg_get(comm, what = "physical")
#>   state community      flow share n_communities
#> 1     a         2 0.1664584   1.0             1
#> 2     b         2 0.1706234   1.0             1
#> 3     c         1 0.1664584   1.0             1
#> 4     d         1 0.1706234   1.0             1
#> 5     h         1 0.1629182   0.5             2
#> 6     h         2 0.1629182   0.5             2
```
