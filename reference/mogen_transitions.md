# Extract Transition Table from a MOGen Model

Returns a data frame of all transitions at a given Markov order, sorted
by count (descending). Each row shows the full path as a readable
sequence of states, along with the observed count and transition
probability.

## Usage

``` r
mogen_transitions(x, order = NULL, min_count = 1L, top = NULL)
```

## Arguments

- x:

  A `net_mogen` object from
  [`build_mogen()`](https://mohsaqr.github.io/hypernets/reference/build_mogen.md).

- order:

  Integer. Which order's transitions to extract. Must be a whole number;
  a non-integer value is an error rather than being silently truncated.
  Defaults to the optimal order selected by the model.

- min_count:

  Integer. Minimum observed count to include (default 1). Use this to
  filter out rare transitions that have unreliable probabilities.

- top:

  Integer or `NULL`. Return only the first `top` rows, applied after any
  filter and after `sort_by`, so `sort_by` and `top` compose. Default
  `NULL` returns every row.

## Value

A data frame with columns:

- path:

  The full state sequence (e.g., "AI -\> FAIL -\> SOLVE").

- count:

  Number of times this transition was observed.

- probability:

  Transition probability P(to \| from).

- from:

  The context / conditioning states (k-gram source node).

- to:

  The predicted next state.

## Details

At order k, each edge in the De Bruijn graph represents a (k+1)-step
path. For example, at order 2, the edge from node "AI -\> FAIL" to node
"FAIL -\> SOLVE" represents the three-step path AI -\> FAIL -\> SOLVE.
The `path` column reconstructs this full sequence for readability.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
mg <- build_mogen(seqs, max_order = 2)
mogen_transitions(mg, order = 1)
#>     path count probability from to
#> 1 B -> C     3      1.0000    B  C
#> 2 A -> B     2      1.0000    A  B
#> 3 C -> D     2      0.6667    C  D
#> 4 C -> A     1      0.3333    C  A
#> 5 D -> A     1      1.0000    D  A

# \donttest{
trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
              c("B","C","D","A"), c("C","D","A","B"))
m <- build_mogen(trajs, max_order = 3)
mogen_transitions(m, order = 1)
#>     path count probability from to
#> 1 A -> B     3      1.0000    A  B
#> 2 C -> D     3      1.0000    C  D
#> 3 D -> A     2      0.6667    D  A
#> 4 B -> C     2      0.6667    B  C
#> 5 D -> C     1      0.3333    D  C
#> 6 B -> D     1      0.3333    B  D
# }
```
