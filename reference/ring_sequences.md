# Simulated sequences with planted memory modules

A process whose communities can only be seen with memory. Four groups of
four actions sit on a ring and neighbouring groups share one action
(`s1`, ..., `s4`); the other two actions of group \\k\\ (`p<k>_1`,
`p<k>_2`) are its own. A walk starts in a random group at a random
action of it and each step moves to another action of the current group,
chosen at random; on reaching a shared action the walk switches to the
neighbouring group with probability 0.05. At a shared action a
first-order network sees both groups equally, while a second-order state
such as `"p1_1 -> s2"` knows it came from group 1. The planted group of
every state is in
[ring_communities](https://mohsaqr.github.io/hypernets/reference/ring_communities.md).

## Usage

``` r
ring_sequences
```

## Format

A list of 200 character vectors, one walk each, 40 actions long, over
the 12 actions `s1`-`s4` and `p1_1`-`p4_2`.

## Source

Simulated with seed 1 by `data-raw/ring_sequences.R`.

## See also

[`hg_communities()`](https://mohsaqr.github.io/hypernets/reference/hg_communities.md),
which recovers the four groups from these sequences.

## Examples

``` r
length(ring_sequences)
#> [1] 200
ring_hon <- hon(ring_sequences, max_order = 2L)
ring_hon
#> Higher-order network: 12 states, 36 nodes, 189 rules (highest order 2, min_freq 1, 200 sequences)
#>                path       from   to count probability from_order to_order
#>          p1_1 -> s1       p1_1   s1   166 0.317399618          1        2
#>          p1_1 -> s2       p1_1   s2   176 0.336520076          1        2
#>        p1_1 -> p1_2       p1_1 p1_2   181 0.346080306          1        1
#>  p1_1 -> s1 -> p1_1 p1_1 -> s1 p1_1    51 0.310975610          2        1
#>  p1_1 -> s1 -> p1_2 p1_1 -> s1 p1_2    56 0.341463415          2        1
#>  p1_1 -> s1 -> p4_1 p1_1 -> s1 p4_1     2 0.012195122          2        1
#>  p1_1 -> s1 -> p4_2 p1_1 -> s1 p4_2     1 0.006097561          2        1
#>    p1_1 -> s1 -> s2 p1_1 -> s1   s2    49 0.298780488          2        2
#>    p1_1 -> s1 -> s4 p1_1 -> s1   s4     5 0.030487805          2        2
#>  p1_1 -> s2 -> p1_1 p1_1 -> s2 p1_1    59 0.341040462          2        1
#> ... 179 more rows
```
