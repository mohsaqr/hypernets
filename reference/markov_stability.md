# Markov Stability of a First- or Higher-Order Chain

Per-state stability of the random walk carried by a transition matrix:
how persistent a state is, how much of the walk's time it holds, how
long the process dwells in it before leaving, how long until it comes
back, and how many steps separate it from every other state.

Given a `net_hon` the states are **higher-order** states — `"c -> b"` is
the state `b` reached from `c` — so the same quantities describe the
higher-order random walk rather than the memoryless one.
[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
returns a row-stochastic matrix over those states, which is exactly the
object this verb consumes.

## Usage

``` r
markov_stability(x, normalize = TRUE)
```

## Arguments

- x:

  A `net_hon` (from
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)),
  `net_hypa`, `net_mogen`, `cograph_network`, `netobject` or `tna`
  object carrying a weight matrix; a row-stochastic numeric transition
  matrix; or sequence data accepted by
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
  (a wide `data.frame`, one trajectory per row, or a list of character
  vectors), in which case the first-order chain
  `build_hon(x, max_order = 1)` is used.

- normalize:

  Logical. Rescale rows that do not sum to 1 (with a
  `hypernets_renormalized` warning)? Default `TRUE`. `FALSE` makes a
  non-stochastic matrix an error instead.

## Value

An object of class `"net_markov_stability"`. Reach it with
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html):

- `as.data.frame(x)`:

  one row per state, columns `state`, `persistence` (\\P\_{ii}\\),
  `stationary_prob` (\\\pi_i\\), `return_time` (\\1/\pi_i\\),
  `sojourn_time` (\\1/(1 - P\_{ii})\\), `avg_time_to_others` (mean first
  passage leaving the state) and `avg_time_from_others` (mean first
  passage arriving at it).

- `as.data.frame(x, what = "passage_time")`:

  one row per ordered pair of states, columns `from`, `to`, `steps`; the
  diagonal carries the mean recurrence time.

- `as.data.frame(x, what = "stationary")`:

  one row per state, columns `state`, `stationary_prob`, `return_time`.

[`summary()`](https://rdrr.io/r/base/summary.html) names the attractor
and the stickiest state,
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
metrics as a faceted bar chart, and
[`print()`](https://rdrr.io/r/base/print.html) shows the table.

## Details

**Persistence** is the self-transition probability \\P\_{ii}\\.

**Sojourn time** is the expected number of consecutive steps spent in a
state before leaving it, \\1/(1 - P\_{ii})\\ — the mean of the geometric
dwell distribution. A state with \\P\_{ii} = 1\\ is absorbing and its
sojourn time is `Inf`.

**Return time** is the mean recurrence time \\1/\pi_i\\, the expected
number of steps between successive visits.

**avg_time_to_others** / **avg_time_from_others** are the off-diagonal
row and column means of the mean-first-passage matrix: how far the state
is from the rest of the chain, and how accessible it is from it (a small
value marks an attractor).

Mean first passage times use the Kemeny-Snell fundamental matrix \$\$Z =
(I - P + \Pi)^{-1}, \qquad M\_{ij} = \frac{Z\_{jj} -
Z\_{ij}}{\pi_j},\$\$ with \\\Pi\_{ij} = \pi_j\\. This requires an
**irreducible** chain: a reducible one has no unique stationary
distribution and \\(I - P + \Pi)\\ is singular. That is not a rare
corner for a higher-order network — pruning a `net_hon` with `min_freq`
can leave an absorbing or unreachable higher-order state — so
irreducibility is checked up front and reported as
`hypernets_not_ergodic` rather than surfacing as a LAPACK singularity.

## Conditions

Raises `hypernets_bad_input` for an unsupported `x`, a non-square or
single-state matrix, missing values, a state with no outgoing
transitions, and for rows that do not sum to 1 under
`normalize = FALSE`. Raises `hypernets_not_ergodic` (which also inherits
`hypernets_bad_input`) when the chain is reducible. Warns
`hypernets_renormalized` when rows are rescaled.

## References

Kemeny, J. G. and Snell, J. L. (1976). *Finite Markov Chains*.
Springer-Verlag.

Rosvall, M., Esquivel, A. V., Lancichinetti, A., West, J. D. and
Lambiotte, R. (2014). Memory in network flows and its effects on
spreading dynamics and community detection. *Nature Communications*, 5,
4630. [doi:10.1038/ncomms5630](https://doi.org/10.1038/ncomms5630)

Scholtes, I., Wider, N., Pfitzner, R., Garas, A., Tessone, C. J. and
Schweitzer, F. (2014). Causality-driven slow-down and speed-up of
diffusion in non-Markovian temporal networks. *Nature Communications*,
5, 5024. [doi:10.1038/ncomms6024](https://doi.org/10.1038/ncomms6024)

## See also

[`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md)
for the higher-order chain,
[`markov_order_test()`](https://mohsaqr.github.io/hypernets/reference/markov_order_test.md)
for how far the memory reaches,
[`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md)
for path-based centralities of the same network.

## Examples

``` r
# A first-order chain, straight from a transition matrix
P <- matrix(c(0.7, 0.2, 0.1,
              0.3, 0.5, 0.2,
              0.2, 0.3, 0.5),
            nrow = 3, byrow = TRUE,
            dimnames = list(c("A", "B", "C"), c("A", "B", "C")))
markov_stability(P)
#> Markov Stability -- 3 states
#> 
#>  state persistence stationary_prob return_time sojourn_time avg_time_to_others
#>      A         0.7          0.4634        2.16         3.33               6.20
#>      B         0.5          0.3171        3.15         2.00               5.18
#>      C         0.5          0.2195        4.56         2.00               4.03
#>  avg_time_from_others
#>                  3.95
#>                  4.23
#>                  7.22
#> 
#> as.data.frame(x) for the full table, as.data.frame(x, what = "passage_time") for first passage times.

# The higher-order walk: states that carry their memory
hon <- build_hon(split(human_long$code, human_long$session_id),
                 max_order = 2L, min_freq = 50L)
ms <- markov_stability(hon)
as.data.frame(ms, sort_by = "stationary_prob", top = 5L)
#>                  state persistence stationary_prob return_time sojourn_time
#> 1   Specify -> Command   0.0000000      0.11167480    8.954571     1.000000
#> 2   Command -> Specify   0.0000000      0.10021296    9.978749     1.000000
#> 3   Specify -> Specify   0.4842520      0.09193677   10.877041     1.938931
#> 4   Command -> Command   0.5784314      0.09161941   10.914718     2.372093
#> 5 Specify -> Interrupt   0.0000000      0.08410419   11.890014     1.000000
#>   avg_time_to_others avg_time_from_others
#> 1           81.10439             8.358104
#> 2           80.59760            10.552172
#> 3           80.48824            19.269256
#> 4           80.43635            23.299006
#> 5           79.03732            11.067830
as.data.frame(ms, what = "passage_time", sort_by = "steps", top = 5L)
#>                 from                 to    steps
#> 1  Specify -> Refine Specify -> Correct 361.1430
#> 2  Refine -> Specify Specify -> Correct 360.1430
#> 3             Refine Specify -> Correct 358.9391
#> 4 Specify -> Specify Specify -> Correct 357.5931
#> 5 Specify -> Inquire Specify -> Correct 357.1221
```
