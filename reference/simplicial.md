# Simplicial complex from sequences, a network or a memory network

Builds a simplicial complex: the sets of actions that occur together in
windows of consecutive actions (`type = "window"`), every clique of a
thresholded weighted network (`type = "clique"`), every higher-order
pathway of a memory network (`type = "pathway"`), or the Vietoris-Rips
complex of a distance matrix (`type = "vr"`). Every face of every
simplex is included.

## Usage

``` r
simplicial(
  x,
  type = "clique",
  threshold = 0,
  max_dim = 10L,
  max_pathways = NULL,
  anomaly = c("all", "over", "under"),
  max_scale = NULL,
  verify = TRUE,
  direction = c("either", "both"),
  window = 3L,
  min_count = 1L,
  validate = FALSE,
  alpha = 0.05,
  null = c("swap", "hypergeometric"),
  n_null = NULL,
  seed = NULL,
  action = NULL,
  actor = NULL,
  time = NULL,
  session = NULL,
  time_threshold = 900,
  timezone = "UTC",
  ...
)

# S3 method for class 'hypernets_simplicial'
plot(x, y, dismantled = FALSE, top = NULL, ...)

# S3 method for class 'hypernets_simplicial'
hg_get(x, ...)

# S3 method for class 'hypernets_simplicial'
print(x, n = 10L, ...)
```

## Arguments

- x:

  Sequences for `type = "window"`, in any form described in
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md);
  a square weighted matrix, a network object carrying one (`netobject`,
  `tna`) or sequences for `type = "clique"`; a memory network
  ([`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md),
  [`hypa()`](https://mohsaqr.github.io/hypernets/reference/hypa.md) on
  sequences,
  [`mogen()`](https://mohsaqr.github.io/hypernets/reference/mogen.md))
  for `type = "pathway"`; or for `type = "vr"` a distance matrix, a
  `dist` object, or points: a table or matrix with one row per point and
  one numeric column per coordinate (not square), or a data frame of
  numeric coordinates, whose one character column, if present, names the
  points. Points are turned into Euclidean distances. A square numeric
  matrix is read as distances.

- type:

  `"clique"` (default), `"window"`, `"pathway"`, or `"vr"` (alias
  `"rips"`).

- threshold:

  For `type = "clique"`: the smallest absolute edge weight kept (default
  `0`, every non-zero edge).

- max_dim:

  Integer. Highest simplex dimension (a k-simplex has k + 1 nodes).
  Default `10`.

- max_pathways:

  For `type = "pathway"`: keep at most this many pathways, ranked by
  count (HON) or ratio (HYPA). `NULL` keeps all.

- anomaly:

  For a HYPA pathway complex: `"all"` (default), `"over"` or `"under"`.

- max_scale:

  For `type = "vr"`: the longest edge admitted. `NULL` uses the largest
  distance.

- verify:

  Logical. Check a clique complex on construction (see Details). Default
  `TRUE`; set `FALSE` to skip the check on very large complexes.

- direction:

  For `type = "clique"` on a directed network: `"either"` (default, a
  pair is joined when the larger of its two weights reaches `threshold`)
  or `"both"` (when both weights do, as in `tna::cliques()`).

- window:

  For `type = "window"`: the number of consecutive actions in a window.
  Default `3`.

- min_count:

  For `type = "window"`: the number of windows in which a set of actions
  must occur to become a simplex. Default `1`. Not used with
  `validate = TRUE`.

- validate:

  For `type = "window"`: `TRUE` keeps only the sets of actions that pass
  the statistical validation described in Details. Default `FALSE`.

- alpha:

  For `validate = TRUE`: the false discovery rate. Default `0.05`.

- null:

  For `validate = TRUE`: `"swap"` (default) or `"hypergeometric"` (see
  Details).

- n_null:

  For `null = "swap"`: the number of shuffled copies. The smallest
  possible p-value is `1 / (n_null + 1)`, and a set can pass on its own
  only if this is at most `alpha` divided by the number of possible sets
  of its size. `NULL` (default) takes the smallest number that meets
  this for every size, and at least 999. A number given here that does
  not meet it raises the warning `hypernets_low_resolution` with the
  number that would.

- seed:

  For `null = "swap"`: a seed for the shuffles, applied locally so the
  global random state is restored. `NULL` (default) uses the current
  random state.

- action, actor, time, session, time_threshold, timezone:

  For `type = "window"`: long-format arguments (see
  [sequence-input](https://mohsaqr.github.io/hypernets/reference/sequence-input.md));
  leave the column names `NULL` for wide or list input.

- ...:

  For `type = "pathway"` on a network object: arguments passed to the
  memory-network builder.

- y:

  Unused.

- dismantled:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): `TRUE`
  draws each maximal simplex in its own panel; `FALSE` (default) draws
  them together on one layout.

- top:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html): the
  number of maximal simplices to draw, in the order described under
  Value. `NULL` (default) draws all of them.

- n:

  Number of rows of the default table to print. Default `10`.

## Value

A `simplicial_complex` object. Read it with
[`hg_get()`](https://mohsaqr.github.io/hypernets/reference/hg_get.md):
one row per simplex (`what = "simplices"`, the default), the face counts
(`"f_vector"`), the per-node simplicial degree (`"degree"`) or, for a
validated window complex, one row per tested set of actions
(`"validation"`), and the Betti numbers (`"betti"`).
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
maximal simplices, those no larger simplex contains, as regions around
their nodes
([`cograph::plot_simplicial()`](https://sonsoles.me/cograph/reference/plot_simplicial.html);
its arguments pass through `...`): the most significant first for a
validated window complex, the most frequent first for a window complex,
the closest first for a Vietoris-Rips complex.
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) returns the
figure. For the [`plot()`](https://rdrr.io/r/graphics/plot.default.html)
method, `x` is the complex.

## Details

`type = "window"` is a co-occurrence complex (Salnikov et al. 2018)
whose unit of co-occurrence is a window of `window` consecutive actions
in a sequence (Ding et al. 2020). The distinct actions of each window
form a set; a set observed in at least `min_count` windows becomes a
simplex, together with all its subsets. A window of `window` actions
yields simplices of dimension at most `window - 1`.

With `validate = TRUE` a window complex needs no count threshold. A set
of actions becomes a simplex when it occurs in more windows than a null
model predicts. With `null = "swap"` (default) the null is swap
randomization (Gionis et al. 2007): the windows and their actions are
shuffled by checkerboard swaps (Gotelli 2000) so that every window keeps
its number of actions and every action its number of windows, and only
which actions share a window changes. The p-value of a set is the share
of `n_null` shuffled copies in which the set occurs at least as often as
observed, counting the observed data (Phipson & Smyth 2010). With
`null = "hypergeometric"` the null is the statistically validated
hypergraph of Musciotto, Battiston & Mantegna (2021): a set of n actions
is tested against the windows with exactly n distinct actions, each
action falls into these windows independently of the others, and the
p-value is exact (their Eq. 3). That null ignores that a window holds
exactly n actions. It suits sparse co-occurrence, where most possible
sets of actions are rare or absent and each is expected in less than one
window, as with many actors who meet in small groups. On dense data,
such as a few actions that all occur together often, it expects more
co-occurrence than the windows can hold and validates little; when the
median tested set of a size is expected in one window or more it raises
the warning `hypernets_dense_cooccurrence`. For either null the p-values
of the sets of each size are corrected with the Benjamini-Hochberg false
discovery rate over all possible sets of that size, and a set is
validated when its adjusted p-value is at most `alpha`. Every action
remains a vertex. The tests are read with
`hg_get(x, what = "validation")`.

`type = "clique"` builds the clique (flag) complex of a network: every
clique of the thresholded network becomes a simplex (Giusti, Ghrist &
Bassett 2016). A directed network is made undirected by `direction`:
with `"either"` (default) a pair is joined when the larger of its two
weights reaches `threshold`, and with `"both"` when both weights do.
`"both"` gives the cliques of `tna::cliques()` with the same threshold
(Tikka, Lopez-Pernas & Saqr 2025): the triangles of the complex are its
cliques of size three, and so on. Given sequences instead of a network,
the clique complex is built on their transition network: the relative
transition probabilities that the tna family estimates from the same
sequences.

A clique complex is checked on construction (`verify = TRUE`): its
simplices must be exactly the cliques igraph enumerates on the same
thresholded graph (when igraph is installed), and its Euler
characteristic from the face counts must equal the alternating sum of
its Betti numbers (the Euler-Poincare identity). A failed check raises
the warning `hypernets_simplicial_unverified`; the complex is returned
either way.

## References

Hatcher, A. (2002). *Algebraic Topology*. Cambridge University Press.

Giusti, C., Ghrist, R., & Bassett, D. S. (2016). Two's company, three
(or more) is a simplex. *Journal of Computational Neuroscience*, 41,
1-14.
[doi:10.1007/s10827-016-0608-6](https://doi.org/10.1007/s10827-016-0608-6)

Tikka, S., Lopez-Pernas, S., & Saqr, M. (2025). tna: An R package for
transition network analysis. *Applied Psychological Measurement*, 49,
326-328.
[doi:10.1177/01466216251348840](https://doi.org/10.1177/01466216251348840)

Salnikov, V., Cassese, D., Lambiotte, R., & Jones, N. S. (2018).
Co-occurrence simplicial complexes in mathematics: identifying the holes
of knowledge. *Applied Network Science*, 3, 37.
[doi:10.1007/s41109-018-0074-3](https://doi.org/10.1007/s41109-018-0074-3)

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification.
*Proceedings of EMNLP 2020*, 4927-4936.
[doi:10.18653/v1/2020.emnlp-main.399](https://doi.org/10.18653/v1/2020.emnlp-main.399)

Gionis, A., Mannila, H., Mielikainen, T., & Tsaparas, P. (2007).
Assessing data mining results via swap randomization. *ACM Transactions
on Knowledge Discovery from Data*, 1(3), 14.
[doi:10.1145/1297332.1297338](https://doi.org/10.1145/1297332.1297338)

Gotelli, N. J. (2000). Null model analysis of species co-occurrence
patterns. *Ecology*, 81(9), 2606-2621.
[doi:10.1890/0012-9658(2000)081\[2606:NMAOSC\]2.0.CO;2](https://doi.org/10.1890/0012-9658%282000%29081%5B2606%3ANMAOSC%5D2.0.CO%3B2)

Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
zero. *Statistical Applications in Genetics and Molecular Biology*,
9(1), 39.
[doi:10.2202/1544-6115.1585](https://doi.org/10.2202/1544-6115.1585)

Musciotto, F., Battiston, F., & Mantegna, R. N. (2021). Detecting
informative higher-order interactions in statistically validated
hypergraphs. *Communications Physics*, 4, 218.
[doi:10.1038/s42005-021-00710-4](https://doi.org/10.1038/s42005-021-00710-4)

Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
interactions: Structure and dynamics. *Physics Reports*, 874, 1–92.
[doi:10.1016/j.physrep.2020.05.004](https://doi.org/10.1016/j.physrep.2020.05.004)

## See also

[`hg_betti()`](https://mohsaqr.github.io/hypernets/reference/hg_betti.md),
[`hg_euler()`](https://mohsaqr.github.io/hypernets/reference/hg_euler.md),
[`hg_qanalysis()`](https://mohsaqr.github.io/hypernets/reference/hg_qanalysis.md),
[`hg_homology()`](https://mohsaqr.github.io/hypernets/reference/hg_homology.md)

## Examples

``` r
adj <- matrix(c(0, 1, 1, 0,
                1, 0, 1, 1,
                1, 1, 0, 1,
                0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
sc <- simplicial(adj)
hg_get(sc)
#>    simplex dimension size members
#> 1        1         0    1       c
#> 2        2         0    1       b
#> 3        3         1    2    b, c
#> 4        4         0    1       d
#> 5        5         1    2    b, d
#> 6        6         2    3 b, c, d
#> 7        7         1    2    c, d
#> 8        8         0    1       a
#> 9        9         1    2    a, b
#> 10      10         2    3 a, b, c
#> 11      11         1    2    a, c
hg_get(sc, what = "f_vector")
#>   dimension count
#> 1         0     4
#> 2         1     5
#> 3         2     2
```
