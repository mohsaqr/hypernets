# Mixed-membership communities of a hypergraph (Hy-MMSBM)

Fits the hypergraph mixed-membership stochastic block model of Ruggeri
et al. (2023), which gives every node a membership vector over `k`
communities from a likelihood model of the whole hypergraph. The weight
\\A_e\\ of a hyperedge \\e\\ is Poisson with mean \\\lambda_e /
\kappa\_{\|e\|}\\, where \$\$\lambda_e = \sum\_{i \< j \in e} u_i^\top
w\\ u_j,\$\$ \\u\\ is the non-negative \\N \times k\\ membership matrix,
\\w\\ the symmetric non-negative \\k \times k\\ affinity between
communities and \\\kappa_n = \frac{n(n-1)}{2}\binom{N-2}{n-2}\\
normalises for the size \\n\\ of the hyperedge. Summed over every
possible hyperedge up to the largest size \\D\\, the normalisers give
the constant \\C = \sum\_{n=2}^{D} 2 / (n(n-1)) = 2(1 - 1/D)\\ and the
log-likelihood (Eq. 5, up to additive constants that do not depend on
\\u\\, \\w\\) \$\$\mathcal{L} = -C \sum\_{i\<j} u_i^\top w\\ u_j +
\sum_e A_e \log \lambda_e.\$\$

## Usage

``` r
hg_mmsbm(
  hg,
  k,
  assortative = FALSE,
  nstart = 10L,
  max_iter = 5000L,
  tol = 1e-05,
  check_every = 10L,
  criterion = c("membership", "parameters"),
  w_prior = 1,
  u_prior = 0,
  max_size = NULL,
  edge_weights = NULL,
  seed = NULL
)

# S3 method for class 'net_hg_mmsbm'
hg_get(
  x,
  what = c("membership", "nodes", "affinity", "restarts"),
  ...,
  top = NULL
)

# S3 method for class 'net_hg_mmsbm'
print(x, n = 10L, ...)

# S3 method for class 'net_hg_mmsbm'
plot(x, type = c("membership", "affinity", "restarts"), labels = NULL, ...)
```

## Arguments

- hg:

  A hypernets `net_hg`
  ([`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md),
  [`group_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/group_hypergraph.md),
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md),
  ...).

- k:

  Number of communities, a whole number between 1 and the number of
  nodes.

- assortative:

  If `TRUE`, `w` is diagonal (communities interact only within
  themselves); default `FALSE`, the full affinity.

- nstart:

  Number of random starts (default `10L`, the authors'
  `--training_rounds`).

- max_iter:

  Maximum EM steps per start (default `5000L`).

- tol:

  Convergence tolerance (default `1e-5`), compared with the change
  between two consecutive steps every `check_every` steps. A start that
  hits `max_iter` first is reported as not converged; if the kept start
  did not converge, a `hypernets_no_converge` warning is raised.

- check_every:

  Steps between convergence checks (default `10L`, as the authors').

- criterion:

  What must stop changing: `"membership"` (default), the largest
  absolute change of any membership; or `"parameters"`, the authors'
  rule, the Frobenius change of `w` divided by `k` and of `u` divided by
  the number of nodes. The likelihood is unchanged when `u` is
  multiplied by \\c\\ and `w` divided by \\c^2\\; with the default
  `w_prior` and no `u_prior` the fit keeps drifting along that direction
  (`w` shrinks, `u` grows) while the memberships and the log-likelihood
  have settled, so `"parameters"` rarely fires there.

- w_prior, u_prior:

  Rates of the exponential priors on `w` and `u` (MAP estimation); `0`
  means none (maximum likelihood). The defaults `w_prior = 1`,
  `u_prior = 0` are the authors' "half-Bayesian" setting that fixes the
  scale of the parameters. As in the authors' code, the `w` prior acts
  on \\C w\\ (the parametrisation used during the iterations).

- max_size:

  `NULL` (default: every hyperedge) or the largest hyperedge size kept;
  it is also the \\D\\ in \\C\\. With `NULL`, \\D\\ is the largest size
  in the data.

- edge_weights:

  `NULL` or positive hyperedge weights \\A_e\\, one per hyperedge or one
  recycled. `NULL` uses the window counts of a
  [`window_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/window_hypergraph.md)
  and 1 otherwise. The incidence is read as binary: the stored cell
  weights (word counts, tf-idf) do not enter the model.

- seed:

  `NULL` (default: the current random stream) or a whole number for the
  random starts; the caller's stream is restored on exit.

- x:

  A `net_hg_mmsbm` object.

- what:

  Which table (see Value).

- ...:

  Unused; for S3 consistency.

- top:

  `NULL` (default, every row) or the number of first rows to return.

- n:

  Number of rows of the default table to print. Default `10`.

- type:

  For [`plot()`](https://rdrr.io/r/graphics/plot.default.html):
  `"membership"` (default) stacks each node's memberships in a bar,
  nodes grouped in panels by hard community; `"affinity"` draws `w` as a
  labelled tile map; `"restarts"` the log-likelihood of every start,
  shaped by convergence.

- labels:

  For `plot(type = "membership")`: show node names on the axis (default:
  when there are at most 40 nodes).

## Value

An object of class `net_hg_mmsbm`, a list read through
`hg_get(x, what = )`:

- `"membership"` (default): one row per node and community, columns
  `node`, `community`, `membership` (sums to one over a node's rows) and
  `membership_weight` (the fitted, unnormalised membership parameter
  `u`);

- `"nodes"`: one row per node, `node`, `community` (the largest
  membership) and `membership` (its value);

- `"affinity"`: one row per community pair with `from <= to`, columns
  `from`, `to`, `affinity` (the fitted `w`);

- `"restarts"`: one row per start, `run`, `log_likelihood`,
  `iterations`, `converged`, `best`, `ari_to_best` (adjusted Rand index
  of the start's hard partition with the kept one).

[`print()`](https://rdrr.io/r/base/print.html) gives the fit,
[`summary()`](https://rdrr.io/r/base/summary.html) one row per community
(`community`, `size` hard-assigned nodes, `mass` total membership,
`mean_membership` of its assigned nodes, `w_within`), and
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) the
memberships, the affinity or the restarts. Raises `hypernets_bad_input`
for a non-hypergraph, an invalid argument or no hyperedge of size two or
more.

## Details

The fit is the authors' expectation-maximisation (Eqs. 8-9, with the
exponential priors of their Appendix A), ported from their reference
code: in each step \\w\\ is updated, then \\u\\; \\C\\ is absorbed into
\\w\\ during the iterations and divided out at the end. Every quantity
goes through \\\lambda_e = (s_e^\top w s_e - \sum\_{i\in e} u_i^\top w
u_i)/2\\ with \\s_e = \sum\_{i \in e} u_i\\ and sparse products with the
incidence matrix, so one step costs \\O(\mathrm{nnz}\\k + E k^2)\\ and
hyperedges of any size (a word used in every document) are handled
exactly. EM finds a local optimum only: the fit is run from `nstart`
random starts and the one with the highest \\\mathcal{L}\\ is kept
(Algorithm 1). The restarts table reports the log-likelihood, the
iterations, the convergence flag of every start, and the adjusted Rand
index of its hard partition with the kept one, so the stability of the
solution is visible.

Hyperedges of size one carry no pair and are outside the model (their
\\\lambda_e\\ is zero); they are left out and counted. `max_size`
optionally leaves out larger hyperedges too, as the authors'
`--max_hye_size`; the default `NULL` keeps them all. A node left with no
hyperedge has \\u_i = 0\\ after the first step and no membership (`NA`),
with a `hypernets_isolated_nodes` warning. A node that does lie in
fitted hyperedges can still have its memberships driven to zero (or to
subnormal remnants) by the multiplicative EM updates; any row whose
total is below the precision of the largest row is treated as having no
membership (`NA`) and reported with a `hypernets_collapsed_membership`
warning, rather than normalising remnants into spurious exact mixtures.

The membership of node \\i\\ in community \\k\\ is \\u\_{ik} / \sum_q
u\_{iq}\\; the hard `community` is its largest entry. The size of
\\u_i\\ reflects how active the node is, its direction the communities
it mixes; at the maximum of the likelihood many memberships sit at or
near 0 and 1, and mixing is read from the nodes that do not. Communities
are numbered by decreasing total membership. On a
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
with documents as nodes, the communities are topics and the memberships
the documents' topic shares.

## References

Ruggeri, N., Contisciani, M., Battiston, F., & De Bacco, C. (2023).
Community detection in large hypergraphs. *Science Advances*, 9(28),
eadg9159.
[doi:10.1126/sciadv.adg9159](https://doi.org/10.1126/sciadv.adg9159)

Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
Classification*, 2, 193–218.
[doi:10.1007/BF01908075](https://doi.org/10.1007/BF01908075)

## See also

[`hg_membership()`](https://mohsaqr.github.io/hypernets/reference/hg_membership.md)
for soft membership read from a hard spectral partition.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots and salt",
  cooking_2 = "this soup recipe needs onions, carrots and salt",
  cooking_3 = "carrots and onions make a sweet soup",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars",
  space_3 = "a distant galaxy full of stars",
  both = "astronomers eat soup with carrots under the stars"
), stop_words = c("the", "with", "and", "a", "this", "at", "of", "under"))
fit <- hg_mmsbm(hg, k = 2, seed = 1)
fit
#> Hy-MMSBM mixed-membership communities (full affinity), k = 2
#>   Nodes: 7 | Hyperedges fitted: 9 (left out: 9 of size 1, 0 above max_size) | D = 4
#>   Log-likelihood: -5.43538 | Starts: 10 (10 converged) | Kept start: 1, 31 steps, converged
#>       node   community    membership membership_weight
#>       both Community 1  9.999990e-01      1.406661e+00
#>       both Community 2  9.776189e-07      1.375179e-06
#>  cooking_1 Community 1  9.889336e-39      7.809211e-38
#>  cooking_1 Community 2  1.000000e+00      7.896598e+00
#>  cooking_2 Community 1  9.296310e-39      7.361566e-38
#>  cooking_2 Community 2  1.000000e+00      7.918805e+00
#>  cooking_3 Community 1  1.643378e-34      4.208042e-38
#>  cooking_3 Community 2  1.000000e+00      2.560604e-04
#>    space_1 Community 1  1.000000e+00      7.461251e+00
#>    space_1 Community 2 7.973192e-100      5.948999e-99
#> ... 4 more rows
hg_get(fit, what = "nodes")
#>        node   community membership
#> 1      both Community 1   0.999999
#> 2 cooking_1 Community 2   1.000000
#> 3 cooking_2 Community 2   1.000000
#> 4 cooking_3 Community 2   1.000000
#> 5   space_1 Community 1   1.000000
#> 6   space_2 Community 1   1.000000
#> 7   space_3 Community 1   1.000000
hg_get(fit, what = "restarts")
#>    run log_likelihood iterations converged  best ari_to_best
#> 1    1      -5.435376         31      TRUE  TRUE   1.0000000
#> 2    2      -9.516172        121      TRUE FALSE  -0.1666667
#> 3    3      -9.516174        121      TRUE FALSE  -0.2162162
#> 4    4      -5.435588         31      TRUE FALSE   1.0000000
#> 5    5      -9.516228         91      TRUE FALSE  -0.1666667
#> 6    6      -5.435951         31      TRUE FALSE   1.0000000
#> 7    7      -9.516205        101      TRUE FALSE  -0.1666667
#> 8    8      -9.516165        121      TRUE FALSE  -0.1666667
#> 9    9      -5.435666         31      TRUE FALSE   1.0000000
#> 10  10      -9.516185        111      TRUE FALSE  -0.1666667
```
