# Relate higher-order structure to an actor-level outcome

Regresses an outcome measured once per actor on per-actor features built
from a fitted higher-order network. It is an assembly and inference
layer: the predictors come from
[`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md)
and
[`build_honem()`](https://mohsaqr.github.io/hypernets/reference/build_honem.md),
and the only new step is the regression.

## Usage

``` r
hon_outcome(
  x,
  outcome,
  sequences,
  by = "actor",
  features = "centrality",
  nested_in = NULL,
  family = c("gaussian", "binomial", "poisson"),
  outcome_col = NULL,
  action = NULL,
  time = NULL,
  centrality_type = c("pagerank", "betweenness", "closeness"),
  honem_dim = 2L,
  weighted = FALSE,
  standardize = TRUE,
  vcov = c("CR3", "CR1"),
  p_adjust = "BH",
  min_clusters = 20L,
  max_paths = 1e+06,
  ...
)
```

## Arguments

- x:

  A `net_hon` from
  [`build_hon()`](https://mohsaqr.github.io/hypernets/reference/build_hon.md).

- outcome:

  One value per actor: a named numeric vector keyed by actor, or a
  data.frame holding the actor key in column `by` and the outcome in
  `outcome_col` (or in its only other column).

- sequences:

  The sequences the network was built from, keyed by actor: a named list
  of character vectors, a long data.frame together with `action` (and
  optionally `time`), or a wide data.frame whose `by` column holds the
  actor and whose remaining columns are time steps. Required, because a
  `net_hon` does not carry the actor identity of the trajectories it was
  built from.

- by:

  Name of the actor key, used for the column in `outcome`, `sequences`
  and the returned features table. Default `"actor"`.

- features:

  Which per-actor predictors to build: any of `"centrality"` and
  `"honem"`, or a data.frame of ready-made per-actor features holding
  the actor key in column `by` and one numeric column per feature.

- nested_in:

  Optional clustering variable for cluster-robust standard errors: the
  name of a column of `outcome`, or a vector named by actor. `NULL`
  (default) uses model-based standard errors.

- family:

  `"gaussian"` (default), `"binomial"` or `"poisson"`; the last two use
  the canonical link.

- outcome_col:

  Name of the outcome column when `outcome` is a data.frame with more
  than one non-key column.

- action, time:

  Long-format column names of `sequences`: `action` holds the state of
  each event and `time` orders events within an actor (row order when
  `NULL`). Leave both `NULL` for a named list or a wide data.frame.

- centrality_type:

  Which centralities to attach to the higher-order nodes, passed to
  [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md).
  Betweenness and closeness enumerate shortest paths and are the
  expensive ones.

- honem_dim:

  Integer. Number of HONEM dimensions to use as features when `"honem"`
  is requested. Default 2.

- weighted:

  Logical. Passed to
  [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md):
  use the higher-order edge weights instead of the binary pattern.
  Default `FALSE`, matching
  [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md).

- standardize:

  Logical. Z-score each feature over the analysis sample. Default
  `TRUE`.

- vcov:

  Cluster covariance estimator, used only when `nested_in` is given.
  `"CR3"` (default) is the cluster jackknife, matching
  `sandwich::vcovCL(type = "HC3")`; `"CR1"` is the Stata estimator,
  `sandwich::vcovCL(type = "HC1", cadjust = TRUE)`. CR1 under-covers at
  moderate cluster counts – measured 0.933 against a nominal 0.95 at 40
  clusters with a cluster-correlated regressor, where CR3 reached 0.953
  – which is why CR3 is the default. CR1 is kept because it is what most
  published results report. CR3 is undefined for a cluster whose
  leverage is 1 (typically a singleton cluster) and raises there.

- p_adjust:

  Multiplicity correction applied across the feature coefficients,
  passed to [`p.adjust`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"`; `"none"` switches it off and is recorded as such.

- min_clusters:

  Integer. Warn below this many clusters. Default 20.

- max_paths:

  Passed to
  [`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md).

- ...:

  Ignored; present so the verb can grow arguments without breaking
  positional calls.

## Value

An object of class `net_outcome`. Reach into it with
[`as.data.frame()`](https://mohsaqr.github.io/hypernets/reference/as.data.frame.net_outcome.md),
which returns one row per feature with `feature`, `estimate`,
`std_error`, `conf_low`, `conf_high`, `statistic`, `p`, `p_adj` and `n`;
`what = "terms"` adds the intercept, `what = "features"` returns the
per-actor design table, and `what = "dropped"` the actors that were
excluded and why. [`summary()`](https://rdrr.io/r/base/summary.html)
reports the model, the covariance and the correction;
[`plot()`](https://rdrr.io/r/graphics/plot.default.html) draws the
coefficient plot with its intervals.

## Details

**How a per-actor feature is built.** Node-level values are attached to
every higher-order node of `x` – a centrality from
`hon_centrality(x, project = FALSE)`, a HONEM coordinate from
`build_honem(x)`, or both. Each actor's own sequences are then decoded
against that same node set with the *longest-suffix* rule: at position
\\i\\ the actor occupies the longest-history node
`"s[i-k+1] -> ... -> s[i]"` that the network contains, which is the rule
BuildHON itself uses when it rewires edges onto higher-order nodes. The
actor's feature is the mean of the node value over that actor's own
visits, pooling every sequence that actor contributed: \$\$f_a = \sum_v
c\_{av} \\ \mathrm{value}(v) / \sum_v c\_{av}\$\$ with \\c\_{av}\\ the
number of times actor \\a\\ occupies node \\v\\. Positions whose state
never entered the network (pruned by `min_freq`, or unseen) match no
node, are excluded from the mean, and are counted in `n_unmatched` in
the features table.

**Standardization.** With `standardize = TRUE` (the default) every
feature is z-scored over the analysis sample, so an estimate is the
change in the outcome per one standard deviation of that feature, and
features on wildly different scales (a PageRank near \\10^{-3}\\, a
HONEM coordinate near 1) are comparable. The outcome is never
transformed.

**Inference.** Gaussian outcomes are fitted by ordinary least squares
(QR); `family = "binomial"` and `"poisson"` by IRLS with the canonical
link. Without `nested_in` the covariance is the usual model-based one,
with \\t\_{n-k}\\ intervals for the gaussian family and normal intervals
otherwise. With `nested_in` the covariance is the cluster-robust
sandwich \$\$V = c \\ B \left(\sum_g u_g u_g'\right) B, \quad u_g =
\sum\_{i \in g} x_i e_i\$\$ with bread \\B = (X'X)^{-1}\\ (gaussian) or
\\(X'WX)^{-1}\\ (GLM, \\W\\ the IRLS working weights), \\e = y -
\hat\mu\\, and the CR1S ("Stata") small-sample correction \\c =
\frac{G}{G-1}\cdot \frac{n-1}{n-k}\\; intervals then use \\t\_{G-1}\\.
This is exactly `sandwich::vcovCL(type = "HC1", cadjust = TRUE)`, and
the test suite checks it against that package where it is installed.
Cluster-robust standard errors are badly biased with few clusters, so
fewer than `min_clusters` clusters raises a `hypernets_few_clusters`
warning.

**Measured coverage.** In the package's own simulation
(`tests/testthat/test-memory_outcome.R`, 2000 replicates with the design
held fixed and only the outcome redrawn) the model-based 95% interval
covered a known effect 95.3% of the time, and the default cluster-robust
interval (`vcov = "CR3"`) covered 95.3% with 40 clusters *when the
predictor is itself cluster-correlated* – the hardest case, and the one
where ignoring the clustering is catastrophic (the model-based interval
covered only 72.0% there).

`vcov = "CR1"` covers 93.3% on the same design. That shortfall is a
property of the CR1 estimator at a moderate number of clusters, not of
this implementation: the identical simulation run through
`sandwich::vcovCL()` gives the same 93.3%, with 94.1% for CR2 and 95.3%
for CR3. CR3 is therefore the default, and CR1 is kept because it is
what most published results report – read a CR1 interval from a few
dozen clusters as slightly optimistic.

**What the intervals do and do not cover.** The higher-order network is
treated as fixed. The features are estimated from the pooled data, which
includes each actor's own sequences, so the intervals are the usual
regression intervals *conditional on the estimated structure*; they do
not propagate the uncertainty in the network itself. Use
[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md)
for that uncertainty.

**Collinearity.** A rank-deficient design raises
`hypernets_rank_deficient` rather than silently dropping terms, and a
scaled condition index above 30 (Belsley, Kuh and Welsch 1980) warns
with the same class.

**Conditions raised.** `hypernets_bad_input` (an input that cannot be
modelled at all: a malformed outcome, duplicate actor keys, a missing
cluster label, fewer actors than terms, fewer clusters than terms);
`hypernets_rank_deficient` (error on an exactly collinear or constant
design, warning on a near-collinear one); `hypernets_dropped_actors`
(warning: events that matched no node, or actors excluded from the fit –
both are listed in the result); `hypernets_few_clusters` (warning);
`hypernets_no_converge` (warning: IRLS hit its iteration bound);
`hypernets_degenerate_fit` (warning: separation in a binomial fit, whose
standard errors are then meaningless). Nothing is dropped or rescued
silently.

## References

Belsley, D. A., Kuh, E., & Welsch, R. E. (1980). *Regression
Diagnostics: Identifying Influential Data and Sources of Collinearity*.
Wiley.

Benjamini, Y., & Hochberg, Y. (1995). Controlling the false discovery
rate. *Journal of the Royal Statistical Society B*, 57(1), 289-300.
[doi:10.1111/j.2517-6161.1995.tb02031.x](https://doi.org/10.1111/j.2517-6161.1995.tb02031.x)

Cameron, A. C., & Miller, D. L. (2015). A practitioner's guide to
cluster-robust inference. *Journal of Human Resources*, 50(2), 317-372.
[doi:10.3368/jhr.50.2.317](https://doi.org/10.3368/jhr.50.2.317)

Liang, K.-Y., & Zeger, S. L. (1986). Longitudinal data analysis using
generalized linear models. *Biometrika*, 73(1), 13-22.
[doi:10.1093/biomet/73.1.13](https://doi.org/10.1093/biomet/73.1.13)

Scholtes, I., Wider, N., & Garas, A. (2016). Higher-order aggregate
networks in the analysis of temporal networks. *European Physical
Journal B*, 89, 61.
[doi:10.1140/epjb/e2016-60663-0](https://doi.org/10.1140/epjb/e2016-60663-0)

## See also

[`hon_centrality()`](https://mohsaqr.github.io/hypernets/reference/hon_centrality.md),
[`build_honem()`](https://mohsaqr.github.io/hypernets/reference/build_honem.md),
[`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md)

## Examples

``` r
# 60 actors, each a different mixture of four work motifs; the motifs
# carry the second-order dependence the network will pick up.
set.seed(42)
motifs <- list(c("plan", "code", "test"), c("debug", "code", "debug"),
               c("read", "plan", "code"), c("test", "read", "plan"))
actors <- sprintf("s%02d", seq_len(60))
mix <- matrix(stats::runif(60 * 4), nrow = 60)
seqs <- stats::setNames(
  lapply(seq_len(60), function(i)
    unlist(motifs[sample(4, 10, replace = TRUE, prob = mix[i, ])],
           use.names = FALSE)),
  actors)
hon <- build_hon(seqs, max_order = 2)

scores <- stats::setNames(60 + 8 * mix[, 1] + stats::rnorm(60, sd = 2),
                          actors)
fit <- hon_outcome(hon, outcome = scores, sequences = seqs)
fit
#> Higher-order structure -> outcome  [gaussian, model-based covariance]
#>   60 actors / 4 terms / features: centrality
#>      feature estimate std_error conf_low conf_high statistic        p p_adj  n
#>  (Intercept)   64.600     0.360  63.8000     65.30   180.000 5.28e-79    NA 60
#>     pagerank   -1.380     1.840  -5.0600      2.30    -0.752 4.55e-01 0.683 60
#>  betweenness    1.480     0.709   0.0588      2.90     2.090 4.15e-02 0.125 60
#>    closeness   -0.562     1.850  -4.2600      3.14    -0.304 7.62e-01 0.762 60
#> 
#>   p_adj: BH across 3 features
as.data.frame(fit, sort_by = "p_adj")
#>       feature  estimate std_error    conf_low conf_high  statistic          p
#> 1 betweenness  1.479653 0.7092566  0.05884117  2.900465  2.0862025 0.04152706
#> 2    pagerank -1.382241 1.8378594 -5.06391571  2.299434 -0.7520928 0.45514490
#> 3   closeness -0.561793 1.8465105 -4.26079793  3.137212 -0.3042458 0.76206839
#>       p_adj  n
#> 1 0.1245812 60
#> 2 0.6827173 60
#> 3 0.7620684 60

# \donttest{
# actors nested in classes: cluster-robust standard errors
classes <- stats::setNames(rep(sprintf("c%02d", seq_len(20)), each = 3),
                           actors)
nested <- hon_outcome(hon, outcome = scores, sequences = seqs,
                      nested_in = classes)
summary(nested)
#> Higher-order outcome model
#>   family              : gaussian (link: identity)
#>   actors in the fit   : 60 (0 dropped)
#>   features            : centrality [pagerank, betweenness, closeness]
#>   predictors          : z-scored (estimate = change per 1 SD)
#>   covariance          : cluster-robust (CR3) over 20 clusters
#>   small-sample factor : CR3: cluster jackknife, (I - H_gg)^-1 per cluster
#>   intervals           : 95%, t(19)
#>   multiplicity        : BH across 3 features
#>   scaled condition    : 11.8
#>   R-squared           : 0.1217 (adjusted 0.0746)
#>   residual SD         : 2.7857
#> 
#>      feature  estimate std_error   conf_low conf_high  statistic          p
#>     pagerank -1.382241 2.2815480 -6.1575757  3.393094 -0.6058346 0.55179401
#>  betweenness  1.479653 0.6490562  0.1211627  2.838143  2.2796992 0.03435012
#>    closeness -0.561793 2.3094454 -5.3955178  4.271932 -0.2432588 0.81041190
#>      p_adj  n
#>  0.8104119 60
#>  0.1030504 60
#>  0.8104119 60
#> 
#>   The network is treated as fixed: these intervals condition on the
#>   estimated higher-order structure and do not propagate its uncertainty.
plot(nested)

# }
```
