testthat::skip_on_cran()

# ---- hon_outcome() -------------------------------------------------------
#
# The dangerous verb: an assembly layer that could quietly produce confident
# nonsense. Everything numerical here is checked against an independent
# implementation -- stats::lm(), stats::glm(), a hand-written sandwich, and
# sandwich::vcovCL() -- plus a coverage simulation that puts a number on
# whether the intervals mean what they say.

# A corpus with genuine second-order structure: each actor is a different
# mixture of four motifs, so per-actor exposure to the higher-order contexts
# varies.
.mo_corpus <- function(n_actors = 60L, n_motifs = 10L, seed = 42L) {
  set.seed(seed)
  motifs <- list(c("plan", "code", "test"), c("debug", "code", "debug"),
                 c("read", "plan", "code"), c("test", "read", "plan"))
  actors <- sprintf("s%02d", seq_len(n_actors))
  mix <- matrix(stats::runif(n_actors * 4L), nrow = n_actors)
  seqs <- stats::setNames(
    lapply(seq_len(n_actors), function(i)
      unlist(motifs[sample(4L, n_motifs, replace = TRUE, prob = mix[i, ])],
             use.names = FALSE)),
    actors)
  list(actors = actors, seqs = seqs, mix = mix,
       hon = build_hon(seqs, max_order = 2L),
       outcome = stats::setNames(60 + 8 * mix[, 1L] +
                                   stats::rnorm(n_actors, sd = 2), actors))
}

# A minimal real net_hon, for the paths where features are supplied ready-made
# and the network itself is never consulted.
.mo_hon0 <- function() {
  build_hon(list(c("A", "B", "A", "B"), c("B", "A", "B", "A")), max_order = 1L)
}

# A fixed per-actor feature frame: the design the inference runs on, with no
# network work in the loop.
.mo_features <- function(n = 200L, seed = 7L) {
  set.seed(seed)
  data.frame(actor = sprintf("a%03d", seq_len(n)),
             f1 = stats::rnorm(n),
             f2 = stats::rnorm(n),
             stringsAsFactors = FALSE)
}

# ---- the feature aggregation --------------------------------------------

test_that("longest-suffix decoding picks the longest history in the network", {
  nodes <- c("A", "B", "C", "A -> B", "A -> B -> C")
  index <- stats::setNames(seq_along(nodes), nodes)
  idx <- hypernets:::.hoo_decode(c("A", "B", "C", "B"), index, max_order = 3L)
  # position 1: "A"; 2: "A -> B"; 3: "A -> B -> C"; 4: "B" (no "B -> C -> B")
  expect_identical(nodes[idx], c("A", "A -> B", "A -> B -> C", "B"))

  # a state the network never saw matches nothing and is counted as unmatched
  idx2 <- hypernets:::.hoo_decode(c("A", "Z", "B"), index, max_order = 3L)
  expect_identical(idx2, c(1L, NA_integer_, 2L))

  # max_order caps the history length that is even tried
  idx3 <- hypernets:::.hoo_decode(c("A", "B", "C"), index, max_order = 2L)
  expect_identical(nodes[idx3], c("A", "A -> B", "C"))
})

test_that("a per-actor feature is the visit-weighted mean of the node value", {
  cp <- .mo_corpus(n_actors = 20L)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     centrality_type = "pagerank", standardize = FALSE)
  feats <- as.data.frame(fit, what = "features")

  # recompute one actor's feature by hand from the published verbs
  cen <- hon_centrality(cp$hon, type = "pagerank", project = FALSE)
  value <- stats::setNames(cen$pagerank, cen$node)
  one <- cp$seqs[["s01"]]
  index <- stats::setNames(seq_along(names(value)), names(value))
  visited <- names(value)[hypernets:::.hoo_decode(one, index, max_order = 2L)]
  expect_equal(feats$pagerank[feats$actor == "s01"],
               mean(value[visited]), tolerance = 1e-12)
  # every event of this corpus lives in the network
  expect_identical(sum(feats$n_unmatched), 0L)
  expect_identical(feats$n_events, vapply(cp$seqs[feats$actor], length,
                                          integer(1L), USE.NAMES = FALSE))
})

test_that("the three sequence input shapes give the same features", {
  cp <- .mo_corpus(n_actors = 20L)
  from_list <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                           centrality_type = "pagerank")

  long <- data.frame(
    actor = rep(names(cp$seqs), lengths(cp$seqs)),
    turn = unlist(lapply(cp$seqs, seq_along), use.names = FALSE),
    code = unlist(cp$seqs, use.names = FALSE),
    stringsAsFactors = FALSE)
  from_long <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = long,
                           action = "code", time = "turn",
                           centrality_type = "pagerank")

  width <- max(lengths(cp$seqs))
  wide <- data.frame(actor = names(cp$seqs), stringsAsFactors = FALSE)
  wide <- cbind(wide, do.call(rbind, lapply(cp$seqs, function(s)
    c(s, rep(NA_character_, width - length(s))))))
  rownames(wide) <- NULL
  from_wide <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = wide,
                           centrality_type = "pagerank")

  expect_equal(as.data.frame(from_long), as.data.frame(from_list))
  expect_equal(as.data.frame(from_wide), as.data.frame(from_list))
})

test_that("a supplied features data.frame reproduces the derived one", {
  cp <- .mo_corpus(n_actors = 40L)
  derived <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                         centrality_type = c("pagerank", "closeness"))
  feats <- as.data.frame(derived, what = "features")
  supplied <- hon_outcome(cp$hon, outcome = cp$outcome,
                          features = feats[, c("actor", "pagerank",
                                               "closeness")])
  expect_equal(as.data.frame(supplied), as.data.frame(derived))
})

# ---- equivalence with stats::lm() and stats::glm() -----------------------

test_that("gaussian estimates and standard errors match stats::lm() to 1e-10", {
  cp <- .mo_corpus(n_actors = 60L)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     standardize = FALSE)
  d <- as.data.frame(fit, what = "features")
  ref <- stats::lm(outcome ~ pagerank + betweenness + closeness, data = d)
  ref_tab <- summary(ref)$coefficients

  ours <- as.data.frame(fit, what = "terms")
  expect_identical(ours$feature, rownames(ref_tab))
  expect_equal(ours$estimate, unname(ref_tab[, "Estimate"]), tolerance = 1e-10)
  expect_equal(ours$std_error, unname(ref_tab[, "Std. Error"]),
               tolerance = 1e-10)
  expect_equal(ours$statistic, unname(ref_tab[, "t value"]), tolerance = 1e-10)
  expect_equal(ours$p, unname(ref_tab[, "Pr(>|t|)"]), tolerance = 1e-10)
  ci <- stats::confint(ref)
  expect_equal(ours$conf_low, unname(ci[, 1L]), tolerance = 1e-10)
  expect_equal(ours$conf_high, unname(ci[, 2L]), tolerance = 1e-10)
  # fit statistics too
  expect_equal(fit$r_squared, summary(ref)$r.squared, tolerance = 1e-10)
  expect_equal(fit$adj_r_squared, summary(ref)$adj.r.squared, tolerance = 1e-10)
  expect_equal(fit$sigma, summary(ref)$sigma, tolerance = 1e-10)
  expect_equal(fit$aic, stats::AIC(ref), tolerance = 1e-10)
})

test_that("standardizing rescales estimates and leaves the test invariant", {
  cp <- .mo_corpus(n_actors = 60L)
  raw <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     standardize = FALSE)
  std <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     standardize = TRUE)
  d <- as.data.frame(raw, what = "features")
  sds <- vapply(d[, c("pagerank", "betweenness", "closeness")], stats::sd,
                numeric(1L))
  # scale equivariance: beta_std = beta_raw * sd(x); t and p are invariant
  expect_equal(as.data.frame(std)$estimate,
               as.data.frame(raw)$estimate * unname(sds), tolerance = 1e-10)
  expect_equal(as.data.frame(std)$statistic, as.data.frame(raw)$statistic,
               tolerance = 1e-10)
  expect_equal(as.data.frame(std)$p, as.data.frame(raw)$p, tolerance = 1e-10)
})

test_that("binomial and poisson match stats::glm() to 1e-10", {
  cp <- .mo_corpus(n_actors = 60L)
  base <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                      centrality_type = "pagerank", standardize = FALSE)
  d <- as.data.frame(base, what = "features")

  set.seed(3L)
  bin <- stats::setNames(stats::rbinom(nrow(d), 1L, 0.5), d$actor)
  fb <- hon_outcome(cp$hon, outcome = bin, sequences = cp$seqs,
                    centrality_type = "pagerank", family = "binomial",
                    standardize = FALSE)
  d$y_bin <- bin[d$actor]
  rb <- summary(stats::glm(y_bin ~ pagerank, data = d,
                           family = stats::binomial()))$coefficients
  ob <- as.data.frame(fb, what = "terms")
  expect_equal(ob$estimate, unname(rb[, "Estimate"]), tolerance = 1e-10)
  expect_equal(ob$std_error, unname(rb[, "Std. Error"]), tolerance = 1e-10)
  expect_equal(ob$statistic, unname(rb[, "z value"]), tolerance = 1e-10)
  expect_equal(ob$p, unname(rb[, "Pr(>|z|)"]), tolerance = 1e-10)
  expect_identical(fb$distribution, "normal")

  cnt <- stats::setNames(stats::rpois(nrow(d), 3), d$actor)
  fp <- hon_outcome(cp$hon, outcome = cnt, sequences = cp$seqs,
                    centrality_type = "pagerank", family = "poisson",
                    standardize = FALSE)
  d$y_cnt <- cnt[d$actor]
  rp <- summary(stats::glm(y_cnt ~ pagerank, data = d,
                           family = stats::poisson()))$coefficients
  op <- as.data.frame(fp, what = "terms")
  expect_equal(op$estimate, unname(rp[, "Estimate"]), tolerance = 1e-10)
  expect_equal(op$std_error, unname(rp[, "Std. Error"]), tolerance = 1e-10)
  expect_equal(op$p, unname(rp[, "Pr(>|z|)"]), tolerance = 1e-10)
})

# ---- the cluster-robust sandwich ----------------------------------------

test_that("cluster-robust SEs match a hand-computed CR1S sandwich to 1e-10", {
  feats <- .mo_features(n = 120L)
  set.seed(11L)
  cl <- rep(sprintf("g%02d", seq_len(24L)), each = 5L)
  y <- 1 + 0.5 * feats$f1 - 0.3 * feats$f2 +
    rep(stats::rnorm(24L, sd = 0.8), each = 5L) + stats::rnorm(120L, sd = 1)
  outcome <- data.frame(actor = feats$actor, score = y,
                        class = cl, stringsAsFactors = FALSE)

  fit <- hon_outcome(.mo_hon0(), outcome = outcome,
                     features = feats, outcome_col = "score",
                     nested_in = "class", vcov = "CR1", standardize = FALSE)

  # hand-written sandwich, from the definition
  X <- cbind(1, feats$f1, feats$f2)
  bread <- solve(crossprod(X))
  beta <- bread %*% crossprod(X, y)
  e <- as.numeric(y - X %*% beta)
  u <- rowsum(X * e, cl)
  n <- nrow(X); k <- ncol(X); g <- nrow(u)
  adj <- (g / (g - 1)) * ((n - 1) / (n - k))
  v <- adj * (bread %*% crossprod(u) %*% bread)

  ours <- as.data.frame(fit, what = "terms")
  expect_equal(ours$estimate, as.numeric(beta), tolerance = 1e-10)
  expect_equal(ours$std_error, sqrt(diag(v)), tolerance = 1e-10)
  # t(G-1), not t(n-k)
  expect_identical(fit$df, g - 1)
  crit <- stats::qt(0.975, df = g - 1)
  expect_equal(ours$conf_high, as.numeric(beta) + crit * sqrt(diag(v)),
               tolerance = 1e-10)
  expect_equal(ours$p, 2 * stats::pt(-abs(as.numeric(beta) / sqrt(diag(v))),
                                     df = g - 1), tolerance = 1e-10)
  expect_identical(fit$vcov_type, "CR1S")
  expect_identical(fit$cluster_n, g)
})

test_that("cluster-robust SEs match sandwich::vcovCL()", {
  skip_if_not_installed("sandwich")
  feats <- .mo_features(n = 150L)
  set.seed(12L)
  cl <- rep(sprintf("g%02d", seq_len(30L)), each = 5L)
  y <- 2 - 0.4 * feats$f1 + 0.6 * feats$f2 +
    rep(stats::rnorm(30L, sd = 1), each = 5L) + stats::rnorm(150L)
  outcome <- data.frame(actor = feats$actor, score = y, class = cl,
                        stringsAsFactors = FALSE)
  dat <- data.frame(y = y, f1 = feats$f1, f2 = feats$f2, cl = cl,
                    stringsAsFactors = FALSE)

  fit <- hon_outcome(.mo_hon0(), outcome = outcome,
                     features = feats, outcome_col = "score",
                     nested_in = "class", vcov = "CR1", standardize = FALSE)
  ref <- stats::lm(y ~ f1 + f2, data = dat)
  v_ref <- sandwich::vcovCL(ref, cluster = dat$cl)  # HC1 + cadjust: CR1S
  expect_equal(as.data.frame(fit, what = "terms")$std_error,
               unname(sqrt(diag(v_ref))), tolerance = 1e-10)
  expect_equal(unname(fit$vcov), unname(v_ref), tolerance = 1e-10)

  # the GLM path too, against the same estimator at type = "HC1".
  # sandwich::estfun.glm multiplies the working residual by glm's STORED
  # working weights, which IRLS computed one step before the final
  # coefficients; because the convergence test is on the deviance, those
  # weights lag the optimum by about the square root of the tolerance
  # (~1e-6 here). Restarting the fit from its own converged coefficients
  # makes the stored weights the converged ones, and the comparison then
  # holds to machine precision. The estimating function itself is the exact
  # one either way: score_i = x_i (y_i - mu_i) at the canonical link.
  ybin <- as.integer(y > stats::median(y))
  outcome_bin <- data.frame(actor = feats$actor, score = ybin, class = cl,
                            stringsAsFactors = FALSE)
  fb <- hon_outcome(.mo_hon0(),
                    outcome = outcome_bin, features = feats,
                    outcome_col = "score", nested_in = "class", vcov = "CR1",
                    family = "binomial", standardize = FALSE)
  dat$ybin <- ybin
  refb <- stats::glm(ybin ~ f1 + f2, data = dat, family = stats::binomial())

  # exact: the hand-written GLM sandwich, built from glm's own bread and the
  # exact canonical-link score x_i (y_i - mu_i)
  Xb <- stats::model.matrix(refb)
  breadb <- summary(refb)$cov.unscaled
  ub <- rowsum(Xb * (ybin - stats::fitted(refb)), cl)
  gb <- nrow(ub); nb <- nrow(Xb); kb <- ncol(Xb)
  vb <- (gb / (gb - 1)) * ((nb - 1) / (nb - kb)) *
    (breadb %*% crossprod(ub) %*% breadb)
  expect_equal(unname(fb$vcov), unname(vb), tolerance = 1e-10)
  expect_equal(as.data.frame(fb, what = "terms")$std_error,
               unname(sqrt(diag(vb))), tolerance = 1e-10)

  # approximate, and the reason it can only be approximate:
  # sandwich::estfun.glm multiplies the working residual by glm's STORED
  # working weights, which IRLS computed one step before the final
  # coefficients. Because glm's convergence test is on the deviance, those
  # weights lag the optimum by about the square root of the tolerance, so
  # sandwich's score differs from the exact one by ~1e-6 relative. The
  # estimator is the same; only glm's bookkeeping differs.
  vb_ref <- sandwich::vcovCL(refb, cluster = dat$cl, type = "HC1")
  expect_equal(unname(fb$vcov), unname(vb_ref), tolerance = 1e-5)
})

# ---- multiplicity --------------------------------------------------------

test_that("BH is applied across features only, and recorded", {
  cp <- .mo_corpus(n_actors = 60L)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs)
  d <- as.data.frame(fit)
  expect_identical(fit$p_adjust, "BH")
  expect_equal(d$p_adj, stats::p.adjust(d$p, method = "BH"))
  # the intercept is not one of the tested features
  terms_tab <- as.data.frame(fit, what = "terms")
  expect_true(is.na(terms_tab$p_adj[terms_tab$feature == "(Intercept)"]))
  expect_identical(nrow(d), nrow(terms_tab) - 1L)

  none <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                      p_adjust = "none")
  expect_identical(none$p_adjust, "none")
  expect_equal(as.data.frame(none)$p_adj, as.data.frame(none)$p)
  expect_output(summary(none), "multiplicity        : none")
})

# ---- the tidy contract ---------------------------------------------------

test_that("as.data.frame returns the promised tidy tables", {
  cp <- .mo_corpus(n_actors = 40L)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs)

  d <- as.data.frame(fit)
  expect_identical(class(d), "data.frame")
  # the original nine columns, unchanged and first; the per-SD
  # interpretation columns are appended (0.5.1, additive only)
  expect_named(d, c("feature", "estimate", "std_error", "conf_low",
                    "conf_high", "statistic", "p", "p_adj", "n",
                    "feature_sd", "effect_per_sd", "per_sd_low",
                    "per_sd_high", "scale"))
  expect_identical(attr(d, "row.names"), seq_len(nrow(d)))
  expect_identical(d$feature, fit$feature_names)
  expect_true(all(d$n == fit$n))
  expect_true(all(d$conf_low < d$estimate & d$estimate < d$conf_high))

  f <- as.data.frame(fit, what = "features")
  expect_identical(class(f), "data.frame")
  expect_identical(nrow(f), fit$n)
  expect_true(all(c("actor", "outcome", "n_events", "n_unmatched") %in%
                    names(f)))

  expect_identical(nrow(as.data.frame(fit, what = "dropped")), 0L)

  # sort_by and top compose, and top applies last
  by_p <- as.data.frame(fit, sort_by = "p")
  expect_false(is.unsorted(by_p$p))
  expect_identical(as.data.frame(fit, sort_by = "p", top = 2L),
                   `rownames<-`(utils::head(by_p, 2L), NULL))
  by_est <- as.data.frame(fit, sort_by = "estimate")
  expect_identical(by_est$estimate,
                   by_est$estimate[order(-abs(by_est$estimate))])

  # significant = TRUE is the argument that replaces a bracket filter
  sig <- as.data.frame(fit, significant = TRUE, alpha = 1 - 1e-9)
  expect_identical(nrow(sig), nrow(d))
  expect_identical(nrow(as.data.frame(fit, significant = TRUE, alpha = 1e-12)),
                   0L)
})

test_that("print, summary and plot work and say what was done", {
  cp <- .mo_corpus(n_actors = 40L)
  classes <- stats::setNames(rep(sprintf("c%02d", seq_len(20L)), each = 2L),
                             cp$actors)
  expect_warning(
    fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                       nested_in = classes, min_clusters = 30L),
    class = "hypernets_few_clusters")


  expect_output(print(fit), "cluster-robust covariance")
  expect_output(print(fit), "p_adj: BH")
  # the shipped default is the cluster jackknife; summary must NAME the
  # estimator rather than assume one, or a CR3 fit reads as "model-based"
  expect_output(summary(fit), "CR3")
  expect_output(summary(fit), "cluster-robust")
  suppressWarnings(
    cr1 <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                       nested_in = classes, min_clusters = 30L, vcov = "CR1")
  )
  expect_output(summary(cr1), "CR1S")
  expect_output(summary(cr1), "cluster-robust")
  expect_output(summary(fit), "treated as fixed")
  expect_s3_class(plot(fit), "ggplot")
  # summary() returns the tidy table invisibly, print() returns its input
  invisible(utils::capture.output(tidy <- summary(fit)))
  expect_identical(tidy, as.data.frame(fit))
  invisible(utils::capture.output(back <- print(fit)))
  expect_identical(back, fit)
})

# ---- error paths, by class ----------------------------------------------

test_that("bad input raises hypernets_bad_input", {
  cp <- .mo_corpus(n_actors = 20L)
  expect_error(hon_outcome(list(), outcome = cp$outcome, sequences = cp$seqs),
               "must be a net_hon")
  expect_error(
    hon_outcome(cp$hon, outcome = cp$outcome, features = "centrality"),
    class = "hypernets_bad_input")
  expect_error(
    hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                features = "nonsense"),
    class = "hypernets_bad_input")
  expect_error(
    hon_outcome(cp$hon, outcome = unname(cp$outcome), sequences = cp$seqs),
    "named vector keyed by actor")

  dup <- data.frame(actor = c(cp$actors, cp$actors[1L]),
                    score = c(cp$outcome, 1), stringsAsFactors = FALSE)
  expect_error(
    hon_outcome(cp$hon, outcome = dup, sequences = cp$seqs),
    class = "hypernets_bad_input")

  many <- data.frame(actor = cp$actors, score = as.numeric(cp$outcome),
                     other = 1, stringsAsFactors = FALSE)
  expect_error(hon_outcome(cp$hon, outcome = many, sequences = cp$seqs),
               class = "hypernets_bad_input")

  # too few actors for the number of terms
  few <- cp$outcome[seq_len(3L)]
  expect_error(
    suppressWarnings(hon_outcome(cp$hon, outcome = few,
                                 sequences = cp$seqs[seq_len(3L)])),
    class = "hypernets_bad_input")

  # wrong outcome shape for the family
  expect_error(
    hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                family = "binomial"),
    class = "hypernets_bad_input")
  expect_error(
    hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                family = "poisson"),
    class = "hypernets_bad_input")
})

test_that("a rank-deficient design raises hypernets_rank_deficient", {
  feats <- .mo_features(n = 60L)
  feats$f3 <- feats$f1 + feats$f2          # exactly collinear
  outcome <- stats::setNames(stats::rnorm(60L), feats$actor)
  expect_error(
    hon_outcome(.mo_hon0(), outcome = outcome,
                features = feats),
    class = "hypernets_rank_deficient")

  const <- .mo_features(n = 60L)
  const$f2 <- 1
  expect_error(
    hon_outcome(.mo_hon0(), outcome = outcome,
                features = const),
    class = "hypernets_rank_deficient")

  # near-collinear warns rather than failing, with the same class
  near <- .mo_features(n = 60L)
  near$f2 <- near$f1 + stats::rnorm(60L, sd = 1e-4)
  expect_warning(
    hon_outcome(.mo_hon0(), outcome = outcome,
                features = near),
    class = "hypernets_rank_deficient")
})

test_that("cluster problems and dropped actors are surfaced, never silent", {
  feats <- .mo_features(n = 40L)
  outcome <- data.frame(actor = feats$actor,
                        score = stats::rnorm(40L),
                        class = rep(c("a", "b"), each = 20L),
                        stringsAsFactors = FALSE)
  # G = 2 clusters, 3 terms: the sandwich meat is singular
  expect_error(
    hon_outcome(.mo_hon0(), outcome = outcome,
                features = feats, outcome_col = "score", nested_in = "class"),
    class = "hypernets_bad_input")

  outcome$class <- rep(sprintf("g%d", seq_len(8L)), each = 5L)
  expect_warning(
    hon_outcome(.mo_hon0(), outcome = outcome,
                features = feats, outcome_col = "score", nested_in = "class"),
    class = "hypernets_few_clusters")

  # an actor with no outcome is dropped, loudly, and recorded
  short <- outcome[seq_len(35L), c("actor", "score")]
  expect_warning(
    fit <- hon_outcome(.mo_hon0(), outcome = short,
                       features = feats, outcome_col = "score"),
    class = "hypernets_dropped_actors")
  expect_identical(fit$n, 35L)
  dropped <- as.data.frame(fit, what = "dropped")
  expect_identical(nrow(dropped), 5L)
  expect_true(all(dropped$reason == "no outcome value"))

  # a missing cluster label is an error, not a quiet listwise deletion
  outcome$class[1L] <- NA
  expect_error(
    hon_outcome(.mo_hon0(), outcome = outcome,
                features = feats, outcome_col = "score", nested_in = "class"),
    class = "hypernets_bad_input")
})

test_that("a separated binomial fit warns rather than reporting its SEs", {
  feats <- .mo_features(n = 60L)
  feats$f2 <- NULL
  # perfect separation: the sign of f1 determines the outcome, so the MLE
  # does not exist and the standard errors are meaningless
  sep <- stats::setNames(as.integer(feats$f1 > 0), feats$actor)
  expect_warning(
    fit <- withCallingHandlers(
      hon_outcome(.mo_hon0(), outcome = sep, features = feats,
                  family = "binomial"),
      # glm's own diagnosed warning about the same thing, muffled so the
      # assertion below is about OUR classed condition. The prefix matters:
      # our message names the same symptom, and must not muffle itself.
      warning = function(w) {
        if (startsWith(conditionMessage(w), "glm.fit:")) {
          invokeRestart("muffleWarning")
        }
      }),
    class = "hypernets_degenerate_fit")
  expect_gt(as.data.frame(fit)$std_error, 100)
})

test_that("accessor arguments are validated", {
  cp <- .mo_corpus(n_actors = 30L)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     centrality_type = "pagerank")
  expect_error(as.data.frame(fit, what = "nope"), "should be one of")
  expect_error(as.data.frame(fit, sort_by = "nope"), "should be one of")
  expect_error(as.data.frame(fit, what = "features", sort_by = "p"),
               class = "hypernets_bad_input")
  expect_error(as.data.frame(fit, alpha = 0), "`alpha` must be")
  expect_error(as.data.frame(fit, top = 0L), "`top` must be")
  expect_error(plot(fit, alpha = 2), "`alpha` must be")
})

# ---- coverage simulation -------------------------------------------------
#
# The point of this file. A verb that reports a 95% interval must produce
# intervals that cover the truth about 95% of the time. The design matrix is
# held fixed (the regression conditions on it) and only the outcome is
# redrawn, which is exactly the claim the intervals make.

.mo_coverage <- function(n_rep, generate, fit_args, truth = 0.5,
                         feature = "f1", seed = 99L) {
  had <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
  if (had) old <- get(".Random.seed", envir = globalenv())
  on.exit(if (had) assign(".Random.seed", old, envir = globalenv()) else
    rm(".Random.seed", envir = globalenv()), add = TRUE, after = FALSE)
  set.seed(seed)
  hon0 <- .mo_hon0()
  hits <- vapply(seq_len(n_rep), function(i) {
    outcome <- generate()
    fit <- do.call(hon_outcome, c(list(hon0, outcome = outcome), fit_args))
    d <- as.data.frame(fit)
    row <- d$feature == feature
    c(covered = d$conf_low[row] <= truth && truth <= d$conf_high[row],
      reject = d$p[row] < 0.05)
  }, logical(2L))
  list(coverage = mean(hits["covered", ]),
       type1 = mean(hits["reject", ]),
       n_rep = n_rep)
}

test_that("model-based intervals cover a known effect at about 95%", {
  feats <- .mo_features(n = 150L, seed = 5L)
  b1 <- 0.5
  gen <- function() {
    stats::setNames(2 + b1 * feats$f1 - 0.3 * feats$f2 +
                      stats::rnorm(nrow(feats), sd = 1.5), feats$actor)
  }
  res <- .mo_coverage(2000L, gen,
                      list(features = feats, standardize = FALSE),
                      truth = b1)
  # binomial SE at 2000 draws is 0.0049, so +/- 3 SE is about 0.015
  expect_gt(res$coverage, 0.93)
  expect_lt(res$coverage, 0.97)
  cat(sprintf("\n  model-based coverage: %.4f over %d replicates\n",
              res$coverage, res$n_rep))
})

test_that("cluster-robust intervals cover a known effect under clustering", {
  # The regressor is itself cluster-correlated -- students in the same class
  # share most of f1 -- which is the case where ignoring the clustering
  # really does break the intervals. (When the regressor varies freely
  # within cluster, the model-based interval is close to right and the
  # cluster-robust one buys little; that is a fact about the estimator, not
  # a defect, and it is why this test is built the way it is.)
  n <- 200L
  n_g <- 40L
  per <- n / n_g
  set.seed(6L)
  cl <- rep(sprintf("g%02d", seq_len(n_g)), each = per)
  feats <- data.frame(
    actor = sprintf("a%03d", seq_len(n)),
    f1 = rep(stats::rnorm(n_g), each = per) + stats::rnorm(n, sd = 0.3),
    f2 = stats::rnorm(n),
    stringsAsFactors = FALSE)
  b1 <- 0.5
  gen <- function() {
    y <- 2 + b1 * feats$f1 - 0.3 * feats$f2 +
      rep(stats::rnorm(n_g, sd = 1.2), each = per) + stats::rnorm(n, sd = 1)
    data.frame(actor = feats$actor, score = y, class = cl,
               stringsAsFactors = FALSE)
  }
  res <- .mo_coverage(2000L, gen,
                      list(features = feats, outcome_col = "score",
                           nested_in = "class", standardize = FALSE),
                      truth = b1)
  cat(sprintf("  cluster-robust coverage: %.4f over %d reps (%d clusters)\n",
              res$coverage, res$n_rep, n_g))
  expect_gt(res$coverage, 0.92)
  expect_lt(res$coverage, 0.98)

  # the model-based interval, which ignores the clustering, under-covers
  # badly here: the reason nested_in exists at all
  naive <- .mo_coverage(500L, gen,
                        list(features = feats, outcome_col = "score",
                             standardize = FALSE),
                        truth = b1, seed = 101L)
  cat(sprintf("  ignoring the clustering: %.4f over %d replicates\n",
              naive$coverage, naive$n_rep))
  expect_lt(naive$coverage, 0.90)
})

# ---- CR3, the shipped default ---------------------------------------------

test_that("CR3 matches sandwich::vcovCL(type = 'HC3') to 1e-10", {
  skip_if_not_installed("sandwich")
  set.seed(7)
  g <- 12L
  cl <- rep(seq_len(g), each = 8L)
  n <- length(cl)
  X <- cbind(1, x1 = stats::rnorm(n), x2 = stats::rnorm(g)[cl])
  y <- as.vector(X %*% c(1, 0.4, -0.3) + stats::rnorm(g)[cl] + stats::rnorm(n))
  fit <- hypernets:::.hoo_fit(X, y, family = "gaussian",
                              cluster = as.character(cl), vcov = "CR3")
  dat <- data.frame(y = y, x1 = X[, 2L], x2 = X[, 3L], cl = factor(cl))
  ref <- sandwich::vcovCL(stats::lm(y ~ x1 + x2, data = dat),
                          cluster = dat$cl, type = "HC3")
  expect_equal(unname(fit$vcov), unname(ref), tolerance = 1e-10)
  expect_identical(fit$vcov_type, "CR3")
  # CR3 carries no G/(G-1) factor -- the leverage term IS the correction, so
  # it must differ from CR1 rather than coincide with it
  cr1 <- hypernets:::.hoo_fit(X, y, family = "gaussian",
                              cluster = as.character(cl), vcov = "CR1")
  expect_false(isTRUE(all.equal(diag(fit$vcov), diag(cr1$vcov))))
  # and it is the more conservative of the two here
  expect_true(all(diag(fit$vcov) > diag(cr1$vcov)))
})

test_that("CR3 refuses a cluster it cannot jackknife", {
  # a singleton cluster carrying its own indicator has leverage 1, so
  # (I - H_gg) is singular and the jackknife is undefined; CR1 still works
  set.seed(11)
  cl <- c(rep("a", 8L), rep("b", 8L), "c")
  n <- length(cl)
  X <- cbind(1, x = stats::rnorm(n), only_c = as.numeric(cl == "c"))
  y <- as.vector(X %*% c(1, 0.5, 2) + stats::rnorm(n))
  expect_error(
    hypernets:::.hoo_fit(X, y, family = "gaussian", cluster = cl,
                         vcov = "CR3"),
    class = "hypernets_rank_deficient"
  )
  expect_true(is.matrix(
    hypernets:::.hoo_fit(X, y, family = "gaussian", cluster = cl,
                         vcov = "CR1")$vcov
  ))
})

test_that("the covariance type is named, never assumed, in every method", {
  # summary/print/plot used to branch on vcov_type == "CR1S", so a CR3 fit
  # reported itself as "model-based" -- the opposite of the truth
  cp <- .mo_corpus(n_actors = 40L)
  classes <- rep(paste0("c", seq_len(8L)), length.out = length(cp$actors))
  names(classes) <- cp$actors
  for (v in c("CR3", "CR1")) {
    suppressWarnings(
      fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                         nested_in = classes, min_clusters = 30L, vcov = v)
    )
    expect_output(print(fit), "cluster-robust covariance")
    expect_output(summary(fit), "cluster-robust")
    expect_match(fit$vcov_type, if (v == "CR3") "CR3" else "CR1S")
    expect_s3_class(plot(fit), "ggplot")
  }
})

# ---- CR3 for GLMs: the leverage carries the working weights ---------------

test_that("CR3 matches sandwich::vcovCL(type = 'HC3') for binomial and poisson", {
  skip_if_not_installed("sandwich")
  set.seed(23)
  g <- 15L
  cl <- rep(seq_len(g), each = 6L)
  n <- length(cl)
  X <- cbind(1, x1 = stats::rnorm(n), x2 = stats::rnorm(g)[cl])
  eta <- as.vector(X %*% c(-0.2, 0.8, -0.5)) + stats::rnorm(g, sd = 0.5)[cl]
  dat <- data.frame(x1 = X[, 2L], x2 = X[, 3L], cl = factor(cl))
  dat$yb <- stats::rbinom(n, 1L, stats::plogis(eta))
  dat$yp <- stats::rpois(n, exp(eta))
  cases <- list(binomial = "yb", poisson = "yp")
  invisible(Map(function(family, resp) {
    fit <- hypernets:::.hoo_fit(X, dat[[resp]], family = family,
                                cluster = as.character(cl), vcov = "CR3")
    ref_fit <- stats::glm(stats::reformulate(c("x1", "x2"), resp),
                          data = dat, family = family,
                          control = stats::glm.control(epsilon = 1e-12,
                                                       maxit = 100L))
    ref <- sandwich::vcovCL(ref_fit, cluster = dat$cl, type = "HC3")
    # before the fix binomial was off by a relative 250x, poisson by 10x
    expect_equal(unname(fit$vcov), unname(ref), tolerance = 1e-6,
                 label = paste("CR3", family))
  }, names(cases), cases))
})

test_that("hon_outcome() refuses an unknown centrality_type", {
  seqs <- data.frame(actor = rep(c("a", "b", "c", "d"), each = 6L),
                     time = rep(seq_len(6L), 4L),
                     action = rep(c("x", "y", "z"), 8L))
  hon <- build_hon(seqs, action = "action", actor = "actor", time = "time",
                   max_order = 2L)
  scores <- data.frame(actor = c("a", "b", "c", "d"), score = c(1, 2, 3, 4))
  expect_error(
    hon_outcome(hon, outcome = scores, sequences = seqs,
                centrality_type = c("pagerank", "in_strength")),
    class = "hypernets_bad_input"
  )
})

# ---- interpretation layer: per-SD effects and percentile contrasts ----------

test_that("the per-SD effect is the coefficient of a refit on z-scored features", {
  cp <- .mo_corpus(n_actors = 60L)
  # raw scale, so the rescaling is actually exercised (factor sd / 1)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     centrality_type = c("pagerank", "closeness"),
                     standardize = FALSE)
  d <- as.data.frame(fit, what = "features")
  z <- data.frame(y = d$outcome, pagerank = as.numeric(scale(d$pagerank)),
                  closeness = as.numeric(scale(d$closeness)))
  ref <- stats::lm(y ~ pagerank + closeness, data = z)
  ours <- as.data.frame(fit)
  expect_equal(ours$effect_per_sd, unname(stats::coef(ref)[-1L]),
               tolerance = 1e-10)
  expect_equal(ours$per_sd_low, unname(stats::confint(ref)[-1L, 1L]),
               tolerance = 1e-10)
  expect_equal(ours$per_sd_high, unname(stats::confint(ref)[-1L, 2L]),
               tolerance = 1e-10)
  expect_equal(ours$feature_sd, c(stats::sd(d$pagerank),
                                  stats::sd(d$closeness)), tolerance = 1e-12)
  expect_identical(unique(ours$scale), "outcome units")

  # the same number whichever scale the features entered the model on
  std <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     centrality_type = c("pagerank", "closeness"))
  expect_equal(as.data.frame(std)$effect_per_sd, ours$effect_per_sd,
               tolerance = 1e-10)
  # standardized: the per-SD effect IS the estimate (invariant)
  expect_equal(as.data.frame(std)$effect_per_sd,
               as.data.frame(std)$estimate, tolerance = 1e-12)
  # the intercept carries no per-SD effect
  terms_tab <- as.data.frame(std, what = "terms")
  expect_true(is.na(terms_tab$effect_per_sd[terms_tab$feature ==
                                              "(Intercept)"]))
})

test_that("per-SD rescaling is exact under a cluster-robust covariance", {
  feats <- .mo_features(n = 120L)
  feats$f2 <- feats$f2 * 250 + 3   # a feature on a very different scale
  set.seed(11)
  y <- stats::setNames(1 + 0.5 * feats$f1 - 0.004 * feats$f2 +
                         stats::rnorm(120L), feats$actor)
  cl <- stats::setNames(rep(sprintf("g%02d", seq_len(30L)), each = 4L),
                        feats$actor)
  raw <- hon_outcome(.mo_hon0(), outcome = y, features = feats,
                     nested_in = cl, standardize = FALSE)
  std <- hon_outcome(.mo_hon0(), outcome = y, features = feats,
                     nested_in = cl, standardize = TRUE)
  cols <- c("effect_per_sd", "per_sd_low", "per_sd_high")
  expect_equal(as.data.frame(raw)[cols], as.data.frame(std)[cols],
               tolerance = 1e-10)
})

test_that("the gaussian contrast is beta * (q75 - q25) with a delta-method SE", {
  cp <- .mo_corpus(n_actors = 60L)
  classes <- stats::setNames(rep(sprintf("c%02d", seq_len(20L)), each = 3L),
                             cp$actors)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     centrality_type = c("pagerank", "betweenness"),
                     nested_in = classes)
  con <- as.data.frame(fit, what = "contrasts")
  co <- as.data.frame(fit)
  feats <- as.data.frame(fit, what = "features")
  q <- vapply(c("pagerank", "betweenness"), function(j)
    unname(diff(stats::quantile(feats[[j]], c(0.25, 0.75)))), numeric(1L))
  # estimate is per design unit = per SD here, so per raw unit is est / sd
  expect_equal(con$difference, co$estimate / co$feature_sd * unname(q),
               tolerance = 1e-10)
  expect_equal(con$difference, con$predicted_high - con$predicted_low,
               tolerance = 1e-12)
  expect_equal(con$low_value, unname(vapply(c("pagerank", "betweenness"),
    function(j) stats::quantile(feats[[j]], 0.25), numeric(1L))))

  # hand-computed delta method: g = c_high - c_low, se = sqrt(g' V g);
  # for a linear predictor g is (q75 - q25) / s_j on feature j alone
  V <- fit$vcov
  g <- stats::setNames(numeric(ncol(V)), colnames(V))
  g[["betweenness"]] <- q[["betweenness"]] / fit$scale[["betweenness"]]
  se_hand <- sqrt(drop(g %*% V %*% g))
  expect_equal(con$std_error[2L], se_hand, tolerance = 1e-12)
  crit <- stats::qt(0.975, df = fit$cluster_n - 1L)
  expect_equal(con$conf_high[2L] - con$conf_low[2L], 2 * crit * se_hand,
               tolerance = 1e-10)

  # the predicted level at the low end: intercept + others at their means
  b <- stats::setNames(as.data.frame(fit, what = "terms")$estimate,
                       colnames(V))
  row_low <- c(1, 0, (con$low_value[2L] - fit$centre[["betweenness"]]) /
                 fit$scale[["betweenness"]])
  expect_equal(con$predicted_low[2L], sum(row_low * b), tolerance = 1e-12)
  expect_equal(con$predicted_low_ci_high[2L] - con$predicted_low[2L],
               crit * sqrt(drop(row_low %*% V %*% row_low)),
               tolerance = 1e-10)
  expect_true(all(is.na(con$ratio)))
})

test_that("binomial and poisson effects are ratios per SD", {
  cp <- .mo_corpus(n_actors = 60L)
  set.seed(3)
  yb <- stats::setNames(stats::rbinom(60L, 1L,
                                      stats::plogis(4 * (cp$mix[, 1L] - 0.5))),
                        cp$actors)
  fb <- hon_outcome(cp$hon, outcome = yb, sequences = cp$seqs,
                    family = "binomial", centrality_type = "pagerank",
                    standardize = FALSE)
  ob <- as.data.frame(fb)
  expect_equal(ob$effect_per_sd, exp(ob$estimate * ob$feature_sd),
               tolerance = 1e-12)
  expect_equal(ob$per_sd_low, exp(ob$conf_low * ob$feature_sd),
               tolerance = 1e-12)
  expect_identical(ob$scale, "odds ratio")
  # the refit on a z-scored feature gives the same odds ratio
  d <- as.data.frame(fb, what = "features")
  ref <- stats::glm(d$outcome ~ as.numeric(scale(d$pagerank)),
                    family = stats::binomial(),
                    control = stats::glm.control(epsilon = 1e-12))
  expect_equal(ob$effect_per_sd, unname(exp(stats::coef(ref)[2L])),
               tolerance = 1e-8)

  con <- as.data.frame(fb, what = "contrasts")
  expect_true(all(con$predicted_low > 0 & con$predicted_low < 1))
  expect_true(con$predicted_low_ci_low >= 0 && con$predicted_low_ci_high <= 1)
  expect_equal(con$ratio, exp(ob$estimate * (con$high_value - con$low_value)),
               tolerance = 1e-10)
  expect_identical(con$scale, "probability")
  # delta method on the probability scale, by hand:
  # g = p_h (1 - p_h) c_h - p_l (1 - p_l) c_l
  b <- as.data.frame(fb, what = "terms")$estimate
  c_l <- c(1, con$low_value)
  c_h <- c(1, con$high_value)
  p_l <- stats::plogis(sum(c_l * b))
  p_h <- stats::plogis(sum(c_h * b))
  g <- p_h * (1 - p_h) * c_h - p_l * (1 - p_l) * c_l
  expect_equal(con$std_error, sqrt(drop(g %*% fb$vcov %*% g)),
               tolerance = 1e-10)
  expect_equal(con$difference, p_h - p_l, tolerance = 1e-12)

  yp <- stats::setNames(stats::rpois(60L, exp(1 + cp$mix[, 1L])), cp$actors)
  fp <- hon_outcome(cp$hon, outcome = yp, sequences = cp$seqs,
                    family = "poisson", centrality_type = "pagerank")
  op <- as.data.frame(fp)
  expect_identical(op$scale, "rate ratio")
  expect_equal(op$effect_per_sd, exp(op$estimate), tolerance = 1e-12)
})

test_that("plot() draws every view and refuses an unknown one by class", {
  cp <- .mo_corpus(n_actors = 40L)
  classes <- stats::setNames(rep(sprintf("c%02d", seq_len(20L)), each = 2L),
                             cp$actors)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     nested_in = classes)
  views <- c("per_sd", "contrasts", "raw")
  plots <- lapply(views, function(w) plot(fit, what = w))
  invisible(lapply(plots, function(p) expect_s3_class(p, "ggplot")))
  # the default is the per-SD view, and it says how it was computed
  default <- plot(fit)
  expect_identical(default$labels$x, plots[[1L]]$labels$x)
  expect_match(default$labels$caption, "CR3 over 20 clusters \\(classes\\)")
  expect_match(default$labels$subtitle, "40 actors in 20 clusters")
  expect_match(default$labels$subtitle, "R-squared")
  expect_match(default$labels$x, "cp\\$outcome|outcome")
  expect_error(plot(fit, what = "nope"), class = "hypernets_bad_input")

  set.seed(5)
  yb <- stats::setNames(stats::rbinom(40L, 1L, 0.5), cp$actors)
  fb <- hon_outcome(cp$hon, outcome = yb, sequences = cp$seqs,
                    family = "binomial", centrality_type = "pagerank")
  pb <- plot(fb)
  expect_match(pb$labels$x, "Odds ratio")
  expect_match(pb$labels$subtitle, "deviance explained")
  expect_s3_class(plot(fb, what = "contrasts"), "ggplot")
})

test_that("contrasts refuse sort_by / significant like the other side tables", {
  cp <- .mo_corpus(n_actors = 30L)
  fit <- hon_outcome(cp$hon, outcome = cp$outcome, sequences = cp$seqs,
                     centrality_type = "pagerank")
  expect_error(as.data.frame(fit, what = "contrasts", sort_by = "p"),
               class = "hypernets_bad_input")
  expect_identical(nrow(as.data.frame(fit, what = "contrasts")), 1L)
})
