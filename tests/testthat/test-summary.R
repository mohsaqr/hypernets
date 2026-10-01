testthat::skip_on_cran()

# summary() on a result returns every hg_get() table as a named list of data
# frames, adds overall figures, and leaves the result unchanged.

.sm_results <- function() {
  set.seed(1)
  seqs <- ring_sequences[1:40]
  hon_fit <- hon(seqs, max_order = 2L)
  by_group <- hon(seqs, group = rep(c("g1", "g2"), 20L), max_order = 2L)
  w <- matrix(stats::runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
  w <- (w + t(w)) / 2
  diag(w) <- 0
  homology <- hg_homology(w, n_steps = 6L, max_dim = 2L)
  list(
    hon = hon_fit,
    honem = honem(hon_fit, dim = 2L),
    mogen = mogen(seqs, max_order = 2L),
    markov_order = markov_order(seqs, max_order = 2L, n_perm = 19L, seed = 1L),
    memory = memory(seqs, order = 2L),
    stability = hg_markov_stability(seqs),
    hypa = hypa(seqs, order = 2L),
    bootstrap = hg_bootstrap(seqs, n_boot = 5L, max_order = 2L, seed = 1L),
    compare = hg_compare(by_group, n_perm = 19L, seed = 1L),
    hon_communities = hg_communities(hon_fit, trials = 2L),
    group = by_group,
    complex = simplicial(w, threshold = 0.5),
    homology = homology,
    landscape = hg_landscape(homology, dimension = 0L),
    qanalysis = hg_qanalysis(simplicial(w, threshold = 0.5))
  )
}

test_that("every summary element is a data frame and equals its hg_get() table", {
  results <- .sm_results()
  lapply(names(results), \(name) {
    fit <- results[[name]]
    s <- summary(fit)
    expect_s3_class(s, "hypernets_summary")
    expect_true(all(vapply(s, is.data.frame, logical(1L))), info = name)
    tables <- intersect(names(s), .ho_what_values(fit))
    expect_gt(length(tables), 0L)
    lapply(tables, \(what) {
      expect_identical(s[[what]], hg_get(fit, what = what),
                       info = paste(name, what))
    })
    # the summary's first table is its default reader
    expect_identical(hg_get(s), s[[1L]], info = name)
  })
})

test_that("summary() leaves the result unchanged and adds overall figures", {
  results <- .sm_results()
  before <- results$memory
  s <- summary(results$memory)
  expect_identical(results$memory, before)
  expect_named(s$overall, c("n_contexts", "n_flips", "kl_weighted",
                            "entropy_drop_weighted"))
  expect_identical(s$overall$n_contexts, nrow(s$contexts))
  expect_named(summary(results$markov_order)$overall,
               c("optimal_order", "aic_order", "bic_order", "n_sequences",
                 "n_permutations"))
})

test_that("a summary prints its tables and names them", {
  s <- summary(.sm_results()$memory)
  out <- capture.output(print(s, n = 2L))
  expect_match(out[1L], "^Memory of order 2")
  expect_identical(out[2L], "Tables: contexts, overall")
  expect_error(hg_get(s, what = "nope"))
})

test_that("a validated complex summarises its tests", {
  # 99 shuffles are fewer than the tests need, which is reported by class
  expect_warning(
    sc <- simplicial(ring_sequences, type = "window", window = 3L,
                     validate = TRUE, n_null = 99L, seed = 1L),
    class = "hypernets_low_resolution")
  s <- summary(sc)
  expect_identical(s$validation, hg_get(sc, what = "validation"))
  expect_true("simplices" %in% names(s))
})
