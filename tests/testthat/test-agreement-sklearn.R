testthat::skip_on_cran()

# Calibration of hg_agreement() against scikit-learn 1.8.0, the reference
# implementation the legal-hypergraphs literature scores partitions with.
# The values below were produced by sklearn.metrics.adjusted_rand_score,
# adjusted_mutual_info_score and normalized_mutual_info_score on exactly the
# partitions .agree_cases() builds; the generating script is
# tmp/agree_run.R. sklearn's AMI and NMI default to average_method =
# "arithmetic", which is the normalisation .thg_ami() / .thg_nmi() implement --
# a switch to "max" or "geometric" would change every value here, so this
# fixture is what pins the convention.
#
# The partitions are built by modular arithmetic rather than sampled, so the
# fixture cannot drift with an RNG change.

.agree_cases <- function() {
  i <- function(n) seq_len(n)
  near <- rep(1:10, each = 20)
  near[c(3, 17, 44, 61, 80, 97, 111, 128, 140, 155, 166, 177, 188, 195, 200)] <-
    c(7, 2, 9, 1, 4, 10, 3, 6, 8, 5, 2, 9, 1, 7, 4)
  list(
    identical      = list(a = rep(1:5, each = 8), b = rep(1:5, each = 8)),
    relabelled     = list(a = rep(1:5, each = 8),
                          b = rep(c(3, 1, 5, 2, 4), each = 8)),
    independent    = list(a = i(200) %% 6, b = (i(200) %/% 5) %% 6),
    unequal_k      = list(a = i(150) %% 3, b = i(150) %% 17),
    many_singleton = list(a = c(rep(1, 40), 42:101), b = i(100) %% 20),
    one_cluster    = list(a = rep(1, 60), b = i(60) %% 4),
    near_identical = list(a = rep(1:10, each = 20), b = near),
    skewed         = list(a = c(rep(1, 300), rep(2, 5), rep(3, 2)),
                          b = (i(307) %/% 11) %% 8),
    large          = list(a = i(3618) %% 50, b = (i(3618) %/% 7) %% 50)
  )
}

.sklearn_reference <- list(
  identical      = c(ari = 1, ami = 1, nmi = 1),
  relabelled     = c(ari = 1, ami = 1, nmi = 1),
  independent    = c(ari = 0.0156624767712796, ami = 0.0679455245479407,
                     nmi = 0.102101678811944),
  unequal_k      = c(ari = -0.0239246627640621, ami = -0.0635306763941731,
                     nmi = 0.00133643825544592),
  many_singleton = c(ari = -0.0251156642432254, ami = -0.0576700607453535,
                     nmi = 0.586885374781399),
  one_cluster    = c(ari = 0, ami = 0, nmi = 0),
  near_identical = c(ari = 0.858156939957116, ami = 0.875424997296611,
                     nmi = 0.888257949991517),
  skewed         = c(ari = -0.000307353301151132, ami = 0.0218473479910435,
                     nmi = 0.0419273306577388),
  large          = c(ari = 0.113499724376495, ami = 0.447878249992796,
                     nmi = 0.502646567943462)
)

.as_partition <- function(labels) {
  data.frame(node = as.character(seq_along(labels)),
             label = as.character(labels), stringsAsFactors = FALSE)
}

test_that("hg_agreement matches scikit-learn on nine partition shapes", {
  cases <- .agree_cases()
  expect_setequal(names(cases), names(.sklearn_reference))
  worst <- vapply(names(cases), function(nm) {
    got <- hg_agreement(.as_partition(cases[[nm]]$a),
                        .as_partition(cases[[nm]]$b),
                        method = c("ari", "ami", "nmi"))
    reference <- .sklearn_reference[[nm]]
    max(abs(c(got$ari - reference[["ari"]], got$ami - reference[["ami"]],
              got$nmi - reference[["nmi"]])))
  }, numeric(1L))
  # measured worst deviation over the 27 comparisons was 1.6e-13
  expect_true(all(worst < 1e-10))
})

test_that("the agreement measures obey the invariants that define them", {
  cases <- .agree_cases()
  # relabelling a partition cannot change a label-permutation-invariant score
  # a permutation of all ten labels: anything shorter would index past the end
  # and silently rewrite the partition with NAs
  relabel <- c(7, 2, 9, 1, 4, 10, 3, 6, 8, 5)
  expect_setequal(relabel, sort(unique(cases$near_identical$a)))
  shuffled <- .as_partition(relabel[cases$near_identical$a])
  straight <- hg_agreement(.as_partition(cases$near_identical$a),
                           .as_partition(cases$near_identical$b),
                           method = c("ari", "ami", "nmi"))
  permuted <- hg_agreement(shuffled, .as_partition(cases$near_identical$b),
                           method = c("ari", "ami", "nmi"))
  expect_equal(permuted$ari, straight$ari)
  expect_equal(permuted$ami, straight$ami)
  expect_equal(permuted$nmi, straight$nmi)
  # symmetric in its arguments
  swapped <- hg_agreement(.as_partition(cases$near_identical$b),
                          .as_partition(cases$near_identical$a),
                          method = c("ari", "ami", "nmi"))
  expect_equal(swapped$ami, straight$ami)
  expect_equal(swapped$ari, straight$ari)
  # a partition against itself is exactly 1, and the adjusted measures sit at
  # 0 in expectation for unrelated partitions, so they may go negative
  self <- hg_agreement(.as_partition(cases$large$a), .as_partition(cases$large$a),
                       method = c("ari", "ami", "nmi"))
  expect_equal(self$ari, 1)
  expect_equal(self$ami, 1)
  expect_equal(self$nmi, 1)
  expect_lt(.sklearn_reference$unequal_k[["ari"]], 0)
})

test_that("the community medoid is the run with the largest summed AMI", {
  skip_if_not_installed("igraph")
  set.seed(11)
  dat <- do.call(rbind, lapply(1:40, function(k) {
    block <- ((k - 1L) %/% 10L) + 1L
    data.frame(member = paste0("n", (block - 1L) * 15L + sample(1:15, 4)),
               event = paste0("e", k), stringsAsFactors = FALSE)
  }))
  fit <- hg_communities(group_hypergraph(dat, "member", "event"),
                        n_runs = 8, trials = 20)
  ami <- fit$similarity$ami
  expect_true(isTRUE(all.equal(ami, t(ami))))
  expect_true(all(diag(ami) == 1))
  # the diagonal is a constant 1 in every row, so including it cannot move the
  # argmax -- the medoid is the run most like all the others either way
  expect_identical(fit$medoid_run, unname(which.max(rowSums(ami))))
  expect_identical(fit$medoid_run, unname(which.max(rowSums(ami) - diag(ami))))
  expect_identical(
    fit$medoid$community,
    fit$partitions$community[fit$partitions$run == fit$medoid_run]
  )
})
