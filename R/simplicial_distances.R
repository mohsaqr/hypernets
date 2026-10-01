# ---- Wasserstein distance between persistence diagrams ----
#
# hg_wasserstein() is hypernets' own verb on the persistent_homology
# objects hg_homology() returns (hg_homology(), hg_bottleneck() and
# hg_landscape() wrap the Nestimate estimators). It needs
# two Nestimate internals, kept as private copies because `:::` is not
# allowed: .ph_as_diagram() and .is_essential(). Their identity with
# Nestimate's is asserted in
# local_testing_and_equivalence/test-identity-nestimate-simplicial.R.

#' @noRd
.ph_as_diagram <- function(x) {
  if (inherits(x, "persistent_homology")) return(x$persistence)
  if (is.data.frame(x)) {
    need <- c("dimension", "birth", "death")
    if (!all(need %in% names(x))) {
      stop("data.frame must have columns: ",
           paste(need, collapse = ", "), call. = FALSE)
    }
    return(x[, c("dimension", "birth", "death")])
  }
  stop("d1/d2 must be persistent_homology or data.frame.", call. = FALSE)
}

#' @noRd
.is_essential <- function(p) {
  is.infinite(p$death) | p$death == 0
}

#' Wasserstein Distance Between Persistence Diagrams
#'
#' Computes the finite-order Wasserstein distance between persistence
#' diagrams. For order \eqn{q}, it minimizes the sum of powered matching
#' costs over bijections between diagram points and copies of the diagonal,
#' then takes the \eqn{q}-th root. Point-to-point costs use an
#' \eqn{L_p} ground metric; `internal_p = Inf` gives the usual
#' \eqn{L_\infty} convention used by GUDHI.
#'
#' Finite points may match the diagonal. Essential classes (`death = Inf`
#' in Vietoris--Rips mode or `death = 0` in clique mode) are matched only to
#' essential classes. A dimension whose essential counts differ has distance
#' `Inf`. The finite assignment is solved exactly with a native Hungarian
#' algorithm, so no optional optimization package is required.
#'
#' @param d1,d2 `persistent_homology` objects from
#'   [hg_homology()], or data.frames with columns
#'   `dimension`, `birth`, `death`.
#' @param dimension Integer vector of dimensions to compare. `NULL`
#'   (default) compares all dimensions appearing in either diagram.
#' @param order Finite Wasserstein order, a number greater than or equal to
#'   one. The common choices are `1` and `2`.
#' @param internal_p Ground-metric order, a number greater than or equal to
#'   one or `Inf` (the default).
#'
#' @return Named numeric vector, one value per requested dimension. Names are
#'   `"dim_<k>"`.
#'
#' @references
#' Kerber, M., Morozov, D., & Nigmetov, A. (2017). Geometry helps to compare
#' persistence diagrams. *Journal of Experimental Algorithmics*, 22, 1--20.
#' \doi{10.1145/3064175}
#'
#' @examples
#' d1 <- data.frame(dimension = 0L, birth = 0, death = 2)
#' d2 <- data.frame(dimension = 0L, birth = 0, death = 3)
#' hg_wasserstein(d1, d2)
#'
#' @export
hg_wasserstein <- function(d1, d2, dimension = NULL, order = 1,
                           internal_p = Inf) {
  df1 <- .ph_as_diagram(d1)
  df2 <- .ph_as_diagram(d2)
  if (length(order) != 1L || !is.finite(order) || order < 1) {
    stop("`order` must be one finite number >= 1.", call. = FALSE)
  }
  if (length(internal_p) != 1L || is.na(internal_p) || internal_p < 1) {
    stop("`internal_p` must be one number >= 1 or Inf.", call. = FALSE)
  }
  dims <- if (is.null(dimension)) {
    sort(unique(c(df1$dimension, df2$dimension)))
  } else {
    if (!is.numeric(dimension) || anyNA(dimension) ||
        any(dimension != as.integer(dimension))) {
      stop("`dimension` must contain integers.", call. = FALSE)
    }
    as.integer(dimension)
  }
  if (length(dims) == 0L) {
    return(stats::setNames(numeric(0), character(0)))
  }
  out <- vapply(dims, function(k) {
    .wasserstein_one_dim(
      df1[df1$dimension == k, , drop = FALSE],
      df2[df2$dimension == k, , drop = FALSE],
      order = order, internal_p = internal_p
    )
  }, numeric(1))
  stats::setNames(out, paste0("dim_", dims))
}

#' @noRd
.wasserstein_one_dim <- function(p1, p2, order, internal_p) {
  ess1 <- .is_essential(p1)
  ess2 <- .is_essential(p2)
  ess_b1 <- sort(p1$birth[ess1])
  ess_b2 <- sort(p2$birth[ess2])
  if (length(ess_b1) != length(ess_b2)) return(Inf)

  essential_power <- if (length(ess_b1)) {
    sum(abs(ess_b1 - ess_b2) ^ order)
  } else {
    0
  }
  fin1 <- as.matrix(p1[!ess1, c("birth", "death"), drop = FALSE])
  fin2 <- as.matrix(p2[!ess2, c("birth", "death"), drop = FALSE])
  finite_power <- .wasserstein_finite_power(fin1, fin2, order, internal_p)
  (essential_power + finite_power) ^ (1 / order)
}

#' @noRd
.wasserstein_finite_power <- function(p1, p2, order, internal_p) {
  n1 <- nrow(p1)
  n2 <- nrow(p2)
  if (n1 + n2 == 0L) return(0)

  diagonal_cost <- function(p) {
    persistence <- abs(p[, "death"] - p[, "birth"])
    if (is.infinite(internal_p)) {
      persistence / 2
    } else {
      persistence / (2 ^ (1 - 1 / internal_p))
    }
  }
  d1_diag <- if (n1) diagonal_cost(p1) else numeric(0)
  d2_diag <- if (n2) diagonal_cost(p2) else numeric(0)

  pair <- matrix(numeric(0), n1, n2)
  if (n1 && n2) {
    db <- abs(outer(p1[, "birth"], p2[, "birth"], "-"))
    dd <- abs(outer(p1[, "death"], p2[, "death"], "-"))
    pair <- if (is.infinite(internal_p)) {
      pmax(db, dd)
    } else {
      (db ^ internal_p + dd ^ internal_p) ^ (1 / internal_p)
    }
  }

  # Augment both diagrams with enough diagonal copies to make the matching
  # square. Diagonal copies are interchangeable: each top-right row has its
  # point's projection cost, each bottom-left column has its point's cost.
  size <- n1 + n2
  cost <- matrix(0, size, size)
  if (n1 && n2) cost[seq_len(n1), seq_len(n2)] <- pair ^ order
  if (n1) {
    cost[seq_len(n1), n2 + seq_len(n1)] <-
      matrix(d1_diag ^ order, nrow = n1, ncol = n1)
  }
  if (n2) {
    cost[n1 + seq_len(n2), seq_len(n2)] <-
      matrix(d2_diag ^ order, nrow = n2, ncol = n2, byrow = TRUE)
  }
  .hungarian_min_cost(cost)
}

# Exact minimum-cost perfect assignment for a finite square cost matrix.
# This potential-based Hungarian implementation is O(n^3) and deliberately
# internal so diagram distances do not acquire a dependency on `clue`.
#' @noRd
.hungarian_min_cost <- function(cost) {
  n <- nrow(cost)
  if (n == 0L) return(0)
  if (ncol(cost) != n || any(!is.finite(cost))) {
    stop("internal assignment cost must be a finite square matrix",
         call. = FALSE)
  }
  u <- numeric(n + 1L)
  v <- numeric(n + 1L)
  p <- integer(n + 1L)
  way <- integer(n + 1L)

  # Sequential augmenting paths update the dual potentials and matching.
  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 1L
    minv <- rep(Inf, n + 1L)
    used <- logical(n + 1L)
    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- 0L
      for (j in 2L:(n + 1L)) {
        if (!used[j]) {
          cur <- cost[i0, j - 1L] - u[i0 + 1L] - v[j]
          if (cur < minv[j]) {
            minv[j] <- cur
            way[j] <- j0
          }
          if (minv[j] < delta) {
            delta <- minv[j]
            j1 <- j
          }
        }
      }
      for (j in seq_len(n + 1L)) {
        if (used[j]) {
          u[p[j] + 1L] <- u[p[j] + 1L] + delta
          v[j] <- v[j] - delta
        } else {
          minv[j] <- minv[j] - delta
        }
      }
      j0 <- j1
      if (p[j0] == 0L) break
    }
    repeat {
      j1 <- way[j0]
      p[j0] <- p[j1]
      j0 <- j1
      if (j0 == 1L) break
    }
  }
  -v[1L]
}

# =========================================================================
# Persistence landscapes (Bubenik 2015)
# =========================================================================

#' Table of a persistence landscape
#'
#' @param x A `persistence_landscape` object from [hg_landscape()].
#' @param what `"landscape"`, the only table.
#' @param ... Additional arguments (ignored).
#' @param k Integer or `NULL`. Keep only this landscape level.
#' @return A data.frame, one row per (level, grid point): `k` (landscape
#'   level), `t` (filtration value), `value` (landscape height).
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @examples
#' set.seed(1)
#' w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
#' w <- (w + t(w)) / 2
#' diag(w) <- 0
#' pl <- hg_landscape(hg_homology(w, n_steps = 10),
#'                             dimension = 0)
#' hg_get(pl, k = 1, top = 5)
#' @export
hg_get.persistence_landscape <- function(x, what = "landscape", ...,
                                         k = NULL, top = NULL) {
  match.arg(what, "landscape")
  out <- x$landscape
  if (!is.null(k)) {
    stopifnot("`k` must be a single integer >= 1" =
                is.numeric(k) && length(k) == 1L && k >= 1)
    out <- out[out$k == as.integer(k), , drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
