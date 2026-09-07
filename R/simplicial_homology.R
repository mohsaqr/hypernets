# ---- Persistent homology (simplicial family) -------------------------------
#
# Betti curves and persistence diagrams over a filtration. The
# filtration and Z/2 reduction machinery it calls lives in
# R/simplicial_filtration.R. Split out of simplicial.R in hypernets
# 0.2.0; the code is unchanged.

#' Persistent Homology
#'
#' @description
#' Computes persistent homology via full boundary-matrix reduction over
#' \eqn{\mathbb{Z}/2} (Edelsbrunner, Letscher & Zomorodian 2000). The
#' returned persistence diagram pairs each k-dimensional homology class
#' to the simplex whose addition creates it (birth) and the simplex whose
#' addition destroys it (death). Essential classes - those never killed -
#' are reported with \code{death = 0} in clique mode (similarity scale,
#' descending) and \code{death = Inf} in VR mode (distance scale, ascending).
#'
#' Two filtration modes are supported:
#' \describe{
#'   \item{\code{type = "clique"}}{Weighted clique filtration. Input is
#'     treated as a similarity matrix; high-weight simplices appear early.
#'     For each k-simplex \eqn{\sigma}, the filtration value is
#'     \eqn{\min_{(i,j) \in \sigma}\,|w(i,j)|}. Thresholds run high to low.}
#'   \item{\code{type = "vr"}}{Vietoris-Rips filtration on a non-negative
#'     distance matrix. For each k-simplex \eqn{\sigma}, the filtration
#'     value is \eqn{\max_{(i,j) \in \sigma}\,d(i,j)}. Thresholds run low
#'     to high. Use \code{max_scale} to cap the filtration diameter.}
#' }
#'
#' @param x A square matrix, \code{tna}, or \code{netobject}. For
#'   \code{type = "vr"}, must be a non-negative distance matrix.
#' @param n_steps Number of grid points for the reported Betti curve
#'   (default 20). The persistence diagram itself is exact - it does not
#'   depend on \code{n_steps}.
#' @param max_dim Maximum simplex dimension to track (default 3).
#' @param type Filtration: \code{"clique"} (default, similarity-weighted)
#'   or \code{"vr"} (Vietoris-Rips on distances).
#' @param max_scale For \code{type = "vr"} only: cap on edge length. Edges
#'   with \code{d(i,j) > max_scale} are excluded. \code{NULL} (default)
#'   uses \code{max(d)}.
#'
#' @return A \code{net_persistent_homology} object with:
#' \describe{
#'   \item{betti_curve}{Data frame: \code{threshold}, \code{dimension},
#'     \code{betti}.}
#'   \item{persistence}{Data frame of birth-death pairs:
#'     \code{dimension}, \code{birth}, \code{death}, \code{persistence}.
#'     Sorted by descending persistence.}
#'   \item{thresholds}{Numeric vector of grid thresholds.}
#'   \item{mode}{Either \code{"clique"} or \code{"vr"}.}
#' }
#'
#' @references
#' Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2000). Topological
#' persistence and simplification. \emph{Discrete & Computational Geometry}
#' \strong{28}, 511-533.
#'
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' ph <- persistent_homology(mat, n_steps = 10)
#' print(ph)
#'
#' @export
persistent_homology <- function(x, n_steps = 20L, max_dim = 3L,
                                type = "clique", max_scale = NULL) {
  stopifnot(
    is.numeric(n_steps), length(n_steps) == 1L,
    !is.na(n_steps), n_steps >= 1L, n_steps == as.integer(n_steps),
    is.numeric(max_dim), length(max_dim) == 1L,
    !is.na(max_dim), max_dim >= 0, max_dim == as.integer(max_dim)
  )
  type <- match.arg(type, c("clique", "vr"))

  # Filtered-complex handoff: build_simplicial(type = "vr") attaches a
  # $filtration vector. Consume it directly instead of rebuilding from a
  # matrix - this is the workflow advertised by the build_simplicial docs.
  if (inherits(x, "net_simplicial") && !is.null(x$filtration)) {
    fc <- .fc_from_filtered_complex(x, max_dim = max_dim,
                                    max_scale = max_scale)
  } else if (type == "vr") {
    d <- if (is.matrix(x)) x else .sc_extract_matrix(x)
    fc <- .filter_vr_complex(d, max_dim = max_dim, max_scale = max_scale)
    if (fc$max_w == 0 && fc$max_filt == 0 &&
        all(fc$dim == 0L)) {
      stop("All distances are zero or excluded; cannot build filtration.",
           call. = FALSE)
    }
  } else {
    mat <- .sc_extract_matrix(x)
    mat <- abs(mat)
    mat <- pmax(mat, t(mat))
    diag(mat) <- 0
    if (max(mat) == 0) {
      stop("All weights are zero; cannot build filtration.", call. = FALSE)
    }
    fc <- .filter_clique_complex(mat, max_dim = max_dim)
  }

  red <- .persistence_pairs_z2(fc)

  # Translate to user-facing scale and assemble persistence table
  if (fc$mode == "clique") {
    max_w <- fc$max_w
    pairs <- red$pairs
    if (nrow(pairs) > 0L) {
      bd_asc <- pairs$birth
      dd_asc <- pairs$death
      pairs$birth <- max_w - bd_asc
      pairs$death <- max_w - dd_asc
      pairs$persistence <- pairs$birth - pairs$death
    }
    ess <- red$essential
    if (nrow(ess) > 0L) {
      ess$birth <- max_w - ess$birth
      ess$death <- 0
      ess$persistence <- ess$birth
    }
    thresholds <- seq(max_w, max_w * 0.01, length.out = n_steps)
  } else {
    pairs <- red$pairs
    ess <- red$essential
    # Keep essential death = Inf so the Betti curve correctly counts them as
    # alive at the final threshold (death > t holds for any finite t). The
    # plot path caps Inf at max_filt for display.
    if (nrow(ess) > 0L) {
      ess$persistence <- Inf
    }
    thresholds <- seq(0, max(fc$max_filt, .Machine$double.eps),
                      length.out = n_steps)
  }

  persistence <- rbind(pairs, ess)
  persistence <- persistence[order(-persistence$persistence), , drop = FALSE]
  rownames(persistence) <- NULL

  # Drop dimensions above max_dim (defensive)
  persistence <- persistence[persistence$dimension <= max_dim, , drop = FALSE]

  betti_curve <- .betti_curve_from_pairs(persistence, thresholds, max_dim,
                                         fc$mode)

  structure(list(
    betti_curve = betti_curve,
    persistence = persistence,
    thresholds = thresholds,
    mode = fc$mode
  ), class = "net_persistent_homology")
}

# =========================================================================
# Simplicial centrality
# =========================================================================

#' Print persistent homology results
#' @param x A \code{net_persistent_homology} object.
#' @param ... Additional arguments (unused).
#' @return The input object, invisibly.
#' @examples
#' mat <- matrix(c(0,.6,.5,.6,0,.4,.5,.4,0), 3, 3)
#' colnames(mat) <- rownames(mat) <- c("A","B","C")
#' ph <- persistent_homology(mat, n_steps = 10)
#' print(ph)
#'
#' @export
print.net_persistent_homology <- function(x, ...) {
  cat("Persistent Homology\n")
  cat(sprintf("  %d filtration steps [%.4f \u2192 %.4f]\n",
              length(x$thresholds),
              max(x$thresholds), min(x$thresholds)))

  if (nrow(x$persistence) > 0L) {
    dims <- sort(unique(x$persistence$dimension))
    parts <- vapply(dims, function(d) {
      sub <- x$persistence[x$persistence$dimension == d, ]
      n_p <- sum(sub$death == 0)
      sprintf("b%d: %d (%d persistent)", d, nrow(sub), n_p)
    }, character(1))
    cat("  Features:", paste(parts, collapse = "  |  "), "\n")

    # Top 3 only
    top <- head(x$persistence[x$persistence$persistence > 0, ], 3)
    if (nrow(top) > 0L) {
      cat("  Longest-lived:\n")
      for (i in seq_len(nrow(top))) {
        cat(sprintf("    b%d: %.4f \u2192 %.4f (life: %.4f)\n",
                    top$dimension[i], top$birth[i], top$death[i],
                    top$persistence[i]))
      }
    }
  }
  invisible(x)
}

#' Plot Persistent Homology
#'
#' Two panels: Betti curve (threshold vs Betti number) and persistence
#' diagram (birth vs death). Persistence pairs come from full boundary-
#' matrix reduction; essential classes are shown at the filtration boundary
#' (\code{death = 0} in clique mode, \code{death = max_scale} in VR mode).
#'
#' @param x A \code{net_persistent_homology} object.
#' @param combined When `TRUE` (default), the two panels are stitched
#'   side-by-side via `gridExtra::arrangeGrob`. When `FALSE`, returns a
#'   named list (`betti_curve`, `persistence`) of ggplots.
#' @param ... Ignored.
#' @return A grid grob (invisibly) when `combined = TRUE`; a named list
#'   of two ggplots when `combined = FALSE`.
#'
#' @examples
#' \donttest{
#' seqs <- data.frame(
#'   V1 = c("A","B","C","A","B"),
#'   V2 = c("B","C","A","B","C"),
#'   V3 = c("C","A","B","C","A")
#' )
#' hon <- build_hon(seqs, max_order = 1)
#' ph  <- persistent_homology(hon)
#' if (requireNamespace("gridExtra", quietly = TRUE)) plot(ph)
#' }
#'
#' @export
plot.net_persistent_homology <- function(x, combined = TRUE, ...) {
  stopifnot(is.logical(combined), length(combined) == 1L)

  filt <- x$betti_curve
  filt$dim_label <- factor(paste0("B", filt$dimension))

  # --- Panel 1: Betti curve ---
  p1 <- ggplot2::ggplot(filt, ggplot2::aes(x = threshold, y = betti,
                                              color = dim_label)) +
    ggplot2::geom_step(linewidth = 1.1, direction = "vh") +
    ggplot2::scale_color_brewer(palette = "Set1") +
    ggplot2::labs(title = "Betti Curve",
                  subtitle = "Betti numbers across inclusive weight thresholds",
                  x = "Weight Threshold", y = "Betti Number",
                  color = NULL) +
    .sc_theme() +
    ggplot2::theme(legend.position = "top")

  # --- Panel 2: Persistence diagram ---
  pers <- x$persistence[x$persistence$persistence > 0, ]
  # Cap both Inf death and Inf persistence for display (essential VR classes).
  # ggplot's continuous size scale drops Inf-valued rows, which would silently
  # erase essential features (e.g., the surviving H_0 in VR mode). Cap to the
  # finite max so every row renders.
  inf_death <- !is.finite(pers$death)
  inf_pers  <- !is.finite(pers$persistence)
  if (any(inf_death) || any(inf_pers)) {
    finite_max <- max(c(pers$birth,
                        pers$death[!inf_death],
                        pers$persistence[!inf_pers],
                        x$thresholds), na.rm = TRUE)
    pers$death[inf_death]      <- finite_max
    pers$persistence[inf_pers] <- finite_max
  }

  if (nrow(pers) > 0L) {
    pers$dim_label <- factor(paste0("B", pers$dimension))
    lim <- max(c(pers$birth, pers$death)) * 1.15

    p2 <- ggplot2::ggplot(pers, ggplot2::aes(x = birth, y = death,
                                                color = dim_label,
                                                size = persistence)) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                            color = "grey60") +
      ggplot2::geom_point(alpha = 0.7) +
      ggplot2::scale_size_continuous(range = c(1.5, 6), guide = "none") +
      ggplot2::scale_color_brewer(palette = "Set1") +
      ggplot2::coord_equal(xlim = c(0, lim), ylim = c(0, lim)) +
      ggplot2::labs(title = "Persistence Diagram",
                    subtitle = "Boundary-matrix reduction over Z/2",
                    x = "Birth", y = "Death", color = NULL) +
      .sc_theme() +
      ggplot2::theme(legend.position = "top")
  } else { # nocov start
    p2 <- ggplot2::ggplot() +
      ggplot2::labs(title = "Persistence Diagram",
                    subtitle = "No features detected") +
      .sc_theme() # nocov end
  }

  panels <- list(betti_curve = p1, persistence = p2)
  if (!combined) return(invisible(panels))
  if (!requireNamespace("gridExtra", quietly = TRUE)) {
    stop("combined = TRUE requires the gridExtra package.", call. = FALSE) # nocov
  }
  combined_grob <- gridExtra::arrangeGrob(p1, p2, ncol = 2,
                                          padding = grid::unit(0.5, "line"))
  grid::grid.newpage()
  grid::grid.draw(combined_grob)
  invisible(combined_grob)
}

#' Coerce a net_persistent_homology to a tidy table
#'
#' @param x A `net_persistent_homology` object.
#' @param row.names Ignored (S3 consistency).
#' @param optional Ignored (S3 consistency).
#' @param ... Additional arguments (ignored).
#' @param what `"persistence"` (default) for the persistence diagram, or
#'   `"betti"` for the Betti curves across the filtration.
#' @param dimension Integer or `NULL`. Keep only this homological dimension.
#' @param sort_by `NULL` (construction order, default) or `"persistence"` -
#'   longest-lived features first. Only for `what = "persistence"`.
#' @return A data.frame. For `what = "persistence"`, one row per feature:
#'   `dimension`, `birth`, `death`, `persistence`. For `what = "betti"`, one
#'   row per (threshold, dimension): `threshold`, `dimension`, `betti`.
#' @param top Integer or `NULL`. Return only the first `top` rows,
#'   applied after any filter and after `sort_by`, so `sort_by` and
#'   `top` compose. Default `NULL` returns every row.
#' @export
as.data.frame.net_persistent_homology <- function(
    x, row.names = NULL, optional = FALSE, ...,
    what = c("persistence", "betti"), dimension = NULL, sort_by = NULL,
    top = NULL) {
  what <- match.arg(what)
  out <- if (what == "betti") x$betti_curve else x$persistence
  if (!is.null(dimension)) {
    stopifnot("`dimension` must be a single integer >= 0" =
                is.numeric(dimension) && length(dimension) == 1L &&
                dimension >= 0)
    out <- out[out$dimension == as.integer(dimension), , drop = FALSE]
  }
  if (!is.null(sort_by)) {
    sort_by <- match.arg(sort_by, "persistence")
    if (what == "betti") {
      stop(errorCondition(
        "`sort_by` applies only to what = \"persistence\"",
        class = "hypernets_bad_input", call = NULL))
    }
    out <- out[order(-out$persistence, out$dimension, out$birth), ,
               drop = FALSE]
  }
  rownames(out) <- NULL
  .ho_top(out, top)
}
