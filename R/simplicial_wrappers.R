# ---- Simplicial-complex constructor and verbs ----
#
# simplicial(), hg_homology(), hg_landscape(), hg_bottleneck(), hg_betti(),
# hg_euler(), hg_qanalysis() and hg_degree() compute nothing new: each calls
# the Nestimate estimator and returns its result unchanged, so every number
# is Nestimate's (identity asserted in tests/testthat/test-simplicial_wrappers.R).
# simplicial() adds one thing: an independent check of the clique complex it
# returns (igraph's clique enumeration and the Euler-Poincare identity),
# surfaced as a classed warning, never a silent pass.

#' Simplicial complex from sequences, a network or a memory network
#'
#' Builds a simplicial complex: the sets of actions that occur together in
#' windows of consecutive actions (`type = "window"`), every clique of a
#' thresholded weighted network (`type = "clique"`), every higher-order
#' pathway of a memory network (`type = "pathway"`), or the Vietoris-Rips
#' complex of a distance matrix (`type = "vr"`). Every face of every simplex
#' is included.
#'
#' `type = "window"` is a co-occurrence complex (Salnikov et al. 2018) whose
#' unit of co-occurrence is a window of `window` consecutive actions in a
#' sequence (Ding et al. 2020). The distinct actions of each window form a
#' set; a set observed in at least `min_count` windows becomes a simplex,
#' together with all its subsets. A window of `window` actions yields
#' simplices of dimension at most `window - 1`.
#'
#' With `validate = TRUE` a window complex needs no count threshold. A set
#' of actions becomes a simplex when it occurs in more windows than a null
#' model predicts. With `null = "swap"` (default) the null is swap
#' randomization (Gionis et al. 2007): the windows and their actions are
#' shuffled by checkerboard swaps (Gotelli 2000) so that every window keeps
#' its number of actions and every action its number of windows, and only
#' which actions share a window changes. The p-value of a set is the share of
#' `n_null` shuffled copies in which the set occurs at least as often as
#' observed, counting the observed data (Phipson & Smyth 2010). With
#' `null = "hypergeometric"` the null is the statistically validated
#' hypergraph of Musciotto, Battiston & Mantegna (2021): a set of n actions
#' is tested against the windows with exactly n distinct actions, each
#' action falls into these windows independently of the others, and the
#' p-value is exact (their Eq. 3). That null ignores that a window holds
#' exactly n actions. It suits sparse co-occurrence, where most possible
#' sets of actions are rare or absent and each is expected in less than one
#' window, as with many actors who meet in small groups. On dense data, such
#' as a few actions that all occur together often, it expects more
#' co-occurrence than the windows can hold and validates little; when the
#' median tested set of a size is expected in one window or more it raises
#' the warning `hypernets_dense_cooccurrence`. For either
#' null the p-values of the sets of each size are corrected with the
#' Benjamini-Hochberg false discovery rate over all possible sets of that
#' size, and a set is validated when its adjusted p-value is at most
#' `alpha`. Every action remains a vertex. The tests are read with
#' `hg_get(x, what = "validation")`.
#'
#' `type = "clique"` builds the clique (flag) complex of a network: every
#' clique of the thresholded network becomes a simplex (Giusti, Ghrist &
#' Bassett 2016). A directed network is made undirected by `direction`: with
#' `"either"` (default) a pair is joined when the larger of its two weights
#' reaches `threshold`, and with `"both"` when both weights do. `"both"` gives
#' the cliques of `tna::cliques()` with the same threshold (Tikka,
#' Lopez-Pernas & Saqr 2025): the triangles of the complex are its cliques of
#' size three, and so on. Given sequences instead of a network, the clique
#' complex is built on their transition network: the relative transition
#' probabilities that the tna family estimates from the same sequences.
#'
#' A clique complex is checked on construction (`verify = TRUE`): its
#' simplices must be exactly the cliques igraph enumerates on the same
#' thresholded graph (when igraph is installed), and its Euler characteristic
#' from the face counts must equal the alternating sum of its Betti numbers
#' (the Euler-Poincare identity). A failed check raises the warning
#' `hypernets_simplicial_unverified`; the complex is returned either way.
#'
#' @param x Sequences for `type = "window"`, in any form described in
#'   [sequence-input]; a square weighted matrix, a network object carrying
#'   one (`netobject`, `tna`) or sequences for `type = "clique"`; a memory
#'   network
#'   ([hon()], [hypa()] on sequences, [mogen()]) for `type = "pathway"`; or
#'   for `type = "vr"` a distance matrix, a `dist` object, or points: a
#'   table or matrix with one row per point and one numeric column per
#'   coordinate (not square), or a data frame of numeric coordinates, whose
#'   one character column, if present, names the points. Points are turned
#'   into Euclidean distances. A square numeric matrix is read as distances.
#' @param type `"clique"` (default), `"window"`, `"pathway"`, or `"vr"`
#'   (alias `"rips"`).
#' @param window For `type = "window"`: the number of consecutive actions in
#'   a window. Default `3`.
#' @param min_count For `type = "window"`: the number of windows in which a
#'   set of actions must occur to become a simplex. Default `1`. Not used
#'   with `validate = TRUE`.
#' @param validate For `type = "window"`: `TRUE` keeps only the sets of
#'   actions that pass the statistical validation described in Details.
#'   Default `FALSE`.
#' @param alpha For `validate = TRUE`: the false discovery rate. Default
#'   `0.05`.
#' @param null For `validate = TRUE`: `"swap"` (default) or
#'   `"hypergeometric"` (see Details).
#' @param n_null For `null = "swap"`: the number of shuffled copies. The
#'   smallest possible p-value is `1 / (n_null + 1)`, and a set can pass on
#'   its own only if this is at most `alpha` divided by the number of
#'   possible sets of its size. `NULL` (default) takes the smallest number
#'   that meets this for every size, and at least 999. A number given here
#'   that does not meet it raises the warning `hypernets_low_resolution`
#'   with the number that would.
#' @param seed For `null = "swap"`: a seed for the shuffles, applied locally
#'   so the global random state is restored. `NULL` (default) uses the
#'   current random state.
#' @param action,actor,time,session,time_threshold,timezone For `type =
#'   "window"`: long-format arguments (see [sequence-input]); leave the
#'   column names `NULL` for wide or list input.
#' @param threshold For `type = "clique"`: the smallest absolute edge
#'   weight kept (default `0`, every non-zero edge).
#' @param direction For `type = "clique"` on a directed network: `"either"`
#'   (default, a pair is joined when the larger of its two weights reaches
#'   `threshold`) or `"both"` (when both weights do, as in `tna::cliques()`).
#' @param max_dim Integer. Highest simplex dimension (a k-simplex has k + 1
#'   nodes). Default `10`.
#' @param max_pathways For `type = "pathway"`: keep at most this many
#'   pathways, ranked by count (HON) or ratio (HYPA). `NULL` keeps all.
#' @param anomaly For a HYPA pathway complex: `"all"` (default), `"over"`
#'   or `"under"`.
#' @param max_scale For `type = "vr"`: the longest edge admitted. `NULL`
#'   uses the largest distance.
#' @param verify Logical. Check a clique complex on construction (see
#'   Details). Default `TRUE`; set `FALSE` to skip the check on very large
#'   complexes.
#' @param ... For `type = "pathway"` on a network object: arguments passed
#'   to the memory-network builder.
#' @return A `simplicial_complex` object. Read it with [hg_get()]: one row
#'   per simplex (`what = "simplices"`, the default), the face counts
#'   (`"f_vector"`), the per-node simplicial degree (`"degree"`) or, for a
#'   validated window complex, one row per tested set of actions
#'   (`"validation"`), and the Betti numbers (`"betti"`). `plot()` draws
#'   the maximal simplices, those no larger simplex contains, as regions
#'   around their nodes (`cograph::plot_simplicial()`; its arguments pass
#'   through `...`): the most significant first for a validated window
#'   complex, the most frequent first for a window complex, the closest
#'   first for a Vietoris-Rips complex. `plot()` returns the figure. For the
#'   `plot()` method, `x` is the complex.
#' @references
#' Hatcher, A. (2002). \emph{Algebraic Topology}. Cambridge University
#' Press.
#'
#' Giusti, C., Ghrist, R., & Bassett, D. S. (2016). Two's company, three (or
#' more) is a simplex. \emph{Journal of Computational Neuroscience}, 41,
#' 1-14. \doi{10.1007/s10827-016-0608-6}
#'
#' Tikka, S., Lopez-Pernas, S., & Saqr, M. (2025). tna: An R package for
#' transition network analysis. \emph{Applied Psychological Measurement}, 49,
#' 326-328. \doi{10.1177/01466216251348840}
#'
#' Salnikov, V., Cassese, D., Lambiotte, R., & Jones, N. S. (2018).
#' Co-occurrence simplicial complexes in mathematics: identifying the holes
#' of knowledge. \emph{Applied Network Science}, 3, 37.
#' \doi{10.1007/s41109-018-0074-3}
#'
#' Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
#' Hypergraph attention networks for inductive text classification.
#' \emph{Proceedings of EMNLP 2020}, 4927-4936.
#' \doi{10.18653/v1/2020.emnlp-main.399}
#'
#' Gionis, A., Mannila, H., Mielikainen, T., & Tsaparas, P. (2007).
#' Assessing data mining results via swap randomization. \emph{ACM
#' Transactions on Knowledge Discovery from Data}, 1(3), 14.
#' \doi{10.1145/1297332.1297338}
#'
#' Gotelli, N. J. (2000). Null model analysis of species co-occurrence
#' patterns. \emph{Ecology}, 81(9), 2606-2621.
#' \doi{10.1890/0012-9658(2000)081[2606:NMAOSC]2.0.CO;2}
#'
#' Phipson, B., & Smyth, G. K. (2010). Permutation p-values should never be
#' zero. \emph{Statistical Applications in Genetics and Molecular Biology},
#' 9(1), 39. \doi{10.2202/1544-6115.1585}
#'
#' Musciotto, F., Battiston, F., & Mantegna, R. N. (2021). Detecting
#' informative higher-order interactions in statistically validated
#' hypergraphs. \emph{Communications Physics}, 4, 218.
#' \doi{10.1038/s42005-021-00710-4}
#'
#' Battiston, F., Cencetti, G., Iacopini, I., Latora, V., Lucas, M.,
#' Patania, A., Young, J.-G., & Petri, G. (2020). Networks beyond pairwise
#' interactions: Structure and dynamics. \emph{Physics Reports}, 874,
#' 1--92. \doi{10.1016/j.physrep.2020.05.004}
#' @seealso [hg_betti()], [hg_euler()], [hg_qanalysis()], [hg_homology()]
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' sc <- simplicial(adj)
#' hg_get(sc)
#' hg_get(sc, what = "f_vector")
#' @export
simplicial <- function(x, type = "clique", threshold = 0, max_dim = 10L,
                       max_pathways = NULL,
                       anomaly = c("all", "over", "under"),
                       max_scale = NULL, verify = TRUE,
                       direction = c("either", "both"),
                       window = 3L, min_count = 1L, validate = FALSE,
                       alpha = 0.05, null = c("swap", "hypergeometric"),
                       n_null = NULL, seed = NULL, action = NULL,
                       actor = NULL, time = NULL, session = NULL,
                       time_threshold = 900, timezone = "UTC", ...) {
  stopifnot(
    "`verify` must be TRUE or FALSE" =
      is.logical(verify) && length(verify) == 1L && !is.na(verify)
  )
  direction <- match.arg(direction)
  if (identical(type, "clique") && .sc_is_sequence_input(x, action)) {
    # sequences: the relative transition network, as the tna family
    # estimates it from the same sequences
    x <- Nestimate::build_network(
      .ho_sequence_input(x, action = action, actor = actor, time = time,
                         session = session, time_threshold = time_threshold,
                         timezone = timezone, lists = "wide",
                         models = "decode"),
      method = "relative")
  }
  if (identical(direction, "both") && identical(type, "clique")) {
    # a pair counts only when both directed weights reach the threshold:
    # the smaller of the two absolute weights must reach it
    w <- abs(.sc_weight_matrix(x))
    x <- pmin(w, t(w))
  }
  if (identical(type, "window")) {
    sc <- .sc_from_windows(x, window = window, min_count = min_count,
                           validate = validate, alpha = alpha,
                           null = match.arg(null), n_null = n_null,
                           seed = seed,
                           max_dim = max_dim, action = action, actor = actor,
                           time = time, session = session,
                           time_threshold = time_threshold,
                           timezone = timezone)
    return(structure(sc, class = c("hypernets_simplicial", class(sc))))
  }
  if (type %in% c("vr", "rips")) x <- .sc_distances(x)
  sc <- Nestimate::build_simplicial(x, type = type, threshold = threshold,
                                    max_dim = max_dim,
                                    max_pathways = max_pathways,
                                    anomaly = anomaly, max_scale = max_scale,
                                    ...)
  if (verify && identical(sc$type, "clique")) {
    .sc_verify(sc, x, threshold = threshold, max_dim = max_dim)
  }
  # one added class, so plot() shows the simplices through cograph
  structure(sc, class = c("hypernets_simplicial", class(sc)))
}

#' @rdname simplicial
#' @param y Unused.
#' @param dismantled For `plot()`: `TRUE` draws each maximal simplex in its
#'   own panel; `FALSE` (default) draws them together on one layout.
#' @param top For `plot()`: the number of maximal simplices to draw, in the
#'   order described under Value. `NULL` (default) draws all of them.
#' @export
plot.hypernets_simplicial <- function(x, y, dismantled = FALSE, top = NULL,
                                      ...) {
  stopifnot(
    "`dismantled` must be TRUE or FALSE" =
      is.logical(dismantled) && length(dismantled) == 1L && !is.na(dismantled),
    "`top` must be NULL or a single whole number >= 1" =
      is.null(top) || (is.numeric(top) && length(top) == 1L &&
                         is.finite(top) && top >= 1 && top == round(top))
  )
  facets <- .sc_facets(x)
  if (length(facets) == 0L) {
    .thg_bad_input("the complex has no simplex of two or more nodes to plot")
  }
  shown <- if (is.null(top)) length(facets) else min(top, length(facets))
  args <- list(
    pathways = facets, ordered = FALSE, max_pathways = NULL,
    pathway_index = seq_len(shown), dismantled = dismantled,
    title = sprintf("%d of %d maximal simplices", shown, length(facets)),
    node_color = "#0072B2", ring_color = "#E69F00",
    blob_colors = rep_len(.thg_okabe_ito, shown))
  # many nodes on one circle sit close together; smaller nodes keep the
  # neighbours apart
  n_drawn <- length(unique(unlist(utils::head(facets, shown))))
  if (n_drawn > 12L) {
    args$node_size <- 12
    args$label_size <- 3.2
  }
  dots <- list(...)
  args[names(dots)] <- dots
  do.call(cograph::plot_simplicial, args)
}

# The maximal simplices (facets) of a complex with two or more nodes, as
# vectors of node names: the simplices no other simplex contains. Order:
# for a validated window complex by the z of its test, for a window complex
# by its count, for a Vietoris-Rips complex by the scale at which it enters
# (closest first), otherwise by size and then name.
.sc_facets <- function(x) {
  simplices <- lapply(x$simplices, sort)
  sizes <- lengths(simplices)
  keys <- vapply(simplices, paste, character(1L), collapse = ",")
  # a simplex is a face of another exactly when the other contains it; a
  # face of one more node suffices, since faces of faces are faces
  covered <- unique(unlist(lapply(simplices[sizes >= 2L], \(s) {
    vapply(seq_along(s), \(i) paste(s[-i], collapse = ","), character(1L))
  })))
  maximal <- which(sizes >= 2L & !keys %in% covered)
  names_of <- lapply(simplices[maximal], \(s) x$nodes[s])
  labels <- vapply(names_of, paste, character(1L), collapse = ", ")
  strength <- if (!is.null(x$validation)) {
    x$validation$z[match(labels, x$validation$members)]
  } else if (!is.null(x$count)) {
    x$count[maximal]
  } else if (identical(x$type, "vr") && !is.null(x$filtration)) {
    -x$filtration[maximal]
  } else {
    rep(0, length(maximal))
  }
  strength[is.na(strength)] <- -Inf
  names_of[order(-sizes[maximal], -strength, labels)]
}

#' @rdname simplicial
#' @export
hg_get.hypernets_simplicial <- function(x, ...) NextMethod()

#' @rdname result-summary
#' @export
summary.hypernets_simplicial <- function(object, ...) .ho_summary(object)

#' @rdname simplicial
#' @param n Number of rows of the default table to print. Default `10`.
#' @export
print.hypernets_simplicial <- function(x, n = 10L, ...) {
  betti <- hg_betti(x)
  kind <- switch(x$type %||% "",
    window = if (is.null(x$validation)) {
      sprintf("Window co-occurrence complex (window = %d, min_count = %d)",
              x$window, x$min_count)
    } else {
      sprintf(paste0("Validated window co-occurrence complex (window = %d, ",
                     "%s, FDR alpha = %s): %d of %d sets of actions ",
                     "validated"),
              x$window,
              if (identical(x$null, "swap")) {
                sprintf("swap null, %d shuffles", x$n_null)
              } else {
                "hypergeometric null"
              },
              format(x$alpha), sum(x$validation$significant),
              nrow(x$validation))
    },
    clique = "Clique complex",
    vr = "Vietoris-Rips complex",
    pathway = "Pathway complex",
    "Simplicial complex")
  cat(kind, "\n", sep = "")
  cat(sprintf("  %d nodes, %d simplices, dimension %d; f-vector (%s); Betti (%s)\n",
              x$n_nodes, x$n_simplices, x$dimension,
              paste(x$f_vector, collapse = ", "),
              paste(betti, collapse = ", ")))
  .ho_print_table(x, n)
  invisible(x)
}

#' Check a clique complex against igraph and Euler-Poincare
#'
#' @param sc The `simplicial_complex` just built.
#' @param x The input it was built from.
#' @param threshold,max_dim The construction arguments.
#' @return `NULL`, invisibly; raises the warning
#'   `hypernets_simplicial_unverified` when a check fails.
#' @noRd
.sc_verify <- function(sc, x, threshold, max_dim) {
  problems <- character(0)
  f <- as.integer(sc$f_vector)
  betti <- as.integer(Nestimate::betti_numbers(sc))
  euler_faces <- sum((-1L)^(seq_along(f) - 1L) * f)
  euler_betti <- sum((-1L)^(seq_along(betti) - 1L) * betti)
  if (euler_faces != euler_betti) {
    problems <- c(problems, sprintf(
      "the Euler characteristic from the face counts (%d) differs from the alternating Betti sum (%d)",
      euler_faces, euler_betti))
  }
  mat <- if (is.matrix(x)) x else x$weights
  if (requireNamespace("igraph", quietly = TRUE) && is.matrix(mat) &&
      is.numeric(mat)) {
    weights <- abs(mat)
    diag(weights) <- 0
    weights <- pmax(weights, t(weights))
    adj <- weights > 0 & weights >= threshold
    g <- igraph::graph_from_adjacency_matrix(adj, mode = "undirected",
                                             diag = FALSE)
    cliques <- igraph::cliques(g, max = as.integer(max_dim) + 1L)
    ours <- sort(vapply(sc$simplices, function(s) {
      paste(sort(as.integer(s)), collapse = ",")
    }, character(1L)))
    theirs <- sort(vapply(cliques, function(cl) {
      paste(sort(as.integer(cl)), collapse = ",")
    }, character(1L)))
    if (!identical(ours, theirs)) {
      problems <- c(problems, sprintf(
        "its %d simplices are not igraph's %d cliques of the thresholded graph",
        length(ours), length(theirs)))
    }
  }
  if (length(problems)) {
    warning(warningCondition(
      paste0("simplicial(): the clique complex failed verification: ",
             paste(problems, collapse = "; "), "."),
      class = "hypernets_simplicial_unverified", call = NULL))
  }
  invisible(NULL)
}

#' Persistent homology of a weighted network or a distance matrix
#'
#' Tracks the homology classes (components, loops, voids) of a filtration:
#' a weighted clique filtration of a similarity matrix (`type = "clique"`,
#' strongest edges first) or a Vietoris-Rips filtration of a distance
#' matrix (`type = "vr"`). Every class is paired with the simplex that
#' creates it (birth) and the one that destroys it (death), by exact
#' boundary-matrix reduction over Z/2.
#'
#' @param x A square matrix, a network object carrying one (`netobject`,
#'   `tna`), points for `type = "vr"` (see [simplicial()]), a
#'   `simplicial_complex` built with `type = "vr"` (its stored
#'   filtration is used as is), or a window complex from
#'   `simplicial(type = "window")`. For a window complex the filtration is
#'   the count filtration (Petri et al. 2013): the complex at count `t` holds
#'   the sets of actions observed in at least `t` windows and their faces, so
#'   a class is born at the count where it appears and dies at a lower count;
#'   `birth`, `death` and the Betti-curve `threshold` are counts, and
#'   `n_steps` is not used (every observed count is a threshold).
#' @param n_steps Grid points of the reported Betti curve. Default `20`; the
#'   persistence diagram itself is exact.
#' @param max_dim Highest simplex dimension tracked. Default `3`.
#' @param type `"clique"` (default) or `"vr"`.
#' @param max_scale For `type = "vr"`: the longest edge admitted. `NULL`
#'   uses the largest distance.
#' @return A `persistent_homology` object. Read it with [hg_get()]: the
#'   persistence diagram (`what = "persistence"`, the default) or the Betti
#'   curves (`"betti"`).
#' @references
#' Edelsbrunner, H., Letscher, D., & Zomorodian, A. (2002). Topological
#' persistence and simplification. \emph{Discrete & Computational Geometry},
#' 28, 511--533. \doi{10.1007/s00454-002-2885-2}
#'
#' Petri, G., Scolamiero, M., Donato, I., & Vaccarino, F. (2013).
#' Topological strata of weighted complex networks. \emph{PLoS ONE}, 8(6),
#' e66506. \doi{10.1371/journal.pone.0066506}
#' @seealso [hg_landscape()], [hg_bottleneck()], [hg_wasserstein()]
#' @examples
#' set.seed(1)
#' w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
#' w <- (w + t(w)) / 2
#' diag(w) <- 0
#' ph <- hg_homology(w, n_steps = 10)
#' hg_get(ph, sort_by = "persistence")
#' @export
hg_homology <- function(x, n_steps = 20L, max_dim = 3L, type = "clique",
                        max_scale = NULL) {
  if (inherits(x, "simplicial_complex") && identical(x$type, "window")) {
    return(.ho_result(.ph_window(x, max_dim = max_dim)))
  }
  if (type %in% c("vr", "rips")) x <- .sc_distances(x)
  .ho_result(Nestimate::persistent_homology(x, n_steps = n_steps,
                                            max_dim = max_dim, type = type,
                                            max_scale = max_scale))
}

#' Persistence landscape of a persistence diagram
#'
#' The persistence landscape: at every filtration value t, the k-th
#' largest of the tent functions \eqn{\max(0, \min(t - b, d - t))} of the
#' diagram's (birth, death) pairs in one homology dimension. Essential
#' classes, which never die, are left out.
#'
#' @param ph A [hg_homology()] result, or a data.frame with columns
#'   `dimension`, `birth`, `death`.
#' @param k_max Highest landscape level. Default `5`.
#' @param dimension Homology dimension. Default `1`.
#' @param t_grid Evaluation points. `NULL` uses 200 points spanning the
#'   pairs.
#' @return A `persistence_landscape` object. Read it with [hg_get()]: one
#'   row per level and grid point (`what = "landscape"`).
#' @references
#' Bubenik, P. (2015). Statistical topological data analysis using
#' persistence landscapes. \emph{Journal of Machine Learning Research}, 16,
#' 77--102.
#' @seealso [hg_homology()]
#' @examples
#' set.seed(1)
#' w <- matrix(runif(36), 6, 6, dimnames = list(letters[1:6], letters[1:6]))
#' w <- (w + t(w)) / 2
#' diag(w) <- 0
#' pl <- hg_landscape(hg_homology(w, n_steps = 10), dimension = 0)
#' hg_get(pl, k = 1, top = 5)
#' @export
hg_landscape <- function(ph, k_max = 5L, dimension = 1L, t_grid = NULL) {
  .ho_result(Nestimate::persistence_landscape(.ho_unresult(ph),
                                              k_max = k_max,
                                              dimension = dimension,
                                              t_grid = t_grid))
}

#' Bottleneck distance between persistence diagrams
#'
#' The largest displacement of the best matching between two persistence
#' diagrams, where a point may also match the diagonal. Essential classes
#' are matched only to essential classes; a dimension whose essential
#' counts differ is `Inf` apart.
#'
#' @param d1,d2 [hg_homology()] results, or data.frames with columns
#'   `dimension`, `birth`, `death`.
#' @param dimension Integer vector of dimensions to compare. `NULL` compares
#'   every dimension present in either diagram.
#' @param tol Tolerance of the binary search. Default
#'   `.Machine$double.eps^0.5`.
#' @return A named numeric vector, one distance per dimension, named
#'   `"dim_<k>"`.
#' @references
#' Cohen-Steiner, D., Edelsbrunner, H., & Harer, J. (2007). Stability of
#' persistence diagrams. \emph{Discrete & Computational Geometry}, 37,
#' 103--120. \doi{10.1007/s00454-006-1276-5}
#'
#' Edelsbrunner, H., & Harer, J. (2010). \emph{Computational Topology: An
#' Introduction}. American Mathematical Society.
#' @seealso [hg_wasserstein()]
#' @examples
#' d1 <- data.frame(dimension = 0L, birth = 0, death = 2)
#' d2 <- data.frame(dimension = 0L, birth = 0, death = 3)
#' hg_bottleneck(d1, d2)
#' @export
hg_bottleneck <- function(d1, d2, dimension = NULL,
                          tol = .Machine$double.eps^0.5) {
  Nestimate::bottleneck_distance(d1, d2, dimension = dimension, tol = tol)
}

#' Betti numbers of a simplicial complex
#'
#' The ranks of the homology groups over Z/2: \eqn{\beta_0} counts connected
#' components, \eqn{\beta_1} independent loops, \eqn{\beta_2} voids, and so
#' on.
#'
#' @param sc A `simplicial_complex` from [simplicial()].
#' @return A named integer vector `c(b0 = , b1 = , ...)`.
#' @references
#' Hatcher, A. (2002). \emph{Algebraic Topology}. Cambridge University
#' Press.
#' @seealso [hg_euler()], [hg_homology()]
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' hg_betti(simplicial(adj))
#' @export
hg_betti <- function(sc) {
  Nestimate::betti_numbers(sc)
}

#' Euler characteristic of a simplicial complex
#'
#' \eqn{\chi = \sum_k (-1)^k f_k}, the alternating sum of the number of
#' k-simplices; by the Euler-Poincare theorem it equals the alternating sum
#' of the Betti numbers.
#'
#' @param sc A `simplicial_complex` from [simplicial()].
#' @return A single integer.
#' @references
#' Hatcher, A. (2002). \emph{Algebraic Topology}. Cambridge University
#' Press.
#' @seealso [hg_betti()]
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' hg_euler(simplicial(adj))
#' @export
hg_euler <- function(sc) {
  Nestimate::euler_characteristic(sc)
}

#' Q-analysis of a simplicial complex
#'
#' Q-connectivity (Atkin 1974): two maximal simplices are q-connected when
#' they share a face of dimension at least q. Reports the number of
#' q-connected components at every level and each node's highest level.
#'
#' @param sc A `simplicial_complex` from [simplicial()].
#' @return A `q_analysis` object. Read it with [hg_get()]: the components
#'   per level (`what = "q_levels"`, the default) or each node's highest
#'   level (`"nodes"`).
#' @references
#' Atkin, R. H. (1974). \emph{Mathematical Structure in Human Affairs}.
#' Heinemann.
#' @seealso [simplicial()]
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' hg_get(hg_qanalysis(simplicial(adj)))
#' @export
hg_qanalysis <- function(sc) {
  .ho_result(Nestimate::q_analysis(sc))
}

#' Simplicial degree of every node
#'
#' How many simplices of each dimension contain each node.
#'
#' @param sc A `simplicial_complex` from [simplicial()].
#' @param normalized Logical. Divide each count by its largest possible
#'   value. Default `FALSE`.
#' @return A base `data.frame`, one row per node, most connected first:
#'   `node`, one column per dimension (`d0`, `d1`, ...), and `total` (the
#'   simplices of dimension one and above). The same table is
#'   `hg_get(sc, what = "degree")`.
#' @references
#' Serrano, D. H., & Gomez, D. S. (2020). Centrality measures in simplicial
#' complexes: Applications of topological data analysis to network science.
#' \emph{Applied Mathematics and Computation}, 382, 125331.
#' \doi{10.1016/j.amc.2020.125331}
#' @seealso [simplicial()]
#' @examples
#' adj <- matrix(c(0, 1, 1, 0,
#'                 1, 0, 1, 1,
#'                 1, 1, 0, 1,
#'                 0, 1, 1, 0), 4, 4, dimnames = list(letters[1:4], letters[1:4]))
#' hg_degree(simplicial(adj))
#' @export
hg_degree <- function(sc, normalized = FALSE) {
  .sc_degree_table(sc, normalized = normalized)
}

#' Simplicial degree table
#'
#' @param sc A `simplicial_complex`.
#' @param normalized Logical.
#' @return data.frame with default row names.
#' @noRd
.sc_degree_table <- function(sc, normalized = FALSE) {
  if (!inherits(sc, "simplicial_complex")) {
    .ho_bad_input("`sc` must be a simplicial_complex from simplicial()")
  }
  out <- Nestimate::simplicial_degree(sc, normalized = normalized)
  rownames(out) <- NULL
  out
}


# Co-occurrence complex of windows of consecutive actions: the distinct
# actions of each window form a set (window_hypergraph() collapses repeats
# and counts identical sets); sets seen in >= min_count windows are simplices,
# closed under subsets.
.sc_from_windows <- function(x, window, min_count, validate, alpha,
                             null = "swap", n_null = NULL, seed = NULL,
                             max_dim, action, actor, time, session,
                             time_threshold = 900, timezone = "UTC") {
  stopifnot(
    "`window` must be a single integer >= 2" =
      is.numeric(window) && length(window) == 1L && is.finite(window) &&
      window >= 2 && window == round(window),
    "`min_count` must be a single count >= 1" =
      is.numeric(min_count) && length(min_count) == 1L &&
      is.finite(min_count) && min_count >= 1,
    "`validate` must be TRUE or FALSE" =
      is.logical(validate) && length(validate) == 1L && !is.na(validate),
    "`alpha` must be a single number in (0, 1)" =
      is.numeric(alpha) && length(alpha) == 1L && is.finite(alpha) &&
      alpha > 0 && alpha < 1,
    "`n_null` must be NULL or a single integer >= 19" =
      is.null(n_null) || (is.numeric(n_null) && length(n_null) == 1L &&
                            is.finite(n_null) && n_null >= 19 &&
                            n_null == round(n_null))
  )
  if (validate && min_count != 1) {
    .thg_bad_input(paste0("`min_count` is a count threshold and ",
                          "`validate = TRUE` replaces it; leave `min_count` ",
                          "at 1"))
  }
  data <- .ho_sequence_input(x, action = action, actor = actor, time = time,
                             session = session,
                             time_threshold = time_threshold,
                             timezone = timezone)
  wh <- window_hypergraph(data, window = as.integer(window))
  validation <- NULL
  if (validate) {
    if (!is.null(seed) && identical(null, "swap")) {
      old_seed <- if (exists(".Random.seed", envir = globalenv())) {
        get(".Random.seed", envir = globalenv())
      }
      on.exit({
        if (is.null(old_seed)) {
          rm(".Random.seed", envir = globalenv())
        } else {
          assign(".Random.seed", old_seed, envir = globalenv())
        }
      }, add = TRUE)
      set.seed(seed)
    }
    validation <- .sc_validate_windows(wh$hyperedges, wh$window_counts,
                                       wh$nodes, alpha, null = null,
                                       n_null = n_null)
    kept <- wh$hyperedges[validation$validated]
  } else {
    kept <- wh$hyperedges[wh$window_counts >= min_count]
  }
  keys <- unique(unlist(lapply(kept, \(set) {
    set <- sort(set)
    sizes <- seq_len(min(length(set), max_dim + 1L))
    unlist(lapply(sizes, \(k) {
      combos <- utils::combn(length(set), k, simplify = FALSE)
      vapply(combos, \(i) paste(set[i], collapse = ","), character(1L))
    }))
  })))
  keys <- as.character(keys)
  simplices <- lapply(strsplit(keys, ",", fixed = TRUE), as.integer)
  dims <- lengths(simplices)
  simplices <- simplices[order(dims, keys)]
  sc <- .sc_make_complex(simplices, wh$nodes, "window")
  sc$window <- as.integer(window)
  sc$min_count <- as.integer(min_count)
  if (validate) {
    sc$validation <- validation$table
    sc$alpha <- alpha
    sc$null <- null
    sc$n_null <- validation$n_null
  }
  # count of each simplex: the largest number of windows of any kept window
  # set that contains it, i.e. the highest min_count at which it is present
  kept_counts <- if (validate) wh$window_counts[validation$validated] else
    wh$window_counts[wh$window_counts >= min_count]
  sc$count <- .sc_face_counts(sc$simplices, kept, kept_counts, max_dim)
  # every action is a vertex of the complex at every min_count, so vertices
  # enter the count filtration at the top
  sc$count[lengths(sc$simplices) == 1L] <- max(sc$count)
  sc
}

.sc_face_counts <- function(simplices, sets, counts, max_dim) {
  if (length(sets) == 0L) return(numeric(length(simplices)))
  face_rows <- lapply(seq_along(sets), \(i) {
    set <- sort(sets[[i]])
    sizes <- seq_len(min(length(set), max_dim + 1L))
    keys <- unlist(lapply(sizes, \(k) {
      vapply(utils::combn(length(set), k, simplify = FALSE),
             \(j) paste(set[j], collapse = ","), character(1L))
    }))
    data.frame(key = keys, count = rep(counts[i], length(keys)),
               stringsAsFactors = FALSE)
  })
  faces <- do.call(rbind, face_rows)
  best <- tapply(faces$count, faces$key, max)
  simplex_keys <- vapply(simplices, \(s) paste(sort(s), collapse = ","),
                         character(1L))
  out <- as.numeric(best[simplex_keys])
  out[is.na(out)] <- 0
  out
}

# Persistent homology of a window complex over its count filtration: the
# complex at threshold t holds the simplices with count >= t, so classes are
# born at high counts and die at lower ones. Nestimate's reduction runs on an
# ascending filtration, so it receives max(count) - count and the result is
# mapped back to counts.
.ph_window <- function(x, max_dim) {
  counts <- x$count
  low <- min(counts)
  asc <- x
  class(asc) <- "simplicial_complex"
  asc$filtration <- max(counts) - counts
  ph <- Nestimate::persistent_homology(asc, max_dim = max_dim)
  pers <- ph$persistence
  pers$birth <- pers$birth + low
  pers$death <- pers$death + low
  pers$persistence <- pers$birth - pers$death
  thresholds <- sort(unique(counts), decreasing = TRUE)
  full_betti <- as.integer(hg_betti(x))
  full_betti <- c(full_betti, rep(0L, max(0L, max_dim + 1L - length(full_betti))))
  curve <- expand.grid(threshold = thresholds, dimension = 0:max_dim)
  curve$betti <- vapply(seq_len(nrow(curve)), \(i) {
    t <- curve$threshold[i]
    d <- curve$dimension[i]
    if (t == low) return(full_betti[d + 1L])
    sum(pers$dimension == d & pers$birth >= t & pers$death < t)
  }, integer(1L))
  structure(list(betti_curve = curve, persistence = pers,
                 thresholds = thresholds, mode = "clique"),
            class = "persistent_homology")
}

# Private copy of Nestimate's internal complex constructor
# (.make_simplicial_complex, deparse-identical; asserted in
# local_testing_and_equivalence/test-identity-nestimate-simplicial.R), so
# window complexes are the same object Nestimate builds.
.sc_make_complex <- function (simplices, nodes, type) 
{
    seen <- new.env(hash = TRUE, parent = emptyenv())
    for (s in simplices) seen[[paste(sort(s), collapse = ",")]] <- TRUE
    for (i in seq_along(nodes)) {
        if (is.null(seen[[as.character(i)]])) {
            simplices[[length(simplices) + 1L]] <- i
        }
    }
    dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))
    max_d <- if (length(dims) > 0L) 
        max(dims)
    else 0L
    f_vec <- vapply(0:max_d, function(d) sum(dims == d), integer(1))
    names(f_vec) <- paste0("dim_", 0:max_d)
    n_simplices <- length(simplices)
    n_nodes <- length(nodes)
    max_possible <- sum(vapply(0:max_d, function(d) {
        choose(n_nodes, d + 1L)
    }, numeric(1)))
    density <- if (max_possible > 0) 
        n_simplices/max_possible
    else 0
    mean_dim <- if (n_simplices > 0L) 
        mean(dims)
    else 0
    structure(list(simplices = simplices, nodes = nodes, n_nodes = n_nodes, 
        n_simplices = n_simplices, dimension = max_d, f_vector = f_vec, 
        density = density, mean_dim = mean_dim, type = type), 
        class = "simplicial_complex")
}



# The weight matrix of a network input (a matrix, or an object carrying one
# in `$weights`: netobject, tna, cograph_network, net_hon).
# Sequences (a long or wide data frame, a list of sequences, or long-format
# arguments) rather than a network: a numeric matrix or a network object
# is a network.
.sc_is_sequence_input <- function(x, action) {
  if (!is.null(action)) return(TRUE)
  if (inherits(x, c("netobject", "netobject_group", "tna", "cograph_network",
                    "net_hon", "net_hypa", "net_mogen"))) {
    return(FALSE)
  }
  if (is.matrix(x)) return(!is.numeric(x))
  is.data.frame(x) || is.list(x)
}

# Distances for a Vietoris-Rips complex: a square numeric matrix is taken
# as distances; a dist object, a non-square numeric matrix or table (rows are
# points), or a data frame of numeric coordinates becomes Euclidean
# distances between its rows.
.sc_distances <- function(x) {
  if (inherits(x, "simplicial_complex")) return(x)
  if (inherits(x, "dist")) return(as.matrix(x))
  if (is.data.frame(x)) {
    numeric_cols <- vapply(x, is.numeric, logical(1L))
    label_cols <- names(x)[vapply(x, \(v) is.character(v) || is.factor(v),
                                  logical(1L))]
    if (!any(numeric_cols) || length(label_cols) > 1L) {
      .thg_bad_input(paste0("points for type = \"vr\" need numeric ",
                            "coordinate columns and at most one label column"))
    }
    points <- as.matrix(x[numeric_cols])
    rownames(points) <- if (length(label_cols)) {
      as.character(x[[label_cols]])
    } else {
      rownames(x)
    }
    return(as.matrix(stats::dist(points)))
  }
  if (is.matrix(x) && is.numeric(x) && nrow(x) != ncol(x)) {
    points <- unclass(x)
    return(as.matrix(stats::dist(matrix(points, nrow(points),
                                        dimnames = dimnames(points)))))
  }
  x
}

.sc_weight_matrix <- function(x) {
  w <- if (is.matrix(x)) x else if (is.list(x)) x$weights
  if (!is.matrix(w) || !is.numeric(w) || nrow(w) != ncol(w)) {
    .thg_bad_input(paste0("`direction = \"both\"` needs a square weighted ",
                          "matrix or a network object carrying one"))
  }
  w
}
