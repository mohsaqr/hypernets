# ---- Filtration + boundary reduction (simplicial family) -------------------
#
# The shared layer under BOTH build_simplicial(type = "vr"), which
# attaches a filtration to the complex it returns, and
# persistent_homology(), which reduces that filtration's boundary
# matrix over Z/2. Split out of simplicial.R in hypernets 0.2.0; the code
# is unchanged.

#' @noRd
.align_filtration <- function(fc, sc) {
  # .make_simplicial_complex may reorder simplices and add isolated vertices;
  # re-key the filtration vector to match sc$simplices order.
  fc_keys <- fc$key
  sc_keys <- vapply(sc$simplices, function(s) paste(sort(s), collapse = ","),
                    character(1))
  out <- numeric(length(sc_keys))
  m <- match(sc_keys, fc_keys)
  out[!is.na(m)] <- fc$filt_asc[m[!is.na(m)]]
  out[is.na(m)] <- 0  # isolated vertices added by .make_simplicial_complex
  out
}

# =========================================================================
# Clique complex - verified against igraph::cliques()
# =========================================================================

#' @noRd
.filter_clique_complex <- function(mat, max_dim = 3L) {
  # Clique filtration in similarity (descending) semantics.
  # filt_asc(sigma) = max_w - min(edge weights in sigma), vertex = 0.
  # So a high-weight simplex has a small filt_asc (enters early in ascending).
  n <- nrow(mat)
  nodes <- rownames(mat) %||% paste0("V", seq_len(n))
  max_w <- max(mat)

  adj <- mat > 0
  diag(adj) <- FALSE

  if (!any(adj)) {
    simps <- as.list(seq_len(n))
    return(list(
      simplices = simps, dim = integer(n), filt_asc = numeric(n),
      key = as.character(seq_len(n)), nodes = nodes,
      max_filt = 0, max_w = max_w, mode = "clique"
    ))
  }

  all_simp <- .find_all_cliques(adj, max_dim)
  dims <- vapply(all_simp, function(s) length(s) - 1L, integer(1))
  filt <- vapply(seq_along(all_simp), function(j) {
    s <- all_simp[[j]]
    if (length(s) == 1L) return(0)
    pairs <- utils::combn(s, 2L)
    max_w - min(mat[cbind(pairs[1L, ], pairs[2L, ])])
  }, numeric(1))

  ord <- order(filt, dims)
  simplices <- all_simp[ord]
  dims <- dims[ord]
  filt <- filt[ord]
  keys <- vapply(simplices, function(s) paste(s, collapse = ","), character(1))

  list(
    simplices = simplices, dim = dims, filt_asc = filt,
    key = keys, nodes = nodes,
    max_filt = if (length(filt) > 0L) max(filt) else 0,
    max_w = max_w, mode = "clique"
  )
}

#' @noRd
.filter_vr_complex <- function(d, max_dim = 3L, max_scale = NULL) {
  # Vietoris-Rips filtration on a non-negative distance matrix.
  # filt(sigma) = max pairwise distance in sigma; vertex = 0.
  stopifnot(is.matrix(d), nrow(d) == ncol(d))
  n <- nrow(d)
  nodes <- rownames(d) %||% paste0("V", seq_len(n))
  d <- pmax(d, t(d))
  diag(d) <- 0
  if (any(d < 0, na.rm = TRUE)) {
    stop("VR filtration requires a non-negative distance matrix.",
         call. = FALSE)
  }
  finite_d <- d
  finite_d[!is.finite(finite_d)] <- NA_real_
  cap <- if (is.null(max_scale)) {
    if (all(is.na(finite_d))) 0 else max(finite_d, na.rm = TRUE)
  } else {
    stopifnot(is.numeric(max_scale), length(max_scale) == 1L,
              !is.na(max_scale), max_scale >= 0)
    max_scale
  }

  # Edges within cap. d(i,j) == 0 for i != j is a valid pseudometric case
  # (duplicate points / equivalence classes), so include zero-distance
  # off-diagonal edges and rely on diag(adj) <- FALSE to exclude self-loops.
  adj <- !is.na(finite_d) & finite_d >= 0 & finite_d <= cap
  diag(adj) <- FALSE

  all_simp <- .find_all_cliques(adj, max_dim)
  dims <- vapply(all_simp, function(s) length(s) - 1L, integer(1))
  filt <- vapply(seq_along(all_simp), function(j) {
    s <- all_simp[[j]]
    if (length(s) == 1L) return(0)
    pairs <- utils::combn(s, 2L)
    max(d[cbind(pairs[1L, ], pairs[2L, ])])
  }, numeric(1))

  ord <- order(filt, dims)
  simplices <- all_simp[ord]
  dims <- dims[ord]
  filt <- filt[ord]
  keys <- vapply(simplices, function(s) paste(s, collapse = ","), character(1))

  list(
    simplices = simplices, dim = dims, filt_asc = filt,
    key = keys, nodes = nodes,
    max_filt = if (length(filt) > 0L) max(filt) else 0,
    max_w = cap, mode = "vr"
  )
}

#' @noRd
.fc_from_filtered_complex <- function(sc, max_dim = 3L, max_scale = NULL) {
  # Convert a simplicial_complex with attached $filtration into the internal
  # filtered-complex shape consumed by .persistence_pairs_z2(). Honors max_dim
  # by dropping simplices above that dimension and re-orders by (filt, dim) so
  # boundary reduction is well-defined.
  stopifnot(inherits(sc, "net_simplicial"),
            !is.null(sc$filtration),
            length(sc$filtration) == length(sc$simplices))
  simplices <- sc$simplices
  filt <- as.numeric(sc$filtration)
  dims <- vapply(simplices, function(s) length(s) - 1L, integer(1))

  # Drop above max_dim
  keep <- dims <= max_dim
  simplices <- simplices[keep]
  filt <- filt[keep]
  dims <- dims[keep]

  # Apply max_scale cap if requested (and the complex is a VR build)
  mode <- if (identical(sc$type, "vr")) "vr" else "clique"
  if (!is.null(max_scale)) {
    stopifnot(is.numeric(max_scale), length(max_scale) == 1L,
              !is.na(max_scale), max_scale >= 0)
    keep <- filt <= max_scale
    simplices <- simplices[keep]
    filt <- filt[keep]
    dims <- dims[keep]
  }

  # Order by (filt asc, dim asc) so faces precede cofaces at the same filt
  ord <- order(filt, dims)
  simplices <- simplices[ord]
  filt <- filt[ord]
  dims <- dims[ord]
  keys <- vapply(simplices, function(s) paste(sort(s), collapse = ","),
                 character(1))

  max_w <- if (mode == "vr") {
    if (!is.null(sc$max_scale)) sc$max_scale
    else if (length(filt) > 0L) max(filt) else 0
  } else {
    if (length(filt) > 0L) max(filt) else 0
  }

  list(
    simplices = simplices, dim = dims, filt_asc = filt,
    key = keys, nodes = sc$nodes,
    max_filt = if (length(filt) > 0L) max(filt) else 0,
    max_w = max_w, mode = mode
  )
}

#' @noRd
.persistence_pairs_z2 <- function(fc) {
  # Standard left-to-right boundary-matrix reduction over Z/2.
  # The j-loop is sequential by construction (column j depends on reduced
  # columns 1..j-1) - this is the package's second documented for-loop
  # exception alongside the permutation loop in sequence_compare.R.
  simplices <- fc$simplices
  dims <- fc$dim
  keys <- fc$key
  filt <- fc$filt_asc
  N <- length(simplices)

  if (N == 0L) {
    empty <- data.frame(dimension = integer(0), birth = numeric(0),
                        death = numeric(0), persistence = numeric(0),
                        stringsAsFactors = FALSE)
    return(list(pairs = empty, essential = empty))
  }

  key_to_idx <- setNames(seq_len(N), keys)

  # Boundary columns: k-simplex (k>=1) has (k+1) (k-1)-faces.
  D <- lapply(seq_len(N), function(j) {
    if (dims[j] < 1L) return(integer(0))
    s <- simplices[[j]]
    face_keys <- vapply(seq_along(s), function(i) {
      paste(s[-i], collapse = ",")
    }, character(1))
    sort.int(as.integer(key_to_idx[face_keys]))
  })

  low_to_col <- integer(N) # low_to_col[r] = column j with low(j)=r, 0 if none
  paired_b <- integer(0L)
  paired_d <- integer(0L)

  for (j in seq_len(N)) {
    col <- D[[j]]
    while (length(col) > 0L) {
      l <- col[length(col)]
      i <- low_to_col[l]
      if (i == 0L) break
      # XOR with D[[i]]: symmetric difference of sorted integer vectors
      col <- sort.int(c(setdiff(col, D[[i]]), setdiff(D[[i]], col)))
    }
    D[[j]] <- col
    if (length(col) > 0L) {
      l <- col[length(col)]
      low_to_col[l] <- j
      paired_b <- c(paired_b, l)
      paired_d <- c(paired_d, j)
    }
  }

  essential_idx <- setdiff(seq_len(N), c(paired_b, paired_d))

  pairs_df <- if (length(paired_b) == 0L) {
    data.frame(dimension = integer(0), birth = numeric(0),
               death = numeric(0), persistence = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    data.frame(
      dimension = dims[paired_b],
      birth = filt[paired_b],
      death = filt[paired_d],
      persistence = filt[paired_d] - filt[paired_b],
      stringsAsFactors = FALSE
    )
  }
  pairs_df <- pairs_df[pairs_df$persistence > 0, , drop = FALSE]

  essential_df <- if (length(essential_idx) == 0L) {
    data.frame(dimension = integer(0), birth = numeric(0),
               death = numeric(0), persistence = numeric(0),
               stringsAsFactors = FALSE)
  } else {
    data.frame(
      dimension = dims[essential_idx],
      birth = filt[essential_idx],
      death = Inf, persistence = Inf,
      stringsAsFactors = FALSE
    )
  }
  list(pairs = pairs_df, essential = essential_df)
}

#' @noRd
.betti_curve_from_pairs <- function(pers, thresholds, max_dim, mode) {
  grid <- expand.grid(threshold = thresholds, dimension = 0:max_dim,
                      KEEP.OUT.ATTRS = FALSE)
  if (nrow(pers) == 0L) {
    grid$betti <- 0L
    return(grid[, c("threshold", "dimension", "betti")])
  }
  # Clique mode: descending thresholds; alive at t iff birth >= t AND death < t.
  # VR mode: ascending thresholds; alive at t iff birth <= t AND death > t.
  alive_per_row <- if (mode == "clique") {
    vapply(seq_len(nrow(grid)), function(k) {
      t <- grid$threshold[k]
      sub <- pers[pers$dimension == grid$dimension[k], , drop = FALSE]
      sum(sub$birth >= t & sub$death < t)
    }, integer(1))
  } else {
    vapply(seq_len(nrow(grid)), function(k) {
      t <- grid$threshold[k]
      sub <- pers[pers$dimension == grid$dimension[k], , drop = FALSE]
      sum(sub$birth <= t & sub$death > t)
    }, integer(1))
  }
  grid$betti <- alive_per_row
  grid[, c("threshold", "dimension", "betti")]
}

