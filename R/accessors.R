# ---------------------------------------------------------------------------
# Shared accessor plumbing
#
# Every hypernets accessor that can return many rows takes `top =`. It is
# applied LAST -- after `what`, after every filter (`order_min`, `min_count`,
# `dim`, `k`, `dimension`, `significant`, ...) and after `sort_by` -- so
# `sort_by` and `top` compose: `top = n` means "the first n rows of the table
# as it would otherwise have been returned". Semantics match
# `path_counts(top =)` and `pathways(top =)`, which shipped first.
#
# This exists so that a caller never has to write head() on the public
# surface. A vignette line that subsets, sorts or filters a returned table is
# a defect; the accessor owns the argument instead.
# ---------------------------------------------------------------------------

#' Truncate a tidy accessor result to its first `top` rows
#'
#' @param x A data.frame, already filtered and ordered.
#' @param top Integer or NULL. NULL returns `x` unchanged.
#' @return `x`, or its first `top` rows with row names reset.
#' @noRd
.ho_top <- function(x, top) {
  if (is.null(top)) return(x)
  stopifnot(
    "`top` must be a single whole number >= 1" =
      is.numeric(top) && length(top) == 1L && is.finite(top) &&
      top == round(top) && top >= 1
  )
  out <- utils::head(x, as.integer(top))
  rownames(out) <- NULL
  out
}
