# Drawing a hypergraph: hypernets computes the layout of the clique projection
# and the member sets, and cograph draws them -- every hyperedge as a smooth
# translucent blob around its members (cograph::plot_simplicial()), the
# convention HypergraphX uses for the decision and tribunal figures of
# Coupette et al. (2024). hypernets adds only what cograph has no notion of:
# the mapping of a hyperedge attribute to blob colour or line type, and its
# legend.

# cograph's blob canvas is fixed at x in [-9, 9], y in [-8.5, 8.5] with a
# blob radius of one unit, so the padding around a hyperedge's members is
# set by how much of the canvas the layout spans.
.thg_blob_canvas <- c(x = 9, y = 8.5)

# Node coordinates in [0, 1]^2 from a layout name or a user table.
.thg_layout <- function(hg, layout, seed) {
  if (is.data.frame(layout)) {
    if (!all(c("node", "x", "y") %in% names(layout))) {
      .thg_bad_input("a `layout` data.frame needs `node`, `x` and `y` columns")
    }
    missing <- setdiff(hg$nodes, as.character(layout$node))
    if (length(missing)) {
      .thg_bad_input(sprintf("`layout` lacks coordinates for %d node(s)",
                             length(missing)))
    }
    idx <- match(hg$nodes, as.character(layout$node))
    xy <- cbind(x = as.numeric(layout$x[idx]), y = as.numeric(layout$y[idx]))
  } else {
    layout <- match.arg(layout, c("spring", "circle"))
    projection <- as.matrix(hg_project(hg, method = "clique",
                                       weighted = FALSE, what = "matrix"))
    net <- cograph::as_cograph(projection, directed = FALSE)
    pos <- if (identical(layout, "spring")) {
      cograph::layout_spring(net, seed = seed)
    } else {
      cograph::layout_circle(net)
    }
    xy <- cbind(x = as.numeric(pos$x), y = as.numeric(pos$y))
  }
  span <- apply(xy, 2L, function(z) diff(range(z)))
  span[span < sqrt(.Machine$double.eps)] <- 1
  xy <- sweep(xy, 2L, apply(xy, 2L, min))
  xy <- sweep(xy, 2L, span, "/")
  data.frame(node = hg$nodes, x = xy[, 1L], y = xy[, 2L],
             stringsAsFactors = FALSE)
}

# One value per hyperedge from a selector: NULL, "size", an edge_data
# column, a named vector keyed by hyperedge, or a vector of length m.
.thg_edge_aesthetic <- function(hg, by, arg) {
  if (is.null(by)) return(NULL)
  edge_names <- colnames(hg$incidence)
  m <- hg$n_hyperedges
  if (is.character(by) && length(by) == 1L) {
    if (identical(by, "size")) return(lengths(hg$hyperedges))
    if (!is.null(hg$edge_data) && by %in% names(hg$edge_data)) {
      values <- hg$edge_data[[by]]
      return(values[match(edge_names, as.character(hg$edge_data$edge))])
    }
    if (m != 1L) {
      .thg_bad_input(sprintf(
        "`%s` must be \"size\", a column of the edge metadata, or one value per hyperedge",
        arg
      ))
    }
  }
  if (!is.null(names(by))) {
    missing <- setdiff(edge_names, names(by))
    if (length(missing)) {
      .thg_bad_input(sprintf("`%s` has no value for %d hyperedge(s)", arg,
                             length(missing)))
    }
    return(unname(by[edge_names]))
  }
  if (length(by) != m) {
    .thg_bad_input(sprintf("`%s` must have one value per hyperedge (%d)", arg, m))
  }
  by
}

# Colours for a fill selector: a named palette over the levels of a
# discrete variable, or a ramp over the range of a numeric one.
.thg_fill_colours <- function(values) {
  if (is.numeric(values)) {
    rng <- range(values)
    unit <- if (diff(rng) > 0) (values - rng[1L]) / diff(rng) else rep(0.5, length(values))
    ramp <- grDevices::colorRamp(.thg_okabe_ito_ramp())
    rgb <- ramp(unit)
    return(list(colours = grDevices::rgb(rgb[, 1L], rgb[, 2L], rgb[, 3L],
                                         maxColorValue = 255),
                palette = NULL))
  }
  levels_fill <- if (is.factor(values)) levels(values) else
    sort(unique(as.character(values)))
  palette <- if (length(levels_fill) <= length(.thg_okabe_ito)) {
    .thg_okabe_ito[seq_along(levels_fill)]
  } else {
    grDevices::colorRampPalette(.thg_okabe_ito_ramp())(length(levels_fill))
  }
  names(palette) <- levels_fill
  list(colours = unname(palette[as.character(values)]), palette = palette)
}

.thg_linetypes <- c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")

#' Plot a hypergraph with hyperedges as blobs
#'
#' Places the nodes with a graph layout of the clique projection and draws
#' every hyperedge as a smooth translucent blob around its members, through
#' [cograph::plot_simplicial()]. Blobs can be coloured by hyperedge size, by
#' a column of the edge metadata, or by any value supplied per hyperedge, and
#' outlined with a line type by the same kinds of selector; hypernets adds the
#' legend for that mapping. This is the drawing convention used for the
#' decision and tribunal figures of Coupette et al. (2024). A hyperedge with
#' a single member is drawn as its node alone.
#'
#' @param x A `net_hypergraph` with at least one hyperedge of two or more
#'   members. Draw a part of a large hypergraph by passing [hg_subset()]
#'   first.
#' @param layout `"spring"` (default; [cograph::layout_spring()] on the
#'   clique projection), `"circle"`, or a data.frame with `node`, `x` and
#'   `y` columns to reuse coordinates across panels.
#' @param seed Seed for the spring layout (default `1`).
#' @param color_by Fill of the blobs: `NULL` (one colour), `"size"`
#'   (hyperedge cardinality, sequential scale), the name of a column in the
#'   edge metadata (`x$edge_data`), a vector named by hyperedge, or a vector
#'   with one value per hyperedge. Character or factor values get the
#'   Okabe-Ito palette; numeric values a sequential scale built from it.
#' @param linetype_by Outline line type of the blobs, resolved like
#'   `color_by`; discrete only.
#' @param labels `TRUE` (default) writes the node names, `FALSE` writes
#'   none, and a character vector named by node replaces the names shown.
#' @param label_size,node_size Text and point sizes.
#' @param alpha Fill transparency of the blobs (default `0.45`).
#' @param padding Blob radius around the member positions as a fraction of
#'   the layout extent (default `0.06`, the radius cograph draws when the
#'   layout fills its canvas; smaller values enlarge the canvas).
#' @param legend_title Legend title for `color_by`; defaults to the
#'   selector's name.
#' @param ... Unused; for S3 consistency.
#' @return A ggplot object.
#' @references Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal
#'   hypergraphs. *Philosophical Transactions of the Royal Society A*,
#'   382(2270), 20230141. \doi{10.1098/rsta.2023.0141}
#' @examples
#' dat <- data.frame(
#'   member = c("a", "b", "c", "b", "c", "d", "d", "e", "f"),
#'   event = c("e1", "e1", "e1", "e2", "e2", "e2", "e3", "e3", "e4")
#' )
#' hg <- group_hypergraph(dat, "member", "event")
#' plot(hg, color_by = "size")
#' @export
plot.net_hypergraph <- function(x, layout = c("spring", "circle"), seed = 1L,
                                color_by = NULL, linetype_by = NULL,
                                labels = TRUE, label_size = 3,
                                node_size = 2.5, alpha = 0.45,
                                padding = 0.06, legend_title = NULL, ...) {
  .thg_check_hg(x)
  stopifnot(
    "`alpha` must be a number in (0, 1]" =
      is.numeric(alpha) && length(alpha) == 1L && alpha > 0 && alpha <= 1,
    "`padding` must be a positive number" =
      is.numeric(padding) && length(padding) == 1L && padding > 0
  )
  if (x$n_nodes == 0L) .thg_bad_input("the hypergraph has no nodes to draw")
  drawable <- lengths(x$hyperedges) >= 2L
  if (!any(drawable)) {
    .thg_bad_input("nothing to draw: no hyperedge has two or more members")
  }
  edge_names <- colnames(x$incidence) %||% paste0("h", seq_len(x$n_hyperedges))
  fill <- .thg_edge_aesthetic(x, color_by, "color_by")
  ltype <- .thg_edge_aesthetic(x, linetype_by, "linetype_by")
  if (!is.null(ltype) && is.numeric(ltype)) ltype <- as.character(ltype)

  # Layout on cograph's canvas: centred, spanning 1 / padding units, so a
  # blob radius of one unit is `padding` of the extent.
  pos <- .thg_layout(x, layout, seed)
  span <- 1 / padding
  coords <- cbind(x = (pos$x - 0.5) * span, y = (pos$y - 0.5) * span)
  rownames(coords) <- pos$node
  sets <- lapply(x$hyperedges[drawable], function(members) x$nodes[members])

  fill_map <- if (is.null(fill)) NULL else .thg_fill_colours(fill)
  blob_colours <- if (is.null(fill_map)) rep(.thg_okabe_ito[[5L]], sum(drawable)) else
    fill_map$colours[drawable]
  ltype_levels <- if (is.null(ltype)) NULL else
    (if (is.factor(ltype)) levels(ltype) else sort(unique(as.character(ltype))))
  ltype_map <- if (is.null(ltype)) NULL else
    stats::setNames(rep_len(.thg_linetypes, length(ltype_levels)), ltype_levels)
  blob_linetypes <- if (is.null(ltype_map)) "solid" else
    unname(ltype_map[as.character(ltype)[drawable]])

  shown <- if (isTRUE(labels)) {
    pos$node
  } else if (isFALSE(labels)) {
    rep("", nrow(pos))
  } else {
    if (is.null(names(labels))) .thg_bad_input("`labels` must be TRUE, FALSE or a vector named by node")
    replacement <- unname(labels[pos$node])
    ifelse(is.na(replacement), pos$node, as.character(replacement))
  }

  projection <- as.matrix(hg_project(x, method = "clique", weighted = FALSE,
                                     what = "matrix"))
  node_fill <- "#E2E0DD"
  # cograph::plot_simplicial() prints as it draws. Its print goes to a null
  # device here, and the legend is added before the caller prints once.
  grDevices::pdf(NULL)
  p <- tryCatch(
    cograph::plot_simplicial(
      projection, pathways = sets, layout = coords,
      blob_colors = blob_colours, blob_linetype = blob_linetypes,
      blob_alpha = alpha, node_color = node_fill, target_color = node_fill,
      ring_color = node_fill, node_size = node_size, label_size = label_size,
      labels = shown, label_color = "black", shadow = FALSE, title = ""
    ),
    finally = grDevices::dev.off()
  )
  p <- p + ggplot2::labs(title = NULL)
  if (span > 2 * min(.thg_blob_canvas)) {
    # A layout wider than cograph's canvas needs wider limits; ggplot2 says
    # so with a message when a coordinate system replaces another, which is
    # exactly what is intended here.
    half <- span / 2 + 1
    p <- suppressMessages(
      p + ggplot2::coord_equal(clip = "off", xlim = c(-half, half),
                               ylim = c(-half, half))
    )
  }

  title <- legend_title %||%
    (if (is.character(color_by) && length(color_by) == 1L) color_by else "value")
  if (!is.null(fill_map)) {
    key <- data.frame(x = coords[1L, 1L], y = coords[1L, 2L],
                      value = if (is.null(fill_map$palette)) fill else
                        factor(as.character(fill), levels = names(fill_map$palette)))
    p <- p + ggplot2::geom_point(
      data = key, ggplot2::aes(x = .data$x, y = .data$y, fill = .data$value),
      shape = 22, size = 0, alpha = 0
    )
    p <- p + if (is.null(fill_map$palette)) {
      ggplot2::scale_fill_gradientn(colours = .thg_okabe_ito_ramp(), name = title,
                                    guide = ggplot2::guide_colourbar())
    } else {
      ggplot2::scale_fill_manual(values = fill_map$palette, name = title,
                                 drop = FALSE,
                                 guide = ggplot2::guide_legend(override.aes = list(
                                   size = 5, alpha = 0.6, colour = NA
                                 )))
    }
  }
  if (!is.null(ltype_map)) {
    # an invisible zero-length path per level carries the line-type legend
    key <- data.frame(x = coords[1L, 1L], y = coords[1L, 2L],
                      value = factor(rep(ltype_levels, each = 2L),
                                     levels = ltype_levels))
    p <- p + ggplot2::geom_path(
      data = key, ggplot2::aes(x = .data$x, y = .data$y, group = .data$value,
                               linetype = .data$value),
      alpha = 0
    ) +
      ggplot2::scale_linetype_manual(
        values = ltype_map, drop = FALSE,
        name = if (is.character(linetype_by) && length(linetype_by) == 1L)
          linetype_by else "line type"
      ) +
      ggplot2::guides(linetype = ggplot2::guide_legend(
        override.aes = list(alpha = 1, colour = "black", linewidth = 0.8)
      ))
  }
  p + ggplot2::theme(legend.position = "right")
}

# Sequential ramp built from Okabe-Ito colours: light yellow through green to
# blue, matching the yellow-to-blue ordering of the paper's cardinality scale.
.thg_okabe_ito_ramp <- function() {
  c("#F0E442", "#009E73", "#0072B2")
}
