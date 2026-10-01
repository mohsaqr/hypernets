# The pathway strings of a memory network are one table of hg_get():
# hg_get(x, what = "pathways"). The pathway tests assert on the character
# vector of that table's single column, so they read it through this helper.
.pathways <- function(x, ...) {
  hg_get(x, what = "pathways", ...)$pathway
}
