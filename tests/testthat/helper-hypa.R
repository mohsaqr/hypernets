# hypa() returns Nestimate's HYPA fit with one added class (hypernets_hypa,
# for printing the selected paths) and the `hypa_view` attribute. Identity
# tests against the Nestimate builder compare everything else.
.unview <- function(x) {
  class(x) <- setdiff(class(x), "hypernets_hypa")
  attr(x, "hypa_view") <- NULL
  x
}

# simplicial() returns Nestimate's complex with one added class
# (hypernets_simplicial, whose plot() uses cograph); compare the rest.
.unclass_sc <- function(x) {
  class(x) <- setdiff(class(x), "hypernets_simplicial")
  x
}
