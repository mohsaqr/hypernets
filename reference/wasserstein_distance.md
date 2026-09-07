# Wasserstein Distance Between Persistence Diagrams

Computes the finite-order Wasserstein distance between persistence
diagrams. For order \\q\\, it minimizes the sum of powered matching
costs over bijections between diagram points and copies of the diagonal,
then takes the \\q\\-th root. Point-to-point costs use an \\L_p\\ ground
metric; `internal_p = Inf` gives the usual \\L\_\infty\\ convention used
by GUDHI.

## Usage

``` r
wasserstein_distance(d1, d2, dimension = NULL, order = 1, internal_p = Inf)
```

## Arguments

- d1, d2:

  `net_persistent_homology` objects, or data.frames with columns
  `dimension`, `birth`, `death`.

- dimension:

  Integer vector of dimensions to compare. `NULL` (default) compares all
  dimensions appearing in either diagram and returns a named numeric
  vector.

- order:

  Finite Wasserstein order, a number greater than or equal to one. The
  common choices are `1` and `2`.

- internal_p:

  Ground-metric order, a number greater than or equal to one or `Inf`
  (the default).

## Value

Named numeric vector, one value per requested dimension. Names are
`"dim_<k>"`.

## Details

Finite points may match the diagonal. Essential classes (`death = Inf`
in Vietoris–Rips mode or `death = 0` in clique mode) are matched only to
essential classes. A dimension whose essential counts differ has
distance `Inf`. The finite assignment is solved exactly with a native
Hungarian algorithm, so no optional optimization package is required.

## References

Kerber, M., Morozov, D., & Nigmetov, A. (2017). Geometry helps to
compare persistence diagrams. *Journal of Experimental Algorithmics*,
22, 1–20. [doi:10.1145/3064175](https://doi.org/10.1145/3064175)

## Examples

``` r
d1 <- data.frame(dimension = 0L, birth = 0, death = 2)
d2 <- data.frame(dimension = 0L, birth = 0, death = 3)
wasserstein_distance(d1, d2)
#> dim_0 
#>     1 
```
