# Count Path Frequencies in Trajectory Data

Counts the frequency of k-step paths (k-grams) across all trajectories.
Useful for understanding which sequences dominate the data before
applying formal models.

## Usage

``` r
path_counts(data, k = 2L, top = NULL)
```

## Arguments

- data:

  A list of character vectors (trajectories) or a data.frame (rows =
  trajectories, columns = time points).

- k:

  Integer. Length of the path / n-gram (default 2). A k of 2 counts
  individual transitions; k of 3 counts two-step paths, etc. Must be a
  whole number; a non-integer value is an error rather than being
  silently truncated.

- top:

  Integer or `NULL`. Return only the first `top` rows of the
  count-sorted table (default `NULL` returns every path).

## Value

A data frame with columns: `path`, `count`, `proportion`.

## Examples

``` r
trajs <- list(c("A","B","C","D"), c("A","B","D","C"))
path_counts(trajs, k = 2)
#>     path count proportion
#> 1 A -> B     2     0.3333
#> 2 B -> C     1     0.1667
#> 3 B -> D     1     0.1667
#> 4 C -> D     1     0.1667
#> 5 D -> C     1     0.1667

# \donttest{
path_counts(trajs, k = 3, top = 10)
#>          path count proportion
#> 1 A -> B -> C     1       0.25
#> 2 A -> B -> D     1       0.25
#> 3 B -> C -> D     1       0.25
#> 4 B -> D -> C     1       0.25
# }
```
