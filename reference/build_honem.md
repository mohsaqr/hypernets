# Build HONEM Embeddings for Higher-Order Networks

Constructs low-dimensional embeddings from a Higher-Order Network (HON)
that preserve higher-order dependencies. Uses exponentially-decaying
matrix powers of the HON transition matrix followed by truncated SVD.

## Usage

``` r
build_honem(hon, dim = 32L, max_power = 10L)
```

## Arguments

- hon:

  A `net_hon` object from
  [`build_hon`](https://mohsaqr.github.io/hypernets/reference/build_hon.md),
  or a square weighted adjacency matrix.

- dim:

  Integer. Embedding dimension (default 32).

- max_power:

  Integer. Maximum walk length for neighborhood computation (default
  10). Higher values capture longer-range structure.

## Value

An object of class `net_honem` with components:

- embeddings:

  Numeric matrix (n_nodes x dim) of node embeddings.

- nodes:

  Character vector of node names.

- singular_values:

  Numeric vector of top singular values.

- explained_variance:

  Proportion of variance explained.

- dim:

  Embedding dimension used.

- max_power:

  Maximum power used.

- n_nodes:

  Number of nodes embedded.

## Details

HONEM is parameter-free and scalable - no random walks, skip-gram, or
hyperparameter tuning required.

## References

Saebi, M., Ciampaglia, G. L., Kaplan, L. M., & Chawla, N. V. (2020).
HONEM: Learning Embedding for Higher Order Networks. *Big Data*, 8(4),
255-269.
[doi:10.1089/big.2019.0169](https://doi.org/10.1089/big.2019.0169)

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon_2 <- build_hon(seqs, max_order = 2)
hem <- build_honem(hon_2, dim = 2)

# \donttest{
trajs <- list(c("A","B","C","D"), c("A","B","D","C"),
              c("B","C","D","A"), c("C","D","A","B"))
hon <- build_hon(trajs, max_order = 2)
emb <- build_honem(hon, dim = 4)
print(emb)
#> HONEM: Higher-Order Network Embedding
#>   Nodes:      4
#>   Dimensions: 3
#>   Max power:  10
#>   Variance explained: 94.6%
plot(emb)

# }
```
