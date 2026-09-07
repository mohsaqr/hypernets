# Build a Higher-Order Network (HON)

Constructs a Higher-Order Network from sequential data, faithfully
implementing the BuildHON algorithm (Xu, Wickramarathne & Chawla, 2016).

The algorithm detects when a first-order Markov model is insufficient to
capture sequential dependencies and automatically creates higher-order
nodes. Uses KL-divergence to determine whether extending a node's
history provides significantly different transition distributions.

## Usage

``` r
build_hon(
  data,
  max_order = 5L,
  min_freq = 1L,
  collapse_repeats = FALSE,
  method = "hon+",
  action = NULL,
  actor = NULL,
  time = NULL
)
```

## Arguments

- data:

  One of:

  - `data.frame`: rows are trajectories, columns are time steps.
    Trailing `NA`s are stripped. All non-NA values are coerced to
    character.

  - `list`: each element is a character (or coercible) vector
    representing one trajectory.

  - `tna`: a tna object with sequence data. Numeric state IDs are
    automatically converted to label names.

  - `netobject`: a netobject with sequence data.

  - long `data.frame`: one event per row, together with `action` (and
    optionally `actor` and `time`).

- max_order:

  Integer. Maximum order of the HON. Default 5. The algorithm may
  produce lower-order nodes if the data do not justify higher orders.

- min_freq:

  Integer. Minimum frequency for a transition to be considered.
  Transitions observed fewer than `min_freq` times are treated as zero.
  Default 1.

- collapse_repeats:

  Logical. If `TRUE`, adjacent duplicate states within each trajectory
  are collapsed before analysis. Default `FALSE`.

- method:

  Character. Algorithm to use: `"hon+"` (default, parameter-free
  BuildHON+ with lazy observation building and MaxDivergence pruning) or
  `"hon"` (original BuildHON with eager observation building).

- action, actor, time:

  Long-format column names: `action` holds the categorical state of each
  event, `actor` groups events into trajectories (one per actor; `NULL`
  treats all rows as one), and `time` orders events within an actor (row
  order when `NULL`). Leave all three `NULL` for wide, list, or model
  input. Same interface as
  [`bootstrap_hon()`](https://mohsaqr.github.io/hypernets/reference/bootstrap_hon.md).

## Value

An S3 object of class `"net_hon"` containing:

- matrix:

  Weighted adjacency matrix (rows = from, cols = to). Rows and columns
  use readable arrow notation (e.g., `"A -> B"`).

- edges:

  Data frame with columns: `path` (full state sequence, e.g., "A -\> B
  -\> C"), `from` (context/conditioning states), `to` (predicted next
  state), `count` (raw frequency), `probability` (transition
  probability), `from_order`, `to_order`.

- nodes:

  data.frame with columns `id`, `label`, `name` (one row per HON node;
  `label`/`name` are the arrow-notation node names). Stored as a
  data.frame for `cograph_network` compatibility.

- n_nodes:

  Number of HON nodes.

- n_edges:

  Number of edges.

- first_order_states:

  Character vector of unique original states.

- max_order_requested:

  The `max_order` parameter used.

- max_order_observed:

  Highest order actually present.

- min_freq:

  The `min_freq` parameter used.

- n_trajectories:

  Number of trajectories after parsing.

- directed:

  Logical. Always `TRUE`.

## Details

**Node naming convention**: Higher-order nodes use readable arrow
notation. A first-order node is simply `"A"`. A second-order node
representing the context "came from A, now at B" is `"A -> B"`.
Third-order: `"A -> B -> C"`, etc.

**Algorithm overview**:

1.  Count all subsequence transitions up to `max_order + 1`.

2.  Build probability distributions, filtering by `min_freq`.

3.  For each first-order source, recursively test whether extending the
    history (adding more context) produces a significantly different
    distribution (via KL-divergence vs. an adaptive threshold).

4.  Build the network from the accepted rules, rewiring edges so
    higher-order nodes are properly connected.

## References

Xu, J., Wickramarathne, T. L., & Chawla, N. V. (2016). Representing
higher-order dependencies in networks. *Science Advances*, 2(5),
e1600028.

Saebi, M., Xu, J., Kaplan, L. M., Ribeiro, B., & Chawla, N. V. (2020).
Efficient modeling of higher-order dependencies in networks: from
algorithm to application for anomaly detection. *EPJ Data Science*,
9(1), 15.

## Examples

``` r
seqs <- list(c("A","B","C","D"), c("A","B","C","A"), c("B","C","D","A"))
hon <- build_hon(seqs, max_order = 2)

# \donttest{
# From list of trajectories
trajs <- list(
  c("A", "B", "C", "D", "A"),
  c("A", "B", "D", "C", "A"),
  c("A", "B", "C", "D", "A")
)
hon <- build_hon(trajs, max_order = 3, min_freq = 1)
print(hon)
#> Higher-Order Network (HON)
#>   Nodes:        4 (4 first-order states)
#>   Edges:        7
#>   Max order:    1 (requested 3)
#>   Min freq:     1
#>   Trajectories: 3
summary(hon)
#> Higher-Order Network (HON) Summary
#>   Nodes: 4 | Edges: 7 | Trajectories: 3
#>   First-order states: A, B, C, D
#>   Max order observed: 1 (requested: 3)
#>   Min frequency: 1
#>   Node order distribution:
#>     Order 1: 4 nodes

# From data.frame (rows = trajectories)
df <- data.frame(T1 = c("A", "A"), T2 = c("B", "B"),
                 T3 = c("C", "D"), T4 = c("D", "C"))
hon <- build_hon(df, max_order = 2)
# }
```
