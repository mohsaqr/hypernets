# Planted communities of the ring sequences

The true community of every node a second-order network built from
[ring_sequences](https://mohsaqr.github.io/hypernets/reference/ring_sequences.md)
can hold. A node `"u -> v"` belongs to the one group holding both `u`
and `v`; the first-order node of an action that is not shared (`"p1_1"`)
belongs to that action's group. First-order nodes of the shared actions
(`"s1"`, ...) have no single group and are not listed.

## Usage

``` r
ring_communities
```

## Format

A data frame with 56 rows (one per node) and 2 columns:

- node:

  Character. A node label as
  [`hon()`](https://mohsaqr.github.io/hypernets/reference/hon.md) writes
  it.

- community:

  Integer. The planted group, 1 to 4.

## Source

Built with
[ring_sequences](https://mohsaqr.github.io/hypernets/reference/ring_sequences.md)
by `data-raw/ring_sequences.R`.

## Examples

``` r
ring_communities
#>            node community
#> 1          p1_1         1
#> 2  p1_1 -> p1_2         1
#> 3    p1_1 -> s1         1
#> 4    p1_1 -> s2         1
#> 5          p1_2         1
#> 6  p1_2 -> p1_1         1
#> 7    p1_2 -> s1         1
#> 8    p1_2 -> s2         1
#> 9    s1 -> p1_1         1
#> 10   s1 -> p1_2         1
#> 11     s1 -> s2         1
#> 12   s2 -> p1_1         1
#> 13   s2 -> p1_2         1
#> 14     s2 -> s1         1
#> 15         p2_1         2
#> 16 p2_1 -> p2_2         2
#> 17   p2_1 -> s2         2
#> 18   p2_1 -> s3         2
#> 19         p2_2         2
#> 20 p2_2 -> p2_1         2
#> 21   p2_2 -> s2         2
#> 22   p2_2 -> s3         2
#> 23   s2 -> p2_1         2
#> 24   s2 -> p2_2         2
#> 25     s2 -> s3         2
#> 26   s3 -> p2_1         2
#> 27   s3 -> p2_2         2
#> 28     s3 -> s2         2
#> 29         p3_1         3
#> 30 p3_1 -> p3_2         3
#> 31   p3_1 -> s3         3
#> 32   p3_1 -> s4         3
#> 33         p3_2         3
#> 34 p3_2 -> p3_1         3
#> 35   p3_2 -> s3         3
#> 36   p3_2 -> s4         3
#> 37   s3 -> p3_1         3
#> 38   s3 -> p3_2         3
#> 39     s3 -> s4         3
#> 40   s4 -> p3_1         3
#> 41   s4 -> p3_2         3
#> 42     s4 -> s3         3
#> 43         p4_1         4
#> 44 p4_1 -> p4_2         4
#> 45   p4_1 -> s1         4
#> 46   p4_1 -> s4         4
#> 47         p4_2         4
#> 48 p4_2 -> p4_1         4
#> 49   p4_2 -> s1         4
#> 50   p4_2 -> s4         4
#> 51   s1 -> p4_1         4
#> 52   s1 -> p4_2         4
#> 53     s1 -> s4         4
#> 54   s4 -> p4_1         4
#> 55   s4 -> p4_2         4
#> 56     s4 -> s1         4
```
