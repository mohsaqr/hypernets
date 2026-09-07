# Citation blocks of the German Federal Constitutional Court

Every citation between the decisions in
[gfcc_decisions](https://mohsaqr.github.io/hypernets/reference/gfcc_decisions.md),
grouped into citation blocks: an uninterrupted run of cited decisions
inside one decision, the hyperedge of Coupette, Hartung and Katz (2024).
One row per citation, dated by both decisions, so backward citations are
the rows with `date_citing > date_cited`.

## Usage

``` r
gfcc_citations
```

## Format

A data frame with 77,284 rows and 5 columns:

- citing:

  Key of the citing decision, the source of the block.

- cited:

  Key of the cited decision, the node.

- block:

  Citation block identifier, `citing:index`, the hyperedge.

- date_citing:

  Date of the citing decision, the time of the block.

- date_cited:

  Date of the cited decision.

## Source

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs (GFCC
dataset). Zenodo.
[doi:10.5281/zenodo.8081511](https://doi.org/10.5281/zenodo.8081511) ,
via the replication archive
[doi:10.5281/zenodo.8081507](https://doi.org/10.5281/zenodo.8081507) .
Licence CC BY-NC 4.0. Built by `data-raw/legal_hypergraphs.R`.

## References

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs.
*Philosophical Transactions of the Royal Society A*, 382(2270),
20230141.
[doi:10.1098/rsta.2023.0141](https://doi.org/10.1098/rsta.2023.0141)

## Examples

``` r
head(gfcc_citations)
#>    citing   cited       block date_citing date_cited
#> 1 001-001 001-014 001-001:000  1951-09-09 1951-10-23
#> 2 001-005 001-004 001-005:000  1951-09-27 1951-09-27
#> 3 001-066 001-014 001-066:000  1951-10-02 1951-10-23
#> 4 001-074 001-322 001-074:000  1951-11-13 1952-05-28
#> 5 001-085 001-117 001-085:000  1951-11-27 1952-02-20
#> 6 001-202 001-184 001-202:000  1952-03-20 1952-03-20
blocks <- temporal_hypergraph(gfcc_citations, actor = "cited",
                              group = "block", time = "date_citing",
                              nodes = gfcc_decisions, sparse = TRUE)
summary(blocks)
#>   n_nodes n_hyperedges n_event_times first_time last_time n_memberships
#> 1    3618        46257          2026          0     25762         77284
#>   mean_edge_size median_edge_size mean_duration  format time_unit
#> 1       1.670753                1            NA contact      days
```
