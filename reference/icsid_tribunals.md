# ICSID arbitration tribunals, one row per seat

The 742 cases registered at the International Centre for Settlement of
Investment Disputes between 1974 and June 2023 for which the three
tribunal members and the dates are known, as released by Coupette,
Hartung and Katz (2024): one row per seat, so a case has three rows.
Arbitrators sharing a case are the co-occurrence data of a temporal
hypergraph whose hyperedges are tribunals active from constitution to
conclusion; the worked analysis is `docs/legal-hypergraphs.Rmd` in the
source repository.

## Usage

``` r
icsid_tribunals
```

## Format

A data frame with 2,226 rows and 10 columns:

- case:

  ICSID case number, the hyperedge.

- seat:

  `"president"`, `"arbitrator_1"` or `"arbitrator_2"`.

- arbitrator:

  Name of the arbitrator holding the seat, the node.

- registered:

  Date the case was registered.

- constituted:

  Date the tribunal was constituted.

- concluded:

  Date the case concluded; `NA` for a pending case, whose tribunal stays
  active through the end of observation (2023-06-15).

- is_concluded:

  Whether the archive marks the case as concluded.

- economic_sector:

  Economic sector of the dispute.

- subject:

  Subject of the dispute.

- respondent:

  Respondent state.

## Source

Coupette, C., Hartung, D., & Katz, D. M. (2024). Legal hypergraphs
(ICSID dataset). Zenodo.
[doi:10.5281/zenodo.8081513](https://doi.org/10.5281/zenodo.8081513) ,
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
head(icsid_tribunals)
#>       case         seat             arbitrator registered constituted
#> 1 ARB/74/1 arbitrator_1 Jacques Michel GROSSEN 1974-03-06  1974-10-07
#> 2 ARB/74/1 arbitrator_2       Dominique PONCET 1974-03-06  1974-10-07
#> 3 ARB/74/1    president           Pierre CAVIN 1974-03-06  1974-10-07
#> 4 ARB/74/2 arbitrator_1           Michael KERR 1974-06-21  1974-12-16
#> 5 ARB/74/2 arbitrator_2           Fuad ROUHANI 1974-06-21  1974-12-16
#> 6 ARB/74/2    president          Jørgen TROLLE 1974-06-21  1974-12-16
#>    concluded is_concluded   economic_sector                           subject
#> 1 1977-08-29         TRUE    Other Industry Production of fibers and textiles
#> 2 1977-08-29         TRUE    Other Industry Production of fibers and textiles
#> 3 1977-08-29         TRUE    Other Industry Production of fibers and textiles
#> 4 1977-02-27         TRUE Oil, Gas & Mining                    Bauxite mining
#> 5 1977-02-27         TRUE Oil, Gas & Mining                    Bauxite mining
#> 6 1977-02-27         TRUE Oil, Gas & Mining                    Bauxite mining
#>      respondent
#> 1 Côte d'Ivoire
#> 2 Côte d'Ivoire
#> 3 Côte d'Ivoire
#> 4       Jamaica
#> 5       Jamaica
#> 6       Jamaica
tribunals <- temporal_hypergraph(icsid_tribunals, actor = "arbitrator",
                                 group = "case", start = "constituted",
                                 end = "concluded")
summary(tribunals)
#>   n_nodes n_hyperedges n_event_times first_time last_time n_memberships
#> 1     441          742          1076          0     17782          2226
#>   mean_edge_size median_edge_size mean_duration   format time_unit
#> 1              3                3      1303.436 interval      days
```
