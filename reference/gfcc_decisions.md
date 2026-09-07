# Decisions of the German Federal Constitutional Court

The 3,618 decisions in volumes 1 to 160 of the official collection of
the German Federal Constitutional Court (BVerfGE), 1951 to 2022, as
released by Coupette, Hartung and Katz (2024). The node table of the
citation-block hypergraph in
[gfcc_citations](https://mohsaqr.github.io/hypernets/reference/gfcc_citations.md):
every decision is a node from its own date, whether or not it is ever
cited.

## Usage

``` r
gfcc_decisions
```

## Format

A data frame with 3,618 rows and 5 columns:

- decision:

  Decision key, `volume-page`.

- date:

  Date of the decision.

- citation:

  Official citation, e.g. `"BVerfGE 1, 1-3"`.

- type:

  Decision type (`"Urteil"` or `"Beschluss"`).

- title:

  Title of the decision, in German.

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
head(gfcc_decisions)
#>   decision       date        citation      type
#> 1  001-001 1951-09-09  BVerfGE 1, 1-3 Beschluss
#> 2  001-003 1951-09-27    BVerfGE 1, 3 Beschluss
#> 3  001-004 1951-09-27  BVerfGE 1, 4-5 Beschluss
#> 4  001-005 1951-09-27  BVerfGE 1, 5-7 Beschluss
#> 5  001-007 1951-10-03  BVerfGE 1, 7-8 Beschluss
#> 6  001-009 1951-10-11 BVerfGE 1, 9-10 Beschluss
#>                                                                                                                                                                                       title
#> 1 Neugliederung im Gebiet der Länder Baden, Württemberg-Baden und Württemberg-Hohenzollern gemäß Art. 118 GG ("Südweststaat" ): Aussetzung der Volksabstimmung durch einstweilige Anordnung
#> 2                                                                                                           Keine Geltendmachung von Schadensersatzansprüchen mittels Verfassungsbeschwerde
#> 3                                                                                                Voraussetzungen der Verfassungsbeschwerde gegen rechtskräftige gerichtliche Entscheidungen
#> 4                                                                                         Keine Verfassungsbeschwerde gegen auf Grund der Entnazifizierungsgesetze ergangene Entscheidungen
#> 5                                                                  Verfassungsbeschwerde gegen rechtskräftige Entscheidungen ordentlicher Gerichte: Umfang der Überprüfung durch das BVerfG
#> 6                                                                   Keine Verfassungsbeschwerde gegen nichtbeschwerdefähige gerichtliche Entscheidungen, die der Urteilsfällung vorausgehen
```
