# Dynet — temporal-network peer

- **Relationship**: separate sibling package, absent from both Imports and
  Suggests. It is not a computational layer under hypernets.
- **Boundary**: Dynet owns pairwise temporal networks. hypernets owns temporal
  hyperedges and snapshots because membership sets remain first-class.
- **Interoperability**: compatible `start`, `end`, `window` and snapshot
  vocabulary is intentional; future adapters may exchange projected graphs.
- **Collision audit 2026-09-02**: no exported name overlaps with hypernets.
- **Version (checked 2026-09-29)**: installed 0.4.10. The current source is
  the sibling checkout `../Dyna` (GitHub `mohsaqr/Dynet`), DESCRIPTION 0.5.1,
  last commit 2026-09-27; r-universe (`mohsaqr.r-universe.dev`) serves 0.5.1.
  Not on CRAN (package page 404). `../temporal`, the path CLAUDE.md names,
  no longer exists; `../temporal-next` is an older checkout at 0.4.7.
  hypernets' vocabulary alignment was written against the 0.4.x API and has
  not been re-checked against 0.5.1.

