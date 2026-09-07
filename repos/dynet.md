# Dynet — temporal-network peer

- **Relationship**: separate sibling package, absent from both Imports and
  Suggests. It is not a computational layer under hypernets.
- **Boundary**: Dynet owns pairwise temporal networks. hypernets owns temporal
  hyperedges and snapshots because membership sets remain first-class.
- **Interoperability**: compatible `start`, `end`, `window` and snapshot
  vocabulary is intentional; future adapters may exchange projected graphs.
- **Collision audit 2026-09-02**: no exported name overlaps with hypernets.
- **Version checked locally**: 0.3.53.

