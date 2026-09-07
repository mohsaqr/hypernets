# cograph — graph and plotting engine

- **Relationship**: explicit hypernets Import. hypernets owns higher-order
  construction and projection; cograph begins when the result is an ordinary
  weighted graph.
- **Used for**: `splot()`, betweenness and closeness on s-line graphs,
  repeated Infomap, community quality and graph rendering.
- **Boundary**: cograph does not construct hypergraphs and hypernets does not
  duplicate its ordinary graph algorithms.
- **Collision audit 2026-09-02**: no exported name overlaps with hypernets.
- **Version checked locally**: 2.4.5.

