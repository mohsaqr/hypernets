# print and summary are stable

    Code
      print(comm)
    Output
      Memory-network communities (map equation)
        Memory nodes: 7 | states: 5 | communities: 3
        (1 community holds only nodes with zero flow, reached by teleportation alone)
        States in more than one community: 1
        Codelength memory:      1.6061 bits (one module 2.2589)
        Codelength first-order: 2.2525 bits (one module 2.2525; communities: 1)
        Runs: 2 (teleportation 0.15, seed 1)
         node state community      flow
            a     a         2 0.1664584
       a -> h     h         2 0.1629182
            b     b         2 0.1706234
            c     c         1 0.1664584
       c -> h     h         1 0.1629182
            d     d         1 0.1706234
            h     h         3 0.0000000

