# print and summary are stable

    Code
      print(comm)
    Output
      Memory-network communities (map equation)
        State nodes: 7 | physical nodes: 5 | modules: 3
        (1 module(s) hold only states with zero flow, reached by teleportation alone)
        Overlapping physical nodes: 1
        Codelength memory:      1.6061 bits (one module 2.2589)
        Codelength first-order: 2.2525 bits (one module 2.2525; modules: 1)
        Trials: 2 (teleportation 0.15, seed 1)
        Tables: as.data.frame(x, what = "states" | "physical" | "modules" | "trials" | "first_order" | "codelength")

---

    Code
      out <- summary(comm)
    Output
      Memory-network communities: summary
        Memory:      1.6061 bits, modules: 3, saves 0.6527 bits (28.9%) vs one module
        First-order: 2.2525 bits, modules: 1, saves 0.0000 bits (0.0%) vs one module
        Memory minus first-order: -0.6463 bits
        Overlap: 1 of 5 physical nodes in more than one module (mean 1.20 modules per node)
        Trials: 2; 2 reached the best codelength; ARI to best: mean 1.000, min 1.000

