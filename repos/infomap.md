# Infomap — map-equation community detection (Python / C++)

- **What**: Rosvall group's implementation of the map equation (Rosvall &
  Bergstrom 2008), with state-node (memory) networks (Rosvall et al. 2014;
  Edler, Bohlin & Rosvall 2017), multilayer networks and hierarchical
  modules.
- **Version**: 2.15.1 is current on PyPI (released 2026-08-04, checked
  2026-09-29) and is the version the oracle test was run with.
- **Role for us**: independent oracle for `hon_communities()`. The HON state
  network (state -> physical node) and a partition are handed to Infomap
  with `directed = True, two_level = True, no_infomap = True`, so Infomap
  only scores the given partition; flow and codelength agree to ~1e-14 bits.
  A second test lets Infomap search and requires hypernets' search to reach
  a codelength within 0.5% of Infomap's (memory and first-order maps); the
  partitions themselves are not required to match.
  `local_testing_and_equivalence/test-equiv-hon-communities-infomap.R`,
  gated by `INFOMAP_PYTHON` (a python that can `import infomap`).
- **Not implemented**: Infomap's coarse-tuning step (Edler et al. 2017,
  Alg. 6). Never a runtime dependency.
- **Links**: https://www.mapequation.org/infomap · https://github.com/mapequation/infomap ·
  `pip install infomap`
