# Co-clustering of nodes and hyperedges (documents and words)

Clusters the nodes and the hyperedges of a hypergraph together, so that
each cluster is a pair: a set of documents and the set of words that
characterises them, from one spectral cut of the document-word bipartite
graph (Dhillon 2001). With \\A\\ the incidence matrix (nodes by
hyperedges, the stored weights: counts or tf-idf on a
[`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)),
\\D_1\\ and \\D_2\\ its row and column sums, the scaled matrix \\A_n =
D_1^{-1/2} A D_2^{-1/2}\\ is decomposed by SVD; its leading singular
pair is trivial, the next \\\ell = \lceil \log_2 k \rceil\\ left and
right singular vectors \\U, V\\ give the embedding \\Z = \[D_1^{-1/2}
U;\\ D_2^{-1/2} V\]\\, one row per node and per hyperedge, and k-means
on the rows of \\Z\\ gives the `k` co-clusters (Dhillon 2001, Algorithm
Multipartition).

## Usage

``` r
hg_cocluster(
  hg,
  k,
  seed = NULL,
  nstart = 25L,
  what = c("clusters", "embedding"),
  role = c("all", "node", "hyperedge")
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hg`) with no empty node or hyperedge.

- k:

  Number of co-clusters, at least 2.

- seed:

  Seed for the k-means starts; the caller's random number stream is
  restored on exit. `NULL` (default) uses the current stream.

- nstart:

  Number of k-means starts (default `25L`).

- what:

  `"clusters"` (default) for the assignment, or `"embedding"` for the
  assignment plus the embedding coordinates `dim1..dimL`.

- role:

  Which rows to return: `"all"` (default), `"node"` (on a document
  hypergraph, the documents) or `"hyperedge"` (the words).

## Value

A base `data.frame`, one row per node and then per hyperedge (filtered
by `role`), with columns `node` (the node or hyperedge name), `role`
(`"node"` or `"hyperedge"`) and `cluster` (`"Cluster 1"`, ..., numbered
by first appearance in that row order); with `what = "embedding"` also
`dim1..dimL`. Raises `hypernets_bad_input` for a `k` below 2 or above
the number of rows, or an empty node or hyperedge (zero degree, where
the scaling is undefined).

## Details

For a binary incidence with unit hyperedge weights, the node half of the
embedding spans the same space as the Zhou et al. (2006) hypergraph
Laplacian eigenvectors that
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
cuts (the squared singular values are one minus the Laplacian
eigenvalues), so co-clustering reads the same cut on both sides of the
bipartite graph; with other weights the two differ, because the Zhou
Laplacian reads only the binary pattern. A cluster may hold nodes and no
hyperedges, or the reverse.

## References

Dhillon, I. S. (2001). Co-clustering documents and words using bipartite
spectral graph partitioning. *Proceedings of the Seventh ACM SIGKDD
International Conference on Knowledge Discovery and Data Mining*,
269–274.
[doi:10.1145/502512.502550](https://doi.org/10.1145/502512.502550)

Zhou, D., Huang, J., & Schölkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *Advances in Neural
Information Processing Systems 19*, 1601–1608.

## Examples

``` r
hg <- text_hypergraph(c(
  cooking_1 = "simmer the soup with onions and carrots",
  cooking_2 = "this soup recipe needs salt on a cold night",
  space_1 = "the telescope revealed a distant galaxy and stars",
  space_2 = "astronomers aimed the telescope at the stars all night"
), stop_words = c("the", "with", "and", "a", "this", "at", "on", "all"))
hg_cocluster(hg, k = 2, seed = 1)
#>           node      role   cluster
#> 1    cooking_1      node Cluster 1
#> 2    cooking_2      node Cluster 1
#> 3      space_1      node Cluster 2
#> 4      space_2      node Cluster 2
#> 5        aimed hyperedge Cluster 2
#> 6  astronomers hyperedge Cluster 2
#> 7      carrots hyperedge Cluster 1
#> 8         cold hyperedge Cluster 1
#> 9      distant hyperedge Cluster 2
#> 10      galaxy hyperedge Cluster 2
#> 11       needs hyperedge Cluster 1
#> 12       night hyperedge Cluster 2
#> 13      onions hyperedge Cluster 1
#> 14      recipe hyperedge Cluster 1
#> 15    revealed hyperedge Cluster 2
#> 16        salt hyperedge Cluster 1
#> 17      simmer hyperedge Cluster 1
#> 18        soup hyperedge Cluster 1
#> 19       stars hyperedge Cluster 2
#> 20   telescope hyperedge Cluster 2
hg_cocluster(hg, k = 2, seed = 1, role = "hyperedge")
#>           node      role   cluster
#> 1        aimed hyperedge Cluster 2
#> 2  astronomers hyperedge Cluster 2
#> 3      carrots hyperedge Cluster 1
#> 4         cold hyperedge Cluster 1
#> 5      distant hyperedge Cluster 2
#> 6       galaxy hyperedge Cluster 2
#> 7        needs hyperedge Cluster 1
#> 8        night hyperedge Cluster 2
#> 9       onions hyperedge Cluster 1
#> 10      recipe hyperedge Cluster 1
#> 11    revealed hyperedge Cluster 2
#> 12        salt hyperedge Cluster 1
#> 13      simmer hyperedge Cluster 1
#> 14        soup hyperedge Cluster 1
#> 15       stars hyperedge Cluster 2
#> 16   telescope hyperedge Cluster 2
```
