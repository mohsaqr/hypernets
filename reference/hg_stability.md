# Stability of a hypergraph clustering across resolutions

Asks whether the partition
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
finds is a property of the data or of the particular sample, one
resolution `k` at a time.

## Usage

``` r
hg_stability(
  hg,
  k,
  type = c("zhou", "random_walk"),
  resample = c("subset", "seeds"),
  n_boot = 100L,
  fraction = 0.5,
  seed = 1L,
  what = c("k", "clusters"),
  seeds = c(1L, 99L),
  nstart = 25L
)
```

## Arguments

- hg:

  A
  [`text_hypergraph()`](https://mohsaqr.github.io/hypernets/reference/text_hypergraph.md)
  (or any hypernets `net_hg`).

- k:

  Vector of cluster counts to test, each at least 2.

- type:

  `"zhou"` or `"random_walk"`, the Laplacian of
  [`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md).

- resample:

  `"subset"` (default) for subsampling stability, or `"seeds"` for the
  two-seed solver check.

- n_boot:

  Number of subsamples (default `100L`).

- fraction:

  Share of the nodes in each subsample, strictly between 0 and 1
  (default `0.5`, as in `fpc::clusterboot()`).

- seed:

  Seed for the subsample draws and the k-means starts of every fit
  (default `1L`).

- what:

  For `resample = "subset"`: `"k"` (default) for one row per resolution,
  or `"clusters"` for one row per resolution and cluster.

- seeds:

  For `resample = "seeds"`: two distinct k-means seeds (default
  `c(1L, 99L)`).

- nstart:

  Number of k-means starts per fit (default `25L`).

## Value

A base `data.frame`.

`resample = "subset"`, `what = "k"`: one row per resolution with `k`,
`n_runs` (subsamples scored), `n_failed` (subsamples that disconnected
the hypergraph), `mean_jaccard` and `min_jaccard` (the mean and the
smallest of the clusters' stabilities), `n_stable` (clusters with
stability at least 0.75), `n_dissolved` (clusters with stability 0.5 or
below) and `eigengap`.

`resample = "subset"`, `what = "clusters"`: one row per resolution and
cluster of the full partition, with `k`, `cluster`, `size` (its nodes),
`jaccard` (its stability, `fpc::clusterboot()`'s `subsetmean`),
`n_dissolved` (subsamples with Jaccard 0.5 or below, `subsetbrd`),
`n_recovered` (subsamples with Jaccard above 0.75, `subsetrecover`) and
`n_runs`.

`resample = "seeds"`: one row per resolution with `k`,
`identical_partition` (are the two partitions literally identical,
labels included) and `ari` (their adjusted Rand index).

Raises `hypernets_bad_input` for a `fraction` outside (0, 1), a
subsample too small to hold `k` clusters, a non-positive `n_boot`, or
`what = "clusters"` with `resample = "seeds"`.

## Details

**`resample = "subset"` (default)** is cluster-wise subsampling
stability (Hennig 2007; the subsetting scheme of Ben-Hur et al. 2002).
The nodes (on a document hypergraph, the documents) are clustered once
in full. Then, `n_boot` times, a subsample of `floor(fraction * n)`
nodes is drawn without replacement, the sub-hypergraph is formed by
restricting every hyperedge to the drawn nodes (hyperedges left empty
are dropped; the incidence weights are kept, so a tf-idf corpus keeps
the full corpus's idf), and it is clustered again with the same `k`,
`type`, `seed` and `nstart`. Each original cluster, restricted to the
drawn nodes, is matched to its most similar subsample cluster by the
Jaccard coefficient \\\|A \cap B\| / \|A \cup B\|\\ (0 when the cluster
has no drawn member). A cluster's stability is its mean Jaccard over the
subsamples. Hennig (2007) reads a mean of 0.5 or below as a dissolved
cluster and 0.75 or above as a stable one. The subsamples are drawn
once, up front, after `set.seed(seed)`, so they are the ones
`fpc::clusterboot(bootmethod = "subset", subtuning = floor(fraction * n), seed = seed)`
draws, and the same subsamples serve every `k`. The caller's random
number stream is restored on exit.

A subsample can disconnect the hypergraph (a word that tied two groups
of documents may not be drawn);
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
cannot cut a disconnected hypergraph, so that run is not scored, is
counted in `n_failed`, and a warning of class
`hypernets_hypergraph_disconnected` reports the count.

The `eigengap` column is the gap \\\lambda\_{k+1} - \lambda_k\\ between
consecutive eigenvalues of the full hypergraph's Laplacian (the `gap`
column of `hg_cluster(what = "eigenvalues")`); the eigengap heuristic
(von Luxburg 2007, Sec. 8.1) prefers the `k` where it is large relative
to the gaps before it. Read it beside the stability: a `k` with a large
gap and stable clusters is supported twice.

**`resample = "seeds"`** keeps the check this verb made before 0.6.0: it
fits the full hypergraph twice per `k` with the two k-means `seeds` and
reports whether the partitions agree. It measures whether the solver is
deterministic on this embedding, not whether the structure survives a
change of sample, and it never resamples the data.

## References

Hennig, C. (2007). Cluster-wise assessment of cluster stability.
*Computational Statistics & Data Analysis*, 52(1), 258–271.
[doi:10.1016/j.csda.2006.11.025](https://doi.org/10.1016/j.csda.2006.11.025)

Ben-Hur, A., Elisseeff, A., & Guyon, I. (2002). A stability based method
for discovering structure in clustered data. *Pacific Symposium on
Biocomputing*, 7, 6–17.
[doi:10.1142/9789812799623_0002](https://doi.org/10.1142/9789812799623_0002)

von Luxburg, U. (2007). A tutorial on spectral clustering. *Statistics
and Computing*, 17(4), 395–416.
[doi:10.1007/s11222-007-9033-z](https://doi.org/10.1007/s11222-007-9033-z)

Hubert, L., & Arabie, P. (1985). Comparing partitions. *Journal of
Classification*, 2, 193–218.

## Examples

``` r
hg <- text_hypergraph(head(covid_abstracts, 40), column = "abstract",
                      id = "doc", stop_words = stop_words_en(),
                      min_count = 3)
#> vocabulary pruned to 360 of 1652 words, retaining 58.59% of tokens
hg_stability(hg, k = 2:4, n_boot = 20)
#>   k n_runs n_failed mean_jaccard min_jaccard n_stable n_dissolved   eigengap
#> 1 2     20        0    0.6010900   0.5806215        0           0 0.01302486
#> 2 3     20        0    0.5524471   0.4963394        0           1 0.01858627
#> 3 4     20        0    0.5950212   0.4436508        0           1 0.01636168
hg_stability(hg, k = 3, n_boot = 20, what = "clusters")
#>   k   cluster size   jaccard n_dissolved n_recovered n_runs
#> 1 3 Cluster 1   16 0.5814574          10           5     20
#> 2 3 Cluster 2   13 0.5795446           9           5     20
#> 3 3 Cluster 3   11 0.4963394          14           3     20
hg_stability(hg, k = 2, resample = "seeds")
#>   k identical_partition ari
#> 1 2                TRUE   1
```
