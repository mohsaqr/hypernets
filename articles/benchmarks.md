# Benchmarks: R8, R52, MR, Ohsumed, 20NG

How far does closed-form hypergraph label spreading go on the five
standard text-classification benchmarks? This article reports test-set
accuracy of `text_hypergraph(sparse = TRUE)` +
[`hg_classify()`](https://mohsaqr.github.io/hypernets/reference/hg_classify.md)
on R8, R52, MR, Ohsumed, and 20-Newsgroups, against a tf-idf
nearest-centroid baseline and against the published accuracy tables in
Ding et al. (2020). Every number below was produced by the harness in
`benchmarks/` of the package repository; nothing is copied from memory.

## Setup

The corpora and train/test splits are the exact files of Yao, Mao & Luo
(2019), which Ding et al. (2020) also use – any other cleaning would
make the numbers incomparable. `benchmarks/download.sh` fetches them.
Each dataset becomes a sparse document–word hypergraph (documents are
vertices, words are hyperedges, tf-idf incidence weights); the training
documents are the labeled seeds; the spreading solution classifies every
test document in one conjugate-gradient solve per class.

``` r

hg <- text_hypergraph(docs, weight = "tfidf", sparse = TRUE)
fit <- hg_classify(hg, labels = train_labels, normalization = "class_mass")
```

Two checks tie the pipeline to the published setup. The split sizes
match Ding et al.’s Table 1 exactly on all five datasets. The vocabulary
sizes match exactly on R8 (7,688 hyperedges) and R52 (8,892); on MR,
Ohsumed and 20NG ours are smaller (18,151 vs 18,764; 13,349 vs 14,157;
36,970 vs 42,757) because the package tokenizer drops digit-containing
tokens. MR has one document unreachable from the rest of the corpus; the
harness runs on the giant component and would score unreachable test
documents as errors (here that document is not in the test set).

## Full-split results

`accuracy` is test-set accuracy with a 1,000-resample percentile
bootstrap CI; `fit_s` is the wall-clock seconds for the transductive
solve itself.

``` r

results <- read.csv("benchmark-results.csv")
knitr::kable(
  subset(results,
         select = c(dataset, method, accuracy, ci_low, ci_high, macro_f1,
                    build_s, fit_s)),
  digits = 4
)
```

| dataset | method | accuracy | ci_low | ci_high | macro_f1 | build_s | fit_s |
|:---|:---|---:|---:|---:|---:|---:|---:|
| 20ng | hypergraph_zhou_cmn | 0.8477 | 0.8399 | 0.8559 | 0.8406 | 41.416 | 1.223 |
| 20ng | hypergraph_random_walk_cmn | 0.7869 | 0.7780 | 0.7963 | 0.7893 | 41.283 | 2.734 |
| 20ng | hypergraph_zhou | 0.1020 | 0.0952 | 0.1091 | 0.0583 | 40.435 | 1.359 |
| 20ng | tfidf_centroid | 0.7796 | 0.7699 | 0.7890 | 0.7725 | 42.593 | 0.076 |
| R8 | hypergraph_zhou_cmn | 0.8451 | 0.8291 | 0.8593 | 0.6730 | 6.541 | 0.167 |
| R8 | hypergraph_random_walk_cmn | 0.8173 | 0.8017 | 0.8342 | 0.6470 | 6.444 | 0.265 |
| R8 | hypergraph_zhou | 0.4947 | 0.4737 | 0.5162 | 0.0827 | 5.921 | 0.104 |
| R8 | tfidf_centroid | 0.9246 | 0.9136 | 0.9361 | 0.8733 | 5.771 | 0.018 |
| R52 | hypergraph_zhou_cmn | 0.5284 | 0.5070 | 0.5471 | 0.3700 | 9.138 | 0.801 |
| R52 | hypergraph_random_walk_cmn | 0.5210 | 0.5004 | 0.5390 | 0.3771 | 5.227 | 1.126 |
| R52 | hypergraph_zhou | 0.4217 | 0.4019 | 0.4412 | 0.0114 | 6.680 | 0.851 |
| R52 | tfidf_centroid | 0.8797 | 0.8668 | 0.8917 | 0.7008 | 7.833 | 0.035 |
| ohsumed | hypergraph_zhou_cmn | 0.4702 | 0.4551 | 0.4860 | 0.4360 | 10.998 | 0.331 |
| ohsumed | hypergraph_random_walk_cmn | 0.5110 | 0.4964 | 0.5278 | 0.4616 | 9.715 | 0.630 |
| ohsumed | hypergraph_zhou | 0.1459 | 0.1348 | 0.1563 | 0.0111 | 9.218 | 0.349 |
| ohsumed | tfidf_centroid | 0.6347 | 0.6196 | 0.6500 | 0.5826 | 8.271 | 0.019 |
| mr | hypergraph_zhou_cmn | 0.7684 | 0.7538 | 0.7811 | 0.7683 | 6.402 | 0.055 |
| mr | hypergraph_random_walk_cmn | 0.7563 | 0.7423 | 0.7698 | 0.7560 | 5.987 | 0.135 |
| mr | hypergraph_zhou | 0.7414 | 0.7273 | 0.7561 | 0.7354 | 5.775 | 0.055 |
| mr | tfidf_centroid | 0.6851 | 0.6708 | 0.7006 | 0.6850 | 3.045 | 0.007 |

Three findings. First, on the two class-balanced datasets the hypergraph
classifier beats the centroid baseline outright: 0.8477 vs 0.7796 on
20NG and 0.7684 vs 0.6851 on MR. Second, on the skewed multi-class
corpora (R8, R52, Ohsumed) the centroid wins; label spreading pays for
class imbalance even after normalization, and R52’s 52 heavily skewed
classes are its worst case (0.5284 vs 0.8797). Third, the
`hypergraph_zhou` rows – the raw Zhou (2006) argmax with no
normalization – collapse onto the majority class (0.4947 on R8 is
exactly the majority-class rate; 0.1020 on 20NG), which is why
`normalization = "class_mass"` (Zhu, Ghahramani & Lafferty 2003) exists
and is used everywhere else in this article. The EDVW walk
(`hypergraph_random_walk_cmn`) is best on Ohsumed (0.5110 vs 0.4702) and
second elsewhere.

The cost side: the largest solve – 20 classes over 18,846 documents x
36,970 hyperedges – fits in 1.2 seconds on a laptop CPU, with no
training, no gradient steps and no GPU. Ding et al.’s Table 3 reports
1,479 MB of GPU memory for TextGCN and 180 MB for HyperGAT on the same
corpus.

## Against the published tables

Published rows are means over 10 runs from Ding et al. (2020), Table 2;
our transductive rows are deterministic given the split.

``` r

published <- read.csv("benchmark-published.csv")
knitr::kable(published, digits = 4)
```

| dataset | model                  | accuracy |     sd | source                    |
|:--------|:-----------------------|---------:|-------:|:--------------------------|
| 20ng    | TextGCN (transductive) |   0.8643 | 0.0009 | Ding et al. 2020, Table 2 |
| R8      | TextGCN (transductive) |   0.9707 | 0.0010 | Ding et al. 2020, Table 2 |
| R52     | TextGCN (transductive) |   0.9356 | 0.0018 | Ding et al. 2020, Table 2 |
| ohsumed | TextGCN (transductive) |   0.6836 | 0.0056 | Ding et al. 2020, Table 2 |
| mr      | TextGCN (transductive) |   0.7674 | 0.0020 | Ding et al. 2020, Table 2 |
| 20ng    | HyperGAT               |   0.8662 | 0.0016 | Ding et al. 2020, Table 2 |
| R8      | HyperGAT               |   0.9797 | 0.0023 | Ding et al. 2020, Table 2 |
| R52     | HyperGAT               |   0.9498 | 0.0027 | Ding et al. 2020, Table 2 |
| ohsumed | HyperGAT               |   0.6990 | 0.0034 | Ding et al. 2020, Table 2 |
| mr      | HyperGAT               |   0.7832 | 0.0027 | Ding et al. 2020, Table 2 |
| 20ng    | fastText               |   0.7938 | 0.0030 | Ding et al. 2020, Table 2 |
| R8      | fastText               |   0.9613 | 0.0021 | Ding et al. 2020, Table 2 |
| R52     | fastText               |   0.9281 | 0.0009 | Ding et al. 2020, Table 2 |
| ohsumed | fastText               |   0.5770 | 0.0049 | Ding et al. 2020, Table 2 |
| mr      | fastText               |   0.7514 | 0.0020 | Ding et al. 2020, Table 2 |

On 20NG the hypergraph classifier (0.8477) is within two points of
transductive TextGCN (0.8643) and HyperGAT (0.8662), and above fastText
(0.7938). On MR (0.7684) it edges out published transductive TextGCN
(0.7674) and trails HyperGAT (0.7832). On R8, R52 and Ohsumed the
trained models are clearly ahead. The honest summary: a closed-form
spectral method with zero trained parameters recovers most of the
accuracy of trained graph neural networks on balanced corpora and is not
competitive on skewed ones.

## The low-label regime

Label spreading is designed for few labels, so the full-split table
above is its worst case. Here the training seeds are subsampled
(stratified, at least one per class, five draws per fraction; the SD
column is over draws) and both methods receive the same seeds.

``` r

lowlabel <- read.csv("benchmark-lowlabel.csv")
knitr::kable(lowlabel, digits = 4)
```

| dataset | fraction | n_seeds | method              | accuracy |     sd | n_draws |
|:--------|---------:|--------:|:--------------------|---------:|-------:|--------:|
| R8      |     0.01 |      59 | hypergraph_zhou_cmn |   0.6965 | 0.0375 |       5 |
| R8      |     0.05 |     278 | hypergraph_zhou_cmn |   0.7523 | 0.0197 |       5 |
| R8      |     0.10 |     552 | hypergraph_zhou_cmn |   0.7869 | 0.0146 |       5 |
| R8      |     0.20 |    1101 | hypergraph_zhou_cmn |   0.8090 | 0.0055 |       5 |
| R8      |     0.01 |      59 | tfidf_centroid      |   0.8917 | 0.0080 |       5 |
| R8      |     0.05 |     278 | tfidf_centroid      |   0.9222 | 0.0093 |       5 |
| R8      |     0.10 |     552 | tfidf_centroid      |   0.9218 | 0.0073 |       5 |
| R8      |     0.20 |    1101 | tfidf_centroid      |   0.9246 | 0.0030 |       5 |
| mr      |     0.01 |      72 | hypergraph_zhou_cmn |   0.5675 | 0.0099 |       5 |
| mr      |     0.05 |     356 | hypergraph_zhou_cmn |   0.6280 | 0.0129 |       5 |
| mr      |     0.10 |     712 | hypergraph_zhou_cmn |   0.6671 | 0.0049 |       5 |
| mr      |     0.20 |    1422 | hypergraph_zhou_cmn |   0.6984 | 0.0069 |       5 |
| mr      |     0.01 |      72 | tfidf_centroid      |   0.5541 | 0.0072 |       5 |
| mr      |     0.05 |     356 | tfidf_centroid      |   0.6201 | 0.0046 |       5 |
| mr      |     0.10 |     712 | tfidf_centroid      |   0.6522 | 0.0115 |       5 |
| mr      |     0.20 |    1422 | tfidf_centroid      |   0.6677 | 0.0025 |       5 |
| ohsumed |     0.01 |      46 | hypergraph_zhou_cmn |   0.1345 | 0.0094 |       5 |
| ohsumed |     0.05 |     178 | hypergraph_zhou_cmn |   0.2120 | 0.0133 |       5 |
| ohsumed |     0.10 |     345 | hypergraph_zhou_cmn |   0.2567 | 0.0120 |       5 |
| ohsumed |     0.20 |     678 | hypergraph_zhou_cmn |   0.3248 | 0.0120 |       5 |
| ohsumed |     0.01 |      46 | tfidf_centroid      |   0.2422 | 0.0350 |       5 |
| ohsumed |     0.05 |     178 | tfidf_centroid      |   0.4094 | 0.0242 |       5 |
| ohsumed |     0.10 |     345 | tfidf_centroid      |   0.4876 | 0.0128 |       5 |
| ohsumed |     0.20 |     678 | tfidf_centroid      |   0.5463 | 0.0026 |       5 |

On MR the hypergraph classifier leads at every fraction, from 72 seeds
(0.5675 vs 0.5541) to 1,422 seeds (0.6984 vs 0.6677): with few labels,
the document–word structure carries information the centroid cannot use.
On R8 and Ohsumed the centroid stays ahead at every fraction – strong
tf-idf class signatures survive subsampling better than spreading does.
Few labels help the case for structure, but they do not overturn the
class balance effect.

## The neural tier: HGNN

[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md)
trains the two-layer hypergraph convolutional network of Feng et
al. (2019) natively in R ({torch}), on the same sparse document–word
hypergraph and tf-idf features. Its propagation matrix is exactly the
Zhou operator the transductive classifier spreads with – verified to
machine precision in the tests, with forward-pass parity against the
official implementation (DHG) at float32 precision. Training is
stochastic, so every number is the mean of three seeds (SD and range
shown); each run holds out 10% of the seeds (stratified) to pick the
best epoch. Two configurations are reported: the original paper’s
(`lr = 0.001`, 200 epochs, tuned for citation graphs with dense
features) and the package default (`lr = 0.01`, 600 epochs, selected by
validation accuracy for high-dimensional sparse text features).

``` r

neural <- read.csv("benchmark-neural.csv")
knitr::kable(
  subset(neural,
         select = c(dataset, epochs, lr, accuracy, sd, acc_min, acc_max,
                    macro_f1, fit_s)),
  digits = 4
)
```

| dataset | epochs |    lr | accuracy |     sd | acc_min | acc_max | macro_f1 |    fit_s |
|:--------|-------:|------:|---------:|-------:|--------:|--------:|---------:|---------:|
| 20ng    |    600 | 0.010 |   0.5355 | 0.0155 |  0.5260 |  0.5534 |   0.4998 | 840.9963 |
| 20ng    |    200 | 0.001 |   0.6347 | 0.0089 |  0.6257 |  0.6435 |   0.6031 | 171.4803 |
| R8      |    600 | 0.010 |   0.9539 | 0.0044 |  0.9511 |  0.9589 |   0.8491 |  79.9967 |
| R8      |    200 | 0.001 |   0.9202 | 0.0049 |  0.9146 |  0.9237 |   0.6366 |  25.7847 |
| R52     |    600 | 0.010 |   0.8440 | 0.0013 |  0.8431 |  0.8454 |   0.2168 | 439.9200 |
| R52     |    200 | 0.001 |   0.7714 | 0.0051 |  0.7667 |  0.7769 |   0.0893 |  33.6647 |
| ohsumed |    600 | 0.010 |   0.4093 | 0.0306 |  0.3740 |  0.4286 |   0.1611 | 121.4023 |
| ohsumed |    200 | 0.001 |   0.3814 | 0.0014 |  0.3799 |  0.3826 |   0.1357 |  40.1323 |
| mr      |    600 | 0.010 |   0.7692 | 0.0045 |  0.7642 |  0.7729 |   0.7684 |  52.3480 |
| mr      |    200 | 0.001 |   0.7671 | 0.0034 |  0.7648 |  0.7710 |   0.7668 |  17.4503 |

The tuned configuration lifts R8 to 0.9539 (above the centroid’s 0.9246;
published TextGCN 0.9707) and R52 to 0.8440 (a 32-point jump over
transduction’s 0.5284, though still under the centroid), and holds MR at
0.7692 – the best HGNN number in this article for that dataset. On
Ohsumed and 20NG, however, HGNN underperforms the closed-form
transduction at every configuration tested (20NG: 0.5355 tuned, 0.6347
at paper defaults, vs 0.8477 for transduction). The loss curves say this
is not an optimization failure alone: the corpus-level document-node
design oversmooths – two rounds of propagation through high-degree word
hyperedges blur exactly the distinctions 20 newsgroups need. The
literature’s answer is document-level hypergraphs with attention over
hyperedges: HyperGAT, next.

## HyperGAT: document-level hypergraphs with dual attention

[`hg_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md)
implements Ding et al. (2020) natively: every document becomes its own
hypergraph (its unique words as vertices, its sentences as hyperedges),
a first attention layer aggregates words into sentence representations,
a second aggregates sentences back into words, and a masked mean pool
feeds the classifier. The implementation mirrors the official code
(verified against its source, with forward-pass parity at float32
precision). `semantic = "none"` uses sentence hyperedges alone;
`semantic = "lda"` adds the paper’s online-LDA topic hyperedges, trained
on labeled documents only, and precomputed official keyword lists can be
passed to replay the upstream preprocessing boundary exactly. The
benchmark table below predates that implementation and remains
explicitly sentence-only, so its closest published comparator is the
paper’s “w/o semantic” ablation. One further honest difference: the
official preprocessing lemmatizes with WordNet; ours does not. The
benchmark runs on the official repository’s own corpus files (the
TextGCN cleaned files used above contain no sentence punctuation at
all), whose R8 split is 5,501/2,190 rather than TextGCN’s 5,485/2,189.
Official training protocol throughout: 10 epochs, batch 8, StepLR,
balanced class weights, 10% validation checkpointing; three torch seeds.

``` r

hypergat <- read.csv("benchmark-hypergat.csv")
knitr::kable(hypergat, digits = 4)
```

| dataset | method | epochs | n_seeds | n_train | n_test | accuracy | sd | acc_min | acc_max | macro_f1 | fit_s |
|:---|:---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| R8 | hypergat | 10 | 3 | 5501 | 2190 | 0.9665 | 0.0034 | 0.9630 | 0.9699 | 0.9046 | 180.9783 |
| R52 | hypergat | 10 | 3 | 6560 | 2570 | 0.9433 | 0.0018 | 0.9412 | 0.9447 | 0.7604 | 227.9630 |

This is the package’s strongest result. R8 reaches 0.9665 (SD 0.0034) –
above every other method in this article, within 1.3 points of the full
published HyperGAT (0.9797) and at published transductive TextGCN level
(0.9707). R52 reaches 0.9433 (SD 0.0018), which **exceeds published
transductive TextGCN (0.9356)** and sits 0.7 points from published
HyperGAT (0.9498) – against our previous R52 best of 0.8797 (the
centroid). Wall-clock: under four minutes per training run on a laptop
CPU. The architecture, not more optimization of the corpus-level design,
was what the skewed news corpora needed.

## Unsupervised UMAP/HDBSCAN baseline

Classification accuracy is not a valid comparator for unsupervised
clustering, so a separate R8 experiment compares
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md)
with fixed TF-IDF/SVD embeddings followed by UMAP and either HDBSCAN or
a matched eight-cluster k-means control. This is explicitly **not** an
execution of the BERTopic package; the actual BERTopic run follows in
the next section.

``` r

unsupervised_effects <- read.csv("benchmark-umap-hdbscan-effects.csv")
knitr::kable(
  subset(unsupervised_effects,
         select = c(comparator, metric, honets_mean, baseline_mean,
                    mean_difference, difference_ci_low,
                    difference_ci_high, paired_dz)),
  digits = 4
)
```

| comparator | metric | honets_mean | baseline_mean | mean_difference | difference_ci_low | difference_ci_high | paired_dz |
|:---|:---|---:|---:|---:|---:|---:|---:|
| umap_hdbscan | ari | 0.4306 | 0.0366 | 0.3940 | 0.3855 | 0.4013 | 29.9417 |
| umap_hdbscan | ami | 0.5416 | 0.3320 | 0.2096 | 0.2024 | 0.2157 | 18.2523 |
| umap_hdbscan | nmi | 0.5425 | 0.3499 | 0.1926 | 0.1860 | 0.1983 | 18.2318 |
| umap_kmeans_k8 | ari | 0.4306 | 0.3934 | 0.0372 | 0.0270 | 0.0478 | 2.0459 |
| umap_kmeans_k8 | ami | 0.5416 | 0.5635 | -0.0219 | -0.0350 | -0.0088 | -0.9807 |
| umap_kmeans_k8 | nmi | 0.5425 | 0.5644 | -0.0219 | -0.0341 | -0.0087 | -0.9825 |

Default HDBSCAN over-splits R8 (142–170 clusters and 17.5–23.0%
outliers), so honets leads all three external metrics. At matched
cluster count, honets leads ARI by 0.0372 (95% paired bootstrap CI
0.0270–0.0478), while the UMAP+k-means baseline leads AMI and NMI by
about 0.0219. Both sides of the result are reported.

## BERTopic baseline (the actual package)

The same R8 experiment was then run through `bertopic.BERTopic` 0.17.4
with its own defaults (`all-MiniLM-L6-v2` sentence embeddings, seeded
UMAP, contrib HDBSCAN, c-TF-IDF), ten UMAP seeds. Three BERTopic
variants are scored: the package as shipped (HDBSCAN picks the topic
count and leaves outliers as -1), the fitted model reduced with
`reduce_topics()` to eight real topics, and that model with
`reduce_outliers(strategy = "embeddings")` so every document is
assigned. Two honets arms cluster with the same k = 8 through
[`hg_cluster()`](https://mohsaqr.github.io/hypernets/reference/hg_cluster.md) +
k-means: the tf-idf document–word hypergraph and a kNN hypergraph
(`text_hypergraph(construction = "knn")`) built on the identical MiniLM
vectors, which isolates the clustering paradigm from the embedding.

``` r

bertopic_seeds <- read.csv("benchmark-bertopic-seeds.csv")
knitr::kable(
  aggregate(cbind(ari, ami, nmi, n_clusters, outlier_fraction) ~ method,
            data = bertopic_seeds, FUN = mean),
  digits = 4, caption = "Means over ten seeds"
)
```

| method              |    ari |    ami |    nmi | n_clusters | outlier_fraction |
|:--------------------|-------:|-------:|-------:|-----------:|-----------------:|
| bertopic_default    | 0.0478 | 0.3326 | 0.3490 |      130.7 |           0.3328 |
| bertopic_k          | 0.2320 | 0.3578 | 0.3597 |        8.0 |           0.3328 |
| bertopic_k_assigned | 0.4099 | 0.5328 | 0.5342 |        8.0 |           0.0000 |
| honets_hypergraph   | 0.4306 | 0.5416 | 0.5425 |        8.0 |           0.0000 |
| honets_knn_sbert    | 0.5124 | 0.6107 | 0.6115 |        8.0 |           0.0000 |

Means over ten seeds {.table}

``` r

bertopic_effects <- read.csv("benchmark-bertopic-effects.csv")
knitr::kable(
  subset(bertopic_effects, comparator == "bertopic_k_assigned",
         select = c(reference, metric, honets_mean, baseline_mean,
                    mean_difference, difference_ci_low,
                    difference_ci_high)),
  digits = 4, caption = "Paired effects against the fully assigned eight-topic BERTopic"
)
```

|  | reference | metric | honets_mean | baseline_mean | mean_difference | difference_ci_low | difference_ci_high |
|:---|:---|:---|---:|---:|---:|---:|---:|
| 5 | honets_hypergraph | ari | 0.4306 | 0.4099 | 0.0207 | 0.0176 | 0.0230 |
| 6 | honets_knn_sbert | ari | 0.5124 | 0.4099 | 0.1026 | 0.0992 | 0.1050 |
| 11 | honets_hypergraph | ami | 0.5416 | 0.5328 | 0.0087 | 0.0011 | 0.0153 |
| 12 | honets_knn_sbert | ami | 0.6107 | 0.5328 | 0.0778 | 0.0703 | 0.0845 |
| 17 | honets_hypergraph | nmi | 0.5425 | 0.5342 | 0.0082 | 0.0001 | 0.0148 |
| 18 | honets_knn_sbert | nmi | 0.6115 | 0.5342 | 0.0772 | 0.0694 | 0.0839 |

Paired effects against the fully assigned eight-topic BERTopic {.table}

BERTopic as shipped finds 121–138 topics and leaves a third of R8 as
outliers, so it trails on every metric. At matched topic count with
every document assigned, the tf-idf hypergraph leads ARI by 0.0207 (95%
paired bootstrap CI 0.0176–0.0230) and AMI/NMI by under 0.01 with
intervals that narrowly exclude zero. On the same sentence embeddings,
the kNN hypergraph leads ARI by 0.1026 (0.0992–0.1050) and AMI/NMI by
about 0.077. The honets arms have essentially no seed variance (the
spectral embedding is deterministic), so the intervals reflect
BERTopic’s UMAP seed variance. Only cluster assignments are scored;
BERTopic’s topic representations are not part of this comparison.

## When to use which

Use `hg_classify(normalization = "class_mass")` when classes are roughly
balanced or labels are scarce, and always over the raw argmax: the raw
rule’s majority-class collapse on skewed seeds is total. Use
`type = "random_walk"` when tf-idf weights should shape the walk itself;
it was best on Ohsumed here. Use the neural tier when plentiful labels
can pay for training, and pick the architecture by corpus shape:
[`hg_hypergat()`](https://mohsaqr.github.io/hypernets/reference/hg_hypergat.md)
for classification of documents with sentence structure – it holds the
package’s best R8 (0.9665) and R52 (0.9433, above published transductive
TextGCN) – and
[`hg_neural()`](https://mohsaqr.github.io/hypernets/reference/hg_neural.md)
(HGNN) when the corpus-level hypergraph with rich vertex features is the
natural object (it holds MR at 0.7692). Neither beats the closed-form
classifier everywhere: on 20NG the transduction’s 0.8477 still leads
everything we trained. For skewed many-class corpora without training, a
tf-idf centroid remains a strong, nearly free baseline; the full
published HyperGAT (with LDA hyperedges) holds the state of the art. The
package’s contribution is that every non-published row above runs
natively in R – the transductive ones in seconds, the neural ones in
minutes on a CPU – each from one verb.

## References

Feng, Y., You, H., Zhang, Z., Ji, R., & Gao, Y. (2019). Hypergraph
neural networks. *AAAI 33*.

Ding, K., Wang, J., Li, J., Li, D., & Liu, H. (2020). Be more with less:
Hypergraph attention networks for inductive text classification. *EMNLP
2020*.

Yao, L., Mao, C., & Luo, Y. (2019). Graph convolutional networks for
text classification. *AAAI 2019*.

Zhou, D., Huang, J., & Scholkopf, B. (2006). Learning with hypergraphs:
Clustering, classification, and embedding. *NeurIPS 19*.

Zhu, X., Ghahramani, Z., & Lafferty, J. (2003). Semi-supervised learning
using Gaussian fields and harmonic functions. *ICML 20*.
