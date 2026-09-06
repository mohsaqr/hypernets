"""Multi-seed benchmark through the actual BERTopic package.

This runner instantiates ``bertopic.BERTopic`` (Grootendorst 2022) with the
package's own defaults: ``all-MiniLM-L6-v2`` sentence embeddings, UMAP
(15 neighbours, 5 components, min_dist 0, cosine), contrib HDBSCAN
(min_cluster_size 10, EOM) and class-based TF-IDF representations. Only the
UMAP ``random_state`` varies across seeds; the embeddings are computed once
and cached so that the R side can build a kNN hypergraph on the identical
vectors.

Three variants are scored per seed against the labelled classes:

* ``bertopic_default`` -- the package as shipped; HDBSCAN picks the topic
  count and leaves outliers labelled -1 (scored as their own cluster, matching
  ``umap_hdbscan_core.py``).
* ``bertopic_k`` -- the fitted model reduced with ``BERTopic.reduce_topics``
  to the number of labelled classes (plus the outlier topic, which BERTopic
  counts toward ``nr_topics``); outliers stay -1.
* ``bertopic_k_assigned`` -- the ``bertopic_k`` model followed by
  ``BERTopic.reduce_outliers(strategy="embeddings")`` so every document is
  assigned, the same completeness the honets spectral path has.
"""

import argparse
import pathlib
import time

import numpy as np
import pandas as pd
from bertopic import BERTopic
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
)
from umap import UMAP

EMBEDDING_MODEL = "all-MiniLM-L6-v2"


def load_textgcn(name: str, root: pathlib.Path):
    meta = pd.read_csv(
        root / f"{name}.txt",
        sep="\t",
        header=None,
        names=["row", "split", "label"],
        dtype=str,
    )
    corpus = (root / "corpus" / f"{name}.clean.txt").read_text(
        encoding="utf-8"
    ).splitlines()
    if len(corpus) != len(meta):
        raise ValueError("corpus and metadata lengths differ")
    return corpus, meta["label"].to_numpy()


def embed(docs, cache: pathlib.Path):
    """Encode once with BERTopic's default sentence model; cache as .npy and
    as a CSV (one row per document, TextGCN order) for the R runner."""
    npy = cache.with_suffix(".npy")
    if npy.exists():
        return np.load(npy), 0.0
    from sentence_transformers import SentenceTransformer

    started = time.perf_counter()
    model = SentenceTransformer(EMBEDDING_MODEL)
    embedding = model.encode(docs, batch_size=64, show_progress_bar=False,
                             convert_to_numpy=True)
    seconds = time.perf_counter() - started
    np.save(npy, embedding)
    pd.DataFrame(embedding).to_csv(cache.with_suffix(".csv"), index=False,
                                   header=False, float_format="%.7g")
    return embedding, seconds


def score(truth, assigned, extra):
    assigned = np.asarray(assigned)
    row = {
        "ari": adjusted_rand_score(truth, assigned),
        "ami": adjusted_mutual_info_score(truth, assigned),
        "nmi": normalized_mutual_info_score(truth, assigned),
        "n_clusters": len(set(assigned.tolist()) - {-1}),
        "outlier_fraction": float(np.mean(assigned == -1)),
    }
    row.update(extra)
    return row


def bertopic_umap(seed):
    # BERTopic's own default UMAP, with the seed threaded through.
    return UMAP(n_neighbors=15, n_components=5, min_dist=0.0,
                metric="cosine", low_memory=False, random_state=seed)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="R8")
    parser.add_argument("--data", default="benchmarks/data")
    parser.add_argument("--output", required=True)
    parser.add_argument("--embedding-cache", required=True,
                        help="path stem; .npy and .csv are written beside it")
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10")
    args = parser.parse_args()

    import bertopic

    seeds = [int(x) for x in args.seeds.split(",")]
    docs, truth = load_textgcn(args.dataset, pathlib.Path(args.data))
    n_classes = len(np.unique(truth))
    embedding, embedding_seconds = embed(docs, pathlib.Path(args.embedding_cache))
    if embedding.shape[0] != len(docs):
        raise ValueError("cached embedding does not match the corpus")

    common = {
        "dataset": args.dataset,
        "bertopic_version": bertopic.__version__,
        "embedding_model": EMBEDDING_MODEL,
        "embedding_s": embedding_seconds,
    }
    rows = []
    for seed in seeds:
        started = time.perf_counter()
        default = BERTopic(umap_model=bertopic_umap(seed), verbose=False)
        topics_default, _ = default.fit_transform(docs, embeddings=embedding)
        rows.append(score(truth, topics_default, {
            **common, "method": "bertopic_default", "seed": seed,
            "fit_s": time.perf_counter() - started}))

        # BERTopic's nr_topics counts the outlier topic (-1) toward the
        # total, so ask for one extra whenever HDBSCAN left outliers; the
        # k arm then has exactly n_classes real topics, like the honets arms.
        started = time.perf_counter()
        outliers = int(-1 in topics_default)
        default.reduce_topics(docs, nr_topics=n_classes + outliers)
        topics_k = list(default.topics_)
        fit_k = time.perf_counter() - started
        rows.append(score(truth, topics_k, {
            **common, "method": "bertopic_k", "seed": seed, "fit_s": fit_k}))

        started = time.perf_counter()
        if -1 in topics_k:
            topics_assigned = default.reduce_outliers(
                docs, topics_k, strategy="embeddings", embeddings=embedding)
        else:
            topics_assigned = topics_k
        rows.append(score(truth, topics_assigned, {
            **common, "method": "bertopic_k_assigned", "seed": seed,
            "fit_s": fit_k + time.perf_counter() - started}))
        # Checkpoint each completed seed so an interrupted run resumes cheaply.
        pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
