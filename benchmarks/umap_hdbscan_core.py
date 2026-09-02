"""Reproducible UMAP/HDBSCAN assignment-stage benchmark.

This runner uses deterministic TF-IDF/SVD document embeddings, UMAP, and
scikit-learn HDBSCAN. It is a useful embedding-clustering baseline, but it is
not an execution of the BERTopic package.
"""

import argparse
import pathlib
import time

import numpy as np
import pandas as pd
from sklearn.cluster import HDBSCAN, KMeans
from sklearn.decomposition import TruncatedSVD
from sklearn.feature_extraction.text import TfidfVectorizer
from sklearn.metrics import (
    adjusted_mutual_info_score,
    adjusted_rand_score,
    normalized_mutual_info_score,
)
from sklearn.preprocessing import normalize
from umap import UMAP


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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", default="R8")
    parser.add_argument("--data", default="benchmarks/data")
    parser.add_argument("--output", required=True)
    parser.add_argument("--seeds", default="1,2,3,4,5,6,7,8,9,10")
    args = parser.parse_args()

    seeds = [int(x) for x in args.seeds.split(",")]
    docs, truth = load_textgcn(args.dataset, pathlib.Path(args.data))

    started = time.perf_counter()
    tfidf = TfidfVectorizer(min_df=2, max_features=20000, sublinear_tf=True)
    x = tfidf.fit_transform(docs)
    dimensions = min(100, x.shape[0] - 1, x.shape[1] - 1)
    embedding = TruncatedSVD(n_components=dimensions, random_state=0).fit_transform(x)
    embedding = normalize(embedding)
    embedding_seconds = time.perf_counter() - started

    rows = []
    for seed in seeds:
        fit_started = time.perf_counter()
        reduced = UMAP(
            n_neighbors=15,
            n_components=5,
            min_dist=0.0,
            metric="cosine",
            random_state=seed,
            n_jobs=1,
        ).fit_transform(embedding)
        # sklearn's min_samples includes the point itself, whereas the
        # contrib-hdbscan backend used by BERTopic does not. Eleven matches
        # BERTopic's effective default of ten neighbors.
        topic = HDBSCAN(
            min_cluster_size=10,
            min_samples=11,
            metric="euclidean",
            cluster_selection_method="eom",
            copy=True,
        ).fit_predict(reduced)
        assignments = {
            "umap_hdbscan": topic,
            # The matched-k control separates embedding/geometry quality
            # from HDBSCAN's natural cluster count.
            "umap_kmeans_k8": KMeans(
                n_clusters=len(np.unique(truth)), n_init=25,
                random_state=seed
            ).fit_predict(reduced),
        }
        for method, assigned in assignments.items():
            rows.append({
                "dataset": args.dataset,
                "method": method,
                "seed": seed,
                "ari": adjusted_rand_score(truth, assigned),
                "ami": adjusted_mutual_info_score(truth, assigned),
                "nmi": normalized_mutual_info_score(truth, assigned),
                "n_clusters": len(set(assigned)) - int(-1 in assigned),
                "outlier_fraction": float(np.mean(assigned == -1)),
                "embedding_s": embedding_seconds,
                "fit_s": time.perf_counter() - fit_started,
            })
        # Checkpoint each completed seed: long UMAP runs remain resumable even
        # if a later seed or downstream summarizer is interrupted.
        pd.DataFrame(rows).to_csv(args.output, index=False)


if __name__ == "__main__":
    main()
