#!/usr/bin/env python3
"""
benchmark_clustering.py
========================
Runs the two Python-based alternative spatial clustering methods
(CellCharter and GraphST) used to benchmark a Louvain-based spatial
meta-clustering pipeline on multi-sample 10x Visium data.

The Louvain reference labels, BayesSpace, BASS, the metrics and the figures
are handled by benchmark_clustering.R. All five methods are compared on the
same 48 samples with the same target number of clusters (k = 11).

Methods
-------
- CellCharter (Varrone et al., Nat. Genet. 2024): a joint scVI latent
  embedding (batch-corrected across samples) is spatially smoothed by
  aggregating each spot's neighborhood, then clustered with a Gaussian
  mixture model.
- GraphST (Long et al., Nat. Commun. 2023): a graph self-supervised
  contrastive embedding is learned independently per sample (GraphST has no
  native multi-sample batch correction), after which the per-sample
  embeddings are concatenated and harmonized jointly (Harmony) before a
  final Gaussian-mixture clustering step. This joint post-hoc harmonization
  is required to make GraphST cluster labels comparable across patients.

Package versions:
    scvi-tools   1.4.2
    squidpy      1.8.2
    cellcharter  0.3.7
    GraphST      1.1.1
    torch        2.5.1+cu121
    scanpy       1.10.3
    harmonypy    2.0.0
    anndata      0.10.9

Usage
-----
    python benchmark_clustering.py --step cellcharter
    python benchmark_clustering.py --step graphst_embed --sample Ech_052_A1
    python benchmark_clustering.py --step graphst_join

Each Visium sample is expected as a standard SpaceRanger output directory
(`<DATA_ROOT>/<sample>/outs/`), i.e. what `scanpy.read_visium()` accepts.
"""

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import anndata as ad

# ----------------------------------------------------------------------
# Configuration -- EDIT THESE PATHS for your own environment
# ----------------------------------------------------------------------
DATA_ROOT = os.environ.get("CRCVISIUM_RAW", "/path/to/visium_spaceranger_outputs")   # one subdir per sample
OUT_DIR = os.path.join(os.environ.get("CRCVISIUM_WORKDIR", "."), "clustering_benchmark")
N_CLUSTERS = 11              # matched to the Louvain reference cluster count
SAMPLES_FILE = os.environ.get("CRCVISIUM_SAMPLES", "data/samples_48.txt")  # one sample name per line

os.makedirs(OUT_DIR, exist_ok=True)


def load_samples():
    with open(SAMPLES_FILE) as fh:
        return [l.strip() for l in fh if l.strip()]


# ========================================================================
# CellCharter
# ========================================================================
def run_cellcharter(n_clusters=N_CLUSTERS, n_latent=15, n_layers=3, seed=42, scvi_seed=0):
    """Joint scVI embedding (batch-corrected) + spatial neighborhood
    aggregation + Gaussian-mixture clustering, across all samples jointly."""
    import scvi
    import squidpy as sq
    import cellcharter as cc

    print(f"scvi-tools {scvi.__version__} | squidpy {sq.__version__} | "
          f"cellcharter {cc.__version__}")

    samples = load_samples()
    adatas = []
    for s in samples:
        a = sq.read.visium(os.path.join(DATA_ROOT, s, "outs"), library_id=s)
        a.var_names_make_unique()
        a.obs["sample"] = s
        adatas.append(a)
    # genes shared by all samples (inner join), raw counts for scVI
    adata = ad.concat(adatas, join="inner", label=None, index_unique="-", uns_merge="unique")
    adata.obs_names_make_unique()
    adata.obs["sample"] = adata.obs["sample"].astype("category")  # squidpy requires a categorical library_key
    adata.layers["counts"] = adata.X.copy()
    sc.pp.filter_genes(adata, min_cells=3)

    # --- joint scVI embedding, batch-corrected on `sample` ---
    scvi.settings.seed = scvi_seed  # scvi-tools default seed
    scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="sample")
    model = scvi.model.SCVI(adata, n_latent=n_latent)
    model.train(max_epochs=400, early_stopping=True)  # default max_epochs is ~40 for 193k spots
    adata.obsm["X_scVI"] = model.get_latent_representation()

    # --- spatial neighborhood aggregation on the scVI latent space ---
    sq.gr.spatial_neighbors(adata, library_key="sample", coord_type="grid", n_neighs=6)
    cc.gr.aggregate_neighbors(adata, n_layers=n_layers, use_rep="X_scVI",
                               out_key="X_cellcharter", sample_key="sample")

    # --- Gaussian-mixture clustering on the aggregated representation ---
    gmm = cc.tl.Cluster(n_clusters=n_clusters, random_state=seed)
    gmm.fit(adata, use_rep="X_cellcharter")
    adata.obs["cellcharter_cluster"] = gmm.predict(adata, use_rep="X_cellcharter")

    out = adata.obs[["sample", "cellcharter_cluster"]].copy()
    out["barcode"] = [i.rsplit("-", 1)[0] for i in adata.obs_names]  # drop the concat suffix
    out = out[["sample", "barcode", "cellcharter_cluster"]]
    out_path = os.path.join(OUT_DIR, "spot_cellcharter_cluster_all48.csv")
    out.to_csv(out_path, index=False)
    print("Wrote:", out_path)


# ========================================================================
# GraphST
# ========================================================================
def run_graphst_embed(sample):
    """Step 1/2: compute the GraphST embedding for a single sample
    (no cross-sample batch correction at this stage -- GraphST embeddings
    are not natively comparable across samples)."""
    import torch
    from GraphST.GraphST import GraphST  # the package exposes a GraphST module; the class lives inside it

    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    adata = sc.read_visium(os.path.join(DATA_ROOT, sample, "outs"))
    adata.var_names_make_unique()

    model = GraphST(adata, device=device, datatype="10X")  # default 600 epochs, random_seed=41
    adata = model.train()

    # raw GraphST embedding (before the per-sample PCA/mclust of GraphST.utils.clustering),
    # so that a single joint PCA + Harmony + clustering can be run across samples
    out_dir = os.path.join(OUT_DIR, "graphst_embeddings")
    os.makedirs(out_dir, exist_ok=True)
    np.save(os.path.join(out_dir, f"{sample}_emb.npy"), adata.obsm["emb"])
    np.save(os.path.join(out_dir, f"{sample}_barcodes.npy"), adata.obs_names.to_numpy())
    print(f"Wrote embedding for {sample} ({adata.obsm['emb'].shape[0]} spots).")


def run_graphst_join(n_clusters=N_CLUSTERS, seed=42):
    """Step 2/2: concatenate all per-sample GraphST embeddings, harmonize
    jointly across samples (Harmony, batch key = sample), and cluster the
    corrected embedding with a Gaussian mixture model. This is the step
    that makes GraphST labels comparable across patients."""
    import harmonypy
    from sklearn.decomposition import PCA
    from sklearn.mixture import GaussianMixture

    samples = load_samples()
    emb_dir = os.path.join(OUT_DIR, "graphst_embeddings")
    embs, metas = [], []
    for s in samples:
        emb = np.load(os.path.join(emb_dir, f"{s}_emb.npy"))
        bc = np.load(os.path.join(emb_dir, f"{s}_barcodes.npy"), allow_pickle=True)
        embs.append(emb)
        metas.append(pd.DataFrame({"sample": s, "barcode": bc}))
    meta = pd.concat(metas, ignore_index=True)
    emb_pca = PCA(n_components=20, random_state=seed).fit_transform(np.vstack(embs))

    ho = harmonypy.run_harmony(emb_pca, meta, ["sample"], random_state=seed)
    emb_harmony = ho.Z_corr
    if emb_harmony.shape[0] != emb_pca.shape[0]:
        emb_harmony = emb_harmony.T

    gmm = GaussianMixture(n_components=n_clusters, covariance_type="full",
                           random_state=seed, n_init=5)
    labels = gmm.fit_predict(emb_harmony)

    out = meta[["barcode", "sample"]].copy()
    out["graphst_cluster_harmony"] = labels
    out_path = os.path.join(OUT_DIR, "spot_graphst_cluster_harmony_all48.csv")
    out.to_csv(out_path, index=False)
    np.save(os.path.join(OUT_DIR, "emb_harmony_all48.npy"), emb_harmony)
    print("Wrote:", out_path)

    # sanity check: each cluster should span multiple patients, not be
    # driven by a single sample (see patient-mixing entropy in the metrics
    # script for the quantitative version of this check)
    check = out.groupby("graphst_cluster_harmony")["sample"].nunique()
    print("n distinct samples per cluster:\n", check)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--step", required=True,
                         choices=["cellcharter", "graphst_embed", "graphst_join"])
    parser.add_argument("--sample", help="required for --step graphst_embed")
    args = parser.parse_args()

    if args.step == "cellcharter":
        run_cellcharter()
    elif args.step == "graphst_embed":
        if not args.sample:
            parser.error("--sample is required for --step graphst_embed")
        run_graphst_embed(args.sample)
    elif args.step == "graphst_join":
        run_graphst_join()
