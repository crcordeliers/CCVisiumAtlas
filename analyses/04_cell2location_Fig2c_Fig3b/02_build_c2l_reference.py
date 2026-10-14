"""
Cell2location reference, step 2: AnnData of the down-sampled 98-category Pelka reference.

Input : <C2L_DIR>/reference/{matrix.mtx.gz, barcodes.tsv.gz, features.tsv.gz, reference_cellIDs_98cat.tsv} (step 1)
Output: <C2L_DIR>/reference/pelka_reference_98cat.h5ad
Usage : python analyses/04_cell2location_Fig2c_Fig3b/02_build_c2l_reference.py
"""
import gzip
import os

import anndata as ad
import pandas as pd
import scipy.io as sio

REF_DIR = os.path.join(os.environ.get("CRCVISIUM_C2L", os.path.expanduser("~/data/cell2location")), "reference")

cell_ids = pd.read_csv(f"{REF_DIR}/reference_cellIDs_98cat.tsv", sep="\t")
with gzip.open(f"{REF_DIR}/barcodes.tsv.gz", "rt") as f:
    all_barcodes = [line.strip() for line in f]
with gzip.open(f"{REF_DIR}/features.tsv.gz", "rt") as f:
    all_features = [line.strip() for line in f]

barcode_to_idx = {b: i for i, b in enumerate(all_barcodes)}
missing = [b for b in cell_ids["cellID"] if b not in barcode_to_idx]
assert not missing, f"{len(missing)} reference cell IDs not found"
select_idx = [barcode_to_idx[b] for b in cell_ids["cellID"]]   # order of reference_cellIDs_98cat.tsv

with gzip.open(f"{REF_DIR}/matrix.mtx.gz", "rt") as f:
    mat = sio.mmread(f).tocsc()                                   # genes x cells
mat_sub = mat[:, select_idx]

obs = pd.DataFrame({"cell_type": cell_ids["cell_type"].values}, index=cell_ids["cellID"].values)
adata = ad.AnnData(X=mat_sub.T.tocsr(), obs=obs, var=pd.DataFrame(index=all_features))   # cells x genes
adata.obs["cell_type"] = adata.obs["cell_type"].astype("category")
assert adata.obs["cell_type"].nunique() == 98
adata.write_h5ad(f"{REF_DIR}/pelka_reference_98cat.h5ad")
print(adata)
