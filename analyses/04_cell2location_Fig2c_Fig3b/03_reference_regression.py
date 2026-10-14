"""
Cell2location reference, step 3: negative-binomial regression of the reference -> cell-state expression signatures.

Input : <C2L_DIR>/reference/pelka_reference_98cat.h5ad (step 2)
Output: <C2L_DIR>/reference/inf_aver_98cat.csv  (genes x 98 categories)
Usage : python analyses/04_cell2location_Fig2c_Fig3b/03_reference_regression.py
Env   : cell2location 0.1.5, scvi-tools 1.3.3 (CPU is sufficient, ~1 h).
"""
import os

import scanpy as sc
from cell2location.models import RegressionModel
from cell2location.utils.filtering import filter_genes

REF_DIR = os.path.join(os.environ.get("CRCVISIUM_C2L", os.path.expanduser("~/data/cell2location")), "reference")

adata_ref = sc.read_h5ad(f"{REF_DIR}/pelka_reference_98cat.h5ad")
# standard cell2location gene filter (default thresholds)
selected = filter_genes(adata_ref, cell_count_cutoff=5, cell_percentage_cutoff2=0.03, nonz_mean_cutoff=1.12)
adata_ref = adata_ref[:, selected].copy()

import scvi
scvi.settings.seed = 42
RegressionModel.setup_anndata(adata_ref, labels_key="cell_type")
mod_ref = RegressionModel(adata_ref)
mod_ref.train(max_epochs=250)
adata_ref = mod_ref.export_posterior(adata_ref, sample_kwargs={"num_samples": 1000, "batch_size": 2500})

inf_aver = adata_ref.varm["means_per_cluster_mu_fg"]
inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
inf_aver.to_csv(f"{REF_DIR}/inf_aver_98cat.csv")
print("signatures:", inf_aver.shape, "| NaN:", int(inf_aver.isna().sum().sum()))
