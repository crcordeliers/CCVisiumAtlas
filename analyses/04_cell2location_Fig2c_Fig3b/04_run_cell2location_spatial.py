"""
Cell2location, step 4: spatial mapping of the 98 reference cell states, one Visium sample per run (30,000 epochs, GPU).

Input : Space Ranger output <CRCVISIUM_RAW>/<SAMPLE>/outs ; <C2L_DIR>/reference/inf_aver_98cat.csv (step 3)
Output: <C2L_DIR>/objects/cell2location_98cat_30000ep_<SAMPLE>.h5ad ; ELBO history in <C2L_DIR>/elbo_history/
Usage : SAMPLE=Ech_074_D1 python analyses/04_cell2location_Fig2c_Fig3b/04_run_cell2location_spatial.py
        (one run per sample of data/samples_48.txt; each sample is fitted independently)
Env   : cell2location 0.1.5, scvi-tools 1.3.3, CUDA 12.1 (L40S / A100).
"""
import os

import pandas as pd
import scanpy as sc
from cell2location.models import Cell2location
from lightning.pytorch.callbacks import Callback

SAMPLE = os.environ["SAMPLE"]
RAW_DIR = os.environ.get("CRCVISIUM_RAW", os.path.expanduser("~/data/processed"))
C2L_DIR = os.environ.get("CRCVISIUM_C2L", os.path.expanduser("~/data/cell2location"))
OUTDIR = os.path.join(C2L_DIR, "objects")
ELBODIR = os.path.join(C2L_DIR, "elbo_history")
EPOCHS = 30000
os.makedirs(OUTDIR, exist_ok=True)
os.makedirs(ELBODIR, exist_ok=True)


class ElboLogger(Callback):
    """Writes the training ELBO every `every` epochs."""
    def __init__(self, path, every=500):
        self.path, self.every, self.rows = path, every, []

    def on_train_epoch_end(self, trainer, pl_module):
        self.rows.append((trainer.current_epoch, float(trainer.callback_metrics.get("elbo_train", float("nan")))))
        if trainer.current_epoch % self.every == 0:
            pd.DataFrame(self.rows, columns=["epoch", "elbo_train"]).to_csv(self.path, index=False)

    def on_train_end(self, trainer, pl_module):
        pd.DataFrame(self.rows, columns=["epoch", "elbo_train"]).to_csv(self.path, index=False)


adata_vis = sc.read_visium(os.path.join(RAW_DIR, SAMPLE, "outs"))
adata_vis.var_names_make_unique()
adata_vis.obs["sample"] = SAMPLE

inf_aver = pd.read_csv(os.path.join(C2L_DIR, "reference", "inf_aver_98cat.csv"), index_col=0)
shared_genes = [g for g in adata_vis.var_names if g in inf_aver.index]
adata_vis = adata_vis[:, shared_genes].copy()

import scvi
scvi.settings.seed = 42
Cell2location.setup_anndata(adata_vis)
mod_vis = Cell2location(adata_vis, cell_state_df=inf_aver.loc[shared_genes, :],
                        N_cells_per_location=30, detection_alpha=20)
mod_vis.train(max_epochs=EPOCHS, batch_size=None, train_size=1,
              callbacks=[ElboLogger(os.path.join(ELBODIR, f"elbo_history_98cat_{SAMPLE}.csv"))])
adata_vis = mod_vis.export_posterior(adata_vis, sample_kwargs={"num_samples": 1000, "batch_size": adata_vis.n_obs})
adata_vis.write_h5ad(os.path.join(OUTDIR, f"cell2location_98cat_30000ep_{SAMPLE}.h5ad"))
print(SAMPLE, adata_vis.n_obs, "spots done")
