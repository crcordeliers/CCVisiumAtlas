"""Shared inputs of the stromal-gradient analyses (Supplementary Fig. S10f,g)."""
import glob
import os

import numpy as np
import pandas as pd

REPO = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
WORK = os.environ.get("CRCVISIUM_WORKDIR", os.path.expanduser("~/data"))
GRADIENT_DIR = os.environ.get("CRCVISIUM_GRADIENT", os.path.join(WORK, "gradients"))
YAPREV_DIR = os.environ.get("CRCVISIUM_YAPREV", os.path.expanduser("~/data/yaprev_spots"))
OUT = os.path.join(WORK, "Figures", "gradients")
os.makedirs(OUT, exist_ok=True)
YAP_THR = 0.40
TUMOR_MC = {"0", "1", "7"}   # tumour meta-clusters used to select stroma-expressed genes (Tumor-fibroblasts excluded)


def relapse_groups():
    clin = pd.read_csv(os.path.join(REPO, "data", "clinical_48.csv"))
    return dict(zip(clin["sample"], np.where(clin.recurr_cens == 1, "Relapse", "No Relapse")))


def load_sample_tables(columns=None):
    """Yields (sample, table) with raw gene counts, nCount, x, y, YAPRev_raw, stroma/tumor flags for each sample."""
    labels = pd.read_csv(os.path.join(WORK, "Tables", "spot_labels_all48.csv"), dtype={"merged_meta_clusters": str})
    for f in sorted(glob.glob(os.path.join(GRADIENT_DIR, "*_genesets.csv.gz"))):
        x = pd.read_csv(f, usecols=columns)
        s = x["sample"].iloc[0]
        xy = pd.read_csv(os.path.join(YAPREV_DIR, "merged_spots", f"{s}_merged.csv"), usecols=["barcode", "x", "y", "YAPRev_raw"])
        lab = labels[labels["sample"] == s][["barcode", "merged_meta_clusters", "stromal_subcluster"]]
        x = x.merge(xy, on="barcode").merge(lab, on="barcode")
        x["stroma"] = x.stromal_subcluster.notna()
        x["tumor"] = x.merged_meta_clusters.isin(TUMOR_MC)
        yield s, x
