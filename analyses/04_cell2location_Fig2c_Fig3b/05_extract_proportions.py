"""
Cell2location, step 5: per-spot cell-state proportions for the 48 samples.

Proportion = q05 cell abundance / sum of q05 abundances of the spot (sums to 1 per spot).
Input : <C2L_DIR>/objects/cell2location_98cat_30000ep_<SAMPLE>.h5ad (step 4)
Output: <C2L_DIR>/tables/c2l_98cat_proportions_<SAMPLE>.csv (barcode x 98 categories) and
        <C2L_DIR>/tables/c2l_98cat_proportions_all48_long.csv
Usage : python analyses/04_cell2location_Fig2c_Fig3b/05_extract_proportions.py
"""
import os

import pandas as pd
import scanpy as sc

C2L_DIR = os.environ.get("CRCVISIUM_C2L", os.path.expanduser("~/data/cell2location"))
SAMPLES_FILE = os.environ.get("CRCVISIUM_SAMPLES", "data/samples_48.txt")
OBJDIR, OUTDIR = os.path.join(C2L_DIR, "objects"), os.path.join(C2L_DIR, "tables")
os.makedirs(OUTDIR, exist_ok=True)

all_long = []
for sample in open(SAMPLES_FILE).read().split():
    adata = sc.read_h5ad(os.path.join(OBJDIR, f"cell2location_98cat_30000ep_{sample}.h5ad"))
    props = adata.obsm["q05_cell_abundance_w_sf"].copy()
    props.columns = [c.replace("q05cell_abundance_w_sf_", "") for c in props.columns]
    props = props.div(props.sum(axis=1), axis=0)
    props.index.name = "barcode"
    props.to_csv(os.path.join(OUTDIR, f"c2l_98cat_proportions_{sample}.csv"))
    long = props.reset_index().melt(id_vars="barcode", var_name="celltype", value_name="proportion")
    long.insert(0, "sample", sample)
    all_long.append(long)

pd.concat(all_long, ignore_index=True).to_csv(os.path.join(OUTDIR, "c2l_98cat_proportions_all48_long.csv"), index=False)
