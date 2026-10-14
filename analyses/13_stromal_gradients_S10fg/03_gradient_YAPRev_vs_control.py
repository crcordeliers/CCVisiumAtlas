"""
Stromal gradients, step 3 - collagen VI (S10f) and TGF-beta signaling (S10g) around YAP-Rev high vs control spots.

Reference spots per patient:
  YAP-Rev high = stromal spots with YAPRev_raw >= 0.40 (as in the CellChat analyses);
  control      = the same number of other stromal spots with the highest EPCAM expression (CP10k) among stromal spots
                 with a YAP-Rev score below the stromal median of the section (stromal spots with YAP-Rev-negative
                 tumour cells).
Score in the remaining stromal spots: collagen VI = log(1 + CP10k of COL6A1 + COL6A2 + COL6A3); TGF-beta = same with
the 33 stroma-expressed genes of HALLMARK_TGF_BETA_SIGNALING (02_*.py). Distance to the nearest reference spot in spot
units (distance / median centre-to-centre spacing), classes 1-10 (x 100 µm); per-patient means of classes with >= 5 spots.
Test: per patient, score <= 2.5 units (<= 250 µm) minus score >= 3.5 units (>= 350 µm), YAP-Rev high vs control,
paired Wilcoxon test.
Outputs: <WORKDIR>/Figures/gradients/gradient_<SIG>_curves_per_patient.csv, gradient_<SIG>_per_patient.csv,
         gradient_<SIG>_summary.csv (mean, SEM per class and reference)
Usage  : python analyses/13_stromal_gradients_S10fg/03_gradient_YAPRev_vs_control.py COL6|TGFB
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

from common import GRADIENT_DIR, OUT, YAP_THR, load_sample_tables, relapse_groups

SIG = sys.argv[1] if len(sys.argv) > 1 else "COL6"
if SIG == "COL6":
    GENES = ["COL6A1", "COL6A2", "COL6A3"]
else:
    sg = pd.read_csv(os.path.join(OUT, "stromal_genes.tsv"), sep="\t")
    GENES = sg.genes[sg.signature == "TGF-beta signaling"].iloc[0].split(",")
    assert len(GENES) == 33
grp = relapse_groups()
EP = pd.read_csv(os.path.join(GRADIENT_DIR, "visium48_EPCAM_KRT20.csv.gz"))

rows, stats = [], []
for s, x in load_sample_tables(columns=["barcode", "sample", "nCount"] + GENES):
    x = x.merge(EP[EP["sample"] == s][["barcode", "EPCAM"]], on="barcode")
    st = x.stroma.values
    yh = st & (x.YAPRev_raw.values >= YAP_THR)
    n = int(yh.sum())
    if n < 3:
        continue
    ep = 1e4 * x.EPCAM.values / x.nCount.values
    cand = np.where(st & ~yh & (x.YAPRev_raw.values < np.median(x.YAPRev_raw.values[st])))[0]
    eh = np.zeros(len(x), bool)
    eh[cand[np.argsort(-ep[cand], kind="stable")[:n]]] = True
    score = np.log1p(1e4 * x[GENES].sum(axis=1).values / x.nCount.values)
    P = x[["x", "y"]].values
    sp = np.median(cKDTree(P).query(P, k=2)[0][:, 1])
    tgt = st & ~yh & ~eh
    r = dict(sample=s, group=grp[s], n_ref=n)
    for ref, m in [("YAP-Rev-high spots", yh), ("EPCAM-high spots", eh)]:
        d = cKDTree(P[m]).query(P[tgt])[0] / sp
        v = score[tgt]
        b = np.clip(np.round(d).astype(int), 1, 10)
        k = d <= 10.5
        bm = pd.DataFrame({"bin": b[k], "v": v[k]}).groupby("bin").v.agg(["mean", "size"])
        for bi, rr in bm[bm["size"] >= 5].iterrows():
            rows.append(dict(sample=s, group=grp[s], reference=ref, bin=bi, mean=rr["mean"]))
        near, rest = v[d <= 2.5], v[d >= 3.5]
        r[f"{ref}|diff"] = near.mean() - rest.mean() if len(near) >= 5 and len(rest) >= 5 else np.nan
    stats.append(r)

C, R = pd.DataFrame(rows), pd.DataFrame(stats)
C.to_csv(os.path.join(OUT, f"gradient_{SIG}_curves_per_patient.csv"), index=False)
R.to_csv(os.path.join(OUT, f"gradient_{SIG}_per_patient.csv"), index=False)
summ = C.groupby(["reference", "bin"])["mean"].agg(["mean", "std", "count"]).reset_index()
summ["sem"] = summ["std"] / np.sqrt(summ["count"])
summ.to_csv(os.path.join(OUT, f"gradient_{SIG}_summary.csv"), index=False)
x = R.dropna(subset=["YAP-Rev-high spots|diff", "EPCAM-high spots|diff"])
pw = wilcoxon(x["YAP-Rev-high spots|diff"], x["EPCAM-high spots|diff"]).pvalue
nyes = int((x["YAP-Rev-high spots|diff"] > x["EPCAM-high spots|diff"]).sum())
print(f"{SIG}: higher around YAP-Rev high in {nyes}/{len(x)} patients; paired Wilcoxon P = {pw:.2g}")
