"""
Stromal gradients, step 2 - gene selection and comparison with random stromal anchors (Methods "Spatial gradient").

1. Signature genes restricted to those more expressed in stromal than in tumour spots (pooled 48 sections; tumour =
   Tumor-cycling, -intermediate, -infiltrated), to avoid epithelial signal. Gives the 33 stroma-expressed genes of
   HALLMARK_TGF_BETA_SIGNALING used in S10g (written to stromal_genes.tsv).
2. Per patient: score of stromal spots <= 2.5 spot units (<= 250 µm) from the nearest YAP-Rev high spot minus that of
   stromal spots >= 3.5 units (>= 350 µm), compared with the same difference around the same number of randomly chosen
   stromal spots (200 draws; paired Wilcoxon test over patients).
Score per spot = log(1 + CP10k of the summed raw counts of the genes).
Inputs : see common.py (01_extract_counts_visium48.R, YAP-Rev spot tables, spot labels, data/clinical_48.csv)
Outputs: <WORKDIR>/Figures/gradients/stromal_genes.tsv, near_vs_rest_random_anchors_per_patient.csv, random_anchor_stats.csv
Usage  : python analyses/13_stromal_gradients_S10fg/02_stromal_genes_and_random_anchors.py
"""
import os

import numpy as np
import pandas as pd
from scipy.spatial import cKDTree
from scipy.stats import wilcoxon

from common import OUT, REPO, YAP_THR, load_sample_tables

gs = pd.read_csv(os.path.join(REPO, "analyses", "13_stromal_gradients_S10fg", "genesets_MSigDB.tsv"), sep="\t",
                 header=None, names=["set", "genes"])
SETS = {"Collagen deposition": gs.genes[gs.set == "NABA_COLLAGENS"].iloc[0].split(","),
        "TGF-beta signaling": gs.genes[gs.set == "HALLMARK_TGF_BETA_SIGNALING"].iloc[0].split(","),
        "Collagen VI": ["COL6A1", "COL6A2", "COL6A3"]}

D = pd.concat([x for _, x in load_sample_tables()], ignore_index=True)

KEEP = {}
for k, g in SETS.items():
    g = [x for x in g if x in D.columns]
    cs = D.loc[D.stroma, g].sum() / D.loc[D.stroma, "nCount"].sum()
    ct = D.loc[D.tumor, g].sum() / D.loc[D.tumor, "nCount"].sum()
    KEEP[k] = [x for x in g if cs[x] > ct[x] and D.loc[D.stroma, x].sum() > 0]
    D[k] = np.log1p(1e4 * D[KEEP[k]].sum(axis=1) / D.nCount)
    print(f"{k}: {len(KEEP[k])}/{len(g)} genes")
pd.DataFrame({"signature": list(KEEP), "genes": [",".join(v) for v in KEEP.values()]}).to_csv(
    os.path.join(OUT, "stromal_genes.tsv"), sep="\t", index=False)

rng = np.random.default_rng(53)
rows = []
for s, x in D.groupby("sample"):
    P = x[["x", "y"]].values
    sp = np.median(cKDTree(P).query(P, k=2)[0][:, 1])            # spot pitch in pixels
    hi = (x.YAPRev_raw.values >= YAP_THR) & x.stroma.values
    st = x.stroma.values & ~hi
    if hi.sum() < 3:
        continue
    dist = cKDTree(P[hi]).query(P)[0] / sp
    r = dict(sample=s, n_yap_high=int(hi.sum()))
    for k in SETS:
        v = x[k].values
        near, rest = st & (dist <= 2.5), st & (dist >= 3.5)
        if near.sum() < 5 or rest.sum() < 5:
            continue
        cand, nul = np.where(st)[0], []
        for _ in range(200):
            a = rng.choice(cand, size=min(hi.sum(), len(cand) - 1), replace=False)
            dd = cKDTree(P[a]).query(P)[0] / sp
            m = st.copy(); m[a] = False
            n2, r2 = m & (dd <= 2.5), m & (dd >= 3.5)
            if n2.sum() >= 5 and r2.sum() >= 5:
                nul.append(v[n2].mean() - v[r2].mean())
        r.update({f"{k}|diff": v[near].mean() - v[rest].mean(), f"{k}|null": np.mean(nul)})
    rows.append(r)
R = pd.DataFrame(rows)
R.to_csv(os.path.join(OUT, "near_vs_rest_random_anchors_per_patient.csv"), index=False)
S = pd.DataFrame([dict(signature=k, n=len(R.dropna(subset=[f"{k}|diff"])),
                       P_vs_random=wilcoxon(R.dropna(subset=[f"{k}|diff"])[f"{k}|diff"],
                                            R.dropna(subset=[f"{k}|diff"])[f"{k}|null"]).pvalue) for k in SETS])
S.to_csv(os.path.join(OUT, "random_anchor_stats.csv"), index=False)
print(S)
