"""
Supplementary Fig. S10f (collagen VI) and S10g (TGF-beta signaling): score of stromal spots by 100-µm class of
distance to the nearest YAP-Rev high stromal spot (red) or control spot (EPCAM-high, YAP-Rev-negative stromal spot,
blue). Lines: mean of per-patient values; ribbons: SEM; shown up to 700 µm. P from 03_gradient_YAPRev_vs_control.py.
Usage: python analyses/13_stromal_gradients_S10fg/04_FigS10fg_gradients.py   (after 03_*.py COL6 and TGFB)
"""
import os

import matplotlib
matplotlib.use("pdf")
import matplotlib.pyplot as plt
import pandas as pd
from scipy.stats import wilcoxon

from common import OUT, WORK

plt.rcParams.update({"font.family": "Arial", "font.size": 9, "pdf.fonttype": 42,
                     "axes.spines.top": False, "axes.spines.right": False})
COL = {"YAP-Rev-high spots": "#D7301F", "EPCAM-high spots": "#2171B5"}
LAB = {"YAP-Rev-high spots": "YAP-Rev-high (stromal) spots", "EPCAM-high spots": "EPCAM-high, YAP-Rev-neg. stromal spots"}
PANELS = {"COL6": ("FigS10f_collagenVI_gradient", "Collagen VI", (2.6, 3.4)),
          "TGFB": ("FigS10g_TGFbeta_gradient", "TGF-β signaling", (3.4, 3.8))}

for sig, (fn, name, ylim) in PANELS.items():
    d = pd.read_csv(os.path.join(OUT, f"gradient_{sig}_summary.csv"))
    d["um"] = d.bin * 100
    d = d[d.um <= 700]
    r = pd.read_csv(os.path.join(OUT, f"gradient_{sig}_per_patient.csv")).dropna(
        subset=["YAP-Rev-high spots|diff", "EPCAM-high spots|diff"])
    P = wilcoxon(r["YAP-Rev-high spots|diff"], r["EPCAM-high spots|diff"]).pvalue
    fig, ax = plt.subplots(figsize=(3.6, 3.0))
    for ref, c in COL.items():
        s = d[d.reference == ref].sort_values("um")
        ax.fill_between(s.um, s["mean"] - s["sem"], s["mean"] + s["sem"], color=c, alpha=.15, lw=0)
        ax.plot(s.um, s["mean"], "-o", ms=3, lw=1.4, color=c, label=LAB[ref])
    ax.set_ylim(*ylim); ax.set_xlim(50, 750); ax.set_xticks(range(100, 701, 100))
    ax.set_xlabel("Distance to nearest reference spot (µm)")
    ax.set_ylabel(f"{name}\n(stromal spots, log CP10k)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper right")
    ax.set_title(f"{name} (48 patients)", fontsize=9)
    ax.text(.02, .02, f"Near minus rest of stroma,\nYAP-Rev-high vs EPCAM-high: P = {P:.2g}", transform=ax.transAxes,
            fontsize=7, va="bottom")
    fig.tight_layout()
    fig.savefig(os.path.join(WORK, "Figures", f"{fn}.pdf"))
    plt.close(fig)
