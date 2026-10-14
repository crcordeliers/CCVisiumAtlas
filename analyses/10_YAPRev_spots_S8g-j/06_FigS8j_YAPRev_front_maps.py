"""
Supplementary Fig. S8i (right panels) / S8j - YAP-Rev high spots (red, YAPRev_raw >= 0.40) and invasive front (orange
line) in a representative No relapse (#21-77_B1) and Relapse (#395_C1) sample. One PDF per sample.

Inputs: <YAPREV_DIR>/merged_spots/<sample>_merged.csv (02_yaprev_pool_export.R), <YAPREV_DIR>/fi_lines/<sample>_lines.json
Output: <WORKDIR>/Figures/FigS8j_YAPRev_front_<sample>.pdf
Usage : python analyses/10_YAPRev_spots_S8g-j/06_FigS8j_YAPRev_front_maps.py [sample ...]
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

YAPREV_DIR = os.environ.get("CRCVISIUM_YAPREV", os.path.expanduser("~/data/yaprev_spots"))
OUT_DIR = os.path.join(os.environ.get("CRCVISIUM_WORKDIR", os.path.expanduser("~/data")), "Figures")
os.makedirs(OUT_DIR, exist_ok=True)
SAMPLES = sys.argv[1:] or ["Ech_21-77_B1", "Ech_395_C1"]
YAP_HIGH_THR = 0.40

for sample in SAMPLES:
    df = pd.read_csv(os.path.join(YAPREV_DIR, "merged_spots", f"{sample}_merged.csv"))
    fig, ax = plt.subplots(figsize=(8, 8))
    ax.scatter(df["x"], df["y"], c="#D9D9D9", s=10, linewidths=0, alpha=0.9)
    high = df[df["YAPRev_raw"] >= YAP_HIGH_THR]
    ax.scatter(high["x"], high["y"], c="red", s=16, linewidths=0, alpha=0.95)
    lines_path = os.path.join(YAPREV_DIR, "fi_lines", f"{sample}_lines.json")
    if os.path.exists(lines_path):
        for ln in json.load(open(lines_path)).get("lines", []):
            ax.plot([p[0] for p in ln["points"]], [p[1] for p in ln["points"]], color="#FF8C00", linewidth=2.5,
                    solid_capstyle="round")
    ax.set_title("#" + sample.replace("Ech_", ""), fontsize=14)
    ax.invert_yaxis()
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    fig.savefig(os.path.join(OUT_DIR, f"FigS8j_YAPRev_front_{sample}.pdf"), bbox_inches="tight")
    plt.close(fig)
