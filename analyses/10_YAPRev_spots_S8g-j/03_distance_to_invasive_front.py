"""
Distance of every Visium spot to the annotated invasive front (Supplementary Fig. S8i, Methods "Spatial proximity").

The invasive front was drawn by the pathologist in QuPath on 31/48 slides (23 No relapse, 8 Relapse) and exported as
one polyline per slide (<YAPREV_DIR>/fi_lines/<sample>_lines.json, QuPath full-resolution pixel space). For each spot,
the Euclidean point-to-segment distance to the nearest segment of the front is computed and converted to µm with the
pixel size of the slide (<YAPREV_DIR>/pixel_sizes.csv: sample, µm per pixel).
Spot coordinates and YAPRev_raw come from 02_yaprev_pool_export.R (<YAPREV_DIR>/merged_spots/<sample>_merged.csv,
same pixel space as the QuPath project).

Output: <YAPREV_DIR>/spot_distances/<sample>.csv (no header: sample, YAPRev_raw, distance_um; rows in merged_spots order)
Usage : python analyses/10_YAPRev_spots_S8g-j/03_distance_to_invasive_front.py
"""
import json
import os

import numpy as np
import pandas as pd

YAPREV_DIR = os.environ.get("CRCVISIUM_YAPREV", os.path.expanduser("~/data/yaprev_spots"))
OUT_DIR = os.environ.get("CRCVISIUM_FRONT_OUT", os.path.join(YAPREV_DIR, "spot_distances"))
os.makedirs(OUT_DIR, exist_ok=True)

px = pd.read_csv(os.path.join(YAPREV_DIR, "pixel_sizes.csv"), header=None, names=["sample", "um_per_px"])
px = dict(zip(px["sample"], px["um_per_px"]))


def point_to_segments(xy, seg_a, seg_b):
    """Minimum distance from each point (n x 2) to a set of segments [a, b] (m x 2 each)."""
    ab = seg_b - seg_a                                   # m x 2
    ab2 = np.maximum((ab ** 2).sum(1), 1e-12)            # m
    best = np.full(len(xy), np.inf)
    for start in range(0, len(xy), 5000):                # chunks to bound memory
        p = xy[start:start + 5000, None, :]              # c x 1 x 2
        t = np.clip(((p - seg_a) * ab).sum(2) / ab2, 0, 1)   # c x m
        proj = seg_a + t[..., None] * ab
        best[start:start + 5000] = np.sqrt(((p - proj) ** 2).sum(2)).min(1)
    return best


for f in sorted(os.listdir(os.path.join(YAPREV_DIR, "fi_lines"))):
    sample = f.replace("_lines.json", "")
    lines = json.load(open(os.path.join(YAPREV_DIR, "fi_lines", f)))["lines"]
    seg_a = np.vstack([np.asarray(l["points"], float)[:-1] for l in lines])
    seg_b = np.vstack([np.asarray(l["points"], float)[1:] for l in lines])
    m = pd.read_csv(os.path.join(YAPREV_DIR, "merged_spots", f"{sample}_merged.csv"))
    d_um = point_to_segments(m[["x", "y"]].to_numpy(float), seg_a, seg_b) * px[sample]
    pd.DataFrame({"sample": sample, "YAPRev_raw": m["YAPRev_raw"], "distance_um": d_um}).to_csv(
        os.path.join(OUT_DIR, f"{sample}.csv"), header=False, index=False)
    print(sample, len(m), "spots")
