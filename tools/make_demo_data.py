#!/usr/bin/env python3
"""Build two SYNTHETIC demo AnnData objects (Tumor, PBMC) so the viewer has
something to show. Cell types, markers and layout are simulated — not real data."""
import os
import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

rng = np.random.default_rng(7)
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "demo_h5ad")
os.makedirs(OUT, exist_ok=True)

FILLER = [f"GENE{i:04d}" for i in range(400)]


def simulate(types, timepoints, tp_shift, n_total, extra_markers):
    # types: list of (name, weight, center(x,y), spread, markers)
    w = np.array([t[1] for t in types], float); w /= w.sum()
    labels = rng.choice(len(types), size=n_total, p=w)
    tps = rng.choice(len(timepoints), size=n_total)
    # timepoint shifts proportions a little
    for k, (tname, shift) in enumerate(tp_shift.items()):
        ti = [t[0] for t in types].index(tname)
        for tp_i, frac in enumerate(shift):
            m = (tps == tp_i) & (rng.random(n_total) < frac)
            labels[m] = ti

    markers = sorted({g for t in types for g in t[4]} | set(extra_markers))
    genes = markers + FILLER
    gidx = {g: i for i, g in enumerate(genes)}
    G = len(genes)

    xy = np.zeros((n_total, 2), np.float32)
    rows, cols, vals = [], [], []
    base = rng.gamma(0.3, 0.08, size=G)  # background rate per gene
    for ti, (name, _, c, s, mk) in enumerate(types):
        idx = np.where(labels == ti)[0]
        n = len(idx)
        # elongated blob with a little sub-structure
        ang = rng.uniform(0, np.pi)
        cov = np.array([[s * 0.9, 0], [0, s * 0.5]])
        R = np.array([[np.cos(ang), -np.sin(ang)], [np.sin(ang), np.cos(ang)]])
        pts = rng.normal(size=(n, 2)) @ (R @ cov)
        pts += np.array(c) + rng.normal(scale=s * 0.15, size=(n, 2))
        xy[idx] = pts
        rate = np.tile(base, (n, 1))
        for g in mk:
            rate[:, gidx[g]] += rng.gamma(3, 1.2) * (1 + 0.6 * np.tanh(pts[:, 0] - c[0]))
        lib = rng.lognormal(0, 0.35, size=(n, 1))
        counts = rng.poisson(rate * lib)
        r, cc = np.nonzero(counts)
        rows.append(idx[r]); cols.append(cc); vals.append(counts[r, cc])
    X = sp.csr_matrix((np.concatenate(vals).astype(np.float32),
                       (np.concatenate(rows), np.concatenate(cols))), shape=(n_total, G))
    # log1p CP10k
    tot = np.asarray(X.sum(1)).ravel(); tot[tot == 0] = 1
    X = sp.diags(1e4 / tot) @ X
    X.data = np.log1p(X.data)
    obs = pd.DataFrame({
        "cell_type": pd.Categorical([types[i][0] for i in labels], categories=[t[0] for t in types]),
        "timepoint": pd.Categorical([timepoints[i] for i in tps], categories=timepoints),
    }, index=[f"cell{i}" for i in range(n_total)])
    a = ad.AnnData(X=X.tocsr(), obs=obs, var=pd.DataFrame(index=genes))
    a.obsm["X_umap"] = xy
    return a


tumor_types = [
    ("Osteoblastic tumor", 30, (-4, 3), 1.3, ["RUNX2", "SP7", "COL1A1", "ALPL", "SPP1", "MKI67"]),
    ("Chondroblastic tumor", 8, (-7, 0), 0.8, ["SOX9", "COL2A1", "ACAN", "COL1A1"]),
    ("Cycling tumor", 7, (-1.5, 5.5), 0.7, ["MKI67", "TOP2A", "RUNX2", "CDK1"]),
    ("CAF", 8, (-5, -3.5), 0.9, ["COL1A1", "COL3A1", "DCN", "LUM", "THBS2", "PDGFRB"]),
    ("Osteoclast", 5, (2, 4), 0.6, ["ACP5", "CTSK", "MMP9", "CD68"]),
    ("TAM (SPP1+)", 10, (4.5, 1.5), 0.9, ["CD68", "CD163", "SPP1", "APOE", "C1QA", "MIF"]),
    ("Monocyte", 4, (6.5, 3), 0.6, ["CD14", "LYZ", "S100A8", "FCN1"]),
    ("CD8 T", 6, (5.5, -3), 0.7, ["CD3E", "CD8A", "GZMB", "NKG7", "PDCD1"]),
    ("CD4 T / Treg", 4, (3.5, -4.5), 0.6, ["CD3E", "CD4", "FOXP3", "IL7R"]),
    ("NK", 3, (7.5, -1.5), 0.5, ["NKG7", "GNLY", "KLRD1"]),
    ("B / plasma", 2, (1, -6.5), 0.5, ["MS4A1", "CD79A", "JCHAIN"]),
    ("Endothelial", 5, (-1, -5), 0.7, ["PECAM1", "VWF", "CDH5", "KDR"]),
]
tumor = simulate(tumor_types, ["T1", "T2", "T3"],
                 {"TAM (SPP1+)": [0.0, 0.06, 0.12], "CD8 T": [0.08, 0.02, 0.0]},
                 30000, ["CD74", "B2M", "ACTB"])
tumor.write_h5ad(os.path.join(OUT, "tumor_demo.h5ad"))

pbmc_types = [
    ("Naive CD4", 20, (-4, 2), 1.0, ["CD3E", "CD4", "CCR7", "SELL", "IL7R"]),
    ("Memory CD4", 15, (-2, 4.5), 0.9, ["CD3E", "CD4", "IL7R", "S100A4"]),
    ("Treg", 5, (-6, 4.5), 0.6, ["CD3E", "CD4", "FOXP3", "IL2RA"]),
    ("Naive CD8", 12, (3, 3), 0.8, ["CD3E", "CD8A", "CCR7", "SELL"]),
    ("Effector CD8", 18, (5, 0), 1.0, ["CD3E", "CD8A", "GZMB", "GZMK", "NKG7"]),
    ("Exhausted-like CD8", 6, (6.5, -3), 0.6, ["CD3E", "CD8A", "PDCD1", "HAVCR2", "TOX"]),
    ("MAIT", 4, (1, -3.5), 0.5, ["CD3E", "SLC4A10", "KLRB1"]),
    ("γδ T", 3, (-1, -5), 0.5, ["CD3E", "TRDC", "TRGC1"]),
    ("Cycling T", 2, (-4, -2.5), 0.5, ["CD3E", "MKI67", "TOP2A"]),
]
pbmc = simulate(pbmc_types, [f"P{i}" for i in range(1, 9)],
                {"Effector CD8": [0, 0, .05, .1, .15, .1, .05, 0]}, 20000, ["B2M", "ACTB"])
pbmc.write_h5ad(os.path.join(OUT, "pbmc_demo.h5ad"))
print("demo h5ad written:", OUT)
