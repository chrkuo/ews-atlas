#!/usr/bin/env python3
"""Turn the folder written by tools/seurat_export.R into an .h5ad.

    python tools/folder_to_h5ad.py export_tumor tumor.h5ad
"""
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.io
import anndata as ad


def main(folder, out):
    X = scipy.io.mmread(os.path.join(folder, "matrix.mtx")).T.tocsr().astype(np.float32)  # cells x genes
    genes = open(os.path.join(folder, "genes.txt")).read().split("\n")
    cells = open(os.path.join(folder, "cells.txt")).read().split("\n")
    genes = [g for g in genes if g != ""][: X.shape[1]]
    cells = [c for c in cells if c != ""][: X.shape[0]]

    meta = pd.read_csv(os.path.join(folder, "metadata.csv"), index_col=0, dtype=str, keep_default_na=False)
    meta = meta.loc[cells]
    levels_path = os.path.join(folder, "levels.json")
    levels = json.load(open(levels_path)) if os.path.exists(levels_path) else {}
    for col in meta.columns:
        cats = levels.get(col) or sorted(meta[col].unique())
        meta[col] = pd.Categorical(meta[col], categories=cats)

    emb = pd.read_csv(os.path.join(folder, "embedding.csv"), index_col=0).loc[cells]

    a = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=pd.Index(genes).astype(str)))
    a.var_names_make_unique()
    a.obsm["X_umap"] = emb.iloc[:, :2].to_numpy(np.float32)

    colors_path = os.path.join(folder, "colors.json")
    if os.path.exists(colors_path):
        for col, mapping in json.load(open(colors_path)).items():
            if col in a.obs:
                cats = list(a.obs[col].cat.categories)
                if all(c in mapping for c in cats):
                    a.uns[f"{col}_colors"] = [mapping[c] for c in cats]
                else:
                    print(f"colors.json: '{col}' doesn't cover every level; using default palette")

    a.write_h5ad(out)
    print(f"Wrote {out}: {a.n_obs:,} cells x {a.n_vars:,} genes; obs = {list(a.obs.columns)}")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
