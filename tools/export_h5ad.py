#!/usr/bin/env python3
"""
Export an AnnData (.h5ad) object into the static files the atlas viewer reads.

Usage
-----
    python tools/export_h5ad.py my_object.h5ad --id tumor --name "Tumor" \
        --obs cell_type,timepoint,sample --default-color cell_type \
        --embedding X_umap --layer counts_norm

Then add/refresh the entry in site/data/datasets.json (the script does this
for you), and the dataset appears in the viewer's dataset switcher.

Output layout (under site/data/<id>/)
-------------------------------------
    meta.json          cell count, obs fields + categories + colors, gene list,
                       embeddings available
    emb_<name>.bin     float32 [x0, y0, x1, y1, ...]
    obs_<field>.bin    uint8 or uint16 category codes, one per cell
    genes/<k>/<i>.bin  one file per gene (i = gene index, k = i // 1000):
                         uint32 nnz | float32 max | uint32[nnz] cell idx |
                         uint8[nnz] value quantized to 1..255 of max

Expression is taken from --layer if given, else .X. It should already be
normalized/log-transformed (e.g. log1p CP10k) — the viewer shows it as is.
Only categorical obs columns are exported (numeric columns are binned into
quantiles if you list them).
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
import scipy.sparse as sp

# Tableau/glasbey-like palette, readable on light and dark backgrounds
PALETTE = [  # scanpy's default 20 (so colors match scanpy plots), then extra distinct colors
    "#1f77b4", "#ff7f0e", "#279e68", "#d62728", "#aa40fc", "#8c564b", "#e377c2",
    "#b5bd61", "#17becf", "#aec7e8", "#ffbb78", "#98df8a", "#ff9896", "#c5b0d5",
    "#c49c94", "#f7b6d2", "#dbdb8d", "#9edae5", "#ad494a", "#8c6d31",
    "#FFFF00", "#1CE6FF", "#FF34FF", "#FF4A46", "#008941", "#006FA6", "#A30059",
    "#7A4900", "#0000A6", "#63FFAC", "#B79762", "#004D43", "#8FB0FF", "#997D87",
    "#5A0007", "#809693", "#1B4400", "#4FC601", "#3B5DFF", "#4A3B53", "#FF2F80",
]


def to_categorical(series, n_bins=8):
    if isinstance(series.dtype, pd.CategoricalDtype):
        return series
    if pd.api.types.is_numeric_dtype(series) and series.nunique() > 30:
        return pd.qcut(series, q=n_bins, duplicates="drop").astype(str).astype("category")
    return series.astype(str).astype("category")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("h5ad")
    ap.add_argument("--id", required=True, help="short URL-safe id, e.g. tumor")
    ap.add_argument("--name", required=True, help="display name in the switcher")
    ap.add_argument("--description", default="")
    ap.add_argument("--obs", required=True, help="comma-separated obs columns to expose")
    ap.add_argument("--default-color", help="obs column shown on load (default: first --obs)")
    ap.add_argument("--embedding", default="X_umap", help="comma-separated obsm keys (first = default)")
    ap.add_argument("--layer", default=None, help="layer with normalized expression (default .X)")
    ap.add_argument("--genes", default=None, help="optional text file of genes to keep (one per line)")
    ap.add_argument("--min-cells", type=int, default=3, help="drop genes expressed in fewer cells")
    root = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
    web = "docs" if os.path.isdir(os.path.join(root, "docs")) else "site"   # GitHub Pages uses docs/
    ap.add_argument("--out", default=os.path.join(root, web, "data"))
    args = ap.parse_args()

    import anndata as ad
    adata = ad.read_h5ad(args.h5ad)
    print(f"Loaded {adata.n_obs:,} cells x {adata.n_vars:,} genes")

    out = os.path.abspath(os.path.join(args.out, args.id))
    os.makedirs(out, exist_ok=True)

    # --- embeddings
    emb_keys = [e.strip() for e in args.embedding.split(",") if e.strip()]
    embeddings = []
    for key in emb_keys:
        if key not in adata.obsm:
            sys.exit(f"obsm['{key}'] not found; available: {list(adata.obsm.keys())}")
        xy = np.asarray(adata.obsm[key])[:, :2].astype(np.float32)
        label = key.replace("X_", "").upper()
        fn = f"emb_{label.lower()}.bin"
        xy.tofile(os.path.join(out, fn))
        embeddings.append({"name": label, "file": fn})

    # --- obs fields
    fields = []
    for col in [c.strip() for c in args.obs.split(",") if c.strip()]:
        if col not in adata.obs:
            sys.exit(f"obs['{col}'] not found; available: {list(adata.obs.columns)}")
        cat = to_categorical(adata.obs[col]).cat.remove_unused_categories()
        cats = [str(c) for c in cat.cat.categories]
        uns_colors = adata.uns.get(f"{col}_colors")
        if uns_colors is not None and len(uns_colors) == len(cats):
            colors = [str(c)[:7] for c in uns_colors]
        else:
            colors = [PALETTE[i % len(PALETTE)] for i in range(len(cats))]
        codes = cat.cat.codes.to_numpy()
        dtype = np.uint8 if len(cats) < 255 else np.uint16
        codes = np.where(codes < 0, len(cats), codes).astype(dtype)  # NaN -> extra code
        fn = f"obs_{col}.bin"
        codes.tofile(os.path.join(out, fn))
        counts = np.bincount(codes, minlength=len(cats))[: len(cats)].tolist()
        fields.append({"name": col, "file": fn, "dtype": "uint8" if dtype == np.uint8 else "uint16",
                       "categories": cats, "colors": colors, "counts": counts})

    # --- expression
    X = adata.layers[args.layer] if args.layer else adata.X
    X = sp.csc_matrix(X) if not sp.isspmatrix_csc(X) else X
    X.eliminate_zeros()
    genes = adata.var_names.astype(str).tolist()
    if args.genes:
        keep = set(open(args.genes).read().split())
        mask = np.array([g in keep for g in genes])
    else:
        mask = np.ones(len(genes), bool)
    nnz_per_gene = np.diff(X.indptr)
    mask &= nnz_per_gene >= args.min_cells
    kept = np.where(mask)[0]

    gene_names, gene_nnz = [], []
    for new_i, j in enumerate(kept):
        lo, hi = X.indptr[j], X.indptr[j + 1]
        idx = X.indices[lo:hi].astype(np.uint32)
        vals = X.data[lo:hi].astype(np.float32)
        order = np.argsort(idx)
        idx, vals = idx[order], vals[order]
        vmax = float(vals.max()) if len(vals) else 0.0
        q = np.clip(np.rint(vals / vmax * 255), 1, 255).astype(np.uint8) if vmax > 0 else np.zeros(0, np.uint8)
        d = os.path.join(out, "genes", str(new_i // 1000))
        os.makedirs(d, exist_ok=True)
        with open(os.path.join(d, f"{new_i}.bin"), "wb") as fh:
            fh.write(np.array([len(idx)], np.uint32).tobytes())
            fh.write(np.array([vmax], np.float32).tobytes())
            fh.write(idx.tobytes())
            fh.write(q.tobytes())
        gene_names.append(genes[j])
        gene_nnz.append(int(len(idx)))
        if new_i % 2000 == 0:
            print(f"  genes written: {new_i:,}/{len(kept):,}")

    default_color = args.default_color or fields[0]["name"]
    meta = {
        "id": args.id, "name": args.name, "description": args.description,
        "n_cells": int(adata.n_obs), "embeddings": embeddings, "fields": fields,
        "default_color": default_color, "genes": gene_names, "gene_nnz": gene_nnz,
    }
    with open(os.path.join(out, "meta.json"), "w") as fh:
        json.dump(meta, fh, separators=(",", ":"))

    # --- register in datasets.json
    reg_path = os.path.join(os.path.dirname(out), "datasets.json")
    reg = json.load(open(reg_path)) if os.path.exists(reg_path) else {"title": "Single-Cell Atlas", "datasets": []}
    reg["datasets"] = [d for d in reg["datasets"] if d["id"] != args.id]
    reg["datasets"].append({"id": args.id, "name": args.name, "description": args.description,
                            "n_cells": int(adata.n_obs)})
    with open(reg_path, "w") as fh:
        json.dump(reg, fh, indent=2)

    print(f"Done: {len(gene_names):,} genes, {len(fields)} obs fields -> {out}")


if __name__ == "__main__":
    main()
