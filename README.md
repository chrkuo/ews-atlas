# Sarcoma Single-Cell Atlas

A static, no-build web viewer for scRNA-seq data, modeled on the osteosarc.com
scRNA-seq explorer: dataset tabs, UMAP colored by metadata or gene, split views,
dot plots and violins. Everything is plain HTML/JS + binary files, so it hosts
anywhere that serves static files (GitHub Pages, Netlify, S3, a lab server).

**The bundled Tumor and CD3+ PBMC datasets are synthetic placeholders.**

## Preview locally

    cd site
    python3 -m http.server 8000
    # open http://localhost:8000/#dataset=tumor

(Opening index.html directly from disk won't work — browsers block fetch() on file://.)

## Add your own data

1. Your AnnData needs: normalized/log expression in `.X` (or a layer), a 2-D
   embedding in `.obsm` (e.g. `X_umap`), and categorical `.obs` columns.
2. Export:

       pip install anndata numpy scipy pandas
       python tools/export_h5ad.py path/to/object.h5ad \
           --id tumor --name "Tumor" \
           --description "Tumor biopsies, T1–T3" \
           --obs cell_type,timepoint,sample \
           --default-color cell_type \
           --embedding X_umap            # comma-separate for more, e.g. X_umap,X_tsne
           # --layer lognorm             # if .X holds raw counts
           # --genes keep.txt            # optional gene whitelist to shrink the site

   Colors come from `adata.uns["<col>_colors"]` when present (scanpy's
   convention), so the site matches your existing figures.
3. Re-running with the same `--id` replaces that dataset. Delete
   `site/data/<id>/` and its entry in `site/data/datasets.json` to remove one.
4. Edit `site/data/datasets.json` for the site title, tab order and footer
   (remove the "synthetic" footer once real data is in).

For Seurat objects, convert first (e.g. `SeuratDisk::SaveH5Seurat` →
`Convert(..., dest = "h5ad")`, or zellkonverter/`sceasy`).

### Size guide
Each gene is one small file holding only the cells that express it (8-bit
quantized), so the browser downloads just the genes someone looks at. A 50k-cell
dataset with ~20k genes is typically 300–800 MB on disk; use `--genes` (e.g.
HVGs + markers) or `--min-cells` to trim. GitHub Pages allows ~1 GB per site;
for bigger atlases use S3/Cloudflare R2 or a university web server.

## Deploy (GitHub Pages)

    git init && git add . && git commit -m "atlas"
    git branch -M main && git remote add origin git@github.com:<you>/<repo>.git
    git push -u origin main
    # Pages only serves the repo root or /docs, so rename site/ → docs/ first
    # (and pass --out docs/data to the exporter), then:
    # Repo → Settings → Pages → Deploy from branch: main, folder /docs

## Shareable links
State lives in the URL hash, e.g.
`#dataset=tumor&gene=SPP1&split=timepoint` or `#dataset=pbmc&color=cell_type`.

## Files
    site/index.html, style.css, app.js   the viewer
    site/data/datasets.json              registry: title, footer, dataset tabs
    site/data/<id>/                       written by the exporter
    tools/export_h5ad.py                  AnnData → site data
    tools/make_demo_data.py               regenerates the synthetic demo

Note: expression values are shown as provided (no further normalization), and
dot-plot colors are scaled 0–1 per gene across groups, as in scanpy's
`standard_scale="var"`.
