#!/usr/bin/env python3
import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp
from scipy.io import mmwrite
import gzip
import shutil
import os


def _ensure_csr(X):
    if sp.issparse(X):
        return X.tocsr()
    # If dense, convert to sparse CSR (may be large!)
    return sp.csr_matrix(X)


def _gzip_file(src_path: Path, dst_path: Path):
    with open(src_path, "rb") as f_in, gzip.open(dst_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    src_path.unlink()  # remove uncompressed


def main():
    p = argparse.ArgumentParser(description="Export AnnData layer to 10X MTX (gz) + obs metadata CSV")
    p.add_argument("--in", dest="inp", required=True, help="Input .h5ad")
    p.add_argument("--outdir", required=True, help="Output folder for 10X MTX")
    p.add_argument("--layer", default="raw_counts", help="Layer to export; if missing uses X")
    p.add_argument("--meta_csv", required=True, help="Output CSV for adata.obs (metadata)")
    args = p.parse_args()

    inp = Path(args.inp).expanduser().resolve()
    outdir = Path(args.outdir).expanduser().resolve()
    outdir.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(inp)

    # choose matrix
    if args.layer in adata.layers:
        X = adata.layers[args.layer]
    else:
        X = adata.X

    X = _ensure_csr(X)

    # 10X expects genes x cells in matrix.mtx
    # AnnData is cells x genes, so transpose
    M = X.T.tocoo()

    # --- barcodes.tsv.gz ---
    barcodes = pd.Series(adata.obs_names.astype(str))
    barcodes_path = outdir / "barcodes.tsv"
    barcodes.to_csv(barcodes_path, index=False, header=False, sep="\t")
    _gzip_file(barcodes_path, outdir / "barcodes.tsv.gz")

    # --- features.tsv.gz ---
    # 10X v3 features.tsv columns: gene_id, gene_name, feature_type
    var = adata.var.copy()
    gene_id = var["ensembl_id"].astype(str) if "ensembl_id" in var.columns else pd.Series(var.index.astype(str))
    gene_name = var["gene_symbol"].astype(str) if "gene_symbol" in var.columns else pd.Series(var.index.astype(str))
    feature_type = pd.Series(["Gene Expression"] * adata.n_vars)

    features = pd.DataFrame({0: gene_id.values, 1: gene_name.values, 2: feature_type.values})
    features_path = outdir / "features.tsv"
    features.to_csv(features_path, index=False, header=False, sep="\t")
    _gzip_file(features_path, outdir / "features.tsv.gz")

    # --- matrix.mtx.gz ---
    mtx_path = outdir / "matrix.mtx"
    mmwrite(str(mtx_path), M)  # writes coordinate format
    _gzip_file(mtx_path, outdir / "matrix.mtx.gz")

    # --- metadata CSV ---
    adata.obs.to_csv(args.meta_csv)

    print("Wrote 10X MTX (gz) to:", outdir)
    print("Wrote metadata CSV to:", args.meta_csv)
    print("Matrix shape (genes x cells):", (adata.n_vars, adata.n_obs))
    print("Nonzeros:", M.nnz)


if __name__ == "__main__":
    main()
