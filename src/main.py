from __future__ import annotations
import argparse
import os
import anndata as ad
import matplotlib.pyplot as plt
import pandas as pd
from .utils import save_dataframe

def parse_args():
    parser = argparse.ArgumentParser(description="Quickly inspect an AnnData file.")
    parser.add_argument("--input", required=True, help="Path to .h5ad file")
    parser.add_argument("--output", required=True, help="Output directory")
    parser.add_argument("--count-columns", nargs="*", default=["sample", "leiden", "cell_type"])
    return parser.parse_args()

def main():
    args = parse_args()
    os.makedirs(args.output, exist_ok=True)

    adata = ad.read_h5ad(args.input)

    summary_lines = [
        f"n_obs: {adata.n_obs}",
        f"n_vars: {adata.n_vars}",
        f"obs columns: {list(adata.obs.columns)}",
        f"var columns: {list(adata.var.columns)}",
        f"obsm keys: {list(adata.obsm.keys())}",
        f"layers: {list(adata.layers.keys())}",
    ]

    with open(os.path.join(args.output, "summary.txt"), "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))

    save_dataframe(pd.DataFrame({"obs_columns": list(adata.obs.columns)}),
                   os.path.join(args.output, "obs_columns.csv"))
    save_dataframe(pd.DataFrame({"var_columns": list(adata.var.columns)}),
                   os.path.join(args.output, "var_columns.csv"))

    for col in args.count_columns:
        if col in adata.obs.columns:
            counts = adata.obs[col].astype(str).value_counts(dropna=False).reset_index()
            counts.columns = [col, "count"]
            save_dataframe(counts, os.path.join(args.output, f"categorical_counts_{col}.csv"))

    if "X_umap" in adata.obsm:
        coords = adata.obsm["X_umap"]
        plt.figure(figsize=(6, 5))
        plt.scatter(coords[:, 0], coords[:, 1], s=3)
        plt.xlabel("UMAP1")
        plt.ylabel("UMAP2")
        plt.title("UMAP preview")
        plt.tight_layout()
        plt.savefig(os.path.join(args.output, "umap_preview.png"), dpi=300)
        plt.close()

if __name__ == "__main__":
    main()
