from pathlib import Path
from typing import Iterable

import anndata as ad
import pandas as pd


def ensure_output_dir(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)


def build_summary_lines(adata: ad.AnnData) -> list[str]:
    summary_lines = [
        "Genetic Codon AnnData Quick Explorer Summary",
        "=" * 45,
        f"n_obs: {adata.n_obs}",
        f"n_vars: {adata.n_vars}",
        f"obs columns: {list(adata.obs.columns)}",
        f"var columns: {list(adata.var.columns)}",
        f"obsm keys: {list(adata.obsm.keys())}",
        f"layers: {list(adata.layers.keys())}",
        f"uns keys: {list(adata.uns.keys())}",
    ]

    if "X_umap" in adata.obsm.keys():
        summary_lines.append("umap_available: True")
    else:
        summary_lines.append("umap_available: False")

    return summary_lines


def save_summary_txt(output_file: Path, lines: Iterable[str]) -> None:
    output_file.write_text("\n".join(lines), encoding="utf-8")


def save_list_csv(output_file: Path, column_name: str, values: list[str]) -> None:
    df = pd.DataFrame({column_name: values})
    df.to_csv(output_file, index=False)
