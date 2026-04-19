import argparse
from pathlib import Path

import anndata as ad

from utils import (
    build_summary_lines,
    ensure_output_dir,
    save_list_csv,
    save_summary_txt,
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Quickly inspect an AnnData (.h5ad) file and export a compact summary."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to input .h5ad file",
    )
    parser.add_argument(
        "--output-dir",
        default="results",
        help="Directory to save summary outputs (default: results)",
    )
    parser.add_argument(
        "--count-columns",
        nargs="*",
        default=[],
        help="Optional obs columns for value counts, e.g. --count-columns sample cell_type",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    output_dir = Path(args.output_dir)

    if not input_path.exists():
        raise FileNotFoundError(f"Input file not found: {input_path}")

    ensure_output_dir(output_dir)

    print(f"[INFO] Reading AnnData file: {input_path}")
    adata = ad.read_h5ad(input_path)

    summary_lines = build_summary_lines(adata)

    save_summary_txt(output_dir / "summary.txt", summary_lines)
    save_list_csv(output_dir / "obs_columns.csv", "obs_column", list(adata.obs.columns))
    save_list_csv(output_dir / "var_columns.csv", "var_column", list(adata.var.columns))
    save_list_csv(output_dir / "obsm_keys.csv", "obsm_key", list(adata.obsm.keys()))
    save_list_csv(output_dir / "layers.csv", "layer", list(adata.layers.keys()))

    if args.count_columns:
        for col in args.count_columns:
            if col in adata.obs.columns:
                counts = adata.obs[col].astype(str).value_counts(dropna=False)
                count_lines = [f"{idx},{val}" for idx, val in counts.items()]
                save_summary_txt(
                    output_dir / f"{col}_value_counts.csv",
                    ["value,count", *count_lines],
                )
                print(f"[INFO] Saved counts for obs column: {col}")
            else:
                print(f"[WARN] Skipping '{col}' because it is not present in adata.obs")

    print(f"[DONE] Summary files saved to: {output_dir.resolve()}")


if __name__ == "__main__":
    main()
