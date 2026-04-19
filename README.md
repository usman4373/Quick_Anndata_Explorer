# Genetic Codon AnnData Quick Explorer

A lightweight tool to inspect `.h5ad` files quickly and generate a compact summary of:
- AnnData shape
- `obs` and `var` columns
- sample / cluster distributions
- quick UMAP preview if available

## Use cases
- sanity-check a newly received AnnData
- understand metadata before analysis
- generate a quick summary for collaborators

## Installation
```bash
pip install -r requirements.txt
```

## CLI usage
```bash
python -m src.main --input data/example.h5ad --output results
```

## Output
- `summary.txt`
- `obs_columns.csv`
- `var_columns.csv`
- `categorical_counts_<column>.csv`
- `umap_preview.png` if `X_umap` is present

## Structure
```text
geneticcodon-anndata-quick-explorer/
├── README.md
├── requirements.txt
├── .gitignore
├── LICENSE
├── data/
├── results/
└── src/
    ├── __init__.py
    ├── main.py
    └── utils.py
```

## Roadmap
- add Streamlit mini-app
- add top expressed genes summary
- add layer inspection
