# Genetic Codon — AnnData Quick Explorer

An interactive **Streamlit** app to explore `.h5ad` (**AnnData**) files quickly and visually.

---

## Overview

Single-cell datasets stored in **AnnData** format can be complex to inspect without running full analysis pipelines.

**AnnData Quick Explorer** provides a fast, interactive interface to:

- Understand dataset structure
- Explore metadata (`obs` / `var`)
- Visualize embeddings (UMAP / PCA)
- Generate quick summaries
- Export results for downstream use

All directly from your browser.

---

## Features

- Dataset summary (cells, genes, metadata counts)
- Inspect `obs` and `var` metadata
- Visualize embeddings (UMAP, PCA, etc.)
- Color embeddings by metadata columns
- Generate value counts for any `obs` column
- Plot expression of selected genes on embeddings
- Download summaries (CSV, TXT, Excel)
- Download gene-expression summary statistics
- Lightweight, fast, and easy to use

---

## Installation

### 1) Clone the repository

```bash
git clone https://github.com/GeneticCodon/Quick_Anndata_Explorer.git
cd Quick_Anndata_Explorer
```

### 2) Install dependencies

```bash
pip install -r requirements.txt
```

---

## Usage

Run the Streamlit app:

```bash
streamlit run app.py
```

Then:

1. Open the browser (Streamlit may auto-launch).
2. Upload your `.h5ad` file.
3. Explore interactively.

---

## Demo dataset (recommended)

For best results, use a processed dataset that already contains embeddings:

```python
import scanpy as sc

adata = sc.datasets.pbmc68k_reduced()
adata.write("demo_pbmc68k_reduced.h5ad")
```

Upload `demo_pbmc68k_reduced.h5ad` in the app to see:

- Embedding plots
- Metadata exploration
- Richer outputs

---

## Project structure

```text
.
├── app.py              # Streamlit application
├── requirements.txt    # Dependencies
├── README.md
├── data/               # Optional input files
├── results/            # Optional outputs
└── src/                # (Optional) CLI version
```

---

## Use cases

- Quick inspection of `.h5ad` datasets
- Teaching and demonstrations
- Data validation before analysis
- Debugging metadata issues
- Rapid exploration for collaborators

---

## Why this tool?

Most workflows require loading full pipelines just to understand a dataset.

This tool provides immediate insight into AnnData objects with minimal overhead.

---

## Requirements

- Python 3.8+
- Streamlit
- Scanpy / AnnData
- Pandas
- Matplotlib


---

## About Genetic Codon

Genetic Codon is a virtual bioinformatics lab focused on:

- Multi-omics analysis
- Spatial transcriptomics
- Reproducible research pipelines
- Training and mentorship

- Website: www.geneticcodon.com

---

## Support

If you find this useful:

- Star the repository
- Share with your network
- Contribute improvements
