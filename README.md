#  Genetic Codon — AnnData Quick Explorer

An interactive Streamlit app to explore `.h5ad` (AnnData) files quickly and visually.

---

## Overview

Single-cell datasets stored in AnnData format can be complex to inspect without running full analysis pipelines.

**AnnData Quick Explorer** provides a fast, interactive interface to:

- Understand dataset structure
- Explore metadata
- Visualize embeddings (UMAP/PCA)
- Generate quick summaries
- Export results for downstream use

All directly from your browser.

---

##  Features

- Dataset summary (cells, genes, metadata counts)
- Inspect `obs` and `var` metadata
- Visualize embeddings (UMAP, PCA, etc.)
- Color embeddings by metadata columns
- Generate value counts for any `obs` column
- Download summaries (CSV, TXT, Excel)
- Lightweight, fast, and easy to use
- Plot expression of selected genes on embeddings
- Download gene expression summary statistics

---

##  Installation

Clone the repository:

```bash
git clone https://github.com/GeneticCodon/Quick_Anndata_Explorer.git
cd Quick_Anndata_Explorer

Install dependencies:

pip install -r requirements.txt
▶️ Usage

Run the Streamlit app:

streamlit run app.py

Then:

Open the browser (auto-launch)
Upload your .h5ad file
Explore interactively
 Demo Dataset (Recommended)

For best results, use a processed dataset with embeddings:

import scanpy as sc

adata = sc.datasets.pbmc68k_reduced()
adata.write("demo_pbmc68k_reduced.h5ad")

Upload this file in the app to see:

embedding plots
metadata exploration
richer outputs
📁 Project Structure
.
├── app.py               # Streamlit application
├── requirements.txt    # Dependencies
├── README.md
├── data/               # Optional input files
├── results/            # Optional outputs
├── src/                # CLI version 
🎯 Use Cases
Quick inspection of .h5ad datasets
Teaching and demonstrations
Data validation before analysis
Debugging metadata issues
Rapid exploration for collaborators
🧠 Why this tool?

Most workflows require loading full pipelines just to understand a dataset.

This tool provides:

Immediate insight into AnnData objects with zero overhead

🔧 Requirements
Python 3.8+
Streamlit
Scanpy / AnnData
Pandas
Matplotlib


Add a screenshot here after running the app for better GitHub visibility.

About Genetic Codon

Genetic Codon is a virtual bioinformatics lab focused on:

Multi-omics analysis
Spatial transcriptomics
Reproducible research pipelines
Training and mentorship

🌐 Website: www.geneticcodon.com

⭐ Support

If you find this useful:

⭐ Star the repository
🔁 Share with your network
🤝 Contribute improvements
