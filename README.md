# 🧬 Genetic Codon — AnnData Quick Explorer

A lightweight tool to quickly inspect `.h5ad` (AnnData) files and generate instant summaries with minimal setup.

---

## 🔍 Overview

Working with single-cell datasets often requires loading full analysis pipelines just to understand the structure of the data.

This tool provides a fast way to:
- Inspect dataset dimensions
- Explore metadata
- Preview structure and distributions
- Visualize embeddings (if available)

---

## ⚡ Features

- Display AnnData shape (cells × genes)
- List `obs` (cell metadata) columns
- List `var` (gene metadata) columns
- Show sample / cluster distributions
- Quick UMAP visualization (if present)
- Lightweight and fast execution

---

## 🚀 Installation

Clone the repository:

```bash
git clone https://github.com/GeneticCodon/Quick_Anndata_Explorer.git
cd Quick_Anndata_Explorer
