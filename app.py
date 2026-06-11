from io import BytesIO
import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st
import os
import math
import gc


st.set_page_config(
    page_title="Genetic Codon Quick AnnData Explorer",
    page_icon="🧬",
    layout="wide",
)

st.title("Genetic Codon — Quick AnnData Explorer")
st.caption("Interactive exploration of AnnData (.h5ad) objects for rapid dataset inspection.")

st.sidebar.header("Settings")
point_size = st.sidebar.slider("Embedding point size", min_value=4, max_value=40, value=10, step=1)
preview_rows = st.sidebar.slider("Preview rows", min_value=5, max_value=50, value=15, step=5)
downsample_frac = st.sidebar.slider(
    "Downsample fraction for plots (reduce for very large datasets)",
    min_value=0.05, max_value=1.0, value=1.0, step=0.05,
    help="Randomly sample this fraction of cells to speed up plotting. Set to 1.0 to use all cells."
)

st.markdown(
    """

Upload an `.h5ad` file to inspect:

- **Dataset overview** – dimensions, metadata columns, embeddings, layers
- **obs / var metadata** – preview any columns
- **Embeddings** – plot UMAP, PCA, etc.  
  - **Color by** any categorical/continuous `obs` column  
  - **Split by** any `obs` column – view the embedding in separate facets  
  - **Save plots** as PNG (automatically saved to `./results/`)
- **Gene expression** – enter a gene name, visualise on any embedding, split by metadata, and download expression summary
- **Quality Control** – automatic calculation of **nCount, nFeature, mitochondrial & ribosomal fractions**  
  - Histograms and embedding plots coloured by QC metrics  
  - Download full QC table (CSV)
- **Value counts** – quick bar plots for categorical `obs` columns
- **Downloadable summaries** – TXT, CSV, or Excel exports
"""
)

def dataframe_download(df: pd.DataFrame, filename: str, label: str) -> None:
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=label,
        data=csv,
        file_name=filename,
        mime="text/csv",
    )


def text_download(text: str, filename: str, label: str) -> None:
    st.download_button(
        label=label,
        data=text.encode("utf-8"),
        file_name=filename,
        mime="text/plain",
    )


def build_summary_text(adata: ad.AnnData) -> str:
    lines = [
        "Genetic Codon AnnData Quick Explorer Summary",
        "=" * 45,
        f"n_obs: {adata.n_obs}",
        f"n_vars: {adata.n_vars}",
        f"obs columns: {list(adata.obs.columns)}",
        f"var columns: {list(adata.var.columns)}",
        f"obsm keys: {list(adata.obsm.keys())}",
        f"layers: {list(adata.layers.keys())}",
        f"uns keys: {list(adata.uns.keys())}",
        f"umap_available: {'X_umap' in adata.obsm.keys()}",
    ]
    return "\n".join(lines)


def dataframe_to_excel_bytes(sheets: dict[str, pd.DataFrame]) -> bytes:
    output = BytesIO()
    with pd.ExcelWriter(output, engine="openpyxl") as writer:
        for sheet_name, df in sheets.items():
            df.to_excel(writer, sheet_name=sheet_name[:31], index=False)
    return output.getvalue()


def extract_gene_expression(adata: ad.AnnData, gene_name: str) -> np.ndarray:
    gene_idx = adata.var_names.get_loc(gene_name)
    expr = adata.X[:, gene_idx]

    if hasattr(expr, "toarray"):
        expr = expr.toarray().flatten()
    else:
        expr = np.asarray(expr).flatten()

    return expr

# Compute QC metrices
def compute_qc_metrics_full(adata: ad.AnnData):
    """
    Compute nCount_RNA, nFeature_RNA, mitochondrial fraction, and ribosomal fraction.
    Returns: (n_counts, n_features, mito_frac, ribo_frac)
    """
    # nCount: total counts per cell
    total_counts = adata.X.sum(axis=1)
    if hasattr(total_counts, "toarray"):
        total_counts = total_counts.toarray().flatten()
    else:
        total_counts = np.asarray(total_counts).flatten()
    
    # nFeature: number of genes with non-zero count per cell
    n_features = (adata.X > 0).sum(axis=1)
    if hasattr(n_features, "toarray"):
        n_features = n_features.toarray().flatten()
    else:
        n_features = np.asarray(n_features).flatten()
    
    # Mitochondrial fraction (human MT- / mouse mt-)
    mito_genes = [g for g in adata.var_names if g.startswith(("MT-", "mt-"))]
    mito_sum = 0
    if mito_genes:
        mito_counts = adata[:, mito_genes].X.sum(axis=1)
        if hasattr(mito_counts, "toarray"):
            mito_sum = mito_counts.toarray().flatten()
        else:
            mito_sum = np.asarray(mito_counts).flatten()
    mito_frac = mito_sum / (total_counts + 1e-8)
    
    # Ribosomal fraction (human RPS/RPL, mouse Rps/Rpl)
    ribo_genes = [g for g in adata.var_names if g.startswith(("RPS", "RPL", "Rps", "Rpl"))]
    ribo_sum = 0
    if ribo_genes:
        ribo_counts = adata[:, ribo_genes].X.sum(axis=1)
        if hasattr(ribo_counts, "toarray"):
            ribo_sum = ribo_counts.toarray().flatten()
        else:
            ribo_sum = np.asarray(ribo_counts).flatten()
    ribo_frac = ribo_sum / (total_counts + 1e-8)
    
    return total_counts, n_features, mito_frac, ribo_frac

# Helper to save matplotlib figures to ./results/
def save_figure(fig, filename, dpi=100):
    """
    Save a matplotlib figure to the './results/' directory.
    Creates the directory if it does not exist.
    """
    save_dir = "results"
    os.makedirs(save_dir, exist_ok=True)
    filepath = os.path.join(save_dir, filename)
    fig.savefig(filepath, dpi=dpi, bbox_inches="tight")
    st.success(f"Figure saved as `{filepath}`")

def downsample_adata_for_plot(adata, fraction, random_state=42):
    """Return a subset of AnnData (views) for plotting only, without modifying original."""
    if fraction >= 1.0:
        return adata
    np.random.seed(random_state)
    n_cells = adata.n_obs
    n_subset = max(1, int(n_cells * fraction))
    subset_idx = np.random.choice(n_cells, size=n_subset, replace=False)
    subset_idx.sort()
    # Return a view (no copy) to save memory
    return adata[subset_idx, :]

def plot_split_embedding(embedding, categories, split_values, color_values, point_size, ncols=3):
    """
    Create a grid of embedding subplots, split by unique values in split_values.
    Each subplot shows the embedding, colored by the corresponding category.
    Returns the matplotlib figure.
    """
    unique_splits = sorted(set(split_values))
    n_splits = len(unique_splits)
    nrows = math.ceil(n_splits / ncols)
    
    fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols, 4 * nrows))
    # flatten axes for easy indexing
    if nrows == 1 and ncols == 1:
        axes = [axes]
    else:
        axes = axes.flatten()
    
    for idx, split_val in enumerate(unique_splits):
        ax = axes[idx]
        mask = (split_values == split_val)
        emb_sub = embedding[mask]
        cat_sub = color_values[mask]
        
        # Use the same coloring logic as before
        unique_cats = sorted(set(cat_sub))
        if len(unique_cats) <= 20:
            for cat in unique_cats:
                idx_cat = (cat_sub == cat)
                ax.scatter(emb_sub[idx_cat, 0], emb_sub[idx_cat, 1],
                           s=point_size, label=cat, alpha=0.75)
            ax.legend(bbox_to_anchor=(1.02, 1), loc='upper left',
                      fontsize=6, frameon=False)
        else:
            ax.scatter(emb_sub[:, 0], emb_sub[:, 1], s=point_size, alpha=0.75)
        
        ax.set_title(f"{split_val}")
        ax.set_xlabel("Dim1")
        ax.set_ylabel("Dim2")
        ax.grid(False)
    
    # Hide any unused subplots
    for idx in range(n_splits, len(axes)):
        axes[idx].axis('off')
    
    plt.tight_layout()
    return fig

uploaded_file = st.file_uploader("Upload AnnData file", type=["h5ad"])

if uploaded_file is None:
    st.info("Upload a `.h5ad` file to begin. For a polished demo, use `pbmc68k_reduced()` from Scanpy.")
    st.stop()

try:
    adata = ad.read_h5ad(uploaded_file)
except Exception as e:
    st.error(f"Failed to load AnnData file: {e}")
    st.stop()

st.success("File loaded successfully.")

# Top summary cards
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Cells", f"{adata.n_obs:,}")
c2.metric("Genes", f"{adata.n_vars:,}")
c3.metric("obs columns", len(adata.obs.columns))
c4.metric("var columns", len(adata.var.columns))
c5.metric("Embeddings", len(adata.obsm.keys()))

summary_text = build_summary_text(adata)

with st.expander("Dataset Summary", expanded=True):
    st.code(summary_text, language="text")

tab1, tab2, tab3, tab4, tab5, tab6, tab7, tab8 = st.tabs(
    ["Overview", "obs Metadata", "var Metadata", "Embeddings", "Gene Expression", "Quality Control", "Counts", "Downloads"]
)

with tab1:
    st.subheader("Overview")

    left, right = st.columns([1, 1])

    with left:
        st.markdown("**obs columns**")
        if len(adata.obs.columns) > 0:
            st.dataframe(pd.DataFrame({"obs_column": adata.obs.columns}), use_container_width=True)
        else:
            st.info("No obs columns available.")

        st.markdown("**obsm keys**")
        if len(adata.obsm.keys()) > 0:
            st.dataframe(pd.DataFrame({"obsm_key": list(adata.obsm.keys())}), use_container_width=True)
        else:
            st.info("No embeddings available.")

    with right:
        st.markdown("**var columns**")
        if len(adata.var.columns) > 0:
            st.dataframe(pd.DataFrame({"var_column": adata.var.columns}), use_container_width=True)
        else:
            st.info("No var columns available.")

        st.markdown("**layers**")
        if len(adata.layers.keys()) > 0:
            st.dataframe(pd.DataFrame({"layer": list(adata.layers.keys())}), use_container_width=True)
        else:
            st.info("No layers available.")

with tab2:
    st.subheader("obs Metadata")

    if len(adata.obs.columns) == 0:
        st.info("No obs metadata found in this AnnData object.")
    else:
        obs_columns_df = pd.DataFrame({"obs_column": adata.obs.columns})
        st.dataframe(obs_columns_df, use_container_width=True)

        selected_obs_cols = st.multiselect(
            "Select obs columns to preview",
            list(adata.obs.columns),
            default=list(adata.obs.columns[: min(5, len(adata.obs.columns))]),
        )

        if selected_obs_cols:
            st.dataframe(adata.obs[selected_obs_cols].head(preview_rows), use_container_width=True)

with tab3:
    st.subheader("var Metadata")

    if len(adata.var.columns) == 0:
        st.info("No var metadata found in this AnnData object.")
    else:
        var_columns_df = pd.DataFrame({"var_column": adata.var.columns})
        st.dataframe(var_columns_df, use_container_width=True)

        selected_var_cols = st.multiselect(
            "Select var columns to preview",
            list(adata.var.columns),
            default=list(adata.var.columns[: min(5, len(adata.var.columns))]),
        )

        if selected_var_cols:
            st.dataframe(adata.var[selected_var_cols].head(preview_rows), use_container_width=True)

with tab4:
    st.subheader("Embeddings")

    obsm_keys = list(adata.obsm.keys())

    if len(obsm_keys) == 0:
        st.info("No embeddings found in `adata.obsm`.")
    else:
        preferred = "X_umap" if "X_umap" in obsm_keys else obsm_keys[0]
        selected_embedding = st.selectbox("Select embedding", obsm_keys, index=obsm_keys.index(preferred))

        # Downsample if needed
        if downsample_frac < 1.0:
            adata_sub = downsample_adata_for_plot(adata, downsample_frac)
            embedding_data = adata_sub.obsm[selected_embedding][:, :2]
            obs_sub = adata_sub.obs
        else:
            embedding_data = adata.obsm[selected_embedding][:, :2]
            obs_sub = adata.obs

        if embedding_data.shape[1] < 2:
            st.warning("Selected embedding has fewer than 2 dimensions.")
        else:
            color_options = ["None"] + list(obs_sub.columns) if len(obs_sub.columns) > 0 else ["None"]
            color_column = st.selectbox("Color by", color_options)

            split_options = ["None"] + list(obs_sub.columns) if len(obs_sub.columns) > 0 else ["None"]
            split_column = st.selectbox("Split by (optional)", split_options, key="split_embed")

            if split_column == "None":
                # Single plot
                fig, ax = plt.subplots(figsize=(7, 5))

                if color_column != "None":
                    color_values = obs_sub[color_column].astype(str).values
                    unique_vals = pd.Series(color_values).unique()
                    if len(unique_vals) <= 20:
                        for val in unique_vals:
                            idx = (color_values == val)
                            ax.scatter(
                                embedding_data[idx, 0],
                                embedding_data[idx, 1],
                                s=point_size,
                                label=val,
                                alpha=0.75,
                            )
                        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0, fontsize=8, frameon=False)
                    else:
                        ax.scatter(embedding_data[:, 0], embedding_data[:, 1], s=point_size, alpha=0.75)
                        st.warning(f"'{color_column}' has too many unique values for a readable legend.")
                else:
                    ax.scatter(embedding_data[:, 0], embedding_data[:, 1], s=point_size, alpha=0.75)

                ax.set_xlabel("Dimension 1")
                ax.set_ylabel("Dimension 2")
                ax.set_title(f"{selected_embedding} preview")
                ax.grid(False)
                st.pyplot(fig)
                gc.collect()

                if color_column != "None":
                    safe_color = color_column.replace("/", "_").replace("\\", "_").strip()
                    embed_filename = f"{selected_embedding}_{safe_color}.png"
                else:
                    embed_filename = f"{selected_embedding}.png"
                st.button(
                    "Save plot",
                    on_click=save_figure,
                    args=(fig, embed_filename),
                    key=f"save_embed_{selected_embedding}_{color_column}"
                )
                gc.collect()
            else:
                # Split plot
                split_values = obs_sub[split_column].astype(str).values
                if color_column != "None":
                    color_values = obs_sub[color_column].astype(str).values
                else:
                    color_values = np.array(["None"] * len(obs_sub))
                fig = plot_split_embedding(embedding_data, color_column, split_values, color_values, point_size, ncols=3)
                st.pyplot(fig)
                gc.collect()

                safe_color = color_column.replace("/", "_").replace("\\", "_").strip() if color_column != "None" else "none"
                safe_split = split_column.replace("/", "_").replace("\\", "_").strip()
                embed_filename = f"{selected_embedding}_{safe_color}_split_by_{safe_split}.png"
                st.button(
                    "Save plot",
                    on_click=save_figure,
                    args=(fig, embed_filename),
                    key=f"save_embed_split_{selected_embedding}_{color_column}_{split_column}"
                )
                gc.collect()

with tab5:
    st.subheader("Gene Expression")

    obsm_keys = list(adata.obsm.keys())

    if len(obsm_keys) == 0:
        st.info("Gene expression plotting requires at least one embedding in `adata.obsm`.")
    else:
        preferred = "X_umap" if "X_umap" in obsm_keys else obsm_keys[0]
        selected_embedding_gene = st.selectbox(
            "Select embedding for gene expression",
            obsm_keys,
            index=obsm_keys.index(preferred),
        )

        gene_query = st.text_input("Enter gene symbol exactly as present in var_names", value="MS4A1")

        if gene_query:
            if gene_query in adata.var_names:
                # Downsample if needed
                if downsample_frac < 1.0:
                    adata_sub = downsample_adata_for_plot(adata, downsample_frac)
                    embedding = adata_sub.obsm[selected_embedding_gene]
                    if embedding.shape[1] < 2:
                        st.warning("Selected embedding has fewer than 2 dimensions.")
                        st.stop()
                    expr = extract_gene_expression(adata_sub, gene_query)
                    obs_sub = adata_sub.obs
                else:
                    embedding = adata.obsm[selected_embedding_gene]
                    if embedding.shape[1] < 2:
                        st.warning("Selected embedding has fewer than 2 dimensions.")
                        st.stop()
                    expr = extract_gene_expression(adata, gene_query)
                    obs_sub = adata.obs

                split_options = ["None"] + list(obs_sub.columns) if len(obs_sub.columns) > 0 else ["None"]
                split_column_expr = st.selectbox("Split by (optional)", split_options, key="split_expr")

                if split_column_expr == "None":
                    # Single plot
                    fig_expr, ax_expr = plt.subplots(figsize=(7, 5))
                    sc = ax_expr.scatter(
                        embedding[:, 0],
                        embedding[:, 1],
                        c=expr,
                        s=point_size,
                        alpha=0.8,
                    )
                    ax_expr.set_xlabel("Dimension 1")
                    ax_expr.set_ylabel("Dimension 2")
                    ax_expr.set_title(f"{gene_query} expression on {selected_embedding_gene}")
                    ax_expr.grid(False)
                    plt.colorbar(sc, ax=ax_expr, label="Expression")
                    st.pyplot(fig_expr)

                    safe_gene = gene_query.replace("/", "_").replace("\\", "_").strip()
                    expr_filename = f"{selected_embedding_gene}_{safe_gene}.png"
                    st.button(
                        "Save plot",
                        on_click=save_figure,
                        args=(fig_expr, expr_filename),
                        key=f"save_expr_{selected_embedding_gene}_{gene_query}"
                    )
                    gc.collect()
                else:
                    # Split plot
                    split_values = obs_sub[split_column_expr].astype(str).values
                    unique_splits = sorted(set(split_values))
                    n_splits = len(unique_splits)
                    nrows = math.ceil(n_splits / 3)
                    fig_expr, axes = plt.subplots(nrows, 3, figsize=(12, 4 * nrows))
                    if nrows == 1 and 3 == 1:
                        axes = [axes]
                    else:
                        axes = axes.flatten()

                    for idx, split_val in enumerate(unique_splits):
                        ax = axes[idx]
                        mask = (split_values == split_val)
                        emb_sub = embedding[mask]
                        expr_sub = expr[mask]
                        sc = ax.scatter(emb_sub[:, 0], emb_sub[:, 1],
                                        c=expr_sub, s=point_size, alpha=0.8)
                        ax.set_title(split_val)
                        ax.set_xlabel("Dim1")
                        ax.set_ylabel("Dim2")
                        ax.grid(False)
                        plt.colorbar(sc, ax=ax, label="Expression")
                    for idx in range(n_splits, len(axes)):
                        axes[idx].axis('off')
                    plt.tight_layout()
                    st.pyplot(fig_expr)

                    safe_gene = gene_query.replace("/", "_").replace("\\", "_").strip()
                    safe_split = split_column_expr.replace("/", "_").replace("\\", "_").strip()
                    expr_filename = f"{selected_embedding_gene}_{safe_gene}_split_by_{safe_split}.png"
                    st.button(
                        "Save plot",
                        on_click=save_figure,
                        args=(fig_expr, expr_filename),
                        key=f"save_expr_split_{selected_embedding_gene}_{gene_query}_{split_column_expr}"
                    )
                    gc.collect()

                st.markdown("**Expression summary**")
                expr_summary = pd.DataFrame(
                    {
                        "metric": ["min", "max", "mean", "median", "nonzero_cells"],
                        "value": [
                            float(np.min(expr)),
                            float(np.max(expr)),
                            float(np.mean(expr)),
                            float(np.median(expr)),
                            int(np.sum(expr > 0)),
                        ],
                    }
                )
                st.dataframe(expr_summary, width='stretch')
                dataframe_download(
                    expr_summary,
                    f"{gene_query}_expression_summary.csv",
                    f"Download {gene_query} expression summary",
                )
            else:
                st.error("Gene not found in `adata.var_names`. Check the exact gene symbol in your dataset.")

with tab6:
    st.subheader("Quality Control Metrics")
    
    # Compute QC metrics once and store in session state to avoid recomputation
    if "qc_metrics" not in st.session_state or st.session_state.qc_adata_n_obs != adata.n_obs:
        with st.spinner("Computing QC metrics (this may take a moment)..."):
            n_counts, n_features, mito_frac, ribo_frac = compute_qc_metrics_full(adata)
            st.session_state.qc_metrics = {
                "n_counts": n_counts,
                "n_features": n_features,
                "mito_frac": mito_frac,
                "ribo_frac": ribo_frac,
            }
            st.session_state.qc_adata_n_obs = adata.n_obs
    else:
        n_counts = st.session_state.qc_metrics["n_counts"]
        n_features = st.session_state.qc_metrics["n_features"]
        mito_frac = st.session_state.qc_metrics["mito_frac"]
        ribo_frac = st.session_state.qc_metrics["ribo_frac"]
    
    # Add to obs for easy reference (optional)
    # Consider whether you really need these in obs; if not, comment out.
    #adata.obs["nCount_RNA"] = n_counts
    #adata.obs["nFeature_RNA"] = n_features
    #adata.obs["percent_mito"] = mito_frac * 100
    #adata.obs["percent_ribo"] = ribo_frac * 100
    
    # Summary cards
    st.markdown("### Summary Statistics")
    col1, col2, col3, col4 = st.columns(4)
    col1.metric("Median nCount", f"{np.median(n_counts):.0f}")
    col2.metric("Median nFeature", f"{np.median(n_features):.0f}")
    col3.metric("Median % Mito", f"{np.median(mito_frac*100):.1f}%")
    col4.metric("Median % Ribo", f"{np.median(ribo_frac*100):.1f}%")
    
    # Distribution plots
    st.markdown("### Distributions (linear & log scales)")
    
    # nCount
    fig1, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(n_counts, bins=50, color='steelblue', alpha=0.7)
    axes[0].set_xlabel("nCount_RNA")
    axes[0].set_ylabel("Cells")
    axes[0].set_title("Linear scale")
    axes[1].hist(np.log10(n_counts+1), bins=50, color='steelblue', alpha=0.7)
    axes[1].set_xlabel("log10(nCount_RNA+1)")
    axes[1].set_title("Log scale")
    st.pyplot(fig1)

    st.button("Save nCount plot", on_click=save_figure, args=(fig1, "nCountRNA.png"), key="save_ncount")
    
    # nFeature
    fig2, axes = plt.subplots(1, 2, figsize=(12, 4))
    axes[0].hist(n_features, bins=50, color='seagreen', alpha=0.7)
    axes[0].set_xlabel("nFeature_RNA")
    axes[0].set_ylabel("Cells")
    axes[0].set_title("Linear scale")
    axes[1].hist(np.log10(n_features+1), bins=50, color='seagreen', alpha=0.7)
    axes[1].set_xlabel("log10(nFeature_RNA+1)")
    axes[1].set_title("Log scale")
    st.pyplot(fig2)

    st.button("Save nFeature plot", on_click=save_figure, args=(fig2, "nFeatureRNA.png"), key="save_nfeature")
    
    # Mito and Ribo distributions (side by side)
    colA, colB = st.columns(2)
    with colA:
        fig3, ax = plt.subplots()
        ax.hist(mito_frac*100, bins=50, color='red', alpha=0.7)
        ax.set_xlabel("% Mitochondrial reads")
        ax.set_ylabel("Cells")
        st.pyplot(fig3)
        st.button("Save % Mitochondrial plot", on_click=save_figure, args=(fig3, "mt_percent.png"), key="save_mito")
    with colB:
        fig4, ax = plt.subplots()
        ax.hist(ribo_frac*100, bins=50, color='green', alpha=0.7)
        ax.set_xlabel("% Ribosomal reads")
        ax.set_ylabel("Cells")
        st.pyplot(fig4)
        st.button("Save % Ribosomal plot", on_click=save_figure, args=(fig4, "rb_percent.png"), key="save_ribo")
    
    # Embedding visualizations (if available)
    if len(adata.obsm.keys()) > 0:
        st.markdown("### QC metrics on embedding")
        obsm_keys = list(adata.obsm.keys())
        preferred = "X_umap" if "X_umap" in obsm_keys else obsm_keys[0]
        embed_qc = st.selectbox("Embedding", obsm_keys, index=obsm_keys.index(preferred), key="qc_embed")
        
        # Downsample QC data if needed
        if downsample_frac < 1.0:
            adata_sub = downsample_adata_for_plot(adata, downsample_frac)
            emb = adata_sub.obsm[embed_qc][:, :2]
            n_counts_sub = compute_qc_metrics_full(adata_sub)[0]
            n_features_sub = compute_qc_metrics_full(adata_sub)[1]
            mito_frac_sub = compute_qc_metrics_full(adata_sub)[2]
            ribo_frac_sub = compute_qc_metrics_full(adata_sub)[3]
            obs_sub = adata_sub.obs
        else:
            emb = adata.obsm[embed_qc][:, :2]
            n_counts_sub = n_counts
            n_features_sub = n_features
            mito_frac_sub = mito_frac
            ribo_frac_sub = ribo_frac
            obs_sub = adata.obs
        
        metric_choice = st.selectbox("Color by", ["nCount_RNA", "nFeature_RNA", "% Mitochondrial", "% Ribosomal"], key="qc_metric")
        
        split_options_qc = ["None"] + list(obs_sub.columns) if len(obs_sub.columns) > 0 else ["None"]
        split_column_qc = st.selectbox("Split by (optional)", split_options_qc, key="split_qc_embed")
        
        if metric_choice == "nCount_RNA":
            color_vals = n_counts_sub
            color_label = "UMI count"
            cmap = "viridis"
        elif metric_choice == "nFeature_RNA":
            color_vals = n_features_sub
            color_label = "Genes detected"
            cmap = "plasma"
        elif metric_choice == "% Mitochondrial":
            color_vals = mito_frac_sub * 100
            color_label = "% Mito"
            cmap = "Reds"
        else:
            color_vals = ribo_frac_sub * 100
            color_label = "% Ribo"
            cmap = "Greens"
        
        if split_column_qc == "None":
            fig5, ax = plt.subplots(figsize=(8, 6))
            sc = ax.scatter(emb[:,0], emb[:,1], c=color_vals, s=point_size, alpha=0.7, cmap=cmap)
            plt.colorbar(sc, label=color_label)
            ax.set_title(f"{metric_choice} on {embed_qc}")
            st.pyplot(fig5)
            
            metric_map = {
                "nCount_RNA": "nCount_RNA",
                "nFeature_RNA": "nFeature_RNA",
                "% Mitochondrial": "percent_mito",
                "% Ribosomal": "percent_ribo"
            }
            safe_metric = metric_map.get(metric_choice, metric_choice.replace(" ", "_"))
            qc_embed_filename = f"{embed_qc}_{safe_metric}.png"
            st.button(
                "Save QC embedding plot",
                on_click=save_figure,
                args=(fig5, qc_embed_filename),
                key="save_qc_embed"
            )
        else:
            split_values = obs_sub[split_column_qc].astype(str).values
            unique_splits = sorted(set(split_values))
            n_splits = len(unique_splits)
            nrows = math.ceil(n_splits / 3)
            fig5, axes = plt.subplots(nrows, 3, figsize=(12, 4 * nrows))
            if nrows == 1 and 3 == 1:
                axes = [axes]
            else:
                axes = axes.flatten()
            
            for idx, split_val in enumerate(unique_splits):
                ax = axes[idx]
                mask = (split_values == split_val)
                emb_sub = emb[mask]
                vals_sub = color_vals[mask]
                sc = ax.scatter(emb_sub[:,0], emb_sub[:,1], c=vals_sub, s=point_size, alpha=0.7, cmap=cmap)
                ax.set_title(split_val)
                ax.set_xlabel("Dim1")
                ax.set_ylabel("Dim2")
                ax.grid(False)
                plt.colorbar(sc, ax=ax, label=color_label)
            for idx in range(n_splits, len(axes)):
                axes[idx].axis('off')
            plt.tight_layout()
            st.pyplot(fig5)
            
            metric_map = {
                "nCount_RNA": "nCount_RNA",
                "nFeature_RNA": "nFeature_RNA",
                "% Mitochondrial": "percent_mito",
                "% Ribosomal": "percent_ribo"
            }
            safe_metric = metric_map.get(metric_choice, metric_choice.replace(" ", "_"))
            safe_split = split_column_qc.replace("/", "_").replace("\\", "_").strip()
            qc_embed_filename = f"{embed_qc}_{safe_metric}_split_by_{safe_split}.png"
            st.button(
                "Save QC embedding plot",
                on_click=save_figure,
                args=(fig5, qc_embed_filename),
                key=f"save_qc_embed_split_{embed_qc}_{metric_choice}_{split_column_qc}"
            )
    
    # Downloadable QC table
    qc_df = pd.DataFrame({
        "nCount_RNA": n_counts,
        "nFeature_RNA": n_features,
        "percent_mito": mito_frac * 100,
        "percent_ribo": ribo_frac * 100,
    })
    dataframe_download(qc_df, "qc_metrics_full.csv", "Download full QC metrics (CSV)")

with tab7:
    st.subheader("Value Counts from obs")

    if len(adata.obs.columns) == 0:
        st.info("No obs columns available for counting.")
    else:
        selected_count_col = st.selectbox(
            "Choose an obs column",
            list(adata.obs.columns),
        )

        counts = (
            adata.obs[selected_count_col]
            .astype(str)
            .value_counts(dropna=False)
            .reset_index()
        )
        counts.columns = ["value", "count"]

        left, right = st.columns([1.1, 1])

        with left:
            st.dataframe(counts, use_container_width=True)

        with right:
            if len(counts) <= 30:
                fig_counts, ax_counts = plt.subplots(figsize=(7, 4))
                ax_counts.bar(counts["value"], counts["count"])
                ax_counts.set_title(f"Counts for {selected_count_col}")
                ax_counts.set_xlabel("Value")
                ax_counts.set_ylabel("Count")
                ax_counts.tick_params(axis="x", rotation=45)
                ax_counts.grid(False)
                st.pyplot(fig_counts)
            else:
                st.info("Too many categories to display as a readable bar plot.")

with tab8:
    st.subheader("Downloads")

    obs_df = pd.DataFrame({"obs_column": adata.obs.columns})
    var_df = pd.DataFrame({"var_column": adata.var.columns})
    obsm_df = pd.DataFrame({"obsm_key": list(adata.obsm.keys())})
    layers_df = pd.DataFrame({"layer": list(adata.layers.keys())})

    text_download(summary_text, "summary.txt", "Download summary.txt")
    dataframe_download(obs_df, "obs_columns.csv", "Download obs_columns.csv")
    dataframe_download(var_df, "var_columns.csv", "Download var_columns.csv")
    dataframe_download(obsm_df, "obsm_keys.csv", "Download obsm_keys.csv")
    dataframe_download(layers_df, "layers.csv", "Download layers.csv")

    excel_bytes = dataframe_to_excel_bytes(
        {
            "obs_columns": obs_df,
            "var_columns": var_df,
            "obsm_keys": obsm_df,
            "layers": layers_df,
        }
    )

    st.download_button(
        label="Download all summaries as Excel",
        data=excel_bytes,
        file_name="anndata_quick_explorer_summary.xlsx",
        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
    )
