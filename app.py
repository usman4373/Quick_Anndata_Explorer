from io import BytesIO

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st


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

st.markdown(
    """
Upload an `.h5ad` file to inspect:

- dataset dimensions
- metadata columns
- embeddings such as UMAP/PCA
- value counts from `obs`
- gene expression on embeddings
- downloadable summaries
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

tab1, tab2, tab3, tab4, tab5, tab6, tab7 = st.tabs(
    ["Overview", "obs Metadata", "var Metadata", "Embeddings", "Gene Expression", "Counts", "Downloads"]
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

        embedding = adata.obsm[selected_embedding]

        if embedding.shape[1] < 2:
            st.warning("Selected embedding has fewer than 2 dimensions.")
        else:
            plot_df = pd.DataFrame(
                {
                    "dim1": embedding[:, 0],
                    "dim2": embedding[:, 1],
                }
            )

            color_options = ["None"] + list(adata.obs.columns) if len(adata.obs.columns) > 0 else ["None"]
            color_column = st.selectbox("Color by", color_options)

            fig, ax = plt.subplots(figsize=(7, 5))

            if color_column != "None":
                color_values = adata.obs[color_column].astype(str).values
                unique_vals = pd.Series(color_values).unique()

                if len(unique_vals) <= 20:
                    for val in unique_vals:
                        idx = color_values == val
                        ax.scatter(
                            plot_df.loc[idx, "dim1"],
                            plot_df.loc[idx, "dim2"],
                            s=point_size,
                            label=val,
                            alpha=0.75,
                        )
                    ax.legend(
                        bbox_to_anchor=(1.02, 1),
                        loc="upper left",
                        borderaxespad=0,
                        fontsize=8,
                        frameon=False,
                    )
                else:
                    ax.scatter(plot_df["dim1"], plot_df["dim2"], s=point_size, alpha=0.75)
                    st.warning(f"'{color_column}' has too many unique values for a readable legend.")
            else:
                ax.scatter(plot_df["dim1"], plot_df["dim2"], s=point_size, alpha=0.75)

            ax.set_xlabel("Dimension 1")
            ax.set_ylabel("Dimension 2")
            ax.set_title(f"{selected_embedding} preview")
            ax.grid(False)
            st.pyplot(fig)

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
                embedding = adata.obsm[selected_embedding_gene]

                if embedding.shape[1] < 2:
                    st.warning("Selected embedding has fewer than 2 dimensions.")
                else:
                    expr = extract_gene_expression(adata, gene_query)

                    expr_df = pd.DataFrame(
                        {
                            "dim1": embedding[:, 0],
                            "dim2": embedding[:, 1],
                            "expression": expr,
                        }
                    )

                    fig_expr, ax_expr = plt.subplots(figsize=(7, 5))
                    sc = ax_expr.scatter(
                        expr_df["dim1"],
                        expr_df["dim2"],
                        c=expr_df["expression"],
                        s=point_size,
                        alpha=0.8,
                    )
                    ax_expr.set_xlabel("Dimension 1")
                    ax_expr.set_ylabel("Dimension 2")
                    ax_expr.set_title(f"{gene_query} expression on {selected_embedding_gene}")
                    ax_expr.grid(False)
                    plt.colorbar(sc, ax=ax_expr, label="Expression")
                    st.pyplot(fig_expr)

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
                    st.dataframe(expr_summary, use_container_width=True)
                    dataframe_download(
                        expr_summary,
                        f"{gene_query}_expression_summary.csv",
                        f"Download {gene_query} expression summary",
                    )
            else:
                st.error("Gene not found in `adata.var_names`. Check the exact gene symbol in your dataset.")

with tab6:
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

with tab7:
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
