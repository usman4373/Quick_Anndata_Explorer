from io import StringIO

import anndata as ad
import pandas as pd
import streamlit as st
import matplotlib.pyplot as plt


st.set_page_config(
    page_title="Genetic Codon AnnData Quick Explorer",
    page_icon="🧬",
    layout="wide",
)

st.title("🧬 Genetic Codon — AnnData Quick Explorer")
st.write(
    "Upload an `.h5ad` file to quickly inspect its structure, metadata, embeddings, "
    "and selected cell-level distributions."
)


def dataframe_download(df: pd.DataFrame, filename: str, label: str) -> None:
    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        label=label,
        data=csv,
        file_name=filename,
        mime="text/csv",
    )


uploaded_file = st.file_uploader("Upload AnnData file", type=["h5ad"])

if uploaded_file is not None:
    try:
        adata = ad.read_h5ad(uploaded_file)

        st.success("File loaded successfully.")

        col1, col2, col3, col4 = st.columns(4)
        col1.metric("Cells (n_obs)", f"{adata.n_obs:,}")
        col2.metric("Genes (n_vars)", f"{adata.n_vars:,}")
        col3.metric("obs columns", len(adata.obs.columns))
        col4.metric("var columns", len(adata.var.columns))

        st.subheader("Dataset Summary")

        summary_lines = [
            f"n_obs: {adata.n_obs}",
            f"n_vars: {adata.n_vars}",
            f"obs columns: {list(adata.obs.columns)}",
            f"var columns: {list(adata.var.columns)}",
            f"obsm keys: {list(adata.obsm.keys())}",
            f"layers: {list(adata.layers.keys())}",
            f"uns keys: {list(adata.uns.keys())}",
            f"umap_available: {'X_umap' in adata.obsm.keys()}",
        ]
        st.code("\n".join(summary_lines), language="text")

        tab1, tab2, tab3, tab4, tab5 = st.tabs(
            ["obs", "var", "Embeddings", "Counts", "Downloads"]
        )

        with tab1:
            st.subheader("obs Metadata")
            if len(adata.obs.columns) > 0:
                obs_df = pd.DataFrame({"obs_column": adata.obs.columns})
                st.dataframe(obs_df, use_container_width=True)
                dataframe_download(obs_df, "obs_columns.csv", "Download obs columns")
                preview_n = min(20, adata.n_obs)
                st.write(f"Preview of first {preview_n} cells")
                st.dataframe(adata.obs.head(preview_n), use_container_width=True)
            else:
                st.info("No obs columns found in this AnnData object.")

        with tab2:
            st.subheader("var Metadata")
            if len(adata.var.columns) > 0:
                var_df = pd.DataFrame({"var_column": adata.var.columns})
                st.dataframe(var_df, use_container_width=True)
                dataframe_download(var_df, "var_columns.csv", "Download var columns")
                preview_n = min(20, adata.n_vars)
                st.write(f"Preview of first {preview_n} genes")
                st.dataframe(adata.var.head(preview_n), use_container_width=True)
            else:
                st.info("No var columns found in this AnnData object.")

        with tab3:
            st.subheader("Embeddings")
            obsm_keys = list(adata.obsm.keys())

            if len(obsm_keys) == 0:
                st.info("No embeddings found in `adata.obsm`.")
            else:
                st.write("Available embedding keys:")
                st.write(obsm_keys)

                selected_embedding = st.selectbox("Select embedding", obsm_keys)

                embedding = adata.obsm[selected_embedding]
                if embedding.shape[1] >= 2:
                    plot_df = pd.DataFrame(
                        {
                            "dim1": embedding[:, 0],
                            "dim2": embedding[:, 1],
                        }
                    )

                    color_column = None
                    if len(adata.obs.columns) > 0:
                        color_column = st.selectbox(
                            "Color by obs column (optional)",
                            ["None"] + list(adata.obs.columns),
                        )

                    fig, ax = plt.subplots(figsize=(7, 5))

                    if color_column and color_column != "None":
                        color_values = adata.obs[color_column].astype(str).values
                        unique_vals = pd.Series(color_values).unique()

                        if len(unique_vals) <= 20:
                            for val in unique_vals:
                                idx = color_values == val
                                ax.scatter(
                                    plot_df.loc[idx, "dim1"],
                                    plot_df.loc[idx, "dim2"],
                                    s=8,
                                    label=val,
                                    alpha=0.7,
                                )
                            ax.legend(
                                bbox_to_anchor=(1.02, 1),
                                loc="upper left",
                                borderaxespad=0,
                                fontsize=8,
                            )
                        else:
                            ax.scatter(plot_df["dim1"], plot_df["dim2"], s=8, alpha=0.7)
                            st.warning(
                                f"'{color_column}' has too many unique values for a readable legend."
                            )
                    else:
                        ax.scatter(plot_df["dim1"], plot_df["dim2"], s=8, alpha=0.7)

                    ax.set_xlabel("Dimension 1")
                    ax.set_ylabel("Dimension 2")
                    ax.set_title(f"{selected_embedding} preview")
                    st.pyplot(fig)
                else:
                    st.warning("Selected embedding has fewer than 2 dimensions.")

        with tab4:
            st.subheader("Value Counts from obs")
            if len(adata.obs.columns) > 0:
                selected_count_col = st.selectbox(
                    "Choose an obs column to summarize",
                    list(adata.obs.columns),
                )
                counts = (
                    adata.obs[selected_count_col]
                    .astype(str)
                    .value_counts(dropna=False)
                    .reset_index()
                )
                counts.columns = ["value", "count"]
                st.dataframe(counts, use_container_width=True)
                dataframe_download(
                    counts,
                    f"{selected_count_col}_value_counts.csv",
                    f"Download {selected_count_col} counts",
                )
            else:
                st.info("No obs columns available for counting.")

        with tab5:
            st.subheader("Download Summary Files")

            summary_text = "\n".join(summary_lines)
            st.download_button(
                label="Download summary.txt",
                data=summary_text,
                file_name="summary.txt",
                mime="text/plain",
            )

            obs_df = pd.DataFrame({"obs_column": adata.obs.columns})
            var_df = pd.DataFrame({"var_column": adata.var.columns})
            obsm_df = pd.DataFrame({"obsm_key": list(adata.obsm.keys())})
            layers_df = pd.DataFrame({"layer": list(adata.layers.keys())})

            dataframe_download(obs_df, "obs_columns.csv", "Download obs_columns.csv")
            dataframe_download(var_df, "var_columns.csv", "Download var_columns.csv")
            dataframe_download(obsm_df, "obsm_keys.csv", "Download obsm_keys.csv")
            dataframe_download(layers_df, "layers.csv", "Download layers.csv")

    except Exception as e:
        st.error(f"Failed to read file: {e}")

else:
    st.info("Upload a `.h5ad` file to begin.")
