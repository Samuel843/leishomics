import streamlit as st
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from adjustText import adjust_text  

# Streamlit Configuration
st.set_page_config(layout="wide", page_title="Volcano Plot")

# File Upload
def load_data():
    uploaded_file = st.file_uploader("Upload your dataset (CSV, TSV, or TXT)", type=["csv", "tsv", "txt"])
    if uploaded_file:
        data = pd.read_csv(uploaded_file)
        data.columns = data.columns.str.strip()  # Clean column names
        return data
    return None

data = load_data()

if data is not None:
    st.title("Volcano Plot")

    # Let user select columns
    logfc_col = st.selectbox("Select Log2 Fold Change column", options=data.columns)
    pvalue_col = st.selectbox("Select p-value column", options=data.columns)
    gene_col = st.selectbox("Select Gene Identifier column", options=data.columns)

    if logfc_col and pvalue_col and gene_col:
        # Convert p-value column to numeric
        data[pvalue_col] = pd.to_numeric(data[pvalue_col], errors="coerce")
        data = data.dropna(subset=[pvalue_col, logfc_col])  # Remove NaNs
        data = data[data[pvalue_col] > 0]  # Ensure p-values are positive

        # Convert p-value to -log10 scale
        data['-log10(pvalue)'] = -np.log10(data[pvalue_col])

        # Thresholds (adjust these values based on your R plot)
        logfc_threshold = st.slider("Set Log2 Fold Change threshold", 0.0, max(abs(data[logfc_col].max()), abs(data[logfc_col].min())), 2.0, step=0.1)
        pvalue_threshold = st.slider("Set p-value threshold (-log10 scale)", 0.0, 10.0, 1.3, step=0.1)

        # Assign categories (match R color scheme)
        data['color_group'] = "Not significant"
        data.loc[(data[logfc_col] <= -logfc_threshold) & (data['-log10(pvalue)'] > pvalue_threshold), 'color_group'] = "Downregulated"
        data.loc[(data[logfc_col] >= logfc_threshold) & (data['-log10(pvalue)'] > pvalue_threshold), 'color_group'] = "Upregulated"

        # Define color mapping
        color_palette = {"Downregulated": "blue", "Upregulated": "red", "Not significant": "gray"}

        # Select Top Genes for Labeling
        top_genes = data.nlargest(30, '-log10(pvalue)')

        # Plot
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.scatterplot(data=data, x=logfc_col, y='-log10(pvalue)', hue='color_group',
                        palette=color_palette, alpha=0.7, ax=ax)

        # Add threshold lines (dashed)
        ax.axhline(y=pvalue_threshold, color='blue', linestyle='dashed', linewidth=1)
        ax.axvline(x=logfc_threshold, color='red', linestyle='dashed', linewidth=1)
        ax.axvline(x=-logfc_threshold, color='red', linestyle='dashed', linewidth=1)

        # Label top genes
        text_annotations = []
        for _, row in top_genes.iterrows():
            text_annotations.append(ax.text(row[logfc_col], row['-log10(pvalue)'], row[gene_col], fontsize=8, ha='right', color='black'))

        # Adjust text to prevent overlap
        adjust_text(text_annotations, arrowprops=dict(arrowstyle="-", color='black', lw=0.5))

        # Customize Plot
        ax.set_xlabel("Log2 Fold Change")
        ax.set_ylabel("-Log10(p-value)")
        ax.set_title("Volcano Plot BC")

        # Display legend only if needed
        if "Downregulated" in data['color_group'].unique() or "Upregulated" in data['color_group'].unique():
            ax.legend(title="color_group")

        st.pyplot(fig)

    else:
        st.error("Please select valid columns.")
else:
    st.write("Please upload a valid data file.")
