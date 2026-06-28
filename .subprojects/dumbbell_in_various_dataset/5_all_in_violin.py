# %%
import scanpy as sc
import numpy as np
import pandas as pd
import anndata
import seaborn as sns
import matplotlib.pyplot as plt
# %%
Amel = sc.read_h5ad("./data/Amel_dumbbell_marked.h5ad")
# %%
sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=[ "log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
Amel.obs["total_counts_scran"] = Amel.layers["scran_normalization"].sum(1)
Amel.obs["total_counts_log1p"] = Amel.layers["log1p_norm"].sum(1)
Amel.obs["total_counts_pearson"] = Amel.layers["analytic_pearson_residuals"].sum(1)
# %%
# %%
#! ---------------------------------------------- All in Violin ---------------------------------------------- 
#? log1p
sc.pl.embedding(Amel, 
                basis="X_umap_harmony_log1p", 
                color=["leiden_harmony_log1p_res0.50"],
                legend_loc="on data")
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
# sc.pl.violin(
#     Amel, 
#     keys='pct_counts_mt', 
#     groupby='leiden_harmony_log1p_res0.50', 
#     rotation=90, 
#     show=False
# )
# plt.tight_layout()
# plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='total_counts_log1p', 
    groupby='leiden_harmony_log1p_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()

# %%
#? scran
sc.pl.embedding(Amel, 
                basis="X_umap_harmony_scran", 
                color=["leiden_harmony_scran_res0.50"],
                legend_loc="on data")
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_scran_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_scran_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='total_counts_scran', 
    groupby='leiden_harmony_scran_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
#? pearson
sc.pl.embedding(Amel, 
                basis="X_umap_harmony_pearson", 
                color=["leiden_harmony_pearson_res0.50"],
                legend_loc="on data")
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_n_genes_by_counts', 
    groupby='leiden_harmony_pearson_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='log1p_total_counts', 
    groupby='leiden_harmony_pearson_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
sc.pl.violin(
    Amel, 
    keys='total_counts_pearson', 
    groupby='leiden_harmony_pearson_res0.50', 
    rotation=90, 
    show=False
)
plt.tight_layout()
plt.show()
# %%
