# %%
import scanpy as sc
import seaborn as sns
import matplotlib.pyplot as plt
import anndata
import numpy as np
# %%
Hsal = sc.read_h5ad("../Sheng_SA_2020_Hsal/data/7_cluster-output/concat.h5ad")
Acer = sc.read_h5ad("../Acer/data/7_cluster-output/concat.h5ad")
Amel = sc.read_h5ad("../Zhang_iScience_2022_Amel/data/7_cluster-output/concat.h5ad")
# %%
# # Hsal
# sc.pl.embedding(Hsal, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# sc.pl.embedding(Hsal, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# sc.pl.embedding(Hsal, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# # %%
# # Acer
# sc.pl.embedding(Acer, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# sc.pl.embedding(Acer, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# sc.pl.embedding(Acer, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# # %%
# # Amel
# sc.pl.embedding(Amel, basis="X_umap_harmony_log1p", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# sc.pl.embedding(Amel, basis="X_umap_harmony_scran", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# sc.pl.embedding(Amel, basis="X_umap_harmony_pearson", color=["log1p_total_counts", "log1p_n_genes_by_counts", "pct_counts_mt", "pct_counts_in_top_20_genes"])
# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(Amel.obs["total_counts"], bins=100, kde=False, ax=axes[0])
sns.histplot(Amel.layers["analytic_pearson_residuals"].sum(1), bins=100, kde=False, ax=axes[1])
# fig
# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(Amel.obs["total_counts"], bins=100, kde=False, ax=axes[0])
sns.histplot(Amel.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
# fig
# %%
fig, axes = plt.subplots(1, 2, figsize=(10, 5))
sns.histplot(Amel.obs["total_counts"], bins=100, kde=False, ax=axes[0])
sns.histplot(Amel.layers["scran_normalization"].sum(1), bins=100, kde=False, ax=axes[1])
# %%
def visual_counts(
    adata: anndata.Anndata
) -> None:
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    # 1. 原始 total_counts (adata.obs 中的 Series 本身是一维的，无需处理)
    sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0, 0], color="skyblue")
    axes[0, 0].set_title("Raw total_counts (Library Size)", fontsize=12)
    axes[0, 0].set_xlabel("UMI Counts")
    
    # 2. log1p_norm (使用 np.ravel 展平二维矩阵)
    sns.histplot(np.ravel(adata.layers["log1p_norm"].sum(axis=1)), bins=100, kde=False, ax=axes[0, 1], color="lightgreen")
    axes[0, 1].set_title("log1p_norm sum", fontsize=12)
    
    # 3. analytic_pearson_residual (使用 np.ravel 展平)
    sns.histplot(np.ravel(adata.layers["analytic_pearson_residuals"].sum(axis=1)), bins=100, kde=False, ax=axes[1, 0], color="salmon")
    axes[1, 0].set_title("Pearson Residuals sum", fontsize=12)
    
    # 4. scran_normalization (使用 np.ravel 展平)
    sns.histplot(np.ravel(adata.layers["scran_normalization"].sum(axis=1)), bins=100, kde=False, ax=axes[1, 1], color="plum")
    axes[1, 1].set_title("Scran sum", fontsize=12)
    
    # 自动调整子图间距，防止标题重叠
    plt.tight_layout()
    plt.show()
# %%
visual_counts(Amel)
# %%
visual_counts(Acer)
# %%
visual_counts(Hsal)
# %%
